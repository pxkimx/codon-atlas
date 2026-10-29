"""Genetic-code inference from nucleotide sequence (the Codetta method, Shulgina & Eddy 2021).

1. Six-frame translation of the genome with the standard code, stop codons as 'X'.
2. Search Pfam profile HMMs against the translation (HMMER via pyhmmer); keep domains with
   independent E-value below the cutoff, and aligned consensus columns with posterior
   probability >= 0.95.
3. Each codon position keeps the single consensus column from its most significant domain.
4. For each codon Z, the aligned columns C_1..C_N give
       log P(C | aa=a) - log P(C | ?) = sum_i [ log e_{C_i}(a) - log pbar(a) ]
   where e is the column's emission probability and pbar(a) its average over every aligned
   column in this genome. The non-specific model '?' scores 0. With a uniform prior the
   decoding probabilities are a softmax over the 21 models; an amino acid is called when its
   probability exceeds the threshold (default 0.9999).
   No single column may contribute more than 1% of a codon's columns (each unique column
   once if the codon has fewer than 100).
"""
from __future__ import annotations

import os
import time
from dataclasses import asdict, dataclass, field

import numpy as np

from .gencode import AA20, CODON_INDEX, CODONS, MODELS, STOPS_STANDARD, revcomp, table
from .genome import Genome
from .pfamdb import CLASSES, ProfileSet

_CODE = {c: ("X" if a == "*" else a) for c, a in table(1).items()}


@dataclass
class Settings:
    profiles: str = "enone"
    evalue: float = 1e-10
    min_pp: float = 0.95
    exclude: list = field(default_factory=lambda: list(CLASSES))
    threshold: float = 0.9999
    column_cap: float = 0.01
    reference_table: int = 11
    fragment: int = 100_000
    overlap: int = 3_000
    cpus: int = 0
    bootstrap: int = 200

    @classmethod
    def from_dict(cls, d: dict | None) -> "Settings":
        s = cls()
        for k, v in (d or {}).items():
            if hasattr(s, k) and v is not None:
                cur = getattr(s, k)
                setattr(s, k, type(cur)(v) if isinstance(cur, (int, float)) and not isinstance(cur, bool) else v)
        return s


# ------------------------------------------------------------------ six-frame translation
@dataclass
class Fragment:
    contig: int
    start: int
    end: int
    strand: int
    frame: int
    codons: np.ndarray  # uint8 codon index per residue, 255 = ambiguous

    def pos(self, i: int) -> int:
        """Genomic coordinate of the 5' base (coding strand) of residue i (0-based)."""
        if self.strand > 0:
            return self.start + self.frame + 3 * i
        return self.end - 1 - (self.frame + 3 * i)


def six_frames(genome: Genome, size: int, overlap: int):
    frags: list[Fragment] = []
    prots: list[str] = []
    for ci, c in enumerate(genome.contigs):
        L = len(c.seq)
        step = max(3, size - overlap)
        s = 0
        while s < L:
            e = min(L, s + size)
            piece = c.seq[s:e]
            for strand, nt in ((1, piece), (-1, revcomp(piece))):
                for f in range(3):
                    sub = nt[f:]
                    n = len(sub) // 3
                    if n < 10:
                        continue
                    cods = [sub[3 * i:3 * i + 3] for i in range(n)]
                    prots.append("".join(_CODE.get(x, "X") for x in cods))
                    frags.append(Fragment(ci, s, e, strand, f,
                                          np.fromiter((CODON_INDEX.get(x, 255) for x in cods), dtype=np.uint8, count=n)))
            if e >= L:
                break
            s += step
    return frags, prots


# ------------------------------------------------------------------ search
OBS_DTYPE = np.dtype([("contig", np.int32), ("strand", np.int8), ("pos", np.int64), ("codon", np.uint8),
                      ("hmm", np.int32), ("k", np.int32), ("evalue", np.float64)])


def search(genome: Genome, prof: ProfileSet, st: Settings, progress=None, should_stop=None) -> np.ndarray:
    """Align profiles to the six-frame translation. Returns one row per codon position (best domain)."""
    import pyhmmer
    from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
    from pyhmmer.plan7 import HMMFile

    aa = Alphabet.amino()
    if progress:
        progress("Translating the genome in six frames", 0.0)
    frags, prots = six_frames(genome, st.fragment, st.overlap)
    block = DigitalSequenceBlock(aa, [TextSequence(name=str(i).encode(), sequence=p).digitize(aa)
                                      for i, p in enumerate(prots)])
    n_res = sum(len(p) for p in prots)
    total = len(prof)
    done = [0]
    t0 = time.time()

    def cb(query, n):
        done[0] += 1
        if should_stop and should_stop():
            raise InterruptedError("stopped")
        if progress and done[0] % 100 == 0:
            frac = done[0] / max(total, 1)
            el = time.time() - t0
            eta = el / frac - el if frac > 0.02 else None
            progress(f"Searching {total:,} Pfam profiles against {n_res:,} translated residues", 0.02 + 0.88 * frac,
                     eta=eta)

    best: dict[tuple, tuple] = {}
    cpus = st.cpus or max(1, (os.cpu_count() or 2) - 1)
    pp_ok = _pp_chars(st.min_pp)
    with HMMFile(prof.profiles_path()) as hf:
        profiles = hf.optimized_profiles()
        for hits in pyhmmer.hmmer.hmmsearch(profiles, block, cpus=cpus, callback=cb, Z=float(total),
                                            E=1e-3, domE=st.evalue):
            qn = hits.query.name
            hi = prof.index.get(qn.decode() if isinstance(qn, bytes) else str(qn))
            if hi is None:
                continue
            for hit in hits:
                fr = frags[int(hit.name)]
                for dom in hit.domains:
                    ev = dom.i_evalue
                    if ev >= st.evalue:
                        continue
                    al = dom.alignment
                    k = al.hmm_from
                    t = al.target_from - 1  # 0-based residue index
                    hs, ts, pp = al.hmm_sequence, al.target_sequence, al.posterior_probabilities
                    for col in range(len(hs)):
                        h = hs[col]
                        if h == ".":
                            t += 1
                            continue
                        if ts[col] == "-":
                            k += 1
                            continue
                        if pp[col] in pp_ok:
                            cod = fr.codons[t]
                            if cod != 255:
                                key = (fr.contig, fr.strand, fr.pos(t))
                                old = best.get(key)
                                if old is None or ev < old[2]:
                                    best[key] = (int(cod), hi, ev, k)
                        k += 1
                        t += 1
    if progress:
        progress("Collecting aligned columns", 0.91)
    obs = np.empty(len(best), dtype=OBS_DTYPE)
    for i, ((c, s, p), (cod, hi, ev, k)) in enumerate(best.items()):
        obs[i] = (c, s, p, cod, hi, k, ev)
    obs.sort(order=["contig", "pos", "strand"])
    return obs


def _pp_chars(min_pp: float) -> set[str]:
    # HMMER posterior characters: '*' >= 0.95, digit d means [d/10 - 0.05, d/10 + 0.05)
    ok = {"*"}
    for d in range(10):
        if d / 10 - 0.05 >= min_pp - 1e-9:
            ok.add(str(d))
    return ok


# ------------------------------------------------------------------ inference
def column_emissions(prof: ProfileSet, hmm: np.ndarray, k: np.ndarray) -> np.ndarray:
    rows = prof.offsets[hmm] + (k - 1)
    uniq, inv = np.unique(rows, return_inverse=True)
    em = np.asarray(prof.emissions[uniq], dtype=np.float64)
    return np.clip(em[inv], 1e-30, 1.0)


def _capped_weights(col_ids: np.ndarray, cap_frac: float) -> tuple[np.ndarray, np.ndarray]:
    """Unique columns and their capped counts (Codetta's 1% rule)."""
    uniq, counts = np.unique(col_ids, return_counts=True)
    n = counts.sum()
    if n < 100:
        w = np.ones_like(counts)
    else:
        cap = max(1, int(cap_frac * n))
        w = np.minimum(counts, cap)
    return uniq, w


def infer(obs: np.ndarray, prof: ProfileSet, st: Settings, seed: int = 7) -> dict:
    excl = prof.excluded_mask(st.exclude)
    use = obs[~excl[obs["hmm"]]] if len(obs) else obs
    ref = table(st.reference_table)
    rng = np.random.default_rng(seed)
    if len(use) == 0:
        pbar = np.full(20, 0.05)
    else:
        pbar = column_emissions(prof, use["hmm"], use["k"]).mean(axis=0)
    log_pbar = np.log(pbar)
    col_id = use["hmm"].astype(np.int64) * 100_000 + use["k"]
    codons = []
    for z, cod in enumerate(CODONS):
        sel = use["codon"] == z
        n_raw = int(sel.sum())
        row = {"codon": cod, "standard": ref[cod], "n_raw": n_raw, "n_used": 0, "n_unique": 0}
        if n_raw == 0:
            post = np.zeros(21)
            post[20] = 1.0
            row.update(_call(post, ref[cod], st.threshold))
            row.update(posterior=post.tolist(), loglik=[0.0] * 21, support=None, consensus={}, families=[])
            codons.append(row)
            continue
        uniq, w = _capped_weights(col_id[sel], st.column_cap)
        hm, kk = (uniq // 100_000).astype(np.int64), (uniq % 100_000).astype(np.int64)
        em = column_emissions(prof, hm, kk)
        logem = np.log(em) - log_pbar  # per column log-odds
        ll = np.concatenate([(w[:, None] * logem).sum(axis=0), [0.0]])
        post = _softmax(ll)
        call = _call(post, ref[cod], st.threshold)
        # bootstrap: resample the capped column set
        support = None
        if st.bootstrap and len(w) > 0:
            p = w / w.sum()
            same = 0
            for _ in range(st.bootstrap):
                draw = rng.multinomial(int(w.sum()), p)
                llb = np.concatenate([(draw[:, None] * logem).sum(axis=0), [0.0]])
                same += int(np.argmax(llb) == np.argmax(ll))
            support = same / st.bootstrap
        cons = np.argmax(em, axis=1)
        cons_counts = np.bincount(cons, weights=w, minlength=20)
        fam_counts = np.bincount(hm, weights=w)
        top = np.argsort(fam_counts)[::-1][:8]
        row.update(call)
        row.update(n_used=int(w.sum()), n_unique=int(len(uniq)), posterior=post.tolist(), loglik=ll.tolist(),
                   support=support,
                   consensus={AA20[i]: float(cons_counts[i]) for i in np.argsort(cons_counts)[::-1][:6] if cons_counts[i] > 0},
                   families=[{"name": prof.names[i], "acc": prof.accessions[i], "desc": prof.descriptions[i],
                              "columns": int(fam_counts[i])} for i in top if fam_counts[i] > 0])
        codons.append(row)
    return {"codons": codons, "pbar": dict(zip(AA20, pbar.tolist())),
            "n_columns": int(len(use)), "n_columns_all": int(len(obs)),
            "n_families": int(len(np.unique(use["hmm"]))) if len(use) else 0,
            "code": "".join(c["inferred"] for c in codons)}


def _softmax(ll: np.ndarray) -> np.ndarray:
    m = ll.max()
    e = np.exp(ll - m)
    return e / e.sum()


def _call(post: np.ndarray, std: str, thr: float) -> dict:
    best = int(np.argmax(post))
    model = MODELS[best]
    inferred = model if (best < 20 and post[best] >= thr) else "?"
    std_i = AA20.find(std) if std in AA20 else -1
    if std_i >= 0:
        # how far the evidence favours the best alternative over the standard amino acid (log10)
        alt = max((i for i in range(20) if i != std_i), key=lambda i: post[i])
        lbf = float((np.log(max(post[alt], 1e-300)) - np.log(max(post[std_i], 1e-300))) / np.log(10))
    else:
        alt, lbf = (best if best < 20 else None), None
    if inferred == "?":
        status = "stop" if std == "*" else "uninferred"
    elif std == "*":
        status = "stop-to-sense"
    elif inferred == std:
        status = "standard"
    else:
        status = "reassigned"
    return {"inferred": inferred, "best": model, "best_prob": float(post[best]), "status": status,
            "standard_prob": float(post[std_i]) if std_i >= 0 else None,
            "alternative": AA20[alt] if alt is not None and alt < 20 else None,
            "log10_alt_vs_standard": lbf}


# ------------------------------------------------------------------ stop codons inside domains
def stop_events(obs: np.ndarray, prof: ProfileSet, genome: Genome, gene_index=None) -> list[dict]:
    stop_ids = [CODON_INDEX[c] for c in STOPS_STANDARD]
    sel = np.isin(obs["codon"], stop_ids)
    rows = obs[sel]
    if len(rows) == 0:
        return []
    em = column_emissions(prof, rows["hmm"].astype(np.int64), rows["k"].astype(np.int64))
    out = []
    for r, e in zip(rows, em):
        order = np.argsort(e)[::-1]
        hi = int(r["hmm"])
        classes = prof.classes[hi]
        cod = CODONS[int(r["codon"])]
        g = gene_index.find(int(r["contig"]), int(r["strand"]), int(r["pos"])) if gene_index else None
        out.append({"contig": genome.contigs[int(r["contig"])].name, "contig_i": int(r["contig"]),
                    "pos": int(r["pos"]) + 1, "strand": "+" if r["strand"] > 0 else "-", "codon": cod,
                    "family": prof.names[hi], "acc": prof.accessions[hi], "desc": prof.descriptions[hi],
                    "column": int(r["k"]), "expected": AA20[order[0]], "expected_p": float(e[order[0]]),
                    "second": AA20[order[1]], "second_p": float(e[order[1]]),
                    "evalue": float(r["evalue"]), "classes": classes,
                    "gene": g, "kind": _stop_kind(cod, classes, e)})
    return out


def _stop_kind(cod: str, classes: list[str], e: np.ndarray) -> str:
    if cod == "TGA" and "seleno" in classes:
        return "Selenocysteine candidate"
    if cod == "TAG" and "pyl" in classes:
        return "Pyrrolysine candidate"
    if "mobile" in classes or "viral" in classes:
        return "Mobile element / phage (often pseudogenes)"
    if e[AA20.index("C")] > 0.5 and cod == "TGA":
        return "UGA at a conserved Cys column (Sec?)"
    return "In-domain stop (readthrough, frameshift, pseudogene or error)"


# ------------------------------------------------------------------ gene lookup
class GeneIndex:
    """Find the annotated/predicted gene (and residue position) covering a codon."""

    def __init__(self, genome: Genome):
        self.genome = genome
        self.by_key: dict[tuple, tuple[int, int]] = {}
        for gi, g in enumerate(genome.genes):
            for ai, pos in genome.codon_map(g):
                self.by_key[(g.contig, g.strand, pos)] = (gi, ai)

    def find(self, contig: int, strand: int, pos: int):
        hit = self.by_key.get((contig, strand, pos))
        if not hit:
            return None
        g = self.genome.genes[hit[0]]
        return {"id": g.ident, "locus_tag": g.locus_tag, "product": g.product, "aa_pos": hit[1] + 1,
                "gene_i": hit[0]}


def settings_dict(st: Settings) -> dict:
    return asdict(st)
