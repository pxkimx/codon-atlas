"""Per-site lists for a decoding hypothesis (codon Z read as amino acid A) and protein databases.

The database package uses <codon>2<aa> naming (e.g. TGA2W) and a layout that Recode Detector's
Database check reads directly:
  <stem>_standard.fasta, <stem>_<TAG>.fasta, <stem>_<TAG>_variants_only.fasta,
  <stem>_combined.fasta, <stem>_<TAG>_mapping.tsv, <stem>_<TAG>_affected_proteins.tsv, README.md
with headers  >record|protein_id|<TAG or CANONICAL> record=.. locus_tag=.. gene=.. product=..
"""
from __future__ import annotations

import io
import json
import re
import zipfile

import numpy as np

from .gencode import AA20, CODON_INDEX, THREE, table, translate
from .genome import Genome, canonical_protein


def tag(codon: str, aa: str) -> str:
    return f"{codon.upper()}2{aa.upper()}"


def site_list(genome: Genome, codon: str, aa: str, obs: np.ndarray | None = None, prof=None,
              table_id: int = 11) -> list[dict]:
    codon = codon.upper().replace("U", "T")
    ref = table(table_id)
    lookup = {}
    if obs is not None and len(obs):
        z = CODON_INDEX[codon]
        sel = obs[obs["codon"] == z]
        for r in sel:
            lookup[(int(r["contig"]), int(r["strand"]), int(r["pos"]))] = r
    out = []
    for gi, g in enumerate(genome.genes):
        seq = genome.contigs[g.contig].seq
        nt = genome.cds_nt(g)
        cmap = genome.codon_map(g)
        n_codons = len(nt) // 3
        for ai, pos in cmap:
            if ai >= n_codons:
                break
            if nt[3 * ai:3 * ai + 3] != codon:
                continue
            if ai == n_codons - 1 and ref.get(codon) == "*":
                continue  # the terminal stop of the gene
            prot_std = ref.get(codon, "X")
            row = {"protein_id": g.protein_id, "locus_tag": g.locus_tag, "gene": g.gene, "product": g.product,
                   "record": genome.contigs[g.contig].name, "strand": "+" if g.strand > 0 else "-",
                   "genome_pos": pos + 1, "aa_pos": ai + 1, "protein_length": n_codons - (1 if nt[-3:] and ref.get(nt[-3:]) == "*" else 0),
                   "codon": codon, "standard": prot_std, "variant": aa, "gene_i": gi,
                   "context": _context(nt, ai, ref, aa), "pfam": "", "column": None,
                   "p_standard": None, "p_variant": None, "column_consensus": "", "log10_variant_vs_standard": None}
            r = lookup.get((g.contig, g.strand, pos))
            if r is not None and prof is not None:
                hi = int(r["hmm"])
                e = np.asarray(prof.emissions[prof.offsets[hi] + int(r["k"]) - 1], dtype=float)
                ps = e[AA20.index(prot_std)] if prot_std in AA20 else None
                pv = e[AA20.index(aa)] if aa in AA20 else None
                row.update(pfam=prof.names[hi], column=int(r["k"]), column_consensus=AA20[int(np.argmax(e))],
                           p_standard=None if ps is None else float(ps), p_variant=None if pv is None else float(pv),
                           log10_variant_vs_standard=None if (ps is None or pv is None) else float(np.log10(max(pv, 1e-30) / max(ps, 1e-30))))
            out.append(row)
    return out


def _context(nt: str, ai: int, ref: dict, aa: str, flank: int = 7) -> str:
    a = max(0, ai - flank)
    left = translate(nt[3 * a:3 * ai], ref)
    right = translate(nt[3 * ai + 3:3 * (ai + flank + 1)], ref)
    return f"{left}[{aa}]{right}"


def database_package(genome: Genome, hypotheses: list[tuple[str, str]], table_id: int = 11,
                     stem: str | None = None, obs=None, prof=None, note: str = "") -> tuple[str, bytes]:
    """ZIP with standard, variant, variants-only and combined FASTA plus mapping tables and README."""
    hyp = [(c.upper().replace("U", "T"), a.upper()) for c, a in hypotheses]
    T = "_".join(tag(c, a) for c, a in hyp)
    stem = stem or re.sub(r"[^A-Za-z0-9_.-]+", "_", (genome.organism or genome.name).strip()) or "genome"
    if genome.contigs and genome.contigs[0].name and genome.contigs[0].name not in stem:
        stem = f"{stem}_{genome.contigs[0].name.split('.')[0]}"
    ref = table(table_id)
    sites_by_gene: dict[int, list[dict]] = {}
    all_sites = []
    for c, a in hyp:
        for s in site_list(genome, c, a, obs, prof, table_id):
            sites_by_gene.setdefault(s["gene_i"], []).append(s)
            all_sites.append(s)
    std_fa, var_fa, only_fa, comb_fa = io.StringIO(), io.StringIO(), io.StringIO(), io.StringIO()
    mapping = io.StringIO()
    cols = ["record", "protein_id", "locus_tag", "gene", "product", "strand", "CDS_location", "protein_length_aa",
            f"{T}_count", f"{T}_amino_acid_positions_1_based", f"canonical_residues_at_{T}", "codons"]
    mapping.write("\t".join(cols) + "\n")
    affected = io.StringIO()
    affected.write("\t".join(cols) + "\n")
    n_aff = 0
    for gi, g in enumerate(genome.genes):
        prot = canonical_protein(genome, g, table_id)
        rec = genome.contigs[g.contig].name
        pid = g.protein_id or g.locus_tag or f"gene{gi + 1}"
        desc = f"record={rec} locus_tag={g.locus_tag or 'NA'} gene={g.gene or 'NA'} product={g.product or 'NA'}"
        std_fa.write(f">{rec}|{pid}|CANONICAL {desc}\n{_wrap(prot)}")
        comb_fa.write(f">{rec}|{pid}|CANONICAL {desc}\n{_wrap(prot)}")
        sites = sorted(sites_by_gene.get(gi, []), key=lambda s: s["aa_pos"])
        var = list(prot)
        canon = []
        for s in sites:
            i = s["aa_pos"] - 1
            if i < len(var):
                canon.append(var[i])
                var[i] = s["variant"]
        var = "".join(var)
        line = [rec, pid, g.locus_tag, g.gene, g.product, "1" if g.strand > 0 else "-1", g.location_text(), str(len(prot)),
                str(len(sites)), ",".join(str(s["aa_pos"]) for s in sites) or "none", ",".join(canon) or "none",
                ",".join(s["codon"] for s in sites) or "none"]
        mapping.write("\t".join(line) + "\n")
        if sites:
            n_aff += 1
            affected.write("\t".join(line) + "\n")
            pos = ",".join(str(s["aa_pos"]) for s in sites)
            vh = f">{rec}|{pid}|{T} {desc} {T}_sites={pos}\n{_wrap(var)}"
            var_fa.write(vh)
            only_fa.write(vh)
            comb_fa.write(vh)
        else:
            var_fa.write(f">{rec}|{pid}|{T} {desc} {T}_sites=none\n{_wrap(prot)}")
    hyp_txt = ", ".join(f"{c} ({THREE.get(ref.get(c, 'X'), '?')}) → {THREE.get(a, a)}" for c, a in hyp)
    readme = f"""# {genome.organism or genome.name}: {T} protein databases

Made by Codon Atlas for the decoding hypothesis {hyp_txt}.

- Source file: {genome.name} (SHA-256 {genome.sha256[:16]}…)
- Records: {len(genome.contigs)}, {genome.length:,} bp
- Genes: {len(genome.genes):,} ({'annotated CDS from the GenBank record' if genome.gene_source == 'annotation' else 'predicted with Prodigal (pyrodigal)'})
- Reference translation: NCBI table {table_id}
- Proteins with at least one site: {n_aff}
- Sites: {len(all_sites)}
{note}
## Files

- `{stem}_standard.fasta`: every protein, standard translation.
- `{stem}_{T}.fasta`: every protein translated under the hypothesis (unchanged where no site).
- `{stem}_{T}_variants_only.fasta`: only the altered proteins.
- `{stem}_combined.fasta`: standard proteins plus the altered variants (unique identifiers),
  ready for a single database search where both versions compete.
- `{stem}_{T}_mapping.tsv`: all proteins with site counts and 1-based positions.
- `{stem}_{T}_affected_proteins.tsv`: the altered proteins only.
- `{stem}_{T}_sites.tsv`: one row per site, with sequence context and, where the site lies in a
  Pfam domain, the domain column's probabilities for the standard and the variant residue.

Only in-frame codons of annotated/predicted genes are changed; the final stop codon is never
changed. Positions are 1-based.
"""
    sites_tsv = io.StringIO()
    keys = ["record", "protein_id", "locus_tag", "gene", "product", "strand", "genome_pos", "aa_pos", "codon",
            "standard", "variant", "context", "pfam", "column", "column_consensus", "p_standard", "p_variant",
            "log10_variant_vs_standard"]
    sites_tsv.write("\t".join(keys) + "\n")
    for s in sorted(all_sites, key=lambda s: (s["record"], s["genome_pos"])):
        sites_tsv.write("\t".join("" if s.get(k) is None else (f"{s[k]:.4g}" if isinstance(s[k], float) else str(s[k])) for k in keys) + "\n")
    buf = io.BytesIO()
    folder = f"{stem}_{T}_database_package"
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{folder}/README.md", readme)
        z.writestr(f"{folder}/{stem}_standard.fasta", std_fa.getvalue())
        z.writestr(f"{folder}/{stem}_{T}.fasta", var_fa.getvalue())
        z.writestr(f"{folder}/{stem}_{T}_variants_only.fasta", only_fa.getvalue())
        z.writestr(f"{folder}/{stem}_combined.fasta", comb_fa.getvalue())
        z.writestr(f"{folder}/{stem}_{T}_mapping.tsv", mapping.getvalue())
        z.writestr(f"{folder}/{stem}_{T}_affected_proteins.tsv", affected.getvalue())
        z.writestr(f"{folder}/{stem}_{T}_sites.tsv", sites_tsv.getvalue())
        z.writestr(f"{folder}/validation.json", json.dumps({
            "source": genome.name, "sha256": genome.sha256, "genes": len(genome.genes), "gene_source": genome.gene_source,
            "hypotheses": [{"codon": c, "standard": ref.get(c), "variant": a} for c, a in hyp],
            "altered_proteins": n_aff, "sites": len(all_sites), "table": table_id}, indent=2))
    return f"{folder}.zip", buf.getvalue()


def _wrap(s: str, w: int = 60) -> str:
    return "\n".join(s[i:i + w] for i in range(0, len(s), w)) + "\n"
