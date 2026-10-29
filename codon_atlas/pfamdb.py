"""The profile database: Pfam profile HMMs, pressed for fast search, plus a compact emission table.

Two profile sets can be prepared:
  * "enone": rebuilt from the Pfam seed alignments with entropy weighting turned off
             (hmmbuild --enone), as Codetta does. Emission probabilities stay close to the
             amino-acid frequencies in the seed alignments, which sharpens decoding evidence.
  * "pfam":  Pfam-A.hmm exactly as distributed (entropy-weighted). Faster to set up; gives the
             same direction of evidence but weaker per column, so fewer codons pass threshold.

Families are tagged by keyword into the five classes Codetta leaves out of inference by
default (mobile elements, viral, selenoprotein, pyrrolysine, mitochondrial-encoded).
"""
from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import time
from pathlib import Path

import numpy as np

HOME = Path(os.environ.get("CA_HOME", Path.home() / "Library" / "Application Support" / "Codon Atlas"))
PFAM_DIR = HOME / "pfam"
SETS = {"enone": "Pfam-A seed, rebuilt without entropy weighting (Codetta default)",
        "pfam": "Pfam-A.hmm as distributed"}

CLASSES = {
    "mobile": "Transposons and other mobile elements",
    "viral": "Viral and phage proteins",
    "seleno": "Selenoproteins (Sec at UGA)",
    "pyl": "Pyrrolysine proteins (Pyl at UAG)",
    "mito": "Mitochondrion-encoded proteins",
}
_RULES = {
    "mobile": re.compile(r"transpos|integrase|\brve\b|retrotrans|reverse transcriptase|\bRVT|insertion (element|sequence)|"
                         r"\bIS\d|\bDDE\b|resolvase|group II intron|maturase|HTH_Tnp|Tnp_|plasmid recombination|"
                         r"conjugal transfer|relaxase|\bMob_|homing endonuclease|LAGLIDADG|GIY-YIG", re.I),
    "viral": re.compile(r"phage|virus|viral|virion|capsid|prophage|\bterminase|portal protein|baseplate|tail fib|"
                        r"tail tube|tail sheath|tail tape|Gag|Env_|polyprotein|coat protein|nucleocapsid", re.I),
    "seleno": re.compile(r"seleno|\bSel[A-Z]?_|\bSelP|\bSelR\b|\bSelW|Sep15|SelK|SelM|SelT|\bRdx\b|\bGRDA\b|\bGRDB\b|"
                         r"glycine reductase|iodothyronine deiodinase|\bGSHPx\b|glutathione peroxidase|\bDIO\b", re.I),
    "pyl": re.compile(r"pyrrolysine|\bMtmB\b|\bMttB\b|\bMtbB\b|monomethylamine methyltransferase|"
                      r"trimethylamine methyltransferase|dimethylamine methyltransferase", re.I),
    "mito": re.compile(r"mitochondria(l)? (encoded|genome)|\bCOX[123]\b|^Cytochrom_B_[NC]|\bCytochrome_B\b|^ATP-synt_8|"
                       r"^Oxidored_q[1-6]$|^NADHdh|^Proton_antipo_N|^NADH_dehy_S2_C|_mt$|\bmtDNA", re.I),
}


def classify(name: str, desc: str) -> list[str]:
    text = f"{name} {desc}"
    out = []
    for k, rx in _RULES.items():
        if rx.search(text) or (k == "mito" and rx.search(name)):
            out.append(k)
    return out


def set_dir(which: str) -> Path:
    return PFAM_DIR / which


def status() -> dict:
    out = {"home": str(HOME), "sets": {}}
    for which, label in SETS.items():
        d = set_dir(which)
        idx = d / "index.json"
        info = {"label": label, "ready": False}
        if idx.exists() and (d / "profiles.h3m").exists() and (d / "emissions.npy").exists():
            meta = json.loads(idx.read_text())
            info.update(ready=True, families=len(meta["names"]), built=meta.get("built"),
                        release=meta.get("release", ""), size_mb=round(sum(f.stat().st_size for f in d.iterdir()) / 1e6))
        out["sets"][which] = info
    out["sources"] = {
        "hmm": _exists("Pfam-A.hmm"),
        "seed": _exists("Pfam-A.seed"),
        "release": _release(),
    }
    out["default"] = "enone" if out["sets"]["enone"]["ready"] else ("pfam" if out["sets"]["pfam"]["ready"] else None)
    return out


def _exists(base: str) -> dict:
    p, gz = PFAM_DIR / base, PFAM_DIR / (base + ".gz")
    if p.exists():
        return {"path": str(p), "mb": round(p.stat().st_size / 1e6)}
    if gz.exists():
        return {"path": str(gz), "mb": round(gz.stat().st_size / 1e6), "gz": True}
    return {}


def _release() -> str:
    v = PFAM_DIR / "Pfam.version"
    if v.exists():
        m = re.search(r"Pfam release\s*:\s*(\S+)", v.read_text())
        if m:
            return m.group(1)
    return ""


def _ungz(base: str, progress=None) -> Path:
    p, gz = PFAM_DIR / base, PFAM_DIR / (base + ".gz")
    if p.exists():
        return p
    if not gz.exists():
        raise FileNotFoundError(f"{base} is not in {PFAM_DIR}. Download it from the Database page first.")
    if progress:
        progress(f"Unpacking {gz.name}", 0.0)
    tmp = p.with_suffix(p.suffix + ".part")
    with gzip.open(gz, "rb") as fi, open(tmp, "wb") as fo:
        shutil.copyfileobj(fi, fo, 16 << 20)
    tmp.rename(p)
    return p


def build_enone_hmms(progress=None, should_stop=None) -> Path:
    """Rebuild every Pfam family from its seed alignment with --enone. Returns the .hmm path."""
    import pyhmmer
    from pyhmmer.easel import Alphabet, MSAFile
    from pyhmmer.plan7 import Background, Builder

    seed = _ungz("Pfam-A.seed", progress)
    out = set_dir("enone")
    out.mkdir(parents=True, exist_ok=True)
    target = out / "Pfam-A_enone.hmm"
    tmp = target.with_suffix(".hmm.part")
    aa = Alphabet.amino()
    builder = Builder(aa, effective_number="none")
    bg = Background(aa)
    total = _count_families(seed)
    n = 0
    with MSAFile(str(seed), digital=True, alphabet=aa, format="stockholm") as msas, open(tmp, "wb") as fo:
        for msa in msas:
            if should_stop and should_stop():
                raise InterruptedError("stopped")
            name = msa.name or f"family{n}".encode()
            msa.name = name
            hmm, _, _ = builder.build_msa(msa, bg)
            hmm.name = name
            if msa.accession:
                hmm.accession = msa.accession
            if msa.description:
                hmm.description = msa.description
            hmm.write(fo)
            n += 1
            if progress and n % 200 == 0:
                progress(f"Building profiles from seed alignments ({n:,} of {total:,})", n / max(total, 1))
    tmp.rename(target)
    return target


def _count_families(seed: Path) -> int:
    n = 0
    with open(seed, "rb") as f:
        for line in f:
            if line.startswith(b"# STOCKHOLM"):
                n += 1
    return n


def prepare(which: str, progress=None, should_stop=None) -> dict:
    """Make a profile set ready: press it for searching and extract the emission table."""
    from pyhmmer.plan7 import HMMFile
    import pyhmmer

    d = set_dir(which)
    d.mkdir(parents=True, exist_ok=True)
    if which == "enone":
        src = d / "Pfam-A_enone.hmm"
        if not src.exists():
            src = build_enone_hmms(progress, should_stop)
    else:
        src = _ungz("Pfam-A.hmm", progress)
    names, accs, descs, lengths, offsets, classes = [], [], [], [], [], []
    rows = []
    off = 0
    if progress:
        progress("Reading and pressing profiles", 0.0)
    for f in d.glob("profiles.h3*"):
        f.unlink()

    def stream():  # one pass: collect emissions while hmmpress writes the search files
        nonlocal off
        with HMMFile(str(src)) as hf:
            for hmm in hf:
                if should_stop and should_stop():
                    raise InterruptedError("stopped")
                nm, ds = _s(hmm.name), _s(hmm.description)
                em = np.asarray(hmm.match_emissions, dtype=np.float32)[1:, :20]
                rows.append(em.copy())
                names.append(nm)
                accs.append(_s(hmm.accession))
                descs.append(ds)
                lengths.append(int(hmm.M))
                offsets.append(off)
                classes.append(classify(nm, ds))
                off += em.shape[0]
                if progress and len(names) % 500 == 0:
                    progress(f"Reading and pressing profiles ({len(names):,})", min(0.95, len(names) / 26000))
                yield hmm

    pyhmmer.hmmer.hmmpress(stream(), str(d / "profiles"))
    np.save(d / "emissions.npy", np.concatenate(rows).astype(np.float32))
    meta = {"names": names, "accessions": accs, "descriptions": descs, "lengths": lengths,
            "offsets": offsets, "classes": classes, "set": which, "release": _release(),
            "built": time.strftime("%Y-%m-%d %H:%M")}
    (d / "index.json").write_text(json.dumps(meta))
    # the pressed files hold everything the search needs; drop the large unpacked text copies
    if which == "enone":
        src.unlink(missing_ok=True)
        if (PFAM_DIR / "Pfam-A.seed.gz").exists():
            (PFAM_DIR / "Pfam-A.seed").unlink(missing_ok=True)
    elif (PFAM_DIR / "Pfam-A.hmm.gz").exists():
        src.unlink(missing_ok=True)
    if progress:
        progress("Ready", 1.0)
    return status()["sets"][which]


def _s(v) -> str:
    if v is None:
        return ""
    return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)


class ProfileSet:
    """A prepared profile set, loaded for one analysis."""

    def __init__(self, which: str):
        d = set_dir(which)
        if not (d / "index.json").exists():
            raise FileNotFoundError(f"The '{which}' profile set is not prepared. Open the Database page.")
        self.which = which
        self.dir = d
        meta = json.loads((d / "index.json").read_text())
        self.meta = meta
        self.names = meta["names"]
        self.index = {n: i for i, n in enumerate(self.names)}
        self.accessions = meta["accessions"]
        self.descriptions = meta["descriptions"]
        self.lengths = np.asarray(meta["lengths"], dtype=np.int32)
        self.offsets = np.asarray(meta["offsets"], dtype=np.int64)
        self.classes = meta["classes"]
        self.release = meta.get("release", "")
        self.emissions = np.load(d / "emissions.npy", mmap_mode="r")

    def __len__(self) -> int:
        return len(self.names)

    def excluded_mask(self, exclude: list[str]) -> np.ndarray:
        ex = set(exclude)
        return np.array([bool(ex.intersection(c)) for c in self.classes], dtype=bool)

    def class_counts(self) -> dict[str, int]:
        out = {k: 0 for k in CLASSES}
        for c in self.classes:
            for k in c:
                out[k] += 1
        return out

    def profiles_path(self) -> str:
        return str(self.dir / "profiles")  # base name of the pressed .h3m/.h3f/.h3p/.h3i files
