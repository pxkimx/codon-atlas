"""tRNA genes (ARAGORN via pyaragorn), which codons each anticodon can read, and identity elements.

Identity calls follow the rules used by Shulgina & Eddy (2021) for arginine-codon reassignments,
applied to tRNAs with arginine anticodons (and Gly-type to UCA anticodons):
  Arg-type  : A20 in the D loop and A/G73 discriminator
  Trp-type  : G73 and not A20 (provisional; A/G1:U72 supports)
  Met-type  : A73 and not A20
  Gln-type  : not A20, A37, A/G73, weak 1:72 pair
  Gly(UCA)  : G1:C72, C2:G71, G3:C70 and U73
Positions 20 and 37 are located by a simple structural heuristic (D-stem pairing, G18-G19),
so these are leads to check on the folded structure, not final assignments.
"""
from __future__ import annotations

from .gencode import CODONS, revcomp, table

PAIR = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("G", "T"), ("T", "G")}
WC = {"A": "T", "T": "A", "G": "C", "C": "G"}

# wobble base (anticodon 34, as written in the gene) -> third codon bases it can read
WOBBLE = {
    "G": (("C", "Watson-Crick"), ("T", "wobble G·U")),
    "C": (("G", "Watson-Crick"),),
    "A": (("T", "inosine I34"), ("C", "inosine I34"), ("A", "inosine I34")),  # A34 is edited to I in bacteria
    "T": (("A", "Watson-Crick"), ("G", "wobble U·G")),
}


def reads(anticodon: str) -> list[tuple[str, str]]:
    """Codons an anticodon (5'->3', DNA letters) can pair with, and how."""
    ac = anticodon.upper().replace("U", "T")
    if len(ac) != 3 or any(b not in "ACGT" for b in ac):
        return []
    first, second = WC[ac[2]], WC[ac[1]]
    return [(first + second + third, how) for third, how in WOBBLE[ac[0]]]


def _discriminator(seq: str) -> str:
    s = seq.upper()
    if s.endswith("CCA") and len(s) > 70:
        return s[-4]
    return s[-1]


def _n20(seq: str) -> str:
    """Base at Sprinzl position 20: the base after G18-G19 in the D loop (heuristic)."""
    s = seq.upper()
    if len(s) < 60:
        return "?"
    stem5 = s[9:13]  # positions 10-13
    best = None
    for j in range(19, 30):  # 0-based start of 3' D-stem (positions 22-25)
        part = s[j:j + 4]
        if len(part) < 4:
            break
        pairs = sum((a, b) in PAIR for a, b in zip(stem5, part[::-1]))
        if pairs >= 3 and (best is None or pairs > best[1]):
            best = (j, pairs)
    if not best:
        return "?"
    loop = s[13:best[0]]
    i = loop.find("GG", 2)
    if i < 0 or i + 2 >= len(loop):
        return "?"
    return loop[i + 2]


def identity(seq: str, ac_offset: int, amino_acid: str = "", anticodon: str = "") -> dict:
    s = seq.upper()
    d73 = _discriminator(s)
    n20 = _n20(s)
    a0 = ac_offset - 1
    n37 = s[a0 + 3] if 0 <= a0 + 3 < len(s) else "?"
    body = s[:-3] if s.endswith("CCA") else s
    p1, p72 = body[0], body[-2] if len(body) > 70 else "?"
    p2, p71 = body[1], body[-3]
    p3, p70 = body[2], body[-4]
    weak_1_72 = (p1, p72) not in {("G", "C"), ("C", "G")}
    calls = []
    arg_family = amino_acid == "Arg" or anticodon in ("CCG", "TCG", "ACG", "GCG", "CCT", "TCT")
    if not arg_family:
        if anticodon == "TCA" and (p1, p72, p2, p71, p3, p70) == ("G", "C", "C", "G", "G", "C") and d73 == "T":
            calls.append("Gly-type (G1:C72 C2:G71 G3:C70, U73)")
        return {"n73": d73, "n20": n20, "n37": n37, "pair_1_72": f"{p1}:{p72}", "identity": calls}
    if n20 == "A" and d73 in "AG":
        calls.append("Arg-type (A20, A/G73)")
    if d73 == "G" and n20 not in ("A", "?"):
        calls.append("Trp-type (G73, no A20)" + (" + A/G1:U72" if p1 in "AG" and p72 == "T" else ""))
    if d73 == "A" and n20 not in ("A", "?"):
        calls.append("Met-type (A73, no A20)")
    if n20 not in ("A", "?") and n37 == "A" and d73 in "AG" and weak_1_72:
        calls.append("Gln-type (no A20, A37, A/G73, weak 1:72)")
    if not calls:
        calls.append("unclassified (check structure)")
    return {"n73": d73, "n20": n20, "n37": n37, "pair_1_72": f"{p1}:{p72}", "identity": calls}


def find_trnas(genome, table_id: int = 11) -> list[dict]:
    import pyaragorn

    finder = pyaragorn.RNAFinder(translation_table=table_id, tmrna=False, trna=True)
    out = []
    for ci, c in enumerate(genome.contigs):
        for t in finder.find_rna(c.seq.encode()):
            if t.type != "tRNA" and not hasattr(t, "anticodon"):
                continue
            ac = (t.anticodon or "").upper().replace("U", "T")
            seq = t.sequence().upper()
            ident = identity(seq, t.anticodon_offset, t.amino_acid, ac)
            out.append({"contig": c.name, "begin": int(t.begin), "end": int(t.end),
                        "strand": "+" if t.strand > 0 else "-", "amino_acid": t.amino_acid,
                        "anticodon": ac, "reads": [{"codon": cd, "how": how} for cd, how in reads(ac)],
                        "score": round(float(t.energy), 1), "length": int(t.length), "sequence": seq, **ident})
    return out


def codon_readers(trnas: list[dict], reference_table: int = 11) -> dict[str, list[dict]]:
    """For each codon, the tRNAs whose anticodon can pair with it."""
    ref = table(reference_table)
    by = {c: [] for c in CODONS}
    for i, t in enumerate(trnas):
        for r in t["reads"]:
            by[r["codon"]].append({"i": i, "anticodon": t["anticodon"], "amino_acid": t["amino_acid"],
                                   "how": r["how"], "identity": t["identity"], "n73": t["n73"], "n20": t["n20"]})
    return by


def summary_for_codon(codon: str, readers: list[dict], std_aa: str) -> str:
    if not readers:
        return "no tRNA with a pairing anticodon found"
    wc = [r for r in readers if r["how"] == "Watson-Crick"]
    names = sorted({f"tRNA-{r['amino_acid']}({r['anticodon']})" for r in (wc or readers)})
    return ("Watson-Crick: " if wc else "wobble only: ") + ", ".join(names)


def anticodon_for(codon: str) -> str:
    return revcomp(codon)
