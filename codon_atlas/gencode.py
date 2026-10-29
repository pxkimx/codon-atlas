"""Codons, amino acids and NCBI translation tables."""
from __future__ import annotations

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]  # NCBI order, TTT .. GGG
CODON_INDEX = {c: i for i, c in enumerate(CODONS)}
AA20 = "ACDEFGHIKLMNPQRSTVWY"  # HMMER alphabet order
AA_INDEX = {a: i for i, a in enumerate(AA20)}
MODELS = list(AA20) + ["?"]
STOPS_STANDARD = ("TAA", "TAG", "TGA")

THREE = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His",
         "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln",
         "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr", "*": "Stop",
         "?": "?", "U": "Sec", "O": "Pyl", "X": "X"}

# NCBI tables as the 64-letter string in TCAG order (from the NCBI gc.prt file).
TABLES = {
    1: ("Standard", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    4: ("Mold, protozoan, Mycoplasma/Spiroplasma", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    11: ("Bacterial, archaeal and plant plastid", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    25: ("Candidate division SR1 and Gracilibacteria", "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    6: ("Ciliate, dasycladacean and Hexamita nuclear", "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    12: ("Alternative yeast nuclear", "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    2: ("Vertebrate mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
}
START_CODONS = {1: {"ATG", "CTG", "TTG"}, 11: {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"},
                4: {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}, 25: {"ATG", "GTG", "TTG"}}

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def table(n: int = 11) -> dict[str, str]:
    return dict(zip(CODONS, TABLES[n][1]))


STANDARD = table(1)


def translate(nt: str, code: dict[str, str] | None = None, stop_as: str = "*") -> str:
    code = code or STANDARD
    out = []
    for i in range(0, len(nt) - 2, 3):
        aa = code.get(nt[i:i + 3].upper(), "X")
        out.append(stop_as if aa == "*" else aa)
    return "".join(out)


def codon_label(c: str) -> str:
    return c.replace("T", "U")


def gene_table_for(inferred: dict[str, str]) -> int:
    """Pick the NCBI table to call genes with, given inferred stop-codon meanings."""
    tga = inferred.get("TGA", "?")
    if tga == "W":
        return 4
    if tga == "G":
        return 25
    return 11
