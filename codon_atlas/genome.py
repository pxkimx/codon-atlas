"""Read genomes (FASTA or GenBank, optionally gzipped) and their protein-coding genes.

GenBank CDS features are used when present; otherwise genes are predicted with Prodigal
(pyrodigal). Coordinates are 0-based, end-exclusive internally.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

from .gencode import revcomp, table, translate

GENOME_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".ffn", ".gb", ".gbk", ".gbff", ".genbank", ".embl")


@dataclass
class Contig:
    name: str
    seq: str  # upper case DNA
    description: str = ""


@dataclass
class Gene:
    contig: int
    parts: list[tuple[int, int]]  # exons in biological order, 0-based [start, end)
    strand: int
    locus_tag: str = ""
    protein_id: str = ""
    gene: str = ""
    product: str = ""
    codon_start: int = 1
    translation: str = ""  # curated translation from the GenBank record, if any
    source: str = "annotation"
    partial: bool = False
    pseudo: bool = False

    @property
    def start(self) -> int:
        return min(p[0] for p in self.parts)

    @property
    def end(self) -> int:
        return max(p[1] for p in self.parts)

    @property
    def ident(self) -> str:
        return self.protein_id or self.locus_tag or f"gene_{self.contig}_{self.start + 1}"

    def location_text(self) -> str:
        s = "+" if self.strand > 0 else "-"
        parts = ",".join(f"[{a}:{b}]" for a, b in sorted(self.parts))
        return f"{parts}({s})"


@dataclass
class Genome:
    name: str
    contigs: list[Contig]
    genes: list[Gene] = field(default_factory=list)
    gene_source: str = "none"  # annotation | prodigal | none
    source_format: str = "fasta"
    sha256: str = ""
    organism: str = ""
    transl_table: int | None = None

    @property
    def length(self) -> int:
        return sum(len(c.seq) for c in self.contigs)

    def gc(self) -> float:
        s = "".join(c.seq for c in self.contigs)
        n = sum(s.count(b) for b in "ACGT")
        return (s.count("G") + s.count("C")) / n if n else 0.0

    def cds_nt(self, g: Gene) -> str:
        seq = self.contigs[g.contig].seq
        if g.strand > 0:
            nt = "".join(seq[a:b] for a, b in g.parts)
        else:
            nt = "".join(revcomp(seq[a:b]) for a, b in g.parts)
        return nt[g.codon_start - 1:]

    def codon_map(self, g: Gene) -> list[tuple[int, int]]:
        """Genomic coordinate of the first base (5' on the coding strand) of each codon."""
        coords: list[int] = []
        for a, b in g.parts:
            coords.extend(range(a, b) if g.strand > 0 else range(b - 1, a - 1, -1))
        coords = coords[g.codon_start - 1:]
        return [(i // 3, coords[i]) for i in range(0, len(coords) - 2, 3)]


def _open_text(raw: bytes) -> str:
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", errors="replace")


def read_genome(path_or_bytes, name: str | None = None) -> Genome:
    if isinstance(path_or_bytes, (str, Path)):
        p = Path(path_or_bytes)
        raw = p.read_bytes()
        name = name or re.sub(r"(\.gz)$", "", p.name)
    else:
        raw = path_or_bytes
        name = name or "genome"
    sha = hashlib.sha256(raw).hexdigest()
    text = _open_text(raw)
    stripped = text.lstrip()
    if stripped.startswith("LOCUS"):
        g = _parse_genbank(text, name)
    elif stripped.startswith(">"):
        g = _parse_fasta(text, name)
    else:
        raise ValueError("Not a FASTA or GenBank file (it should start with '>' or 'LOCUS').")
    if not g.contigs or g.length == 0:
        raise ValueError("No DNA sequence found in the file.")
    acgt = sum(sum(c.seq.count(b) for b in "ACGT") for c in g.contigs)
    if acgt < 0.8 * g.length:
        raise ValueError("This does not look like a DNA sequence (fewer than 80% A/C/G/T). "
                         "Codon Atlas needs nucleotide input, not protein.")
    g.sha256 = sha
    return g


def _clean(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")


def _parse_fasta(text: str, name: str) -> Genome:
    contigs: list[Contig] = []
    head, buf = None, []
    for line in io.StringIO(text):
        if line.startswith(">"):
            if head is not None:
                contigs.append(Contig(head.split()[0] if head else f"seq{len(contigs) + 1}", _clean("".join(buf)), head))
            head, buf = line[1:].strip(), []
        else:
            buf.append(line)
    if head is not None:
        contigs.append(Contig(head.split()[0] if head else f"seq{len(contigs) + 1}", _clean("".join(buf)), head))
    organism = ""
    if contigs and contigs[0].description:
        organism = " ".join(contigs[0].description.split()[1:4])
    return Genome(name, [c for c in contigs if c.seq], source_format="fasta", organism=organism)


# ---------------------------------------------------------------- GenBank
_LOC_RANGE = re.compile(r"<?(\d+)\.\.>?(\d+)|<?(\d+)(?!\.)")


def parse_location(loc: str) -> tuple[list[tuple[int, int]], int, bool]:
    """GenBank location -> (parts in biological order, strand, partial)."""
    loc = loc.replace(" ", "")
    strand = 1
    partial = "<" in loc or ">" in loc
    if loc.startswith("complement(") and loc.endswith(")"):
        strand = -1
        loc = loc[len("complement("):-1]
    parts = []
    inner_complement = "complement(" in loc
    for m in re.finditer(r"(complement\()?<?(\d+)(?:\.\.>?(\d+))?", loc):
        a = int(m.group(2))
        b = int(m.group(3) or a)
        parts.append((a - 1, b))
    if inner_complement:  # join(complement(..),complement(..)): already listed 5'->3'
        strand = -1
    if strand < 0 and not inner_complement:
        parts = parts[::-1]
    return parts, strand, partial


def _parse_genbank(text: str, name: str) -> Genome:
    contigs: list[Contig] = []
    genes: list[Gene] = []
    organism = ""
    tables: set[int] = set()
    for rec in re.split(r"\n//\s*(?:\n|$)", text):
        if "LOCUS" not in rec:
            continue
        m = re.search(r"^LOCUS\s+(\S+)", rec, re.M)
        acc = re.search(r"^VERSION\s+(\S+)", rec, re.M)
        cname = acc.group(1) if acc else (m.group(1) if m else f"record{len(contigs) + 1}")
        om = re.search(r"^\s+ORGANISM\s+(.+)$", rec, re.M)
        if om and not organism:
            organism = om.group(1).strip()
        defn = re.search(r"^DEFINITION\s+(.+)$", rec, re.M)
        if "\nORIGIN" not in rec:
            continue
        head, origin = rec.split("\nORIGIN", 1)
        seq = _clean(re.sub(r"\d", "", origin))
        if not seq:
            continue
        ci = len(contigs)
        contigs.append(Contig(cname, seq, defn.group(1).strip() if defn else ""))
        if "\nFEATURES" not in head:
            continue
        feats = head.split("\nFEATURES", 1)[1].split("\n", 1)[1] if "\n" in head.split("\nFEATURES", 1)[1] else ""
        for key, loc, quals in _features(feats):
            if key != "CDS":
                continue
            if "pseudo" in quals or "pseudogene" in quals:
                continue
            try:
                parts, strand, partial = parse_location(loc)
            except Exception:
                continue
            if not parts:
                continue
            tt = quals.get("transl_table")
            if tt and tt.isdigit():
                tables.add(int(tt))
            genes.append(Gene(ci, parts, strand, quals.get("locus_tag", ""), quals.get("protein_id", ""),
                              quals.get("gene", ""), quals.get("product", ""),
                              int(quals.get("codon_start", "1") or 1),
                              quals.get("translation", "").replace(" ", ""), "annotation", partial))
    g = Genome(name, contigs, genes, "annotation" if genes else "none", "genbank", organism=organism)
    g.transl_table = min(tables) if tables else None
    return g


def _features(block: str):
    key = loc = None
    quals: dict[str, str] = {}
    cur_q = None
    for line in block.splitlines():
        if len(line) < 21:
            continue
        if line[5] != " " and not line.startswith(" " * 21):
            if key:
                yield key, loc, quals
            key = line[5:21].strip()
            loc = line[21:].strip()
            quals, cur_q = {}, None
            continue
        body = line[21:]
        if body.startswith("/"):
            cur_q = None
            if "=" not in body:
                quals[body[1:].strip()] = ""
                continue
            q, v = body[1:].split("=", 1)
            if v.startswith('"'):
                if len(v) > 1 and v.endswith('"'):
                    quals[q] = v[1:-1]
                else:
                    quals[q] = v[1:]
                    cur_q = q
            else:
                quals[q] = v
        elif cur_q:
            joiner = "" if cur_q == "translation" else " "
            v = body.strip()
            quals[cur_q] += joiner + v.rstrip('"')
            if v.endswith('"'):
                cur_q = None
        elif quals == {} and key:
            loc += body.strip()
    if key:
        yield key, loc, quals


# ---------------------------------------------------------------- gene prediction
def predict_genes(genome: Genome, table_id: int = 11) -> list[Gene]:
    import pyrodigal

    total = genome.length
    meta = total < 100_000
    finder = pyrodigal.GeneFinder(meta=meta, mask=True)
    if not meta:
        finder.train(*[c.seq.encode() for c in genome.contigs if len(c.seq) >= 1000][:200] or
                     [genome.contigs[0].seq.encode()], translation_table=table_id)
    genes = []
    n = 0
    for ci, c in enumerate(genome.contigs):
        for pg in finder.find_genes(c.seq.encode()):
            n += 1
            a, b = pg.begin - 1, pg.end
            genes.append(Gene(ci, [(a, b)], pg.strand, f"{c.name}_{n}", "", "", "predicted protein",
                              1, "", "prodigal", bool(pg.partial_begin or pg.partial_end)))
    return genes


def ensure_genes(genome: Genome, table_id: int = 11) -> None:
    if genome.genes:
        return
    genome.genes = predict_genes(genome, table_id)
    genome.gene_source = "prodigal"


def canonical_protein(genome: Genome, g: Gene, table_id: int = 11) -> str:
    """Protein under the given table; the GenBank translation is preferred when it matches."""
    if g.translation:
        return g.translation
    nt = genome.cds_nt(g)
    prot = translate(nt, table(table_id))
    if prot.endswith("*"):
        prot = prot[:-1]
    if prot and not g.partial:
        prot = "M" + prot[1:]
    return prot
