import io
import json
import zipfile

import numpy as np
import pytest

from codon_atlas import engine, sites, trna
from codon_atlas.gencode import AA20, CODON_INDEX, CODONS, TABLES, revcomp, table, translate
from codon_atlas.genome import parse_location, read_genome
from codon_atlas.lifecycle import Lifecycle

GB = """LOCUS       TEST1                    60 bp    DNA     linear   BCT 01-JAN-2026
DEFINITION  Test organism chromosome.
VERSION     TEST1.1
SOURCE      Testus testi
  ORGANISM  Testus testi
FEATURES             Location/Qualifiers
     CDS             1..15
                     /locus_tag="T_0001"
                     /protein_id="P1.1"
                     /product="first protein"
                     /translation="MRWK"
     CDS             complement(31..48)
                     /locus_tag="T_0002"
                     /product="second protein with a long
                     product name"
ORIGIN
        1 atgcggtggaaatagcccgg gtttaaaccc ttacgtcatt ccgcatcggc ttaacccggg
//
"""


def test_tables_have_64_codons():
    for n, (_, s) in TABLES.items():
        assert len(s) == 64, n
    assert table(4)["TGA"] == "W" and table(11)["TGA"] == "*" and table(25)["TGA"] == "G"


def test_genbank_parsing_and_cds(tmp_path):
    p = tmp_path / "t.gb"
    p.write_text(GB)
    g = read_genome(p)
    assert g.organism == "Testus testi" and len(g.contigs) == 1 and len(g.genes) == 2
    a, b = g.genes
    assert a.translation == "MRWK" and translate(g.cds_nt(a)) == "MRWK*"
    assert b.strand == -1 and b.product == "second protein with a long product name"
    assert g.cds_nt(b) == revcomp(g.contigs[0].seq[30:48])
    # codon map gives the 5' base of each codon on the coding strand
    assert [c for _, c in g.codon_map(a)][:3] == [0, 3, 6]
    assert [c for _, c in g.codon_map(b)][:2] == [47, 44]


def test_location_forms():
    assert parse_location("complement(join(1..10,20..30))") == ([(19, 30), (0, 10)], -1, False)
    assert parse_location("join(complement(20..30),complement(1..10))") == ([(19, 30), (0, 10)], -1, False)
    assert parse_location("<1..>90")[2] is True


def test_fasta_rejects_protein(tmp_path):
    p = tmp_path / "p.fa"
    p.write_text(">x\nMKLVFFEEQQRRSSTTWWYY\n")
    with pytest.raises(ValueError):
        read_genome(p)


def test_six_frames_positions():
    seq = "ATGAAATGACCCGGGTTTAAACCCGGG" * 3
    g = type("G", (), {})()
    g.contigs = [type("C", (), {"seq": seq, "name": "c"})()]
    frags, prots = engine.six_frames(g, 30, 6)
    for fr, p in zip(frags, prots):
        for i in range(len(p)):
            pos = fr.pos(i)
            cod = seq[pos:pos + 3] if fr.strand > 0 else revcomp(seq[pos - 2:pos + 1])
            assert CODONS[fr.codons[i]] == cod


class FakeProf:
    """FamA: columns 1-50 strongly W, 51-100 strongly R. FamB (a mobile element): 3 columns."""

    def __init__(self):
        em = np.full((103, 20), 0.01, dtype=np.float32)
        em[0:50, AA20.index("W")] = 0.81
        em[50:100, AA20.index("R")] = 0.81
        em[100:, :] = 0.05
        self.emissions = em
        self.offsets = np.array([0, 100])
        self.names = ["FamA", "FamB"]
        self.accessions = ["PF0", "PF1"]
        self.descriptions = ["a", "b"]
        self.classes = [[], ["mobile"]]

    def excluded_mask(self, exclude):
        return np.array([bool(set(exclude) & set(c)) for c in self.classes])


def _obs(rows):
    o = np.zeros(len(rows), dtype=engine.OBS_DTYPE)
    for i, (cod, h, k) in enumerate(rows):
        o[i] = (0, 1, i * 3, CODON_INDEX[cod], h, k, 1e-20)
    return o


def test_inference_calls_reassignment_and_respects_threshold():
    prof = FakeProf()
    rows = ([("CGG", 0, k) for k in range(1, 31)] * 2 + [("AGA", 0, k) for k in range(51, 91)]
            + [("TGG", 0, k) for k in range(1, 41)] + [("CGT", 0, 60)] * 30 + [("CGA", 0, 1)])
    res = engine.infer(_obs(rows), prof, engine.Settings(exclude=["mobile"], bootstrap=50))
    by = {c["codon"]: c for c in res["codons"]}
    assert by["CGG"]["inferred"] == "W" and by["CGG"]["status"] == "reassigned" and by["CGG"]["support"] == 1.0
    assert by["AGA"]["status"] == "standard" and by["TGG"]["status"] == "standard"
    assert by["TAA"]["status"] == "stop" and by["TAA"]["n_used"] == 0
    # below 100 columns each unique column counts once: CGG has 60 aligned positions, 30 columns
    assert by["CGG"]["n_used"] == 30 and by["CGG"]["n_raw"] == 60
    # one repeated column or a single position is not enough to call anything
    assert by["CGT"]["inferred"] == "?" and by["CGA"]["inferred"] == "?"
    assert by["CGA"]["alternative"] == "W" and by["CGA"]["log10_alt_vs_standard"] > 0


def test_excluded_families_are_ignored():
    prof = FakeProf()
    rows = [("CGG", 1, k) for k in (1, 2, 3)] * 3
    res = engine.infer(_obs(rows), prof, engine.Settings(exclude=["mobile"], bootstrap=0))
    assert {c["codon"]: c for c in res["codons"]}["CGG"]["n_raw"] == 0


def test_capped_weights():
    ids = np.array([1] * 500 + list(range(2, 502)))
    uniq, w = engine._capped_weights(ids, 0.01)
    assert w[uniq == 1][0] == 10 and w.sum() == 10 + 500


def test_trna_reading_rules_and_identity():
    assert [c for c, _ in trna.reads("CCG")] == ["CGG"]
    assert [c for c, _ in trna.reads("ACG")] == ["CGT", "CGC", "CGA"]
    assert [c for c, _ in trna.reads("TCA")] == ["TGA", "TGG"]
    # a tRNA-Arg(CCG) with G73 and C20 -> Trp-type
    s = "GGGTCTGTAGTTTAGTGGCAAAACAAGGGTGTCCGAAGCCAAAGATGGGGGTTCGATTCCCTCCAGGCCTGCCA"
    idt = trna.identity(s, 33, "Arg", "CCG")
    assert idt["n73"] == "G" and idt["n20"] == "C" and any(x.startswith("Trp-type") for x in idt["identity"])
    s2 = "GCGTCATTAGCTCAGCTGGATAGAGCGTCTGGCTACGGACCAGAAGGCCTGGGGTTCGAATCCTCAATGACGCGCCA"
    assert trna.identity(s2, 35, "Arg", "ACG")["identity"][0].startswith("Arg-type")


def test_database_package(tmp_path):
    p = tmp_path / "t.gb"
    p.write_text(GB)
    g = read_genome(p)
    name, data = sites.database_package(g, [("CGG", "Q")])
    z = zipfile.ZipFile(io.BytesIO(data))
    files = {n.split("/")[-1]: n for n in z.namelist()}
    only = z.read(files[[f for f in files if f.endswith("CGG2Q_variants_only.fasta")][0]]).decode()
    assert ">TEST1.1|P1.1|CGG2Q" in only and "CGG2Q_sites=2" in only and "MQWK" in only
    v = json.loads(z.read(files["validation.json"]))
    assert v["sites"] == 1 and v["altered_proteins"] == 1


def test_lifecycle_grace():
    l = Lifecycle(on_quit=lambda: None, grace=1.0, never_opened=5.0)
    assert not l.should_quit(l.started + 1)
    assert l.should_quit(l.started + 6)
    l.window_opened()
    assert not l.should_quit(l.started + 100)
    l.window_closed()
    assert not l.should_quit(l.last_zero + 0.5) and l.should_quit(l.last_zero + 1.5)
