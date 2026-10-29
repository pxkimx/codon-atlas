"""One genome from file to result: search, inference, stop codons, tRNAs, genes, sites."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import __version__, engine, trna
from .gencode import CODONS, STOPS_STANDARD, THREE, gene_table_for, table
from .genome import ensure_genes, read_genome
from .pfamdb import CLASSES, ProfileSet

MAX_STOP_EVENTS = 3000


def analyze(genome_path: str, settings: dict, outdir: Path, progress=None, prof: ProfileSet | None = None) -> dict:
    t0 = time.time()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    st = engine.Settings.from_dict(settings)
    say = progress or (lambda *a, **k: None)
    say("Reading the genome", 0.0)
    genome = read_genome(genome_path)
    prof = prof or ProfileSet(st.profiles)
    obs = engine.search(genome, prof, st, progress=lambda m, f, **k: say(m, 0.02 + 0.80 * f, **k))
    np.save(outdir / "columns.npy", obs)
    say("Inferring the decoding of each codon", 0.83)
    inf = engine.infer(obs, prof, st)
    inferred = {c["codon"]: c["inferred"] for c in inf["codons"]}

    say("Finding protein-coding genes", 0.88)
    gene_table = genome.transl_table or gene_table_for(inferred)
    ensure_genes(genome, gene_table)
    gidx = engine.GeneIndex(genome)

    say("Looking for stop codons inside conserved domains", 0.91)
    stops = engine.stop_events(obs, prof, genome, gidx)
    for s in stops:  # a stop codon read as an amino acid genome-wide explains its in-domain positions
        aa = inferred.get(s["codon"], "?")
        if aa != "?" and s["kind"] not in ("Mobile element / phage (often pseudogenes)",):
            s["kind"] = f"Reassigned stop codon (read as {THREE.get(aa, aa)})"
    stop_summary = {}
    for c in STOPS_STANDARD:
        ev = [s for s in stops if s["codon"] == c]
        kinds = {}
        for s in ev:
            kinds[s["kind"]] = kinds.get(s["kind"], 0) + 1
        fams = {}
        for s in ev:
            fams[s["family"]] = fams.get(s["family"], 0) + 1
        exp = {}
        for s in ev:
            exp[s["expected"]] = exp.get(s["expected"], 0) + 1
        stop_summary[c] = {"events": len(ev), "kinds": kinds, "genes": len({s["gene"]["id"] for s in ev if s.get("gene")}),
                           "families": dict(sorted(fams.items(), key=lambda kv: -kv[1])[:10]),
                           "expected": dict(sorted(exp.items(), key=lambda kv: -kv[1])[:6])}

    say("Finding tRNA genes (ARAGORN)", 0.94)
    try:
        trnas = trna.find_trnas(genome, 11)
        trna_error = ""
    except Exception as e:  # never lose the whole run over the tRNA step
        trnas, trna_error = [], f"{type(e).__name__}: {e}"
    readers = trna.codon_readers(trnas)
    ref = table(st.reference_table)
    for c in inf["codons"]:
        rd = readers[c["codon"]]
        c["trnas"] = rd
        c["trna_summary"] = trna.summary_for_codon(c["codon"], rd, ref[c["codon"]])

    # codon usage inside genes
    usage = {c: 0 for c in CODONS}
    for g in genome.genes:
        nt = genome.cds_nt(g)
        for i in range(0, len(nt) - 2, 3):
            k = nt[i:i + 3]
            if k in usage:
                usage[k] += 1
    for c in inf["codons"]:
        c["usage"] = usage[c["codon"]]

    hypotheses = [{"codon": c["codon"], "aa": c["inferred"], "standard": c["standard"], "prob": c["best_prob"],
                   "status": c["status"]} for c in inf["codons"] if c["status"] in ("reassigned", "stop-to-sense")]
    leads = [{"codon": c["codon"], "aa": c["alternative"], "standard": c["standard"], "n": c["n_used"],
              "log10": c["log10_alt_vs_standard"], "prob": c["posterior"][engine.MODELS.index(c["alternative"])]}
             for c in inf["codons"]
             if c["status"] in ("uninferred",) and c.get("log10_alt_vs_standard") is not None
             and c["log10_alt_vs_standard"] > 1 and c["n_used"] >= 3]
    leads.sort(key=lambda x: -x["log10"])
    trna_flags = [t for t in trnas if any(not s.startswith(("Arg-type", "unclassified")) for s in t["identity"])]

    result = {
        "version": __version__,
        "created": time.time(),
        "genome": {"name": genome.name, "organism": genome.organism, "length": genome.length,
                   "contigs": len(genome.contigs), "gc": genome.gc(), "genes": len(genome.genes),
                   "gene_source": genome.gene_source, "gene_table": gene_table, "format": genome.source_format,
                   "sha256": genome.sha256, "annotated_table": genome.transl_table,
                   "contig_names": [c.name for c in genome.contigs[:50]]},
        "settings": engine.settings_dict(st),
        "profiles": {"set": prof.which, "release": prof.release, "families": len(prof),
                     "excluded": {k: v for k, v in prof.class_counts().items() if k in st.exclude}},
        "inference": inf,
        "hypotheses": hypotheses,
        "leads": leads,
        "stops": {"summary": stop_summary, "events": stops[:MAX_STOP_EVENTS], "total": len(stops)},
        "trnas": trnas,
        "trna_error": trna_error,
        "trna_flags": [{"anticodon": t["anticodon"], "amino_acid": t["amino_acid"], "identity": t["identity"],
                        "reads": [r["codon"] for r in t["reads"]], "contig": t["contig"], "begin": t["begin"]}
                       for t in trna_flags],
        "seconds": round(time.time() - t0, 1),
        "classes": CLASSES,
    }
    (outdir / "result.json").write_text(json.dumps(result, allow_nan=False, default=_jsonable))
    say("Done", 1.0)
    return result


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o).__name__)


def summary_row(result: dict) -> dict:
    """One line for the batch table."""
    g = result["genome"]
    inf = result["inference"]
    calls = [f"{h['codon']}→{h['aa']}" for h in result["hypotheses"]]
    leads = [f"{l['codon']}→{l['aa']} (10^{l['log10']:.0f}, n={l['n']})" for l in result["leads"][:3]]
    best_lead = max([l["log10"] for l in result["leads"]] or [0])
    score = 1000 * len(result["hypotheses"]) + best_lead + 5 * len(result["trna_flags"])
    return {"name": g["name"], "organism": g["organism"], "length": g["length"], "gc": g["gc"],
            "genes": g["genes"], "columns": inf["n_columns"], "code": inf["code"],
            "reassignments": calls, "leads": leads, "trna_flags": [f"tRNA-{t['amino_acid']}({t['anticodon']}): {'; '.join(t['identity'])}" for t in result["trna_flags"]],
            "in_domain_stops": {c: s["events"] for c, s in result["stops"]["summary"].items()},
            "uninferred_sense": sum(1 for c in inf["codons"] if c["status"] == "uninferred"),
            "score": score, "seconds": result["seconds"]}


def describe(h: dict) -> str:
    return f"{h['codon']} ({THREE.get(h['standard'], '?')}) → {THREE.get(h['aa'], h['aa'])}"
