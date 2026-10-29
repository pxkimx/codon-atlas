"""Tables (TSV) and a self-contained HTML report for one genome, and the batch table."""
from __future__ import annotations

import html as H
import re
import time

from .gencode import AA20, BASES, MODELS, THREE


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (list, tuple)):
        return "; ".join(_fmt(x) for x in v)
    if isinstance(v, dict):
        return "; ".join(f"{k}:{_fmt(x)}" for k, x in v.items())
    return str(v).replace("\t", " ").replace("\n", " ")


def _tsv(rows: list[dict], keys: list[str]) -> str:
    return "\t".join(keys) + "\n" + "".join("\t".join(_fmt(r.get(k)) for k in keys) + "\n" for r in rows)


def tables(r: dict) -> dict[str, str]:
    codons = r["inference"]["codons"]
    code_rows = [{**c, "inferred_name": THREE.get(c["inferred"], c["inferred"]), "standard_name": THREE.get(c["standard"], ""),
                  "families": [f["name"] for f in c.get("families", [])][:5]} for c in codons]
    out = {
        "genetic_code.tsv": _tsv(code_rows, ["codon", "standard", "inferred", "status", "best", "best_prob", "standard_prob",
                                             "alternative", "log10_alt_vs_standard", "n_raw", "n_used", "n_unique",
                                             "support", "usage", "trna_summary", "families"]),
        "decoding_probabilities.tsv": "codon\t" + "\t".join(MODELS) + "\n" + "".join(
            c["codon"] + "\t" + "\t".join(f"{p:.4g}" for p in c["posterior"]) + "\n" for c in codons),
        "in_domain_stop_codons.tsv": _tsv([{**s, "gene_id": (s.get("gene") or {}).get("id"),
                                            "gene_product": (s.get("gene") or {}).get("product"),
                                            "gene_aa_pos": (s.get("gene") or {}).get("aa_pos")} for s in r["stops"]["events"]],
                                          ["contig", "pos", "strand", "codon", "kind", "family", "acc", "desc", "column",
                                           "expected", "expected_p", "second", "second_p", "evalue", "gene_id", "gene_product", "gene_aa_pos"]),
        "trnas.tsv": _tsv([{**t, "reads": [x["codon"] for x in t["reads"]]} for t in r["trnas"]],
                          ["contig", "begin", "end", "strand", "amino_acid", "anticodon", "reads", "score", "n20", "n37", "n73",
                           "pair_1_72", "identity", "sequence"]),
    }
    return out


AA_CLASS = {**{a: "np" for a in "AVLIMFWPG"}, **{a: "po" for a in "STCYNQ"}, **{a: "ac" for a in "DE"}, **{a: "ba" for a in "KRH"}}


def html(r: dict) -> str:
    g = r["genome"]
    inf = r["inference"]
    by = {c["codon"]: c for c in inf["codons"]}
    name = H.escape(r.get("meta", {}).get("name", g["name"]))
    rows = []
    for a in BASES:
        for c3 in BASES:
            cells = []
            for b in BASES:
                cod = a + b + c3
                c = by[cod]
                cls = {"reassigned": "re", "stop-to-sense": "re", "uninferred": "un", "stop": "st"}.get(c["status"], "ok")
                cells.append(f'<td class="{cls}"><b>{cod.replace("T", "U")}</b> {H.escape(c["inferred"])}'
                             f'<small> {c["standard"]} · n={c["n_used"]}</small></td>')
            rows.append("<tr>" + "".join(cells) + "</tr>")
    hyp = "".join(f"<li><b>{h['codon'].replace('T', 'U')}</b> {THREE.get(h['standard'], '?')} → <b>{THREE.get(h['aa'], h['aa'])}</b> "
                  f"(decoding probability {h['prob']:.6f})</li>" for h in r["hypotheses"]) or "<li>None: every inferred codon matches the reference code.</li>"
    leads = "".join(f"<li>{l['codon'].replace('T', 'U')} {THREE.get(l['standard'], '?')} → {THREE.get(l['aa'], l['aa'])}: "
                    f"evidence 10<sup>{l['log10']:.1f}</sup> over the standard residue from {l['n']} columns, below the call threshold</li>"
                    for l in r["leads"]) or "<li>None.</li>"
    stops = "".join(f"<tr><td>{c.replace('T', 'U')}</td><td>{s['events']}</td><td>{s['genes']}</td>"
                    f"<td>{H.escape(', '.join(f'{k} {v}' for k, v in s['expected'].items()))}</td>"
                    f"<td>{H.escape(', '.join(f'{k} ({v})' for k, v in s['kinds'].items()))}</td></tr>"
                    for c, s in r["stops"]["summary"].items())
    trs = "".join(f"<tr><td>tRNA-{H.escape(t['amino_acid'])}</td><td>{t['anticodon']}</td><td>{', '.join(x['codon'] for x in t['reads'])}</td>"
                  f"<td>{t['n20']}</td><td>{t['n73']}</td><td>{H.escape('; '.join(t['identity']))}</td></tr>"
                  for t in r["trnas"] if t["identity"])
    st = r["settings"]
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>Codon Atlas report · {name}</title>
<style>
body{{font:14px/1.5 -apple-system,BlinkMacSystemFont,"Helvetica Neue",Arial,sans-serif;color:#1d1d1f;margin:32px auto;max-width:980px;padding:0 20px}}
h1{{font-size:28px;margin:0 0 4px}} h2{{font-size:18px;margin:28px 0 8px}} .sub{{color:#6e6e73}}
table{{border-collapse:collapse;width:100%;margin:8px 0}} td,th{{border:1px solid #e5e5ea;padding:5px 7px;text-align:left;font-size:12.5px}}
th{{background:#f5f5f7}} td.re{{background:#fff1e6;outline:2px solid #ff9f0a}} td.un{{color:#8e8e93;background:#fafafa}} td.st{{color:#8e8e93}}
td small{{color:#8e8e93}} .kv span{{display:inline-block;margin-right:22px}} .note{{color:#6e6e73;font-size:12.5px}}
</style></head><body>
<h1>{name}</h1><div class="sub">{H.escape(g['organism'] or '')} · Codon Atlas {r['version']} · {time.strftime('%Y-%m-%d %H:%M', time.localtime(r['created']))}</div>
<h2>Genome</h2><div class="kv"><span>{g['length']:,} bp</span><span>{g['contigs']} record(s)</span><span>GC {100 * g['gc']:.1f}%</span>
<span>{g['genes']:,} genes ({'annotated' if g['gene_source'] == 'annotation' else 'predicted'})</span><span>{inf['n_columns']:,} aligned Pfam columns from {inf['n_families']:,} families</span></div>
<h2>Inferred genetic code</h2><p><code>{inf['code']}</code></p>
<table>{''.join(rows)}</table>
<p class="note">Orange: inferred meaning differs from NCBI table {st['reference_table']}. Grey: not inferred (‘?’; expected for stop codons and rare codons). n = Pfam consensus columns used after the 1% cap.</p>
<h2>Called reassignments</h2><ul>{hyp}</ul>
<h2>Leads below threshold</h2><ul>{leads}</ul>
<h2>Stop codons inside conserved domains</h2><table><tr><th>Codon</th><th>Aligned positions</th><th>Genes</th><th>Expected residue at column</th><th>Kind</th></tr>{stops}</table>
<h2>tRNAs reading arginine codons (identity elements)</h2><table><tr><th>tRNA</th><th>Anticodon</th><th>Reads</th><th>N20</th><th>N73</th><th>Identity call</th></tr>{trs or '<tr><td colspan=6>None found</td></tr>'}</table>
<h2>Method</h2><p class="note">Six-frame standard-code translation (stops as X) searched with Pfam {H.escape(r['profiles']['release'])} profiles
({'rebuilt from seed alignments without entropy weighting' if r['profiles']['set'] == 'enone' else 'as distributed'}) using HMMER (pyhmmer).
Domains with independent E-value &lt; {st['evalue']:g}; aligned columns with posterior probability ≥ {st['min_pp']}; one column per codon position (best domain);
no column above {100 * st['column_cap']:.0f}% of a codon's columns. Excluded family classes: {H.escape(', '.join(st['exclude']) or 'none')}.
Decoding probabilities from the Codetta model (Shulgina &amp; Eddy, eLife 2021) with threshold {st['threshold']}. tRNAs by ARAGORN (pyaragorn); identity
elements by heuristic position finding. Candidate codes should be checked independently, ideally experimentally.</p>
</body></html>"""


def batch_tsv(b: dict) -> tuple[str, bytes]:
    keys = ["index", "name", "organism", "status", "length", "gc", "genes", "columns", "reassignments", "leads", "trna_flags",
            "in_domain_stops", "uninferred_sense", "code", "score", "seconds", "error"]
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", b.get("meta", {}).get("name", "batch")) + "_codon_atlas.tsv"
    return name, _tsv(b["rows"], keys).encode()
