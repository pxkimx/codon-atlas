# Codon Atlas

Find alternative decoding in a genome: which amino acid each of the 64 codons is read as,
stop codons inside conserved domains, the tRNAs that could read each codon, and every
affected site in the genes, exported as a protein database for proteomics.

The inference is a re-implementation of **Codetta** by Yekaterina Shulgina and Sean R. Eddy
(*eLife* 2021, [doi:10.7554/eLife.71402](https://doi.org/10.7554/eLife.71402); Codetta 2.0:
[github.com/kshulgina/codetta](https://github.com/kshulgina/codetta)). Please cite their paper
if you use the predictions. Codon Atlas is an independent program and is not affiliated with them.

More, with a screenshot and an example report: [paul-h-kim.com/software/codon-atlas](https://paul-h-kim.com/software/codon-atlas).

Requirements: macOS, Python 3.9+, internet once (Python packages and Pfam).

## Use

Double-click **Codon Atlas.app** (a copy is in `~/Desktop/Toolbox`). The first start creates a
private Python environment in `~/Library/Application Support/Codon Atlas/venv` (needs internet
once). Closing the window quits the program and anything it was running.

1. **Database** (once): download Pfam-A.seed (≈195 MB) and prepare the recommended profile set
   (≈10 minutes). Profiles are rebuilt from the seed alignments without entropy weighting, as
   Codetta does. Pfam-A.hmm can be prepared as an alternative set.
2. **Analyze**: drop a FASTA or GenBank genome (gzip is fine) and press *Analyze genome*. A
   2 Mb bacterial genome takes about a minute.
3. **Genetic code**: the codon wheel and table; click a codon for its decoding probabilities,
   bootstrap support, the residues its aligned columns expect, contributing families and the
   tRNAs that can read it.
4. **Stop codons**, **tRNAs**, **Sites & export**: in-domain stop codons, tRNA genes with
   identity elements, per-site lists, and the downloadable protein database package
   (standard, alternative and combined FASTA plus site tables, named like `TGA2W`, readable by
   Recode Detector's Database check), the tables and an HTML report.
5. **Batch**: many genomes at once, ranked by called reassignments, leads and tRNA flags.
6. **Runs**: every analysis is saved and can be reopened.

## Method in brief

Six-frame standard-code translation (stops as X, 100 kb pieces with 3 kb overlap) → Pfam
search with HMMER (pyhmmer), domains with i-Evalue < 1e-10 (E-values scaled to the number of
profiles, as in hmmscan), columns with posterior ≥ 0.95, one column per codon position →
for each codon, `Σ w_i [log e_Ci(a) − log ē(a)]` for 20 amino-acid models against a
non-specific model scored 0, uniform prior, call at probability ≥ 0.9999, no column above 1%
of a codon's columns. Mobile-element, viral, selenoprotein, pyrrolysine and mitochondrial
families are left out of inference by default (keyword classes on Pfam names and
descriptions). Bootstrap support and below-threshold "leads" are added. tRNAs come from
ARAGORN (pyaragorn); arginine-anticodon tRNAs get the Arg/Trp/Met/Gln identity calls used by
Shulgina & Eddy (A20, N73). Genes come from the GenBank annotation or Prodigal (pyrodigal).

Differences from Codetta to be aware of: HMMER is run as `hmmsearch` (profiles against the
translation) with Z set to the number of profiles instead of `hmmscan`; families are excluded
by keyword rather than Codetta's curated lists; the latest Pfam is used.

## Example

The public *Mycoplasma genitalium* G37 genome (NCBI NC_000908.2) is called as UGA → Trp
from 308 aligned Pfam positions with 100% bootstrap support, which is its known genetic code
(NCBI table 4). Fetch it from NCBI and drop it on the Analyze page to reproduce this.

## Third-party software and data

HMMER via pyhmmer (Larralde & Zeller), Prodigal via pyrodigal, ARAGORN via pyaragorn, and
Pfam (Mistry et al. 2021; CC0) are downloaded at first use under their own licences; none
are included here.

## Licence

MIT, see `LICENSE`.

## Developer notes

- `codon_atlas/engine.py` search + inference, `trna.py`, `sites.py` (site lists and database
  package), `pipeline.py` (one genome end to end), `worker.py` (background jobs in their own
  process group), `server.py` (standard-library HTTP server), `web/` (page), `pfamdb.py`
  (profile sets in `~/Library/Application Support/Codon Atlas/pfam`; `CA_HOME` overrides).
- Port 8784 (`CA_PORT`), tests/previews 8785. Tests: `python -m pytest tests`.
- `macos/build_app.sh` rebuilds the app (the code is copied into it).
