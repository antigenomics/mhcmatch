<!-- The SVGs are transparent (the `_bg` variants carry a full-canvas white/dark fill and are not
     used here); the PNG is the fallback. GitHub honours <source> and gets the SVG in both colour
     schemes; PyPI ignores <picture>/<source> and does not render SVG at all, so it falls through
     to the <img> and keeps the PNG it already renders well. Changing the <img> to an SVG would
     blank the logo on the PyPI project page. -->
<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/mhcmatch/master/assets/mhcmatch_dark.svg">
    <source srcset="https://raw.githubusercontent.com/antigenomics/mhcmatch/master/assets/mhcmatch_light.svg">
    <img alt="mhcmatch" src="https://raw.githubusercontent.com/antigenomics/mhcmatch/master/assets/mhcmatch_light.png" width="340">
  </picture>
</p>

<h1 align="center">mhcmatch — which neoantigens are presented, and which ones a T cell will see</h1>

<p align="center">
  <a href="https://pypi.org/project/mhcmatch/"><img alt="PyPI" src="https://img.shields.io/pypi/v/mhcmatch"></a>
  <a href="https://github.com/antigenomics/mhcmatch/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/mhcmatch/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://antigenomics.github.io/mhcmatch/"><img alt="docs" src="https://github.com/antigenomics/mhcmatch/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

Pure Python, no compiled extension beyond the [`seqtree`](https://github.com/antigenomics/seqtree)
search core. MHC-I and MHC-II, human and mouse. Every reference dataset is fetched from
[`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data) on first use, so a fresh
`pip install` runs every example here with no manual downloads.

```bash
pip install mhcmatch

# rank a donor's neoantigen candidates end to end
mhcmatch rank fasta candidates.fasta --alleles donor.alleles --cls mhc1 --tumor SKCM --out ranked.tsv
```

```python
import mhcmatch

store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")  # auto-fetched, cached
store.restriction("NLVPMVATV", calibrated=True)   # ranked alleles + %rank / P(present) / band
store.binder_score("NLVPMVATV")                   # the single-number binder index
store.scan_protein(my_protein, cls="mhc1")        # which windows are presented
```

**📖 [Full documentation](https://antigenomics.github.io/mhcmatch/)** ·
[Getting started](https://antigenomics.github.io/mhcmatch/getting-started.html) ·
[Command reference](https://antigenomics.github.io/mhcmatch/cli.html) ·
[The EPIC scorer](https://antigenomics.github.io/mhcmatch/neoantigen.html) ·
[Notebooks](notebooks/README.md)

## Pick your entry point

The fifteen things people actually come here to do. The
[CLI page](https://antigenomics.github.io/mhcmatch/cli.html) has every command with every flag.

| your question | command | Python |
|---|---|---|
| Which of these peptides does an allele present? | `mhcmatch predict f.fasta --alleles 'HLA-A*02:01'` | `predict.predict_fasta` |
| Which allele presents this peptide? | `mhcmatch restriction PEP --calibrated` | `store.restriction` |
| Is it a binder at all, one number? | `mhcmatch binder PEP` | `store.binder_score` |
| What is the IC50, and vs its wild type? | `mhcmatch affinity PEP --allele A --wt WTPEP` | `store.affinity_model` |
| Will a T cell respond to it? | `mhcmatch complement --peptides p.txt` | `complement.score` |
| **Rank neoantigen candidates for a donor** | `mhcmatch rank fasta ...` | `rank.rank_fasta` |
| Why did *this* candidate rank there? | `mhcmatch explain PEP --allele A` | — |
| Which peptides in this protein are presented? | `mhcmatch scan p.fasta --correction bh` | `store.scan_protein` |
| What self / viral peptide does it mimic? | `mhcmatch mimics --peptides p.txt` | `mimics.neighbours` |
| Where in the proteome does it come from? | `mhcmatch source --peptides p.txt` | `Proteome.find_sources` |
| Is this gene on in a normal tissue? | `mhcmatch expression GENE --tissue T --safety` | `expression.lookup` |
| What does this allele's motif look like? | `mhcmatch logo 'HLA-A*02:01'` | `logo.motif` |
| **Which *k* candidates go in the cassette?** | `mhcmatch cassette select --candidates pool.tsv -k 20` | `cassette.select` |
| **What is this cassette worth, against another?** | `mhcmatch cassette score --cassettes c.tsv --pool pool.tsv` | `cassette.score` / `.lam` |
| Assemble and order the chosen units | `mhcmatch cassette build --candidates units.tsv --screen` | `vector.select` / `vector.order` |

Four things that surprise people, and each is a link away from its detail:

- **`predict` and `rank fasta` drop nothing by default.** `--rank-threshold` takes `sb` / `wb` /
  `none` / a percentage, and the tiers are **class-aware**, because a bare number cannot be: `2.0`
  is the weak cut for class I and the *strong* cut for class II. On class II a flat `2.0` kept
  **0 of 56** scored pairs in testing — an empty table, returncode 0.
- **Two whitelists, because they make two different claims.** `--keep-genes 'TP53,KRAS'` keeps a
  candidate for its gene; `--keep-epitopes builtin` keeps it for being one of 23,299 assay-validated
  immunogenic peptides. Matched rows carry `keep_reason`, because a row kept for its gene is not
  evidence about its peptide. → [CLI](https://antigenomics.github.io/mhcmatch/cli/commands.html)
- **`cassette` and `vector` are two jobs, not two names.** `cassette` *chooses* units and scores a
  finished construct; `vector` *assembles* one already chosen. One `mhcmatch cassette` verb hides
  the split on the CLI; the Python column above does not.
  → [Designing a cassette](https://antigenomics.github.io/mhcmatch/cassette.html)
- **`predict` and `restriction` are different axes.** `predict` asks *is it presented at all*
  (the NetMHCpan `%Rank_EL` analogue); `restriction` asks *which allele*. A peptide can top one and
  not the other — `NLVPMVATV` is unambiguously A\*02:01-restricted yet bands mid-pack against
  A\*02:01's own ligands.

## The shipped model

**E**xpression, **P**resentation, **I**mmunogenic **C**omplementarity — four blocks, listed here
in the order that spells the name:

| letter | block | columns |
|---|---|---|
| `E` | expression | `expr_lvl`, `expr_norm` |
| `P` | presentation | `binder`, `log10a` |
| `I` | immunogenic — physchem | `C_phys_buried`, `C_phys_charge` |
| `C` | complementarity — corpus | `C_corpus_thymus`, `C_corpus_self`, `C_corpus_viral` |

The **fit** enters them presentation first — `presentation, expression, physchem, corpus`, which is
the `blocks` list on the artifact itself — so a later block's coefficient is what that term is
worth *after* the earlier ones. Read a coefficient against that order, not against the acronym.

One artifact per `(cls, species, mode)`, and **no fallback** — asking for a cell that was never
fitted raises rather than scoring it with a neighbour's coefficients.

<!-- BEGIN shipped-models (generated by mhcmatch._modeldoc) -->
| `model_id` | model version | release | terms | rows | positives | intercepts | AUROC | how that AUROC is measured |
|---|--:|---|--:|--:|--:|--:|--:|---|
| `mhc1.human.neoantigen` | **12** | 1.20.0 | 9 | 339,595 | 594 | 7 per screen | **0.7094** | leave-one-screen-out, mean |
| `mhc1.mouse.neoantigen` | **6** | 1.20.0 | 9 | 921 | 379 | 61 per reference | **0.6335** | in-sample, within reference |
| `mhc2.human.neoantigen` | **2** | 1.20.0 | 6 | 1,112 | 656 | 157 per reference | **0.6020** | in-sample, within reference |
| `mhc2.mouse.neoantigen` | **4** | 1.20.0 | 6 | 468 | 177 | 30 per reference | **0.5741** | in-sample, within reference |
| `mhc1.human.pathogen` | **2** | 1.20.0 | 2 | 16,790 | 7,002 | 1 per corpus, global | **0.5988** | in-sample, pooled off the logit |
| `mhc1.mouse.pathogen` | **2** | 1.20.0 | 2 | 10,404 | 2,196 | 1 per corpus, global | **0.5561** | in-sample, pooled off the logit |
| `mhc2.human.pathogen` | **2** | 1.20.0 | 3 | 7,946 | 5,148 | 1 per corpus, global | **0.5824** | in-sample, pooled off the logit |
| `mhc2.mouse.pathogen` | **2** | 1.20.0 | 3 | 11,725 | 3,324 | 1 per corpus, global | **0.6446** | in-sample, pooled off the logit |
<!-- END shipped-models -->

**Human class I is the cell this library is built around.** Human class II is fitted on far fewer
rows from far fewer assays and its AUROC says so; the mouse cells are thinner again, and mouse
class II is the thinnest of all. Read a class-II or mouse number as a transfer that was checked,
not as a second headline.

**Two things the table does not show.**

*A ninth cell ships that is not a default.* `mhc2.human.neoantigen` has a second fit carrying the
three corpus terms, reached as `rank.aggregate("mhc2", "human", variant="corpus")`. It exists
because the corpus block became computable at class II in 1.20.0 — the HLA Ligand Atlas gives
132,818 class-II peptides over 29 benign tissues, against the 27,987 of the class-II thymic
fraction alone. It is **not** the default because computing it answered in the negative: on the
same 1,112 rows all three coefficients span zero, `C_corpus_thymus` at +0.0170 (95 % CI −0.4944 to
+0.3995, sign held in 53 % of resamples). **At class II the corpus block adds nothing**, and the
six-term fit is what a class-II score should use. The same corpus rebuild is what made the mouse
class-I cell fittable, where the block does carry.

*`mode="pathogen"` is a supported surface, not an internal.* Four cells are fitted on IEDB T-cell
corpora whose negative class is "assayed and did not respond", which is a different question from
the neoantigen cells' "nominated and did not respond" — so the two are not comparable term by term,
and a pathogen cell drops expression outright, there being no host gene to measure. Reach them with
`--epitope pathogen`. They also back the physicochemistry work, which is why the pathogen corpora are
maintained alongside the neoantigen ones.

**The AUROC column is three different protocols — do not read it down, and do not average it.**
Coefficients are deliberately not written here: they move with every refit, and this table quoted a
superseded set for a full release each time it was maintained by hand. The record is
**[docs/models.rst](https://antigenomics.github.io/mhcmatch/models.html)** — every coefficient with
its bootstrap interval, what each fit was trained on, how its AUROC was measured, and the caveats
that come with each one, all generated from the artifacts on every docs build. Or ask the artifact:

```bash
mhcmatch models --all                 # which cells ship
mhcmatch rank --coefficients          # every term, its block, its coefficient
mhcmatch rank --holdout               # held-out AUROC, the grouped CVs, the corpus
```

## Python

```python
import mhcmatch
from mhcmatch import complement, expression, mimics, rank

store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")
store.decompose("NLVPMVATV")                         # anchor / TCR-facing split

aff = store.affinity_model("mhc1")
aff.predict_ic50("NLVPMVATV", "HLA-A*02:01")         # 18.9 nM (shortlist tier)
aff.amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")   # Kd_WT/Kd_MT (Łuksza eq. 9)

complement.score(peptides)                           # vectorised: pass the list, not a loop
complement.posterior(peptides, prior=4.2e-4)         # the log-odds carries NO prior; supply yours

rank.models()                                        # every shipped fit and its provenance
rank.aggregate("mhc1", "mouse")                      # the artifact itself

expression.gene_level("Trp53", species="mouse")      # FANTOM5 tissues + syngeneic models
mimics.neighbours(peptides, ref_sets, threads=0)     # threaded C++ neighbour search

pm = mhcmatch.Proteome.from_hf("human")
pm.find_sources(peptides, max_subs=1, threads=0)     # batch; find_source() is the single-query form
pm.wildtype("NLVPMVATV")                             # the WT counterpart, for agretopicity
```

Every one of these takes a **list** and returns one per input. Passing a list beats looping the
single-query form by orders of magnitude — the batch call reaches C++ once with the GIL released.
→ [Batch and threads](README_EXT.md#batch-and-threads--read-this-before-scripting-a-loop)

**Advanced entry points** — `mimicry` (signed per-component mimicry risk), `portfolio` (the
over-dispersion and support machinery under `cassette score`), `luksza` (the published `R` term),
`recognition` (the ESM head), `calibrate`, `known`, `ligand`, `pseudoseq`.
→ [API reference](https://antigenomics.github.io/mhcmatch/api.html)

## Install extras

The base install is `seqtree`, `numpy` and `huggingface_hub`. **Every model that ships by default
runs on it** — `torch` is not required to score recognition, and asking for a head without its
extra raises a named error rather than silently degrading.

| extra | pulls | needed for |
|---|---|---|
| `mhcmatch[notebooks]` | `marimo`, `polars` | the worked examples in `notebooks/` |
| `mhcmatch[logo]` | `logomaker`, `matplotlib`, `pandas` | **drawing** a motif logo. `logo.motif` returns the matrix on the base install |
| `mhcmatch[stats]` | `scipy` | your own over-dispersion (`portfolio.betabinom_rho`) and the exact LP behind `portfolio.linearly_supported` |
| `mhcmatch[esm]` | `torch`, `transformers` | **only** the `esm64_glm` recognition head (~2.4 GB checkpoint on first use) |
| `mhcmatch[structure]` | `tcren` | the structure-based ΔΔG head |
| `mhcmatch[precursor]` | `vdjmatch` | precursor-frequency estimates |
| `mhcmatch[docs]` | `sphinx`, `pydata-sphinx-theme` | building the documentation |

`mhcmatch bootstrap` optionally pre-fetches the ligand panel (~16 MB) — it only decides *when*
reference data is fetched, which matters on a compute node with no outbound network.

## Deeper reading

Everything answering "why does it do that" rather than "how do I run it" is in
**[README_EXT.md](README_EXT.md)**, so this page stays the length of a page:

| | |
|---|---|
| [Batch and threads](README_EXT.md#batch-and-threads--read-this-before-scripting-a-loop) | read before scripting a loop |
| [Caching calibration across jobs](README_EXT.md#caching-calibration-across-jobs) | the one environment variable a cluster run wants |
| [Composition is not ranking](README_EXT.md#composition-is-not-ranking) | why a cassette is a set problem, and what `lam` compares |
| [What `rank` costs](README_EXT.md#what-rank-costs) | where the time goes, and which stages ship off |
| [The two axes](README_EXT.md#the-two-axes) | presentation and recognition, and what each carries |
| [Presentation and affinity are not the same term](README_EXT.md#presentation-and-affinity-are-not-the-same-term) | why both are fitted and why `binder` is neither |
| [Data](README_EXT.md#data) | the four staging tiers, and what each fetches |
| [Deployment](README_EXT.md#deployment) | Nextflow, Snakemake, an overlay for a pipeline you already run, containers and SLURM |
| [Benchmarks](README_EXT.md#benchmarks) | what is measured where |

## Development

```bash
bash setup.sh            # repo-local .venv + editable install (uses a sibling ../seqtree if present)
bash setup.sh --tests    # + pytest
pytest -q
```

Theory and derivations are in the manuscript repo (`latex_sn/`, Methods and Supplementary Note 1);
what is planned and in flight is in [`ROADMAP.md`](ROADMAP.md).
