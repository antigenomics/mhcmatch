# mhcmatch notebooks

Worked examples of the public API, one topic per notebook. These are
[marimo](https://marimo.io) notebooks — **plain Python files** with `@app.cell` decorators, so they
diff and review like source rather than like JSON.

Every notebook bootstraps its own data from public HuggingFace datasets
([`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data),
[`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark)) and caches it, so a
fresh `pip install mhcmatch` is enough to run them. No local paths, no pre-staged files.

## Index

| notebook | module | demonstrates | runtime |
|---|---|---|--:|
| [`01_presentation_and_binder_score.py`](01_presentation_and_binder_score.py) | `Store`, `mhcmatch.predict` | the core workflow — decompose, restriction, `binder_score`; every output column; why `p_binder` is pool-invariant and a within-list percentile is not | ~70 s |
| [`02_immunogenicity_features.py`](02_immunogenicity_features.py) | `mhcmatch.immuno` | the 141 physicochemical features; the three anchor schemes plus the contact scheme; run statistics as contiguity; the unsupervised contact profile | ~1 s |
| [`03_precursor_frequency.py`](03_precursor_frequency.py) | `mhcmatch.precursor` | `observed_mass` / `ball_mass` / `motif_mass` on one epitope's real TCRs; ball union vs naive sum; the CDR3-vs-junction trap | ~20 s |
| [`04_mimicry_and_self.py`](04_mimicry_and_self.py) | `mhcmatch.mimics` | thymus / viral / neoag reference scanning, and the self-identity leakage trap any mimicry feature has to exclude | ~20 s |
| [`05_position_role_bayes.py`](05_position_role_bayes.py) | `mhcmatch.posbayes` | anchor vs TCR-facing amino-acid evidence, the residues whose two roles carry opposite signs, and why the score carries no prior | ~2 s |
| [`06_complementarity.py`](06_complementarity.py) | `mhcmatch.complement` | the six feature blocks of the recognition axis; arrangement vs composition; that the `aa` block *is* `posbayes`; scoring a whole 511k-row corpus in one call | ~30 s |
| [`07_mimicry_risk.py`](07_mimicry_risk.py) | `mhcmatch.mimicry` | notebook 4's scan in fitted form — three references × two channels as signed log-odds; that the channels partition the peptide; pooled vs within-screen AUROC; why the tested-neoantigen database is an annotation and never a fitted term | ~5 s |
| [`08_ranking_and_cassette.py`](08_ranking_and_cassette.py) | `mhcmatch.rank`, `mhcmatch.vector` | the applied pipeline on a **mock scored table** — recompute rather than re-sort; why a unit is the 27-mer and not the 9-mer; safety exclusion before capacity; the spacer as a result; junction-spanning epitopes that belong to no gene | ~18 s |
| [`09_cassette_composition.py`](09_cassette_composition.py) | `mhcmatch.portfolio` | composition as a set problem, on synthetic data checkable by hand — that two cassettes with identical expected responders differ twofold in `P(at least one works)`; that a weighted sum reaches only the upper convex hull; that Chebyshev fixes the second and nothing fixes the first | ~1 s |
| [`10_cassette_select_and_score.py`](10_cassette_select_and_score.py) | `mhcmatch.cassette` | the two operations on a whole published corpus: `select` on the 46 NCI GI patients held out of the EPIC fit, `score` on the 1,631 units two trials manufactured. That `select` wins on `H` (0.219 vs 0.127) and a sort wins on `yield` (1.044 vs 0.768), because they maximise different things; that one offset per donor pins every pool mean to 0.060162602 (sd 2.4e-17) and collapses the `yield` spread from 0.0789 to exactly 0; that `lam` needs no shared calibration at all | ~90 s |
| [`11_linkers_and_mrna.py`](11_linkers_and_mrna.py) | `mhcmatch.vector` | the last step — the named linker presets and why the table refuses to rank itself (two published mechanisms at different positions, pointing opposite ways); assembling with the format already fixed; `mrna()`'s parts map tiling the molecule exactly; that back-translating the whole reading frame in one pass is what repairs the seams the linker created | <1 s |
| [`12_hla_loss_and_coverage.py`](12_hla_loss_and_coverage.py) | `mhcmatch.cassette`, `mhcmatch.portfolio` | what allotype coverage is and why the denominator is the donor's genotype; pricing HLA loss as the exact covariance a lost allele implies, and `captured_loh` on the six TESLA donors; why tumour selectivity is a stated preference and not a refit, given that the shipped model fits normal-tissue expression at its largest positive coefficient | ~60 s |

Runtimes are warm-cache, single core on an M-series Mac. The first run of each notebook also
downloads its reference data (~4 MB for the pmhc panel, ~10 MB for VDJdb, ~30 MB for the three
mimicry references).

## Running them

```bash
pip install 'mhcmatch[notebooks]'          # adds marimo
marimo edit notebooks/01_presentation_and_binder_score.py    # interactive
marimo run  notebooks/01_presentation_and_binder_score.py    # read-only app
python      notebooks/01_presentation_and_binder_score.py    # plain script; prints, no UI
```

Notebook 3 additionally needs the `[precursor]` extra:

```bash
pip install 'mhcmatch[notebooks,precursor]'
```

## Conventions

- Each notebook opens with a markdown cell stating what it demonstrates and what to conclude.
- Every number shown is computed in the notebook. Nothing is transcribed.
- Notebooks demonstrate the **API**. Benchmarks, head-to-head comparisons and result tables live in
  the separate benchmark repository, not here.
- `marimo check notebooks/*.py` is the lint gate; `python notebooks/<name>.py` is the execution gate.
