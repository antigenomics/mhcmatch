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
