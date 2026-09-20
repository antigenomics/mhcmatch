# mhcmatch notebooks

Six worked examples, **one per section of the paper**. These are [marimo](https://marimo.io)
notebooks — plain Python files with `@app.cell` decorators, so they diff and review like source
rather than like JSON.

Every notebook bootstraps its own data from the public HuggingFace dataset
[`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data) and caches it, so a fresh
`pip install mhcmatch` is enough to run them. No local paths, no pre-staged files.

## Index

| notebook | paper section | demonstrates | wall | peak RSS |
|---|---|---|--:|--:|
| [`01_binding_prediction.py`](01_binding_prediction.py) | presentation | decompose, restriction, `binder_score`; why `p_binder` is pool-invariant and a within-list percentile is not | 9 s | 0.44 GB |
| [`02_physchem_and_recognition.py`](02_physchem_and_recognition.py) | biophysics | the TCR-facing strip on a published corpus: burial and charge, the position-role evidence whose two readings carry opposite signs, and the six-block recognition axis | 1 s | 0.17 GB |
| [`03_corpus_overlap_and_similarity.py`](03_corpus_overlap_and_similarity.py) | corpus | what the self / thymic / viral corpora share at three levels — 0 whole peptides, 1 face, 97.8 % of masked 3-mers — and why that makes the channels a density rather than a lookup | 1 s | 0.15 GB |
| [`04_epic_across_datasets.py`](04_epic_across_datasets.py) | EPIC | the whole nine-term aggregate from the command line; the per-term decomposition; `p_response` as a prior you own; and the held-out record read off the artifact | 36 s | 0.63 GB |
| [`05_cassette_on_cohorts.py`](05_cassette_on_cohorts.py) | cassette | `rank pairs → cassette select → cassette score` on the 125 manufactured units of a real trial; why the cohort offset matters and what `lam` does and does not guarantee | 12 s | 0.63 GB |
| [`06_assembly_and_safety.py`](06_assembly_and_safety.py) | assembly | the linker presets and why the table refuses to rank itself; `mrna()`'s parts map tiling the molecule exactly | 1 s | 0.10 GB |

Measured warm-cache, single core, on an M-series Mac, running each as a plain script. The **first**
run of a notebook also downloads its reference data — at most ~5 MB per notebook, and the panel in
notebook 1 is shared with the rest.

## Running them

```bash
pip install 'mhcmatch[notebooks]'          # adds marimo and polars
marimo edit notebooks/01_binding_prediction.py    # interactive
marimo run  notebooks/01_binding_prediction.py    # read-only app
python      notebooks/01_binding_prediction.py    # plain script; prints, no UI
```

## Conventions

- **One notebook per section of the paper**, in the paper's order. A reader who has the manuscript
  open should be able to find the notebook for what they are reading without a lookup table.
- **No rival tool is ever run.** These show a user solving their own problem with `mhcmatch` out of
  the box. Head-to-head comparisons and result tables live in the separate benchmark repository, and
  where a dataset here happens to carry another tool's scores those columns are not read.
- **Library where the point is the method, CLI where the point is a run.** Notebooks 1, 2, 3 and 6
  call the Python API because the subject is the algorithm; 4 and 5 shell out to `mhcmatch`, because
  the subject is what you would actually type.
- Each notebook opens with a markdown cell stating what it demonstrates and what to conclude.
- **Every number shown is computed in the notebook.** Nothing is transcribed, which is also why a
  claim in the prose that the output contradicts is a bug — three were caught that way while these
  were written, and each one is now stated as the data has it.
- `marimo check notebooks/*.py` is the lint gate and runs in CI; `python notebooks/<name>.py` is the
  execution gate.

## What happened to the other twelve

The previous set was organised by **library module** — one notebook per `mhcmatch.*` — which meant a
reader of the paper had no way to find the one they wanted, and two things the paper leads on (the
whole EPIC aggregate, and corpus overlap) had no notebook at all. The twelve are folded into these
six: 02/05/06 into notebook 2, 04/07 into notebook 3, 08/09/10/12 into notebook 5, 11 into notebook
6. Notebook 3 of the old set (`precursor_frequency`) mapped to no section of the paper, needed an
extra install, and is retired; `mhcmatch.precursor` is still shipped and documented.
