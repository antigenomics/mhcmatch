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
search core, MHC-I and MHC-II, human and mouse. Every reference dataset is fetched from
[`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data) on first use, so a fresh
`pip install` runs every example in this file with no manual downloads.

```bash
pip install mhcmatch
mhcmatch bootstrap                                   # optional: pre-fetch the ligand panel (~16 MB)
```

Nothing else has to be downloaded by hand — every reference table is fetched on first use. `bootstrap`
only decides *when*, which matters on a compute node with no outbound network; the four staging
tiers are under [Data](#data).

The library examples below run from a plain `pip install mhcmatch`. The Nextflow pipeline pins its
own version; see [Deployment](#deployment).

**Optional extras.** The base install is `seqtree`, `numpy` and `huggingface_hub` — nothing heavy,
and every model that ships by default runs on it.

| extra | pulls | needed for |
|---|---|---|
| `mhcmatch[esm]` | `torch`, `transformers` | **only** the `esm64_glm` recognition head. Downloads a ~2.4 GB ESM2 checkpoint on first use |
| `mhcmatch[structure]` | `tcren` | the structure-based ΔΔG head |
| `mhcmatch[precursor]` | `vdjmatch` | precursor-frequency estimates |
| `mhcmatch[notebooks]` | `marimo`, `polars` | the worked examples in `notebooks/` |
| `mhcmatch[logo]` | `logomaker`, `matplotlib` | **drawing** a motif logo (`logo.render`). `logo.motif` returns the matrix on the base install |

`torch` is **not** required to score recognition. The default head is the six-block `complement`
score, which is pure numpy, so a user who never installs `[esm]` gets a complete fitted model rather
than a degraded one — `mhcmatch.recognition.score()` just works. Asking for the ESM head without the extra raises a
named error telling you which extra to install; it never silently drops features and returns a
number that looks fine.

```python
from mhcmatch import recognition as rec
rec.default_head("human")            # 'complement' -- no torch involved
rec.score(peps)                      # works on the base install
rec.score(peps, head="esm64_glm")    # ImportError unless mhcmatch[esm] is installed
```

```bash
# rank a donor's neoantigen candidates end to end
mhcmatch rank fasta candidates.fasta --alleles donor.alleles --cls mhc1 --tumor SKCM --out ranked.tsv
```

## Pick your entry point

| your question | command | Python |
|---|---|---|
| Which of these peptides does an allele present? | `mhcmatch predict f.fasta --alleles 'HLA-A*02:01' --cls mhc1` | `predict.predict_fasta` |
| …keeping only conventional binders | `mhcmatch predict ... --rank-threshold wb` | `predict.resolve_rank_threshold` |
| …never dropping candidates in a driver gene | `mhcmatch predict ... --keep-genes 'TP53,KRAS'` | `predict.Keep` |
| …never dropping a validated immunogenic epitope | `mhcmatch predict ... --keep-epitopes builtin` | `predict.Keep` |
| Which allele presents this peptide? | `mhcmatch restriction PEP --calibrated` | `store.restriction` |
| Is it a binder at all, one number? | `mhcmatch binder PEP` | `store.binder_score` |
| What is the IC50, and vs its wild type? | `mhcmatch affinity PEP --allele A --wt WTPEP` | `store.affinity_model` |
| Will a T cell respond to it? | `mhcmatch complement --peptides p.txt` | `complement.score` |
| Rank neoantigen candidates for a donor | `mhcmatch rank fasta ...` | `rank.rank_fasta` |
| How many of the donor's own allotypes present it? | (a `rank` column) | `predict.Prediction.n_alleles_presenting` |
| …with mimicry risk and what each one resembles | `mhcmatch rank ... --extended --annotate` | `mimicry.score` |
| Why did *this* candidate rank there? | `mhcmatch explain PEP --allele A` | — |
| What self / viral / bacterial peptide does it mimic? | `mhcmatch mimics --peptides p.txt` | `mimics.neighbours` |
| Does that mimicry raise or lower the risk, and why? | `mhcmatch mimicry --peptides p.txt` | `mimicry.score` |
| Has this, or something near it, already been tested? | `mhcmatch neoag --peptides p.txt` | `mimicry.annotate` |
| Where in the proteome does it come from? | `mhcmatch source --peptides p.txt --proteome human` | `Proteome.find_sources` |
| Which gene is it from, when the deposit did not say? | `mhcmatch genes pairs.tsv --out annotated.tsv` | `Proteome.assign_genes` |
| Is this gene on in a normal tissue, and where else? | `mhcmatch expression GENE --tissue TISSUE --safety` | `expression.lookup` |
| What does this allele's motif look like? | `mhcmatch logo 'HLA-A*02:01'` | `logo.motif` |
| Which peptides in this protein are presented? | `mhcmatch scan p.fasta --correction bh` | `store.scan_protein` |
| What is the full MHC-II ligand around this core? | `mhcmatch span CORE --protein p.fasta` | `ligand.presented_span` |
| **Which *k* of this donor's candidates go in the cassette?** | `mhcmatch cassette select --candidates pool.tsv -k 20 --tol 3` | `cassette.select` |
| **What is this cassette worth, against one from another donor of another size?** | `mhcmatch cassette score --cassettes c.tsv --pool pool.tsv` | `cassette.score` / `cassette.lam` |
| Build the cassette from ranked candidates | `mhcmatch cassette build --candidates units.tsv --n0 8 --screen` | `vector.select` / `vector.order` |
| …spread over allotype **and** mechanism, not just allotype | — | `vector.select(block=…)` |
| Order units I have already chosen, and pick the spacer | `mhcmatch cassette order --candidates chosen.tsv` | `vector.order` |
| …with the linker already decided, not swept | `mhcmatch cassette order ... --linker GS10` | `vector.order(linker=)` / `vector.assemble` |
| Which linkers are there, and what is each for? | `mhcmatch cassette linkers` | `vector.LINKERS` |
| Turn the finished cassette into an mRNA | `mhcmatch cassette build ... --linker GS10 --mrna c.fa` | `vector.mrna` |
| …and a map of it a viewer can draw | `mhcmatch cassette build ... --map c.tsv --map-json c.json` | `vector.epitope_map` |
| How many *independent* shots is it worth? | (a `cassette score` column) | `portfolio.p_at_least` / `n_effective` |
| Are my own response counts over-dispersed? | — | `portfolio.betabinom_rho` |
| Which candidates can no weighted score ever pick? | — | `portfolio.linearly_supported` |
| Strip frameshift-prone motifs from the CDS | `mhcmatch cassette deslip cassette.fa` | `vector.slippery_sites` |
| Split a peptide into anchor / TCR-facing parts | `mhcmatch decompose PEP` | `store.decompose` |
| How viral-like is it, as a soft sum not a cutoff? | — | `luksza.viral_r` |

Full command reference, grouped by task: [the CLI page](https://antigenomics.github.io/mhcmatch/cli.html).

**`predict` and `rank fasta` drop nothing by default.** `--rank-threshold` takes `sb` / `wb` /
`none` / a percentage, and the tiers are **class-aware** because a bare number cannot be: `2.0` is
the weak cut for class I and the *strong* cut for class II. Measured on class II, a flat `2.0`
keeps **0 of 56** scored pairs, discarding the best window at `%rank 2.364` — an empty table,
returncode 0.

**Two whitelists, because they make two different claims.** `--keep-genes 'TP53,KRAS'` keeps every
candidate in a driver gene; `--keep-epitopes builtin` keeps every candidate that *is* one of the
23,299 peptides an assay has called immunogenic; `--keep-mismatch 1` widens that to one
substitution. Neither is dropped by any `--rank-threshold`, and matched rows carry `keep = 1` plus
`keep_reason` (`gene` / `epitope` / `epitope~1`) — because a row kept for its gene is not evidence
about its peptide, and one flag cannot say which claim held. `builtin` is a **pre-built `seqtree`
index**: it reloads in ~1 ms, so a thousand-sample run pays a read rather than a build and has no
cache to race on.

**Two module names for the cassette, because they are two jobs.** `cassette` *chooses* the units and
scores a finished construct (`select`, `score`, `lam`); `vector` *assembles* one that has been chosen
(`select` over allotype slots, `order`, `assemble`, `mrna`, `epitope_map`, `LINKERS`). The CLI hides
the split behind one `mhcmatch cassette` verb; the Python column above does not, which is why both
names appear in it.

**Chaining `rank` into a cassette takes two flags** — `--prefix` renames the column the next step
looks for, and `rank`'s peptide is the minimal epitope where the assembly wants the long window.
[The tested four-command chain is on the CLI page](https://antigenomics.github.io/mhcmatch/cli/commands.html#rerank-chain).

`predict` is the presentation axis (**is it presented at all**, the NetMHCpan `%Rank_EL` analogue);
`restriction` is the specificity axis (**which allele**). They answer different questions and a
peptide can top one and not the other — `NLVPMVATV` is unambiguously A\*02:01-restricted yet bands
mid-pack against A\*02:01's own ligands.

## Composition is not ranking

Top-*m* by a score maximises the expected **number** of responding units. A vaccine needs the
probability that **at least one** works in *this* donor, and the two agree only if the units respond
independently. They do not: on the adjuvant TNBC mRNA vaccine trial of Sahin et al.
(*Nature* 2026;651:1088–1096) the intra-patient correlation is ρ = 0.124 (p = 1.0×10⁻³), 3.45× the
binomial variance.

So `mhcmatch cassette select` maximises **mean minus variance** of the responding-unit count rather
than the mean, and the objective is derived from that goal rather than fitted to an outcome cohort:

```
H(S) = sum_i [ p_i - (gamma/2) s_i^2 ]  -  gamma sum_{i<j} rho_ij s_i s_j,   s_i = sqrt(p_i(1-p_i))
```

Three inputs, none of them an outcome cohort: `p_i` is the calibrated response probability, `rho` is
one number measured on published per-unit assays, `gamma` is a stated preference (1.0 — one unit of
variance traded for one expected unit, **per unit of the cassette**: a correlated count's mean is
linear in `k` and its variance quadratic, so `gamma` is divided by the design effect `1 + rho(k-1)`
to mean the same trade at every size). `rho_ij` spreads `rho` over pairs by how much two units share
a way of failing.

**`--rule v2` poses the same problem the other way round, and it is the one to read first.** `p_i`
is a probability, so the number of units that respond is a *random variable* and many size-*k* sets
are indistinguishable in it. A sort already maximises the expected count; the sets it cannot tell
apart are not a nuisance, they are the design freedom. So:

> **mhcmatch returns, among all cassettes that are — with stated probability — no worse than the
> ranked list, the one whose units share the fewest ways of failing.**

Four ways of failing, all smooth except one: the restricting **allotype** (discrete, because HLA
is), the source gene's **expression** across tissues, TCR-facing **chemistry**, and BLOSUM-graded
**sequence** similarity of the TCR face. `--not-worse 1.0` returns the sort exactly; lower values
buy diversity and say how often you are willing to be wrong. It is a **per-donor** guarantee — a
cohort-level count needs a tighter floor than intuition suggests.

```bash
mhcmatch cassette select --candidates pool.tsv -k 20 --rule v2 --not-worse 0.7
```

```bash
mhcmatch cassette select --candidates pool.tsv -k 20 --tol 3 --out cassette.tsv
mhcmatch cassette score  --cassettes cassette.tsv --pool pool.tsv
```

```python
from mhcmatch import cassette as CA

c = CA.select(scores, peptides, alleles, k=20, tol=3)     # scores = rank.aggregate_score, WHOLE pool
s = CA.score(scores, peptides, alleles, chosen=c.index,
             pool_scores=scores, pool_peptides=peptides, offset=c.offset)
s["yield"], s["p_at_least"], s["lam"], s["n_effective"]
```

Greedy plus a bounded swap pass, `O(kN)` — and it reaches the brute-force optimum on every pool small
enough to enumerate, which is a test rather than a claim.

**Give it the whole candidate pool, not a shortlist.** `binder` and `expr_lvl` are the two largest
positive coefficients in the shipped model and `expr_norm` is positive too (`mhcmatch rank
--coefficients` prints the sizes, which move at every refit), so a pool already cut on binding and
expression has no range left along them. Measured: on the 46-patient half of the NCI gastrointestinal
screen held out of the EPIC fit — an exhaustive exome screen responding at **0.0144** per mutation —
selection lifts captured responses to **3.92× the base rate** at *k* = 5 (13 of 58 positives against
3.3 expected). On TESLA's *nominated* list, which responds at **0.0612**, every rule sits at the base
rate: the selection had already been done to it.

Full treatment, including why a gradient-boosted score fixes the geometry but not the objective:
[cassette design](https://antigenomics.github.io/mhcmatch/cassette.html) and
[the composition page](https://antigenomics.github.io/mhcmatch/portfolio.html).

### `lam` is what compares two cassettes

`sum p` is a **level** and it is comparable only if every cassette was calibrated together —
`rank.probability` anchors the mean of *the batch it is handed*, so calling it once per donor pins
every donor's pool mean to the declared prevalence. Measured on 7,261 TCGA donors with pools of 1 to
5,221 candidates: every per-donor-anchored mean lands on **0.060163**, sd **2.75 × 10⁻¹⁷**. Read as a
probability, that number is not one.

`lam` needs no shared calibration at all. It is `H(S)` minus the exact log partition function over
every size-*k* subset of that donor's **own** pool, plus `log C(N, k)` — so zero is a uniformly random
subset of the same pool, and both pool depth and *k* divide out. On 3,064 TCGA donors, a cassette
built by sorting the candidate list scores a median **−0.539 nats** — below a random subset — against
**+3.417** for the greedy argmax, a gain of **+4.083**.

## What `rank` costs

`rank --score aggregate` computes **every one** of the model's features before scoring: a model emits
the features it used and refuses to run without them. All three corpus channels — `thymus`, `self`,
`viral` — are affordable because the term contracts a k-mer table rather than searching an index, so
the host-proteome reference index (~7.5 GB, 6 min 15 s) is off the ranking path entirely and
`--no-self` is still allowed with `--score aggregate`. That index is what `--extended` and
`--annotate` cost, because they report *which* reference peptide was hit.

`--score gate` uses the two-term noisy-AND and stays cheaper still.

**The cost that remains is the per-allele `%rank` background, and it is cached.** Building one
allele's calibration background is a 10,000-peptide draw scored under that allele's model, ~0.95 s,
and it is a pure function of `(allele, model, background, footprint, seed, library version)`. The
on-disk cache defaults on at `~/.cache/mhcmatch/calibration`; set `MHCMATCH_CALIBRATION_CACHE` to
relocate or share it, or to `off` to disable. Measured on a 363,324-pair, 2,093-allele build:
**1,788 s → 15 s**, with the two outputs identical in all 40 columns. Budget ~240 kB per
(scorer, allele).

**There is nothing to cache on the corpus path.** `C_corpus` does not search: it contracts a k-mer
frequency table, which is the **exact** Łuksza sum rather than a radius-2 truncation of it, and
costs a 64 KB table per reference deposit instead of a 7.5 GB trie. The tables are memoised per
process and need no lock — see `docs/corpus.rst`.

The indexed search is still there for what genuinely needs it: `features()`, `annotate()` and the
self-mimicry safety scan report *which* reference peptide was hit and from what protein, which a
weighted sum cannot.

## Caching calibration across jobs

Scoring a peptide needs the allele's random-peptide background: 10,000 draws scored through the
model, plus an isotonic fit. That costs ~0.2-3 s **per allele**, once, after which peptides score
at ~80,000/s. In one process it is amortised automatically. Across a SLURM array or a Nextflow
run, every task otherwise repeats it.

```bash
export MHCMATCH_CALIBRATION_CACHE=/shared/scratch/mhcmatch-cal
```

Tasks then share the work: measured **15x** on a 25-allele sweep (13.3 s to 0.9 s). Entries are
written to a temporary file in the same directory and moved into place with `os.replace`, which is
atomic on POSIX and on POSIX-compliant network mounts, so a concurrent reader never sees a partial
file. Two tasks that compute the same allele simultaneously both write and the second rename wins,
which is safe rather than merely tolerated: the payload is a deterministic function of the cache
key, so the racing writers produce identical bytes. There is no lock — a lock would serialise the
fleet to buy nothing.

The key covers the library version, class, background, footprint, head, panel size, draw count,
seed and the positives feeding the isotonic fit. A cache keyed on less than that would be worse
than none, because it would serve a background drawn against a different model as though it were
this one. Set it to `off` and nothing is cached; leave it unset and entries go to `~/.cache/mhcmatch/calibration`. Either directory is disposable -- and use a fresh one for any run meant to establish a number rather than iterate.

## Batch and threads — read this before scripting a loop

**Pass `--peptides FILE` to any peptide-keyed command.** The expensive part of most of them is setup
that a per-peptide invocation pays again every time: the presentation and affinity calibrators are
~5 s, the binder calibrator ~45 s, a human-proteome length index 64.6 s. All of it is cached for the
life of the process, so one process over a whole list is the difference between 49 s per peptide and
thousands per second. The index outlives the process too — it is cached on disk and can be fetched
prebuilt in 3.08 s (`mhcmatch bootstrap --index`), so it is paid once per machine. Measured, both ways, in
`bench/cli/` in [`2026-mhcmatch-code`](https://github.com/repseq/2026-mhcmatch-code) (private; released to reviewers).

```bash
mhcmatch binder     --peptides peptides.txt --alleles "$ALLELES" --top 1 --out binders.tsv
mhcmatch complement --peptides peptides.txt --prior 4.2e-4       --out recognition.tsv
mhcmatch affinity   --peptides pairs.tsv --allele 'HLA-A*02:01'  --out affinity.tsv   # pairs.tsv has
                                                                     # peptide + wt_peptide columns
mhcmatch source     --peptides peptides.txt --proteome human --threads 0 --out sources.tsv
mhcmatch mimics     --peptides peptides.txt --categories thymus,viral,bacterial --threads 0
mhcmatch mimicry    --peptides candidates.tsv --annotate --out risk.tsv   # + what was hit
mhcmatch neoag      --peptides candidates.tsv --out annotated.tsv  # keeps every original column
cut -f1 table.tsv | mhcmatch complement --peptides -              # `-` reads stdin
```

The input is one peptide per line, or a TSV with a `peptide` column (`.gz` fine); the output is TSV
with a header on stdout or `--out`. `--threads` is offered **only** on `source`, `mimics` and
`genes`, whose neighbour search runs in C++ with the GIL released; everywhere else the per-peptide
work is a small numpy product and a thread pool would buy nothing, so the flag is absent rather
than ignored.

## The two axes

Presentation is necessary and not sufficient: most presented peptides are ignored. mhcmatch keeps
the two questions apart and scores them with the fitted **`EPIC`** aggregate, whose `C_phys_*` and
`C_corpus_*` terms are the recognition axis and whose `binder` / `log10a` terms are the presentation
one. It is **hierarchical**: nine columns in four blocks — presentation, expression, physchem,
corpus — entered in pipeline order, so a recognition coefficient is what that term is worth *after*
presentation and expression rather than in competition with them. None of the recognition terms is
fitted on immunogenicity labels.

**Expression enters as two free terms, not one and not a ratio.** `expr_lvl` is what this candidate
is transcribed at and `expr_norm` is the same gene's median in the tumour's matched normal tissue,
both `log2(1 + TPM/c)` on the floor `c` that the tumour type's own transcriptome sets — the 25th
percentile of its non-zero gene medians, 0.1400 to 0.2400 TPM across 35 cancer types. Entering them
separately lets a tumour-versus-normal ratio be *found* rather than imposed, and it is not found: a
difference of logs requires equal and opposite coefficients, and both come back positive.

The unit does not have to be TPM, because `c` is a quantile of the same column and the two cancel —
but only while they are the same column. Where a submitted abundance is on some other scale,
`expression.batch_scale` estimates the factor by median-of-ratios against the reference and
**refuses** unless the input covers half the context's expressed genes. A candidate list cannot
clear that gate, and should not: a mutation reaches one only where the gene was seen in RNA, so the
ratio would measure that conditioning rather than the library. A candidate whose gene is unknown
scores on the terms it does have, flagged, never dropped.

**Both terms are keyed on a gene symbol, and most deposits do not ship one.** Over the neoantigen
corpus the symbol is missing on **356,387 of 695,811 rows (51.2%)** and on **5,205 of 5,833
immunogenic candidates (89.2%)**, and every row without one takes the same mean-imputed value — on
the VACCIMEL screen that left `expr_norm` at standard deviation **exactly 0.0000** and AUROC
**exactly 0.5000**. `mhcmatch genes` recovers the symbol from the peptide, because a neoantigen is a
near-copy of a self peptide: near-exact proteome search, each parent named by its UniProt `GN=`
field. Coverage over that corpus goes to **692,349 of 695,811 rows (99.5%)** and **4,511 of the
5,833 positives** gain a symbol, which takes `expr_norm`'s standard deviation on VACCIMEL from
**0.0000** to **2.520** (`bench/results/epic_gene_repair.md`).

```bash
mhcmatch genes pairs.tsv --species human --out annotated.tsv     # + a `gene` column
mhcmatch rank pairs mouse.tsv --species mouse --tumor B16F10     # mouse: mouse scorer, mouse
                                                                 #   expression, HUMAN corpus
mhcmatch expression Trp53 --species mouse --tissue thymus        # FANTOM5, not GTEx
mhcmatch rank pairs t.tsv --score features --tissue skin         # every fitted column, no score:
                                                                 #   what a refit needs first
mhcmatch rank pairs annotated.tsv --tumor SKCM --out ranked.tsv  # reads it: no join, no rename
```

Three semantics decide whether the answer is the right one. **Only the nearest shell votes** — a
radius-2 shell is ~85× the radius-1 shell inside it, so pooling the two lets a coincidence outvote
a real single-substitution parent. **Ties come back in full, one row each**, because which of
several equally-near parents to score under is a question expression answers and a search cannot;
take the best score per peptide, and expect them — **22,172 of 70,485 8-mers (31.5%)** name more
than one nearest gene, against **7,448 of 97,995 11-mers (7.6%)**. **The radius is 2 because a
neoantigen can carry more than one mutation**, and the second shell is not rounding: **3,004 of the
695,811 rows** find their parent only there.

A **gate** — a product of sigmoids rather than a sum, so a candidate failing either axis cannot be
rescued by the other — is reachable as `mhcmatch rank --score gate`.

**Presentation** — per-allele %rank / `P(present)` / band from a learned anchor model with
cross-allele **pseudosequence diffusion** (rare alleles borrow from groove-similar frequent ones), a
K=3 motif mixture and per-allele register EM for class II; plus a pan-allele **Potts affinity head**
(IC50 nM, Łuksza amplitude `A = Kd_WT/Kd_MT`, DAI). Their calibrated combination is the
**generalized binder score** (`binder_rank`), the recommended single-number binder index.

**Recognition** — `mhcmatch.complement`, a prior-free log-odds over six blocks: physicochemistry and
length; the same components split **MHC-facing vs TCR-facing**; MJ1996 on the anchors and **TCRen
marginalised over 28 M real CDR3 loops** on the TCR-facing side; contiguous-hydrophobic-run motifs;
per-role **residue log-odds**, with per-length and position-zone tables (class I bins 8/9/10/11+ by
relative third of the TCR face; class II bins 14/16/19 by register zone, via `cls="mhc2"`); and
adjacent TCR-facing dipeptides. Fitted per species and never pooled across hosts. Vectorised — a
whole published corpus scores in seconds, so pass a list. `mhcmatch.posbayes` is a strict special
case of it and ships alongside for comparison.

**MJ1996 and TCRen sit on opposite faces because they are different physics**, and the split is
measured rather than assumed: MJ1996 is **96.4% one-body**, a hydrophobicity axis and the right
object for burial in a pocket; TCRen, inverted from 374 TCR:pMHC crystals, is **3.29% one-body**,
below its own composition-matched shuffle floor, with no per-residue scale to extract — which is why
the receptor side is integrated out instead of read off. The potentials and the contact maps are our
own upstream work, [`antigenomics/tcren`](https://github.com/antigenomics/tcren); the tables are
**vendored here**, so `tcren` is a runtime dependency of the optional `[structure]` extra alone.

**The recognition axis reduces to one published scale, and the reduction is measured.** Scoring all
**576** candidate columns (every vendored residue vector × {anchor, TCR} × {sum, mean}) by ΔBIC
*inside* the general model keeps exactly one: **the Rose burial propensity over the TCR face**,
which ships as `C_phys_buried`. Rose's scale is not a hydrophobicity scale — it is the mean fraction
of solvent-accessible area a residue loses on folding ([Rose et al., *Science*
1985](https://doi.org/10.1126/science.4023714)) — so over the exposed face it scores the area a
receptor *could* bury, and because its basis is imported rather than fitted it cannot memorise the
corpus's cysteine artefact. The second column is `C_phys_charge`, Atchley 2005's electrostatic
factor ([Atchley et al., *PNAS* 2005](https://doi.org/10.1073/pnas.0408677102)), selected on its
**residual against Rose** rather than on its own AUROC.

Every number behind those two paragraphs — the shuffle floors, the per-face AUROCs, the 576-column
sweep, the per-column coefficients — is in [docs/burial.rst](docs/burial.rst) and
[docs/complementarity.rst](docs/complementarity.rst), regenerable from `bench/immuno/`.
`mhcmatch rank --coefficients` prints the shipped values; this file deliberately does not transcribe
them, having once carried a superseded one across two refits.

**Evidence that outranks a model.** `mhcmatch.known` carries five reference sets built from the
public deposits — confirmed tumour neoantigens, peptides the screens tested and found
non-immunogenic, IEDB-immunogenic epitopes, the thymic self-immunopeptidome, the viral ligandome. An
exact match is stronger evidence than any score, so `rank` reports it as a flag and floats those
candidates into a tier of their own instead of folding it into the number.

**Pick your tumour type.** `mhcmatch expression --list-contexts` prints all 19 TCGA↔GTEx pairings;
`expression.matched_tissues('BRCA')` gives the matched normal and `expression.lookup(gene,
tumor='BRCA')` the tumour value. **Pass your own tumour type**: it sets the floor both expression
terms are divided by, and a tumour's floor is roughly half its matched normal's, so the pooled
fallback is not a neutral choice. If the origin arrives as free text, `expression.resolve_context`
maps it — `"liver"`, `"LIHC"` and `"hepatocellular"` all resolve, and an unrecognised string raises
rather than quietly returning a number from the wrong distribution.

**Expression, and which normal tissue to compare against.** `--tumor` takes a **TCGA study
abbreviation** (`SKCM`, `LUAD`, …; `CRC` merges TCGA's `COAD` and `READ`) and `--tissue` a **GTEx
`SMTSD`** name (`Skin - Sun Exposed (Lower leg)`). Neither is a clinical coding system — not
ICD-O-3, SNOMED CT or OncoTree — so a pipeline needing one brings its own crosswalk.
`expression.matched_tissues("SKCM")` gives a tumour type's matched normal, which is what makes the
safety read askable without knowing the pairing by heart; `mhcmatch expression --list-contexts`
prints all 19 pairings and, of 104 GTEx tissues in total, the **82** that are no tumour type's
matched normal and are for the safety read only. `HNSC` is marked approximate — it maps to Minor
Salivary Gland / Esophagus - Mucosa, because GTEx has no head-and-neck mucosa.

**Cross-reactivity.** `mhcmatch.mimics` reports near-identical reference peptides per category, and
never sums them, because a hit in each argues something different: **thymus** (presented during
negative selection — tolerance, and autoimmune risk for a vaccine), **self** (encoded but not known
to be presented), **viral** / **bacterial** (a pre-existing repertoire may cross-react, raising
immunogenicity), **neoag** (already tested somewhere).

**Mimicry as risk.** `mhcmatch.mimicry` is the fitted form of that scan: `viral`, `self` and
`thymus`, each split into an **anchor** and a **TCR-facing** channel that partition the peptide, as
six signed log-odds contributions and their sum. A single whole-peptide distance is the wrong
feature, and the sparsity that suggests otherwise is a search artifact. Scores are log-odds;
`probability()` needs a *named* corpus, because a base rate is a property of the pool.

**`mimicry` is a scoring term, not a safety screen.** Flagging candidates by "resembles a
tolerance-side reference" fires on almost everything — influenza `GILGFVFTL` drew 14
essential-tissue hits — because anchor-masked similarity to a *presented* reference is presentation,
not recognition. Exclusion goes through `vector.self_origin_risk`.

**Building the cassette.** `mhcmatch.vector` is the step after ranking: withdraw on safety, order
the units, pick a linker, emit amino acids and a codon-optimised CDS. `LINKERS` is a table of named
presets and `deslip` removes the m1-pseudouridine +1-frameshift motifs a concatemer can create.

```zsh
mhcmatch rank fasta windows.fasta --alleles "$HLA" --out ranked.tsv
mhcmatch cassette build --candidates ranked.tsv --context windows.fasta --n0 8 --screen \
    --fasta cassette.faa --fasta-nt cassette.fna
```

The selection rule, the safety screen and what it does *not* catch, the linker presets and the
cohort calibration are [docs/cassette.rst](docs/cassette.rst), [docs/safety.rst](docs/safety.rst)
and [docs/portfolio.rst](docs/portfolio.rst).

## Model names

### `EPIC` — the shipped model, one letter per *block*

**E**xpression, **P**resentation, **I**mmunogenic **C**omplementarity. Four letters, four blocks,
entered in the pipeline order the table below gives — presentation first, then expression — so a later block's coefficient is what that term is worth *after*
the earlier ones, not in competition with them. Ridge with an unpenalised per-screen intercept at
`tau = 0.25`; `sd`, `z`, `p` and the 95 % CI are a 400-resample cluster bootstrap over
(patient, screen).

| letter | block | columns |
|---|---|---|
| `P` | presentation | `binder`, `log10a` |
| `E` | expression | `expr_lvl`, `expr_norm` |
| `I` | immunogenic — physchem | `C_phys_buried`, `C_phys_charge` |
| `C` | complementarity — corpus | `C_corpus_thymus`, `C_corpus_self`, `C_corpus_viral` |

**The coefficients are not written down here.** They moved with every refit and this table went on
quoting a superseded set for a full release each time. They are printed in full, with bootstrap
intervals, in [docs/models.rst](docs/models.rst) — *generated from the artifacts on every docs
build*, which is the only form of that table that has ever stayed correct. Or ask the artifact
directly, which is the record:

```bash
mhcmatch rank --coefficients                         # every term, its block, its coefficient
mhcmatch rank --holdout                              # held-out AUROC, the grouped CVs, the corpus
mhcmatch rank --coefficients --species mouse         # the mouse class-I fit, not the human one
mhcmatch rank --coefficients --cls mhc2 --species mouse
```

```python
from mhcmatch import rank

rank.models()                       # every shipped fit: model_id, version, release, rows, positives
a = rank.aggregate("mhc1", "human")  # the artifact itself
a["model_id"], a["version"], a["release"], a["features"], a["coef"], a["fit"]["rows"]
```

### Four shipped fits, and a missing one refuses

There is one artifact per `(cls, species, mode)` and **no fallback**: asking for a combination that
was never fitted raises rather than scoring it with another fit's coefficients.

<!-- BEGIN shipped-models (generated by mhcmatch._modeldoc) -->
| `model_id` | model version | release | terms | rows | positives | intercepts | AUROC | how that AUROC is measured |
|---|--:|---|--:|--:|--:|--:|--:|---|
| `mhc1.human.neoantigen` | **11** | 1.6.1 | 9 | 339,599 | 597 | 7 per screen | **0.7102** | leave-one-screen-out, mean |
| `mhc1.mouse.neoantigen` | **5** | 1.13.0 | 9 | 921 | 379 | 61 per reference | **0.6335** | in-sample, within reference |
| `mhc2.human.neoantigen` | **1** | 1.12.0 | 6 | 1,112 | 656 | 157 per reference | **0.6020** | in-sample, within reference |
| `mhc2.mouse.neoantigen` | **3** | 1.12.0 | 6 | 468 | 177 | 30 per reference | **0.5741** | in-sample, within reference |
| `mhc1.human.pathogen` | **1** | 1.13.0 | 5 | 38,106 | 2,634 | 1 per corpus, global | **0.5926** | in-sample, pooled off the logit |
<!-- END shipped-models -->

**The AUROC column is two protocols, not one.** Human class I spans seven independent screens, so
it holds one out whole and is scored on it. The other three are single-deposit fits with no second
screen to hold out, so what they record is an **in-sample within-reference** figure — the slope
term alone, with the fitted per-reference intercepts excluded from the score, macro-averaged over
the references carrying at least three of each class. Do not read the column down, and do not
average it.

**[docs/models.rst](docs/models.rst) is the full record** — every coefficient with its bootstrap
interval, the per-screen held-out table, what each model was fitted on, and the caveats that come
with each one. It is generated from the artifacts on every docs build; this table is generated by
the same code and pinned by `tests/test_modeldoc.py`.

**All four are `neoantigen` fits.** `mhcmatch rank --epitope pathogen` selects a second
immunological mode for a peptide the host does not encode — a viral or bacterial epitope. It drops
the expression block, which is *undefined* rather than missing with no host transcript, and reads
which corpus channels it carries off the artifact's own `features` list rather than off the mode.
No pathogen artifact ships yet: `--score features` computes every column the mode admits,
`--score aggregate` refuses by name until one is installed. See [docs/models.rst](docs/models.rst).

**`release` is not the running library version.** It is the package version the fit was *accepted*
in, and it is stored rather than derived, because a manuscript pins a fit while the library keeps
moving underneath it: `mhc1.human.neoantigen v11 (release 1.6.1)` is a citation and
`mhcmatch 1.13.0` is not. **That fit is pinned and does not get regenerated** — its coefficients,
bootstrap, `loo` and both grouped CVs are what the manuscript cites.

**All four `(cls, species)` cells are fitted from 1.12.0.** An unregistered species still refuses
rather than being served a neighbour's coefficients, which is the mistake the lookup exists to
prevent. `mode` is `neoantigen` on every artifact so far; `pathogen` is a registered spelling with
no fit, because a tumour neoantigen and a pathogen epitope are two mechanisms rather than two
values of one covariate.

**A mouse run reads the human corpus references — all three components, both classes.**
`mimicry.reference_species(species, component)` routes all three mouse components to human, so a
mouse query is matched against the identical `mhc1|{thymus,self,viral}|human|3` tables the human
artifact scores against. The mouse deposits are too small and too groove-skewed to be a reference:
the thymic one is one haplotype (every one of its 2,663 annotated class-I peptides is `H-2Db` or
`H-2Kb`), the viral one samples 9 allotypes against human's 129, and `self` agrees across species
at r = 0.9990 anyway. **Nothing is trained on human data** — a corpus channel is a k-mer density
lookup, and all nine coefficients in `aggregate_mhc1_mouse.json` are fitted on mouse neoantigens.
Human paths are unchanged.
**Expression is not covered by this and must not be**: human and mouse organs and tumours are
different tissues, so `expression.py` stays species-keyed at every rung. See `docs/corpus.rst`.

**What "mouse" means, component by component.** Only presentation/binding and expression rest on
mouse-derived models *and* mouse references. Physicochemistry is species-free by construction, the
whole corpus block reads human tables, and the affinity head is pseudosequence-conditioned rather
than per-species — the coefficients are fitted on mouse throughout, but most of what they index is
not. [docs/models.rst](docs/models.rst) carries the table, one row per component, and it is the
single place that statement lives.

**The two class-II fits carry six terms and no corpus block.** A `C_corpus_*` channel is a density
over a reference set of peptides — thymic, self, viral — and all three deposited sets are class I;
contracting a 15-mer class-II register against a 9-mer density asks the wrong question rather than
answering it weakly. So the block leaves the design, `blocks` lists three entries, and the
corpus-geometry keys are absent rather than declared-and-unused.

**The human class-II fit is a CD4 response model over human self proteins.** 143 of its 1,112 rows
are a cancer and 260 are healthy donors; the largest single disease is type 1 diabetes at 364 rows.
The composition ships inside the artifact as `fit.population`. Ranking class-II *tumour*
neoantigens with it is an extrapolation from that population.

`rank.AGGREGATE_BLOCKS` is the same block structure at runtime. The letters are a mnemonic for the
blocks, **not** the fitting order — presentation enters before expression, and every conditional
coefficient is reported against that order.

### Presentation and affinity are not the same term

**Affinity is not a second presentation term.** Both end up as a `%rank` against the same kind of
background, so the mechanism doesn't separate them — the training data and the target do. The Potts
affinity head is fitted on **measured IEDB IC50**, targeting `Kd`: the biophysics of the groove. The
`AnchorModel` behind `pres` is fitted on the **observed ligand panel**, targeting how *ligand-like* a
peptide is, which carries processing, transport and abundance signal that binding alone does not.
That is the field's binding-affinity vs eluted-ligand split, and the two are measurably not
redundant — on TESLA-608 affinity scores 0.757 AUROC, presentation 0.763, and their Fisher
combination (`binder_rank`) **0.786**. A combination cannot beat both parents by that margin on the
same measurement twice.

**`pres` is a rank, not a similarity search** — worth stating because "presentation" invites the
other reading. It is a score against a random-peptide background (10,000 peptides matched to the
corpus's amino-acid and length distribution). Nothing is retrieved: no reference peptide is looked
up and no anchor-matched protein is searched for. Affinity and `binder_rank` are the same, so the
whole presentation side is **scoring, not retrieval**. The searches are `restriction` (the epitope
panel, anchor-masked), `mimicry` (thymic / viral / proteome windows) and the viral ligandome behind
foreignness.

## Python

```python
import mhcmatch
from mhcmatch import complement, known, mimics

store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")   # auto-fetched from HF, cached

store.restriction("NLVPMVATV", calibrated=True)      # ranked alleles + %rank / P(present) / band
store.binder_score("NLVPMVATV")                      # the single-number binder index
store.scan_protein(my_protein, cls="mhc1")
store.decompose("NLVPMVATV")                         # anchor / TCR-facing split, with X masks

aff = store.affinity_model("mhc1")
aff.predict_ic50("NLVPMVATV", "HLA-A*02:01")             # 18.9 nM (shortlist tier)
aff.amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")   # Kd_WT/Kd_MT (Łuksza eq. 9)

complement.score(peptides)                           # vectorised: pass the list, not a loop
complement.posterior(peptides, prior=4.2e-4)         # the log-odds carries NO prior; supply yours
complement.score(peptides, species="mouse")          # separate table; the hosts are never pooled

from mhcmatch import expression, rank
expression.gene_level("Trp53", species="mouse")      # FANTOM5 tissues + the syngeneic models
rank.aggregate("mhc1", "mouse")                      # the mouse fit; ("mhc2","human") raises

known.lookup("GILGFVFTL")                            # -> 'neoantigen': the FIRST set in
                                                     # SET_NAMES containing it, not the only one
mimics.neighbours(peptides, ref_sets, threads=0)     # threaded C++ neighbour search

pm = mhcmatch.Proteome.from_hf("human")
pm.find_sources(peptides, max_subs=1, threads=0)     # batch; find_source() is the single-query form
pm.assign_genes(peptides)                            # {peptide: [gene, ...]}, ties kept in full
pm.wildtype("NLVPMVATV")                             # the WT counterpart, for agretopicity
```

Full API: [antigenomics.github.io/mhcmatch](https://antigenomics.github.io/mhcmatch/). Twelve
[marimo](https://marimo.io) notebooks in [`notebooks/`](notebooks/README.md) run the workflows end to
end on whole published deposits (`pip install 'mhcmatch[notebooks]'`).

## Data

Everything is fetched on demand from [`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data)
and cached by `huggingface_hub`; `$MHCMATCH_PMHC_DIR` points at a local mirror instead.

**Four staging tiers, each a superset of the need above it.** Everything works with none of them —
they trade disk now for network later, which is the trade a compute node with no outbound route
needs made in advance.

```bash
mhcmatch bootstrap                              # the ligand panel, both tiers         ~16 MB
mhcmatch bootstrap --tier shortlist             # ... one tier only
mhcmatch bootstrap --proteome human,mouse       # + reference proteomes                 51 MB
mhcmatch bootstrap --reference                  # + corpora, known epitopes, mimicry
                                                #   references, expression tables      ~115 MB
mhcmatch bootstrap --index "human:9"            # + a PREBUILT whole-proteome window
                                                #   index, per (species, length)    1.2-2.8 GB
```

`--reference` is the one a cluster wants: it is everything `rank`, `neoag` and `mimicry` read, in one
call. `--index` is separate because it is GB-scale and only the **safety screen** and the **mimicry
annotation** need it — both ship off (see [Deployment](#deployment)). Fetching one takes 3.08 s where
building it locally takes 64.6 s (human, L=9); ask for the lengths you will use, class I being 8-11,
and all eight published indexes together are 17 GB. Anything not published falls back to building
locally, so the command never fails for want of an upload.

**Where an index lands is not where you would guess.** A *published* one is read in place from
`$MHCMATCH_PMHC_DIR` if a mirror has it and downloaded into `$HF_HOME` if not; only one this machine
*built* goes to `$MHCMATCH_CALIBRATION_CACHE/proteome_index/`. Measured, mouse `L=9`: 1.49 GB, 0.0 s
off a mirror against 47.0 s from the hub. Point all three at shared, roomy paths on a cluster —
`$HF_HOME` especially, or ~250 MB of ordinary reference data arrives in your home quota. `bootstrap
--index` prints the size, the wall clock and the directory for every length, so you can see which
route you got.

Pseudosequences (34-mer grooves) and the fitted model parameters are vendored in
`src/mhcmatch/data/` with their `PROVENANCE.md`. Nothing is refitted at import.

## Deployment

`integrations/nextflow/mhcmatch/` is a self-contained nf-core-style module — **nine processes**
(`MHCMATCH_ALLELES`, `_PREDICT`, `_RANK`, `_RERANK`, `_NEOAG`, `_MIMICRY`, `_CASSETTE_SELECT`,
`_CASSETTE`, `_CASSETTE_SCORE`), **two arms** that chain them, and **`pipeline.nf`**, a runnable
entry point over a directory of files:

```bash
nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
    --indir donor_files --outdir results --mode both \
    --mhcmatch_cassette_k 20 --mhcmatch_tumor SKCM
```

| `--mode` | in | the deliverable is |
|---|---|---|
| `rerank` | your candidate table (+ the window FASTA it came from) | **your** table — every column intact, in your order — plus an `mm_` block, re-sorted by the aggregate |
| `denovo` | your mutation-window FASTA | **our** table: binding called from scratch, ranked, annotated |
| `both` | both | both, independently; each arm builds its own cassette |

Both arms end in a cassette, published as six files per donor and arm:

| file | what |
|---|---|
| `<id>.<arm>.vaccine.units.tsv` | the *k* selected **epitopes** (default **k = 20**) — **your own table filtered to what the cassette carries, with nothing removed from the row**: every one of your columns, every `mm_` column, plus the selection's own |
| `<id>.<arm>.cassette.faa` | the assembled construct, with whichever spacer the junction sweep chose |
| `<id>.<arm>.cassette.fna` | its CDS, deslipped |
| `<id>.<arm>.cassette.map.tsv` / `.map.json` | unit / linker / epitope in 1-based coordinates. Both carry the feature rows; **only the JSON carries the per-unit summary**, which is where `self_help` is — `summary.n_units_with_self_help` |
| `<id>.<arm>.cassette.tsv` | the assembly **report** (`section, i, key, value, detail`) — where the safety screen records what it withdrew and why, **when it is enabled**; with the shipped default the `withdrawn` section is empty because no screen ran. **Not** the epitope table |
| `cohort.<arm>.cassette_score.tsv` | one per run and per arm, because `rank` anchors `p_response` on the batch it is handed |

The construct carries no more *units* than *k* and usually fewer — several epitopes can share one
27-mer window, and the safety screen withdraws some **when it is enabled** — so read `units=` from
the FASTA header rather than assuming *k*.

> ### The safety screen is OFF by default. Turn it on before you manufacture anything.
>
> `--mhcmatch_vector_screen true`. With it off — the shipped default — **no unit is
> withdrawn for essential-tissue self-origin** and the cassette carries whatever it was handed.
> Every `MHCMATCH_CASSETTE` task prints that no screen ran, so the absence is never silent, but a
> log line is not a substitute for reading this.
>
> It is off because it cost hours rather than because it is optional: the whole-proteome index it
> needs was rebuilt inside every task. That is fixed — the index is published and loads in **3.08 s**
> against 64.6 s to build — so turning it on is now cheap. Stage it first with
> `mhcmatch bootstrap --index "human:8|9|10|11"`.
>
> The **mimicry annotation** (`--mhcmatch_mimicry`) is off for the same reason and is a different
> case: it is annotation only, and **scores are identical either way**.

> **The pipeline pins its own mhcmatch version.** `environment.yml`, the `Dockerfile`,
> `params.mhcmatch_container` and `templates/setup.sbatch` all name it, and `setup.sbatch` asserts
> what it installed rather than letting the run discover a mismatch inside a task log.

On a cluster, start from **`integrations/nextflow/mhcmatch/templates/`** — four SLURM scripts
(`setup.sbatch` once, then `run_human.sbatch` / `run_mouse.sbatch` for a few samples or
`run_slurm_head.sbatch` for a cohort), each with a single `EDIT THESE` block at the top and nothing
cluster-specific below it.

**`pipeline.nf` is the easy entry point, not the integration surface.** A pipeline that already does
variant calling, HLA typing and expression quantification should `include` the processes in
`main.nf` or the arms in `subworkflows/` into its own channel topology.

The per-process input/output contract, the column rules (`--passthrough` refuses a collision,
`cassette select` renames and says so), every parameter, the SLURM templates and the Docker build
are [`integrations/nextflow/mhcmatch/README.md`](integrations/nextflow/mhcmatch/README.md) and
[docs/pipeline.rst](docs/pipeline.rst). They are the contract; this section is the orientation.

## Benchmarks

> Harness and result tables live in
> [`2026-mhcmatch-code`](https://github.com/repseq/2026-mhcmatch-code) (private; released to reviewers). Paths like
> `bench/results/...` resolve there.

Head-to-head against **NetMHCpan-4.2b** / **NetMHCIIpan-4.3i** on the same per-(peptide, allele)
task, stratified by allele rarity, with bootstrap CIs and paired significance
(`bench/compare/SOURCES.md` for provenance and caveats):

- **Immunogenicity ranking on TESLA-608** (608 candidates, 37 T-cell-validated; predictor-agnostic,
  every tool scores it independently) — mhcmatch's `binder_score` **0.786** AUROC vs NetMHCpan
  **0.747**; each single head also beats it (affinity 0.757, presentation 0.763).
  `bench/results/immuno_binder_score.md`.
- **Allele specificity, MHC-I** — mhcmatch beats NetMHCpan on medium and frequent alleles on AUROC,
  AUPRC and PPV@P (all p < 0.001; frequent AUPRC 0.852 vs 0.769). Rare is a wash (p = 0.41).
- **Presented-vs-random screening** — mhcmatch wins MHC-I frequent (AUPRC 0.881 vs 0.846, p = 0.001);
  medium and rare sit inside the CI. Both tools are ≥ 0.97 here.
- **MHC-II** — mhcmatch wins the **rare** stratum on both tasks; NetMHCIIpan leads medium and
  frequent. That gap is **one locus, not the class**: DP averages −0.305 AUPRC while **DR is at
  parity or better (+0.010)**, and the mechanism is a register-EM convergence failure on
  DPA1\*02:01 that `register_em="converge"` closes 28 % of. `bench/results/register_em_convergence_dp.md`.
- **Mouse MHC-II** — mhcmatch wins all nine cells on specificity.
- **Speed** — MHC-I ~195k peptide-allele scores/s (~68× NetMHCpan); MHC-II ~19k/s (~6.6×
  NetMHCIIpan), heavier because of 3 mixture components × ~7 register frames.
- **Recognition** — the complementarity score beats the shipped `posbayes` sum on all four corpus
  arms and both hosts under peptide-grouped CV; the per-length and relative-position role tables add
  +0.007 to +0.021 AUROC on top, with paired bootstrap CIs excluding zero on every arm.
  `bench/results/complementarity.md`, `bench/results/length_roles.md`.

Read the class-II numbers with `compare/SOURCES.md` in hand: NetMHCIIpan trained on essentially all
public IEDB eluted-ligand data, so its in-corpus medium/frequent strata are contaminated in its
favour and the rare / zero-shot axis is the fair comparison.

### Naming a class-II restriction

Class-II alleles are reported in NetMHCIIpan's own form — `DRB1_0101` for DR (DRA is monomorphic, so
the beta chain names the molecule) and `HLA-DQA10501-DQB10301` for the DP/DQ heterodimers. The two
forms **do not lead with the same chain**: DR leads with its beta, DP and DQ with their alpha. Code
that compares two callers by pulling the leading gene out of the key is therefore matching DR's beta
against DP/DQ's alpha, and splitting DR against itself whenever the DRB gene differs (`DRB1` vs
`DRB3`). Both mistakes are easy to make and neither announces itself.

`--mhc2-report` picks the granularity, on every command that *chooses* an allele (`restriction`,
`binder`, `scan`, `predict`, `rank`):

```bash
mhcmatch restriction PKYVKQNTLKLAT --cls mhc2 --mhc2-report isotype
```

| mode | `DRB1_0101` becomes | `HLA-DQA10501-DQB10301` becomes | use it for |
|---|---|---|---|
| `pair` (default) | `DRB1_0101` | `HLA-DQA10501-DQB10301` | reporting, and string-comparing against NetMHCIIpan |
| `beta` | `DRB1*01:01` | `DQB1*03:01` | comparing alleles across isotypes on the same chain |
| `isotype` | `DR` | `DQ` | "did the two callers pick the same molecule family" |

`mhcmatch.pseudoseq.class2_report(key, mode)` is the same reduction from Python. Commands that are
*handed* an allele (`affinity`, `explain`, `logo`) echo back what the caller typed.

## Development

```bash
bash setup.sh            # repo-local .venv + editable install (uses a sibling ../seqtree if present)
bash setup.sh --tests    # + pytest
pytest -q
```

Theory and derivations are in the manuscript repo (`appendix/mhcmatch.tex`); what is planned and what
is in flight is in [`ROADMAP.md`](ROADMAP.md) and [`CHANGELOG.md`](CHANGELOG.md).
