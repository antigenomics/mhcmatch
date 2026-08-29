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
mhcmatch bootstrap                                   # pre-fetch the panel (optional; ~16 MB)
```

**Optional extras.** The base install is `seqtree`, `numpy` and `huggingface_hub` — nothing heavy,
and every model that ships by default runs on it.

| extra | pulls | needed for |
|---|---|---|
| `mhcmatch[esm]` | `torch`, `transformers` | **only** the `esm64_glm` recognition head. Downloads a ~2.4 GB ESM2 checkpoint on first use |
| `mhcmatch[structure]` | `tcren` | the structure-based ΔΔG head |
| `mhcmatch[precursor]` | `vdjmatch` | precursor-frequency estimates |
| `mhcmatch[notebooks]` | `marimo` | the worked examples in `notebooks/` |
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
| Which of these peptides does an allele present? | `mhcmatch predict f.fasta --cls mhc1` | `predict.predict_fasta` |
| Which allele presents this peptide? | `mhcmatch restriction PEP --calibrated` | `store.restriction` |
| Is it a binder at all, one number? | `mhcmatch binder PEP` | `store.binder_score` |
| What is the IC50, and vs its wild type? | `mhcmatch affinity PEP --wt WTPEP` | `store.affinity_model` |
| Will a T cell respond to it? | `mhcmatch complement --peptides p.txt` | `complement.score` |
| Rank neoantigen candidates for a donor | `mhcmatch rank fasta ...` | `rank.rank_fasta` |
| How many of the donor's own allotypes present it? | (a `rank` column) | `predict.Prediction.n_alleles_presenting` |
| …with mimicry risk and what each one resembles | `mhcmatch rank ... --extended --annotate` | `mimicry.score` |
| Why did *this* candidate rank there? | `mhcmatch explain PEP --allele A` | — |
| What self / viral / bacterial peptide does it mimic? | `mhcmatch mimics --peptides p.txt` | `mimics.neighbours` |
| Does that mimicry raise or lower the risk, and why? | `mhcmatch mimicry --peptides p.txt` | `mimicry.score` |
| Has this, or something near it, already been tested? | `mhcmatch neoag --peptides p.txt` | `mimicry.annotate` |
| Where in the proteome does it come from? | `mhcmatch source --peptides p.txt --proteome human` | `Proteome.find_sources` |
| Is the gene on in the tumour, and in normal tissue? | `mhcmatch expression --list-contexts` | `expression.lookup` |
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
| How many *independent* shots is it worth? | (a `cassette score` column) | `portfolio.p_at_least` / `n_effective` |
| Are my own response counts over-dispersed? | — | `portfolio.betabinom_rho` |
| Which candidates can no weighted score ever pick? | — | `portfolio.linearly_supported` |
| …and a map of it a viewer can draw | `mhcmatch cassette build ... --map c.tsv --map-json c.json` | `vector.epitope_map` |
| Strip frameshift-prone motifs from the CDS | `mhcmatch cassette deslip cassette.fa` | `vector.slippery_sites` |
| Split a peptide into anchor / TCR-facing parts | `mhcmatch decompose PEP` | `store.decompose` |
| How viral-like is it, as a soft sum not a cutoff? | — | `luksza.viral_r` |

Full command reference, grouped by task: [the CLI page](https://antigenomics.github.io/mhcmatch/cli.html).

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
a way of failing: the same allotype, the same 3-mers, the same place on the dominance axis.

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
key, so the racing writers produce identical bytes. There is no lock -- a lock would serialise the
fleet to buy nothing.

The key covers the library version, class, background, footprint, head, panel size, draw count,
seed and the positives feeding the isotonic fit. A cache keyed on less than that would be worse
than none, because it would serve a background drawn against a different model as though it were
this one. Unset the variable and nothing is cached; the directory is disposable.

## Batch and threads — read this before scripting a loop

**Pass `--peptides FILE` to any peptide-keyed command.** The expensive part of most of them is setup
that a per-peptide invocation pays again every time: the presentation and affinity calibrators are
~5 s, the binder calibrator ~45 s, a human-proteome length index ~70 s. All of it is cached for the
life of the process, so one process over a whole list is the difference between 49 s per peptide and
thousands per second. Measured, both ways, in
[`bench/cli/`](https://github.com/antigenomics/2026-mhcmatch-benchmark).

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
with a header on stdout or `--out`. `--threads` is offered **only** on `source` and `mimics`, whose
neighbour search runs in C++ with the GIL released; everywhere else the per-peptide work is a small
numpy product and a thread pool would buy nothing, so the flag is absent rather than ignored.

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

**MJ1996 and TCRen sit on opposite faces because they are different physics, and the difference is
measured rather than assumed.** A generic contact potential is essentially additive: MJ1996 is
**96.4% one-body**, its two leading modes correlating with Kyte–Doolittle at exactly **±0.851** — a
hydrophobicity axis, which is the right object for burial in a pocket. TCRen, inverted from 374
TCR:pMHC crystals, is **3.29% one-body**, *below* its own composition-matched shuffle floor
(**9.68 ± 2.16%**, 500 shuffles), and neither leading mode is a hydropathy axis (**+0.353** on the
receptor side, **−0.295** on the peptide side). There is no per-residue TCRen scale to extract,
which is why the receptor side is integrated out instead of read off. The face assignment is then
checked against coordinates: predicting *does this side chain reach the groove floor* over 3,875
(structure, position) rows from the same 374 crystals, MJ1996 alone reaches **AUROC 0.5818** and
TCRen alone **0.4801** — below chance — which is the ordering the block assumes. On class II
neither separates alone (**0.5171** MJ, **0.5368** TCRen over 94 structures) and that ordering is
not reproduced, so it is not claimed there.

Both are lossy summaries, and they are carried for the distinction they draw rather than for the
ranking they move: on top of the geometric prior the two scalars add **+0.0012** AUROC against a
free 20-way one-hot's **+0.0062** on the same task, and in the shipped complement head their four
columns are small next to the fitted-identity ones — per column standard deviation, `mj_anchor`
**+0.0354**, `mj_tcr` **−0.0927**, `para_tcr` **−0.0299** and `para_sd_tcr` **+0.0052**, against
`aa_tcr`'s **−0.8796** and `kmer_llr`'s **+0.6544**. The potentials, the 374-crystal contact maps
and the spectral analysis are our own upstream work,
[`antigenomics/tcren`](https://github.com/antigenomics/tcren); the tables are **vendored here**, so
`tcren` is a runtime dependency of the optional `[structure]` extra alone.

**The recognition axis reduces to one published scale, and the reduction is measured.** Split into
its chemistry half and its fitted-identity half — exact partial sums via `score(blocks=...)` — the
two behave oppositely: the identity half wins in-corpus and the chemistry half transfers. Scoring
all **576** candidate columns (every vendored residue vector × {anchor, TCR} × {sum, mean}) by ΔBIC
*inside* the general model keeps exactly one: **the Rose burial propensity summed over the TCR
face**, at z **+4.57** — the second-largest coefficient of the ten-term model it was selected in,
behind expression alone — against the sixteen-column chemistry block's +0.18 in the same slot. In
the shipped eight-term model it is `C_phys_buried`, the larger of the block's two fitted terms, and
smaller than at selection because the corpus block now carries three coefficients it previously
shared one with (`mhcmatch rank --coefficients` prints the current value; this file no longer
transcribes it, having carried a superseded one across two refits). Rose's scale is not a hydrophobicity scale
— it is the mean fraction of solvent-accessible area a residue loses on folding ([Rose et al.,
*Science* 1985](https://doi.org/10.1126/science.4023714)) — so summed over the exposed face it
scores the area a receptor *could* bury. Because its basis is imported rather than fitted it cannot
memorise the corpus's cysteine artefact: correlation with per-peptide cysteine count is **+0.108**
against the shipped score's **+0.688**. Full derivation in [docs/burial.rst](docs/burial.rst) and
§11 of the theory appendix; all of it regenerable from `bench/immuno/` in the benchmark repo.

**The second chemistry column is charge, and what it buys is burial's stability.** `C_phys_charge`
is Atchley 2005's fifth factor — electrostatic charge ([Atchley et al., *PNAS*
2005](https://doi.org/10.1073/pnas.0408677102)), read as the mean over the TCR face — and it was
selected on its **residual against Rose**, not on its own AUROC. That is the point of it: all 39
transfer scales swept correlate **0.74 to 0.95** with burial over the twenty residues, so a
hydropathy scale is burial measured a second way and the two are not identified. v3 paired Rose
with Kidera KF4 at *r* = **−0.837** per peptide; AF5 sits at **+0.008**. Swapping the partner
leaves burial's coefficient *smaller* and its evidence *stronger*, which is what dropping a
collinear term does: bootstrap sd **0.0874 → 0.0487**, *z* **+1.71 → +2.34**, *p* **0.088 →
0.020**, sign stability **96.5% → 100%**, over 400 cluster bootstraps on the 354,909 rows and 958
immunogenic peptides the corpus held when the swap was measured. AF5 also carries the **lowest cysteine loading of the 141 complete
residue scales swept, −0.0028**, which matters because the Chowell family runs a 12.5×
mass-spectrometry cysteine enrichment a fitted basis can learn. Its own coefficient is the smaller of the
two — the column being fixed here is burial, not the one swapped in. `mhcmatch.complement.PHYS_SCALE_CHARGE` names the scale;
[docs/burial.rst](docs/burial.rst) carries the arms table.

Four opt-in parameters came out of that work, all defaulting to the shipped behaviour so no recorded
number moves: `blocks=` (score a subset of blocks, exact partial sum), `mask_cys=` (zero cysteine in
the fitted tables, as `posbayes` does by construction), `positions="profile"` (read the chemistry
over the crystallographic per-position TCR-contact profile rather than a binary anchor mask), and
`paratope="contact"` (marginalise TCRen over the receptor residues that actually contact, rather
than over the whole CDR3 loop).

The other half of that axis is not chemistry and is not fitted on labels either. `C_aa` — the
residue log-odds — is estimated on Chowell, which separates **foreign** from **self and presented**:
a statement about *passing thymic selection*, not about whether a T cell responds to a somatic
mutation. `mimicry.corpus_R` replaces it with a label-free neighbour density against three
references, split by *when a T cell meets them*:

| channel | what it is | reads as | sign in the shipped fit |
|---|---|---|:--:|
| `thymus` | thymic immunopeptidome — **the only one that enters selection** | danger | + |
| `self` | host proteome, met in the periphery | tolerance | − |
| `viral` | foreign ligandome — **never seen during selection** | reference | + |

(`mhcmatch rank --coefficients` for the magnitudes. Signs are what the argument rests on, and all
three have held at 100 % of 400 cluster bootstrap resamples across every refit.)

The thymic channel is positive because the thymus is not a random sample of self: mTECs
promiscuously express tissue-restricted antigens under *Aire* and *Fezf2* precisely to purge the
clones that would cause autoimmunity, so thymic display is enriched for the self **worth tolerising
against**. Measured on the burial axis, thymic ligands sit above *non-thymic presented* self at
Cohen's d = +0.1650 (p = 1.0×10⁻⁸⁰) — presentation held constant, so the effect is thymus-specific.
`thymus` and `self` are both similarity to *self* sets, and their opposite signs are what no
single-mechanism account gives. See [docs/corpus.rst](docs/corpus.rst).

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
feature and the sparsity that suggests otherwise is a search artifact — whole-peptide radius-2
thymic coverage is 1.63 %, while the TCR face at radius 1 reaches **53.4 %**. Signs follow the
reference, as designed: `viral` positive on both channels, `self` negative on both, `thymus`
positive on its anchor. `MimicryScore.nearest` carries *which* peptide was hit and what protein it
came from, so `mimicry.safety()` reaches `expression.safety_profile` — a bare distance cannot.
**Scores are log-odds**; `probability()` needs a *named* corpus, because the seven screens behind
the calibration run from 0.048 % to 46.8 % positive. Report the **within-screen** AUROC (0.596), not
the pooled one (0.849). The tested-neoantigen database is `mimicry.annotate` / `mhcmatch neoag` —
prior evidence, and deliberately never a fitted term, since every labelled screen we hold sits
inside it.

**Building the cassette.** `mhcmatch.vector` is the step after ranking: **withdraw on safety**, then
how many units each allotype carries, in what order, joined by what. `screen()` **excludes** rather
than down-ranks — the second-best cassette is cheap, myocarditis is not — and it screens *every*
register of a 27-mer unit, not the mutated one, against near-exact self origin joined to tissue
expression. `select()` grows each allotype while the next candidate beats that allotype's own
expected yield per slot, so diversification falls out of the arithmetic instead of a quota;
`order()` tries **no spacer first** and picks the layout minimising the strongest predicted binder
spanning each junction; `slippery_sites()`/`deslip()` remove the m1Ψ +1-frameshift motif, which
matters more for a concatemer than for a natural ORF.

**Choosing a linker, and building the mRNA.** `LINKERS` is a table of named presets — GS-rich
flexible, class-I favouring, class-II oriented, minimal, protease-cleavable, rigid — each carrying
its family, the class it is *intended* for and where it comes from. `mhcmatch cassette linkers`
prints it and `--linker NAME` pins one instead of sweeping. **The table does not rank itself:** the
two mechanisms that would settle a class-I ranking act at different positions — Gly and Pro are
abundant in the C-terminal regions class-I ligands are cleaved from (PMID 30645615), yet the same
residues immediately flanking a class-I epitope inhibit recognition of the epitope on their
amino-terminal side and move the ratio of two responses from one construct up to fifty-fold
(PMID 8871618) — so `order()` measures each candidate
against the recipient's own allotypes and that is what selects. `GS10` is a **reconstruction**: the
manufactured pentatope format is described as a "non-immunogenic 10-mer glycine/serine linker" and
not published residue by residue.

`mrna()` then assembles the molecule and returns a parts map that tiles it exactly, in nucleotides:
5' UTR, start, leader, units, linkers, trailer, stop, 3' UTR, poly(A). The **coding sequence is
generated for the whole reading frame in one pass**, so homopolymer avoidance and the m1Ψ
+1-frameshift repair act across the seams the designer created rather than inside each unit — which
is where a concatemer's problems actually are. The **backbone is not supplied and defaults to
nothing**: a UTR belongs to a particular vector and a plausible invented one is worse than none.
`checks` reports numbers rather than a verdict, and `translates` — the coding sequence read back in
the frame the construct sets it in — is the one that must hold.

Both ends join to the rest of the library. `--context windows.fasta` takes `rank`'s **minimal
epitopes** and rebuilds them as long units against the FASTA they were called on, one per variant
rather than one per register — a minimal peptide loads onto any cell without costimulation and is
the tolerising configuration, so the reader will not take one. `--fasta-nt` writes the epitope
cassette's coding sequence alone: highest-usage human codon per residue, backed off to shorten
homopolymers, then deslipped. It is **not a codon optimiser** — it fixes the two things that break a
concatemer specifically and leaves GC, structure and CpG to the manufacturer's tooling.

```zsh
mhcmatch rank fasta windows.fasta --alleles "$HLA" --out ranked.tsv
mhcmatch cassette build --candidates ranked.tsv --context windows.fasta --n0 8 --screen \
    --fasta cassette.faa --fasta-nt cassette.fna
```

⚠️ **`mimicry` is a scoring term, not a safety screen.** Flagging candidates by "resembles a
tolerance-side reference" fires on almost everything — influenza `GILGFVFTL` drew 14
essential-tissue hits — because anchor-masked similarity to a *presented* reference is presentation,
not recognition. Exclusion goes through `vector.self_origin_risk`
(`bench/results/vector_safety_screen.md`).

## Model names

### `EPIC` — the shipped model, one letter per *block*

**E**xpression, **P**resentation, **I**mmunogenic **C**omplementarity. Four letters, four blocks,
entered in that pipeline order — so a later block's coefficient is what that term is worth *after*
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
quoting a superseded set for a full release each time. Ask the artifact, which is the record:

```bash
mhcmatch rank --coefficients     # every term, its block, its coefficient
mhcmatch rank --holdout          # per-screen AUROC, the two grouped CVs, the fit's own corpus
```

```python
import json, importlib.resources as R
d = json.loads(R.files("mhcmatch.data").joinpath("aggregate_mhc1.json").read_text())
d["model"], d["version"], d["features"], d["coef"], d["fit"]["rows"], d["fit"]["screens"]
```

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
aff.predict_ic50("NLVPMVATV", "HLA-A*02:01")             # 52.5 nM (shortlist tier)
aff.amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")   # Kd_WT/Kd_MT (Łuksza eq. 9)

complement.score(peptides)                           # vectorised: pass the list, not a loop
complement.posterior(peptides, prior=4.2e-4)         # the log-odds carries NO prior; supply yours
complement.score(peptides, species="mouse")          # separate table; the hosts are never pooled

known.lookup("GILGFVFTL")                            # -> 'neoantigen': the FIRST set in
                                                     # SET_NAMES containing it, not the only one
mimics.neighbours(peptides, ref_sets, threads=0)     # threaded C++ neighbour search

pm = mhcmatch.Proteome.from_hf("human")
pm.find_sources(peptides, max_subs=1, threads=0)     # batch; find_source() is the single-query form
pm.wildtype("NLVPMVATV")                             # the WT counterpart, for agretopicity
```

Full API: [antigenomics.github.io/mhcmatch](https://antigenomics.github.io/mhcmatch/). Nine
[marimo](https://marimo.io) notebooks in [`notebooks/`](notebooks/README.md) run the workflows end to
end on whole published deposits (`pip install 'mhcmatch[notebooks]'`).

## Data

Everything is fetched on demand from [`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data)
and cached by `huggingface_hub`; `$MHCMATCH_PMHC_DIR` points at a local mirror instead.

```bash
mhcmatch bootstrap                              # the reference ligand panel, both tiers (~16 MB)
mhcmatch bootstrap --proteome human,mouse       # + reference proteomes
mhcmatch bootstrap --reference                  # + corpora, known-epitope, mimicry, expression (~115 MB)
```

Pseudosequences (34-mer grooves) and the fitted model parameters are vendored in
`src/mhcmatch/data/` with their `PROVENANCE.md`. Nothing is refitted at import.

## Deployment

`integrations/nextflow/mhcmatch/` is a self-contained nf-core-style module — **six processes**
(`MHCMATCH_PREDICT`, `_RANK`, `_NEOAG`, `_MIMICRY`, `_CASSETTE`, `_CASSETTE_SCORE`) plus a
`subworkflows/mhcmatch.nf`
that chains them, with `nextflow.config`, a `slurm.config` executor profile, `environment.yml` and a
`Dockerfile`. `PREDICT` drops in for MHCflurry class I and the class-II binding subworkflow,
consuming the same `(meta, peptide.fasta, alleles)` channel and emitting a pipeline-compatible
57-column `.scored.csv`; the other five cover ranking, prior evidence, safety, cassette assembly
and cassette scoring,
which have no incumbent. No stub types a column header — each asks the installed library for its own
schema, so `-stub-run` cannot drift from the real shape. The image bootstraps its panel at **build**
time, so compute nodes need no network.

`slurm.config` sizes each process to what it measurably consumes and points every task at one shared
reference and calibration directory — without it a 200-sample run re-derives the same per-allele
background 200 times. `integrations/nextflow/mhcmatch/README.md` is the full contract, including the
two things a cluster gets wrong first: the partition name has no safe default, and the wheel needs
Python ≥ 3.10 where a cluster's system `python3` is often older.

## Benchmarks

> Harness and result tables live in
> [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark). Paths like
> `bench/results/...` resolve there.

Head-to-head against **NetMHCpan-4.2b** / **NetMHCIIpan-4.3i** on the same per-(peptide, allele)
task, stratified by allele rarity, with bootstrap CIs and paired significance
(`bench/compare/SOURCES.md` for provenance and caveats):

- **Immunogenicity ranking on TESLA-608** (608 candidates, 37 T-cell-validated; predictor-agnostic,
  every tool scores it independently) — mhcmatch's `binder_score` **0.786** AUROC vs NetMHCpan
  **0.747**; each single head also beats it (affinity 0.757, presentation 0.763).
  `bench/results/immuno_binder_score.md`.
- **Allele specificity, MHC-I** — mhcmatch beats NetMHCpan on medium and frequent alleles on AUROC,
  AUPRC and PPV@k (all p < 0.001; frequent AUPRC 0.850 vs 0.769). Rare is a wash (p = 0.41).
- **Presented-vs-random screening** — mhcmatch wins MHC-I frequent (AUPRC 0.881 vs 0.846, p = 0.001);
  medium and rare sit inside the CI. Both tools are ≥ 0.97 here.
- **MHC-II** — mhcmatch wins the **rare** stratum on both tasks; NetMHCIIpan leads medium and
  frequent. That gap is **one locus, not the class**: DP averages −0.305 AUPRC while **DR is at
  parity or better (+0.010)**, and the mechanism is a register-EM convergence failure on
  DPA1\*02:01 that `register_em="converge"` closes 28 % of. `bench/results/register_em_convergence_dp.md`.
- **Mouse MHC-II** — mhcmatch wins all nine cells on specificity.
- **Speed** — MHC-I ~195k–260k peptide-allele scores/s (~68× NetMHCpan); MHC-II ~19k/s (~6.6×
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
