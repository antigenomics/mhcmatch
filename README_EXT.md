<!-- Split out of README.md, which had grown to 846 lines. This file holds the sections that
     answer "why does it do that" rather than "how do I run it"; the README holds the second kind.
     Nothing here was rewritten in the split -- these are the same sections, moved. -->

# mhcmatch — the extended README

The [README](README.md) is what you need to run the thing. This is the reasoning behind it: why a
cassette is a set problem and not a ranking, where the time goes, what the two scoring axes each
carry, and what has been measured against what.

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

**Two ways, where an earlier design had three.** The channels are the shared **allotype** and
the shared **sequence** — both mechanisms. A third, *score dominance*, coupled two units for scoring
alike, which is not a mechanism, and it never abstained: over every within-donor pair of the TESLA
and HiTIDE pools it is zero on 0.03% and 0.01% of pairs against 97.5% and 96.5% for the sequence
channel, so it supplied **71–79% of the total channel mass** and the allotype channel — the only
mechanism among the three — reached `H` at a third weight. Dropping it buys allotype entropy
0.9629 → 0.9883 of maximum and Gini 0.1745 → 0.0889 for 4.466 → 4.368 expected responding units on a
20-unit cassette, across 19 donors. It is off by default; `--dominance` restores it, and that is
what reproduces the three-channel form.

A restriction cell naming a whole genotype (`HLA-A*01:01,HLA-A*03:01`) is resolved to the presented
*set* rather than compared as a string, so a unit whose restriction was never resolved no longer
reads as an allotype of its own — it used to be able to hold a slot beside both of its own
constituents, and to count as a seventh allotype on a donor carrying six.

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

### Reading one patient against a cohort

```bash
mhcmatch cassette report --cassettes design.tsv --pool pool.tsv \
    --reference cohort.tsv --reference-column lam --reference-outcome os_months \
    --out report.html
```

One self-contained HTML page: the units, what each contributes, and — with `--reference` — a risk
band placing this patient's `lam` as a percentile of that cohort, with the band's observed outcome
beside it. **The cohort is yours to pass; mhcmatch ships none.** A cohort is data, and one vendored
in a wheel would be stale by the next release.


## What `rank` costs

`rank --score aggregate` computes **every one** of the model's features before scoring: a model emits
the features it used and refuses to run without them — and nothing else, so a `C_corpus_*` column the
fitted model does not carry is absent from the header rather than present and empty. For
`mhc1.human.neoantigen` that is all three corpus channels — `thymus`, `self`, `viral` — and they are
affordable because the term contracts a k-mer table rather than searching an index, so
the host-proteome reference index (~7.5 GB, 6 min 15 s) is off the ranking path entirely and
`--no-self` is still allowed with `--score aggregate`. That index is what `--extended` and
`--annotate` cost, because they report *which* reference peptide was hit.

`--score gate` uses the two-term noisy-AND and stays cheaper still.

**The cost that remains is the per-allele `%rank` background, and it is cached.** Building one
allele's calibration background is a 10,000-peptide draw scored under that allele's model, ~0.95 s,
and it is a pure function of `(allele, model, background, footprint, seed, library version,
scoring-source digest)`. The
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

The key covers the library version, a digest of the scoring source itself (`predict`,
`calibrate`, `diffusion`, `affinity` -- so a change to the scoring code invalidates the cache
without anyone remembering to bump a counter), class, background, footprint, head, panel size,
draw count, seed and the positives feeding the isotonic fit. A cache keyed on less than that would be worse
than none, because it would serve a background drawn against a different model as though it were
this one. Set it to `off` and nothing is cached; leave it unset and entries go to `~/.cache/mhcmatch/calibration`. Either directory is disposable -- and use a fresh one for any run meant to establish a number rather than iterate.


## Batch and threads — read this before scripting a loop

**Pass `--peptides FILE` to any peptide-keyed command.** The expensive part of most of them is setup
that a per-peptide invocation pays again every time: the presentation and affinity calibrators are
~5 s, the binder calibrator ~45 s, the human-proteome text index 0.7 s. All of it is cached for the
life of the process, so one process over a whole list is the difference between 49 s per peptide and
thousands per second. Measured, both ways, in
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
one. It is **hierarchical**: for `mhc1.*.neoantigen`, nine columns in four blocks — presentation, expression, physchem,
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


## Presentation and affinity are not the same term

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
```

`--reference` is the one a cluster wants: it is everything `rank`, `neoag` and `mimicry` read, in one
call.

**There is no index to stage, because there is nothing left worth staging.** Peptide origin search
runs on one `seqtree.TextIndex` per proteome, and a proteome is the only input it needs: measured on
the human proteome (147,506 records, 69,578,135 residues), **0.7 s to build and 0.6 GB peak RSS, for
every peptide length and every substitution radius at once**. What that replaced was one index per
query *length* — ~65 s and 12.6 GB for the first, ~5.5 GB on disk apiece, and ~82 GB for the fifteen
lengths class II admits — which is why there used to be a `bootstrap --index` flag, a
`proteome_index/` cache directory, an `O_EXCL` build marker and a `$MHCMATCH_INDEX_WAIT` timeout.
All four are gone. **The race-free design is the one with nothing to race on**: at 0.7 s, rebuilding
costs less than agreeing about who rebuilds.

Pseudosequences (34-mer grooves) and the fitted model parameters are vendored in
`src/mhcmatch/data/` with their `PROVENANCE.md`. Nothing is refitted at import.


## Deployment

Three integrations, all under `integrations/`, all calling the same CLI, all **local-only** — no
scheduler profile and no sbatch template ships, and `-profile conda` / `-profile docker`
is the whole deployment story. What each was run against is recorded in
the module READMEs under `integrations/`.

| | for |
|---|---|
| [`nextflow/mhcmatch/`](integrations/nextflow/mhcmatch/) | an nf-core-style module, nine processes and a runnable entry point |
| [`snakemake/mhcmatch/`](integrations/snakemake/mhcmatch/) | the same two arms as a Snakemake 8 module, with a config schema that rejects a typo'd key before the DAG is built |
| [`nextflow/overlay/`](integrations/nextflow/overlay/) | attaching mhcmatch to a neoantigen pipeline you already run, in five modes, of which the default schedules nothing |

The two engines agree to the last digit on the same input — same offset, same yields — which is the
cross-check worth having when one library has two front ends.

`integrations/nextflow/mhcmatch/` is a self-contained nf-core-style module — **nine processes**
(`MHCMATCH_ALLELES`, `_PREDICT`, `_RANK`, `_RERANK`, `_NEOAG`, `_MIMICRY`, `_CASSETTE_SELECT`,
`_CASSETTE`, `_CASSETTE_SCORE`), **two arms** that chain them, and **`pipeline.nf`**, a runnable
entry point over a samplesheet:

```bash
nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
    --input samplesheet.csv --outdir results --mode both \
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

> ### The safety screen: ON through either engine, OFF in the bare CLI
>
> **Which layer you call decides it, so the default is worth knowing per layer.** Both engines ship
> it on — `params.mhcmatch_vector_screen = true`, `vector.screen: true` — so a pipeline run
> withdraws units on essential-tissue self-origin. `mhcmatch cassette build` / `order` invoked
> directly is the other way round: `--screen` is a flag you pass, and without it **no safety check
> runs at all** and the cassette carries whatever it was handed. Every `MHCMATCH_CASSETTE` task
> prints a line when no screen ran, so the absence is never silent.
>
> It was off in both engines while it was expensive, and what changed is the **cost, not the
> judgement**: the whole-proteome index was rebuilt inside every task, once per register length. One
> text index now answers every length and builds in **0.7 s**, so there is nothing to stage and no
> longer a reason to turn it off.
>
> The **mimicry annotation** (`--mhcmatch_mimicry`) is still off, and is a different case: it is
> annotation only, and **scores are identical either way**.

> **The pipeline pins its own mhcmatch version.** `environment.yml`, the `Dockerfile` and
> `params.mhcmatch_container` all name it, and a test checks each against `pyproject.toml`. The
> module calls the CLI by subcommand and flag name, so a checkout ahead of the installed release
> passes flags that release has never heard of and the failure is a bare argparse error deep in a
> task log.

**`pipeline.nf` is the easy entry point, not the integration surface.** A pipeline that already does
variant calling, HLA typing and expression quantification should `include` the processes in
`main.nf` or the arms in `subworkflows/` into its own channel topology.

The per-process input/output contract, the column rules (`--passthrough` refuses a collision,
`cassette select` renames and says so), every parameter and the Docker build are [`integrations/nextflow/mhcmatch/README.md`](integrations/nextflow/mhcmatch/README.md) and
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
