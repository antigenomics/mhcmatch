<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/mhcmatch/master/assets/mhcmatch_dark.png">
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
| Build a vaccine cassette from ranked candidates | `mhcmatch vector --candidates units.tsv --n0 8 --screen` | `vector.select` / `vector.order` |
| …and a map of it a viewer can draw | `mhcmatch vector ... --map c.tsv --map-json c.json` | `vector.epitope_map` |
| Strip frameshift-prone motifs from the CDS | `mhcmatch deslip cassette.fa` | `vector.slippery_sites` |
| Split a peptide into anchor / TCR-facing parts | `mhcmatch decompose PEP` | `store.decompose` |
| How viral-like is it, as a soft sum not a cutoff? | — | `luksza.viral_r` |

Full command reference, grouped by task: [the CLI page](https://antigenomics.github.io/mhcmatch/cli.html).

`predict` is the presentation axis (**is it presented at all**, the NetMHCpan `%Rank_EL` analogue);
`restriction` is the specificity axis (**which allele**). They answer different questions and a
peptide can top one and not the other — `NLVPMVATV` is unambiguously A\*02:01-restricted yet bands
mid-pack against A\*02:01's own ligands.

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
the two questions apart and combines them with a **gate** (a product of sigmoids), not a sum, so a
candidate that fails either one cannot be rescued by the other.

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
whole published corpus scores in seconds, so pass a list. `mhcmatch.posbayes` and `mhcmatch.ipred`
are strict special cases of it and ship alongside for comparison.

**Evidence that outranks a model.** `mhcmatch.known` carries five reference sets built from the
public deposits — confirmed tumour neoantigens, peptides the screens tested and found
non-immunogenic, IEDB-immunogenic epitopes, the thymic self-immunopeptidome, the viral ligandome. An
exact match is stronger evidence than any score, so `rank` reports it as a flag and floats those
candidates into a tier of their own instead of folding it into the number.

**Pick your tumour type.** `mhcmatch expression --list-contexts` prints all 19 TCGA↔GTEx pairings;
`expression.matched_tissues('BRCA')` gives the matched normal and `expression.lookup(gene,
tumor='BRCA')` the tumour value. The benchmark's own expression term is a GTEx **cross-tissue
median** so that it is identical on fit and holdout — that is a comparability requirement, not a
recommendation, and a real ranking run should pass its own tumour type.

**Expression, and which normal tissue to compare against.** `--tumor` takes a **TCGA study
abbreviation** (`SKCM`, `LUAD`, …; `CRC` merges TCGA's `COAD` and `READ`) and `--tissue` a **GTEx
`SMTSD`** name (`Skin - Sun Exposed (Lower leg)`). Neither is a clinical coding system — not
ICD-O-3, SNOMED CT or OncoTree — so a pipeline needing one brings its own crosswalk.
`expression.matched_tissues("SKCM")` gives a tumour type's matched normal, which is what makes the
safety read askable without knowing the pairing by heart; `mhcmatch expression --list-contexts`
prints all 19 pairings and the 31 tissues that are safety-read-only. `HNSC` is marked approximate —
GTEx has no head-and-neck mucosa.

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

Both ends join to the rest of the library. `--context windows.fasta` takes `rank`'s **minimal
epitopes** and rebuilds them as long units against the FASTA they were called on, one per variant
rather than one per register — a minimal peptide loads onto any cell without costimulation and is
the tolerising configuration, so the reader will not take one. `--fasta-nt` writes the coding
sequence: highest-usage human codon per residue, backed off to shorten homopolymers, then deslipped.
It is **not a codon optimiser** — it fixes the two things that break a concatemer specifically and
leaves GC, structure and CpG to the manufacturer's tooling.

```zsh
mhcmatch rank fasta windows.fasta --alleles "$HLA" --out ranked.tsv
mhcmatch vector --candidates ranked.tsv --context windows.fasta --n0 8 --screen \
    --fasta cassette.faa --fasta-nt cassette.fna
```

⚠️ **`mimicry` is a scoring term, not a safety screen.** Flagging candidates by "resembles a
tolerance-side reference" fires on almost everything — influenza `GILGFVFTL` drew 14
essential-tissue hits — because anchor-masked similarity to a *presented* reference is presentation,
not recognition. Exclusion goes through `vector.self_origin_risk`
(`bench/results/vector_safety_screen.md`).

## Model names

Every fitted model is named by the **acronym of its parameters** — `aggregate5` and "the full model"
said nothing about what was in them, and two designs were once both "the neoantigen model". One
letter per parameter, in a fixed canonical order:

| letter | parameter | from | what it is |
|---|---|---|---|
| `P` | presentation | `AnchorModel` | `-log10` of the per-allele `%rank`; fitted on **observed ligands** |
| `B` | binder score | `predict.binder_score` | `-log10` of the calibrated combined `%rank` (Fisher of `P` and `A`) |
| `A` | affinity | `PottsAffinity` | `-log10` of the Potts IC50 `%rank`; fitted on **measured IC50** |
| `D` | differential agretopicity | `PottsAffinity.dai` | `log10(Kd_WT / Kd_MT)` vs the recovered wild type |
| `E` | expression | `mhcmatch.expression` | `log1p(TPM)`, observed or reference-imputed |
| `V` | vanilla physicochemistry | `mhcmatch.ipred` | the 13-parameter calibrated log-odds |
| `C` | complementarity | `mhcmatch.complement` | the six-block recognition log-odds |
| `R` | Łuksza recognition | `mhcmatch.luksza` | `Z/(1+Z)`, a soft sum over near-matches rather than a distance cut |
| `F` | foreignness | viral IEDB ligandome | distance to the nearest viral epitope |
| `M` | mimicry | `mhcmatch.mimicry` | the six-channel signed aggregate |

So `PADEC` is presentation + affinity + agretopicity + expression + complementarity, and `PADECM`
adds mimicry. Suffixes are fitting choices rather than parameters: `-bal` (every screen weighted to
the same total mass), `-scr` (screen indicators as nuisance columns).

**`V` is "vanilla", not "ipred".** `ipred` is the *old* recognition term and `complement` is what
replaced it — the same axis at two generations, with `ipred` a strict special case of `complement`.
Naming the letter after the generation rather than the module makes `BDEVF` legible as "the old
model" at a glance.

**`P` is not a second affinity term.** Both end up as a `%rank` against the same kind of background,
so the mechanism doesn't separate them — the training data and the target do. `A` is a Potts model
fitted on **measured IEDB IC50**, targeting `Kd`: the biophysics of the groove. `P` is the
AnchorModel fitted on the **observed ligand panel**, targeting how *ligand-like* a peptide is, which
carries processing, transport and abundance signal that binding alone does not. That is the field's
binding-affinity vs eluted-ligand split, and the two are measurably not redundant — on TESLA-608
`A` scores 0.757 AUROC, `P` 0.763, and their Fisher combination `B` **0.786**. A combination cannot
beat both parents by that margin on the same measurement twice.

**`P` is a rank, not a similarity search** — worth stating because the name invites the other
reading. It is a score against a random-peptide background (10,000 peptides matched to the corpus's
amino-acid and length distribution). Nothing is retrieved: no reference peptide is looked up and no
anchor-matched protein is searched for. `A` and `B` are the same, so the whole presentation side is
**scoring, not retrieval**. The searches are `restriction` (the epitope panel, anchor-masked — not in
any acronym), `M` (thymic/viral/proteome windows) and `F` (the viral ligandome).

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
aff.predict_ic50("NLVPMVATV", "HLA-A*02:01")             # ~64 nM
aff.amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")   # Kd_WT/Kd_MT (Łuksza eq. 9)

complement.score(peptides)                           # vectorised: pass the list, not a loop
complement.posterior(peptides, prior=4.2e-4)         # the log-odds carries NO prior; supply yours
complement.score(peptides, species="mouse")          # separate table; the hosts are never pooled

known.lookup("GILGFVFTL")                            # -> 'viral'
mimics.neighbours(peptides, ref_sets, threads=0)     # threaded C++ neighbour search

pm = mhcmatch.Proteome.from_hf("human")
pm.find_sources(peptides, max_subs=1, threads=0)     # batch; find_source() is the single-query form
pm.wildtype("NLVPMVATV")                             # the WT counterpart, for agretopicity
```

Full API: [antigenomics.github.io/mhcmatch](https://antigenomics.github.io/mhcmatch/). Eight
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

`integrations/nextflow/mhcmatch/` is a self-contained nf-core-style module (`main.nf`,
`nextflow.config`, `environment.yml`, `Dockerfile`) that drops in for MHCflurry class I and the
class-II binding subworkflow, consuming the same `(meta, peptide.fasta, alleles)` channel and
emitting a pipeline-compatible `.scored.csv`. The image bootstraps its panel at **build** time, so
compute nodes need no network.

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
