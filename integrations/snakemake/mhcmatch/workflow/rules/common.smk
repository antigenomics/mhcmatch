"""Samplesheet parsing, the arm predicates, and the option strings every rule shares."""
import csv
import os

OUT = config.get("outdir", "results")
MODE = config["mode"]
ARMS = {"rerank": ["rerank"], "denovo": ["denovo"], "both": ["rerank", "denovo"]}[MODE]
CLASSES = ["mhc1", "mhc2"]

#: **`config['input']` -- a samplesheet CSV -- is the entire input contract**, and it replaced a
#: filename convention for one reason: a convention is a claim about someone else's pipeline. The
#: columns are `sample,class,candidates,windows,hla`; one row per sample per class.
#:
#: What the columns are expected to hold, for the community-standard upstream:
#:
#: * `candidates` -- the rerank arm's input: any table with a peptide column and an allele column.
#:   pVACtools (nf-core/sarek -> VEP -> pVACseq) writes one per donor.
#: * `windows` -- the de novo arm's input and the rerank arm's `--context`: the peptide-window
#:   FASTA, which is what `pvacseq generate_protein_fasta` produces.
#: * `hla` -- a typing file for `mhcmatch alleles` (OptiType for class I, arcasHLA / HLA-LA for
#:   class II). Omit it and the `alleles` / `alleles_mhc2` config literal must be given instead.
_SHEET_COLUMNS = ("sample", "class", "candidates", "windows", "hla")


def _samplesheet():
    """`{sample: {kind: path}}` from `config['input']`, read ONCE at DAG-build time.

    Read once and stored, never re-read inside a rule and never discovered by globbing a directory.
    A rule whose input globs a directory is evaluated against whatever happens to exist when that
    rule is considered, which for the cohort step below would silently fit the offset over a subset
    of the donors -- the exact "every donor's mean equals the declared prevalence" defect the cohort
    step exists to avoid.

    The keys are `table_<cls>` (candidates), `fasta_<cls>` (windows) and `typing`.
    """
    path = config.get("input")
    if not path:
        return {}
    home = os.path.dirname(os.path.abspath(path))

    def _resolve(cell):
        """A cell's path, **resolved against the samplesheet's own directory** when relative.

        Not against the working directory. A samplesheet is written beside the files it names and
        then run from somewhere else -- `--directory`, or a `module` include -- and resolving
        against the cwd turns every relative cell into a file that does not exist. Snakemake then
        reports a missing *input*, which reads as a broken upstream rather than as a samplesheet
        read from the wrong place.
        """
        cell = (cell or "").strip()
        if not cell:
            return None
        return cell if os.path.isabs(cell) else os.path.normpath(os.path.join(home, cell))

    found, seen = {}, {}
    # `utf-8-sig`: a samplesheet that has been through Excel carries a BOM, and the first column
    # would otherwise read as `﻿sample` -- i.e. every row would have no sample id.
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        header = [(c or "").strip() for c in (reader.fieldnames or [])]
        for req in ("sample", "class"):
            if req not in header:
                raise WorkflowError(
                    f"{path}: no `{req}` column (header: {header}). The samplesheet columns are "
                    f"{','.join(_SHEET_COLUMNS)} -- see README.md")
        for lineno, row in enumerate(reader, start=2):
            sid = (row.get("sample") or "").strip()
            cls = (row.get("class") or "").strip()
            if not sid:
                raise WorkflowError(f"{path}:{lineno}: empty `sample`")
            if cls not in CLASSES:
                raise WorkflowError(f"{path}:{lineno}: sample {sid!r} has class {cls!r}, which is "
                                    f"not one of {CLASSES}")
            table, fasta = _resolve(row.get("candidates")), _resolve(row.get("windows"))
            if not table and not fasta:
                raise WorkflowError(
                    f"{path}:{lineno}: sample {sid!r} ({cls}) gives neither `candidates` nor "
                    "`windows`, so there is nothing for either arm to read -- `candidates` is the "
                    "rerank arm's input and `windows` is the de novo arm's")
            # **One row per (sample, class).** Taking the last row silently made the second
            # row's paths win: measured, a sheet with two `S1,mhc1` rows scored S2's
            # candidates and published them under S1's name, exit 0, no warning. Nextflow
            # refuses the same sheet by name, so this is also the two engines agreeing.
            if (sid, cls) in seen:
                raise WorkflowError(
                    f"{path}:{lineno}: sample {sid!r} already has a {cls!r} row at line "
                    f"{seen[(sid, cls)]}. One row per (sample, class) -- a second row's "
                    f"`candidates`/`windows` would silently replace the first's.")
            seen[(sid, cls)] = lineno
            rec = found.setdefault(sid, {})
            if table:
                rec[f"table_{cls}"] = table
            if fasta:
                rec[f"fasta_{cls}"] = fasta
            hla = _resolve(row.get("hla"))
            if hla:
                # One sample's two rows naming two DIFFERENT typing files is a samplesheet error.
                # Taking the last one silently would score one class against the other donor's
                # panel, which is the failure shape `mhcmatch alleles` exists to make loud.
                if rec.get("typing", hla) != hla:
                    raise WorkflowError(f"{path}:{lineno}: sample {sid!r} already carries typing "
                                        f"{rec['typing']!r} and this row says {hla!r}")
                rec["typing"] = hla
    # **No existence check here, on purpose.** Under `module` / `use rule * from mhcmatch` a listed
    # path may be an output the caller's own rules have yet to produce, and refusing it at parse
    # time would break that spelling. Snakemake names a genuinely missing file when it builds the
    # DAG.
    return found


SAMPLES_MAP = _samplesheet()
SAMPLES = sorted(SAMPLES_MAP)
# **An unset `input` used to reach `Nothing to be done` and exit 0**, which reads as a
# successful run. `pipeline.nf` refuses the same case by name; this is the two engines
# agreeing that a workflow with no declared input is a mistake and not an empty one.
if not config.get("input"):
    raise WorkflowError(
        "no `input`: point it at a samplesheet CSV, either in `config/config.yaml` or on the "
        "command line (`--config input=samplesheet.csv`). Columns: "
        + ",".join(_SHEET_COLUMNS) + " -- see README.md")
if not SAMPLES:
    raise WorkflowError(f"{config['input']}: no sample rows -- see README.md for the samplesheet")

# The cassette and the cohort offset are class I, and `cassette score` fits ONE offset over EVERY
# sample in the run, so an arm cannot be run with one sample missing from it. Said here, naming the
# sample and the column, rather than left to a `KeyError: 'table_mhc1'` inside an input lambda.
for _s in SAMPLES:
    for _arm in ARMS:
        _need, _col = (("table_mhc1", "candidates") if _arm == "rerank"
                       else ("fasta_mhc1", "windows"))
        if _need not in SAMPLES_MAP[_s]:
            raise WorkflowError(
                f"sample {_s!r} has no `class: mhc1` row with a `{_col}` path, which the {_arm} arm "
                f"needs: the cassette is class I and its cohort offset is fitted over every sample "
                f"at once. Add the row, drop the sample, or choose another `mode`")


def pool_for(arm, sample, cls="mhc1"):
    """The scored table an arm's cassette is chosen from."""
    return (f"{OUT}/rerank/{sample}.{cls}.epitopes.mhcmatch.tsv" if arm == "rerank"
            else f"{OUT}/denovo/{sample}.{cls}.mhcmatch.ranked.tsv")


def score_column(arm):
    """Which column of the pool holds the aggregate.

    **Do not leave this to the fallback on the rerank arm.** A caller's candidate table HAS a
    `score` column -- theirs -- so an unqualified fallback selects on the upstream tool's ranking
    while looking as though it selected on ours.
    """
    if config.get("cassette", {}).get("score_column"):
        return config["cassette"]["score_column"]
    return f"{config['rerank_prefix']}score" if arm == "rerank" else ""


def opt(flag, value):
    """`--flag value` when the value is set, and nothing at all when it is not."""
    return f"{flag} {value}" if value not in (None, "", False) else ""


def flag(name, value):
    """A bare `--flag` when true. Accepts the strings a CLI `--config` produces."""
    return name if str(value).lower() not in ("false", "0", "no", "none", "") else ""


def collect_targets():
    """What `rule all` asks for: one cassette per sample per arm, plus one cohort file per arm."""
    if not SAMPLES:
        return []
    out = []
    for arm in ARMS:
        out += [f"{OUT}/{arm}/{s}.cassette.faa" for s in SAMPLES]
        out += [f"{OUT}/{arm}/cohort.cassette_score.tsv",
                f"{OUT}/{arm}/cohort.cassette_report.html"]
        # The scored table for EVERY row of the samplesheet, both classes. The cassette is class I,
        # so without these a `class: mhc2` row would drive no job at all and the class-II path
        # would be unreachable in a dry run.
        kind = "table" if arm == "rerank" else "fasta"
        for s in SAMPLES:
            for c in CLASSES:
                if f"{kind}_{c}" in SAMPLES_MAP[s]:
                    out.append(pool_for(arm, s, c))
                    if arm == "denovo":
                        out.append(f"{OUT}/denovo/{s}.{c}.mhcmatch.native.tsv")
    return out
