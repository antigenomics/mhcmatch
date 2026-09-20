"""The `integrations/` tree, checked against the CLI it actually calls.

Both engines invoke `mhcmatch` **by flag name**, from a shell string, so a flag that was renamed --
or one that never existed -- is not a type error, an import error or anything a unit test sees. It
is an `exit 2` inside a task log, hours into a fan-out, or worse: a flag whose value happens to be
unset, so the wrapper emits nothing and the mistake stays latent until somebody sets it.

That is not hypothetical. `integrations/snakemake/.../denovo.smk` passed `--threshold` to
`mhcmatch rank`, which has never accepted it, from the day the rule was written until 1.21.0. It
survived because the default is `none` and the `opt()` helper emits nothing for an unset value --
and it would have exited 2 on every task in the arm the first time anyone set the documented
`rank_threshold` key. `--threshold` is a real flag, just not on that command: `cassette
build`/`order` take one. So a check that asks "is this a flag mhcmatch knows anywhere" passes it.
The check below asks argparse about the SUBCOMMAND the line actually runs.
"""
import argparse
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
INTEGRATIONS = ROOT / "integrations"

# sdist/wheel checkouts do not carry integrations/
pytestmark = pytest.mark.skipif(not INTEGRATIONS.is_dir(), reason="no integrations/ in this tree")


# --- asking argparse what a subcommand accepts --------------------------------------------------

class _Captured(Exception):
    def __init__(self, parser):
        self.parser = parser


def _root_parser():
    """The `mhcmatch` parser, captured on its way into `parse_args`.

    `cli.main` builds the parser inline and parses in the same call, so there is nothing to import.
    Intercepting `parse_args` is the direct route and keeps this test reading the ONE definition
    that ships -- a second copy of the flag list here would be a second thing to keep current, and
    would go stale in exactly the way it exists to catch.
    """
    original = argparse.ArgumentParser.parse_args

    def spy(self, *a, **k):
        raise _Captured(self)

    argparse.ArgumentParser.parse_args = spy
    try:
        from mhcmatch import cli
        cli.main([])
    except _Captured as got:
        return got.parser
    finally:
        argparse.ArgumentParser.parse_args = original
    raise AssertionError("mhcmatch.cli.main no longer calls parse_args; this capture is broken")


def _subparsers(parser):
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices
    return {}


def _accepted_flags(parser):
    return {s for action in parser._actions for s in action.option_strings if s.startswith("--")}


def _resolve(words):
    """`['cassette', 'order']` -> the flags that subcommand accepts.

    A non-literal second word -- `mhcmatch cassette ${verb}` in main.nf, where `verb` is `build` or
    `order` -- resolves to the UNION of the parent's subcommands. That loses nothing here, because
    `build` and `order` are the same parser with one sizing rule skipped and declare identical
    flags; it would lose something the day two siblings diverge, which is when this should be made
    to read the Groovy.
    """
    parser = _root_parser()
    for i, word in enumerate(words):
        children = _subparsers(parser)
        if not children:
            break
        if word in children:
            parser = children[word]
            continue
        if i == 0:
            return None                       # not an `mhcmatch` subcommand at all
        return set().union(*(_accepted_flags(c) for c in children.values()))
    return _accepted_flags(parser)


# --- reading the integration files --------------------------------------------------------------

#: A `--flag` as the CLI spells one. **`_` is the discriminator and it is not a coincidence**: every
#: mhcmatch CLI flag is kebab-case and every Nextflow / Snakemake *parameter* is snake_case
#: (`--mhcmatch_cassette_screen`, `--publish_dir_mode`), so the two vocabularies cannot collide.
_FLAG = re.compile(r"--[a-z][a-z0-9]*(?:-[a-z0-9]+)*\b")
_INVOKE = re.compile(r"\bmhcmatch\s+([a-z][a-z0-9-]*)(?:\s+(\S+))?")

#: Flags belonging to the ENGINE, not to mhcmatch. They reach a `WorkflowError` message telling an
#: operator how to re-run (``--config input=samplesheet.csv``), which is neither a comment nor a
#: docstring, so `_strip_prose` cannot see it and nothing else can tell it apart from a helper
#: building a real argument. A short explicit list is the honest answer; each name here is a
#: snakemake or nextflow CLI flag and none of them is a mhcmatch one.
_ENGINE_FLAGS = {"--config", "--configfile", "--directory", "--cores", "--sdm", "--profile",
                 "--use-conda", "--dry-run", "--snakefile", "--input", "--outdir", "--mode",
                 "--stub-run", "--resume"}


def _strip_prose(path):
    """The file's text with the parts that only TALK about flags blanked out.

    Comments in both syntaxes, and -- for the Python-shaped files only -- docstrings. A `.nf`
    file's triple-quoted block IS the shell command, so the same rule applied there would blank
    every command in the module.

    Blanked rather than deleted, so `^rule` / `^process` still land at the same column and the
    splitter below does not have to care.
    """
    text = path.read_text()
    if path.suffix == ".smk" or path.name == "Snakefile":
        quotes = ('"' * 3, "'" * 3)
        out, i = [], 0
        while True:
            hits = [j for j in (text.find(q, i) for q in quotes) if j != -1]
            if not hits:
                out.append(text[i:])
                break
            start = min(hits)
            quote = text[start:start + 3]
            end = text.find(quote, start + 3)
            if end == -1:
                out.append(text[i:])
                break
            out.append(text[i:start])
            out.append("\n" * text.count("\n", start, end + 3))
            i = end + 3
        text = "".join(out)
    return "\n".join("" if line.strip().startswith(("//", "#")) else line
                     for line in text.splitlines())


def _blocks(path):
    """`{label: text}` for each process / rule in a file, plus `""` for everything shared.

    Per block, not per file. A whole-file union would have passed the `--threshold` defect above,
    because the same file's cassette rule legitimately takes one.

    A flag named in prose is not a flag passed, so comments and (in the Python-shaped files)
    docstrings are blanked first -- see `_strip_prose`. Without that, `cassette select`'s own
    docstring saying it accepts neither `--species` nor `--tier` reads as it passing both.

    Anything at column 0 that is not a `process`/`rule` header -- a helper `def`, a constant --
    closes the current block and goes to the shared bucket, which is checked against the union of
    every subcommand the file runs. `keepArgs()` in main.nf is exactly that: one helper whose flags
    belong to two of the nine processes.
    """
    blocks, label, buf = {}, "", []
    for line in _strip_prose(path).splitlines():
        header = re.match(r"(?:process\s+(\w+)\s*\{|rule\s+(\w+)\s*:)", line)
        if header:
            blocks[label] = blocks.get(label, "") + "\n".join(buf)
            label, buf = header.group(1) or header.group(2), []
            continue
        if line[:1] not in ("", " ", "\t", "}") and label:
            blocks[label] = blocks.get(label, "") + "\n".join(buf)
            label, buf = "", []
        buf.append(line)
    blocks[label] = blocks.get(label, "") + "\n".join(buf)
    return blocks


def _integration_files():
    engines = (sorted(INTEGRATIONS.rglob("*.nf")) + sorted(INTEGRATIONS.rglob("*.smk"))
               + sorted(INTEGRATIONS.rglob("Snakefile")))
    return [p for p in engines if p.is_file()]


def test_every_flag_the_integrations_pass_is_one_that_subcommand_accepts():
    """The check that would have caught `mhcmatch rank --threshold` on the day it was written."""
    files = _integration_files()
    assert files, "found no .nf/.smk/Snakefile under integrations/ -- this guard has gone vacuous"

    checked, bad = 0, []
    for path in files:
        blocks = _blocks(path)
        per_block = {}
        for label, text in blocks.items():
            commands = {(m.group(1), m.group(2)) for m in _INVOKE.finditer(text)}
            accepted, names = set(), []
            for first, second in commands:
                got = _resolve([w for w in (first, second) if w])
                if got is not None:
                    accepted |= got
                    names.append(" ".join(w for w in (first, second) if w))
            per_block[label] = (accepted, names, text)

        # The shared bucket's flags belong to whichever process defined the helper, so they are
        # checked against everything this file runs.
        shared = set().union(*(a for a, _, _ in per_block.values())) or set()
        for label, (accepted, names, text) in per_block.items():
            allowed = shared if label == "" else accepted
            if not allowed:
                continue                      # a block that runs no mhcmatch command
            for flag in sorted(set(_FLAG.findall(text)) - _ENGINE_FLAGS):
                checked += 1
                if flag not in allowed:
                    bad.append(f"{path.relative_to(ROOT)}:{label or '<shared>'}: {flag} "
                               f"is not accepted by `mhcmatch {' / '.join(sorted(names))}`")
    assert checked > 100, f"only {checked} flag(s) reached the check -- the extractor stopped working"
    assert not bad, "flags the CLI does not accept:\n  " + "\n  ".join(bad)


# --- the version pins ---------------------------------------------------------------------------

def _declared_version():
    return re.search(r'^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M).group(1)


def test_nextflow_pins_match_pyproject():
    """The container pins must name the version this checkout builds -- unless it is a dev version.

    A ``.devN`` suffix means no wheel has been published and no image pushed, so there is nothing
    for `mhcmatch==<version>` to resolve to. On a dev version the pins are allowed to lag; they are
    checked again at release, when the suffix is dropped.

    **Match the pin, not any version-shaped string.** This scan used to be
    ``re.findall(r"\\b0\\.\\d+\\.\\d+\\b", ...)``, which went VACUOUS the day 1.0.0 shipped. It also
    never opened the two files whose pins actually drifted: ``nextflow.config`` is not ``*.nf``, so
    ``params.mhcmatch_container`` went unchecked. Anchoring on the spellings of a *mhcmatch* pin
    also keeps Nextflow's own ``22.10.0`` from reading as a stale one.

    A ``README.md`` is in the glob because it is the file a collaborator actually follows, and it
    was the last one left out: at 1.9.0 it still said ``pip install "mhcmatch==1.8.0"`` in twelve
    places while every machine-read pin beside it had moved.

    Rooted at ``integrations/`` and globbed by KIND, so a new subdirectory is covered the day it
    lands rather than the release after someone remembers -- the overlay sat on 1.17.0 under a scan
    rooted one directory too deep, and ``preflight.nf`` compares its pin to ``mhcmatch --version``
    under ``errorStrategy = 'terminate'``, so every overlay mode except ``off`` died at its FIRST
    process.
    """
    want = _declared_version()
    if ".dev" in want:
        pytest.skip("dev version: no wheel published for the pins to name")

    pins = re.compile(r"(?:mhcmatch==|mhcmatch:|mhcmatch-|MHCMATCH_VERSION=|VERSION=|"
                      r"--branch v|tag=\"v|/v|pins |require_version[^\n']*')"
                      r"(\d+\.\d+\.\d+)", re.M)
    scanned, stale = [], {}
    for path in sorted(set(INTEGRATIONS.rglob("*.nf")) | set(INTEGRATIONS.rglob("*.config"))
                       | set(INTEGRATIONS.rglob("*.smk")) | set(INTEGRATIONS.rglob("Snakefile"))
                       | set(INTEGRATIONS.rglob("*.sbatch")) | set(INTEGRATIONS.rglob("README.md"))
                       | set(INTEGRATIONS.rglob("Dockerfile"))
                       | set(INTEGRATIONS.rglob("*.yaml")) | set(INTEGRATIONS.rglob("*.yml"))):
        if not path.is_file():
            continue
        for found in set(pins.findall(path.read_text())):
            scanned.append(str(path.relative_to(ROOT)))
            if found != want:
                stale.setdefault(str(path.relative_to(ROOT)), set()).add(found)
    assert not stale, f"version pins behind pyproject {want}: {stale}"
    # A guard that matches nothing is the failure mode this test has already had. Fail loudly.
    assert scanned, ("found no mhcmatch version pin at all under integrations/ -- the pin spelling "
                     "changed and this guard has gone vacuous again")


# --- the one rule whose shape matters ------------------------------------------------------------

def test_snakemake_cohort_rule_expands_over_every_sample():
    """**The property the whole cohort design rests on**, and it is checkable by reading the rule.

    `cassette score` fits ONE offset over every donor in the run. Its input must be `expand()` over
    the module-level `SAMPLES`, resolved once from the samplesheet -- never a glob over the OUTPUT
    directory, which is evaluated against whatever exists when the DAG is built and would silently
    fit the offset over a subset. That is the 'every donor's mean equals the declared prevalence'
    defect the cohort step exists to prevent, and it would reappear looking fixed.
    """
    smk = (INTEGRATIONS / "snakemake" / "mhcmatch" / "Snakefile").read_text()
    assert "def _all_units(w):" in smk and "def _all_pools(w):" in smk
    units = smk[smk.index("def _all_units(w):"):smk.index("def _all_pools(w):")]
    assert "expand(" in units and "sample=SAMPLES" in units, units
    assert "glob(" not in smk, "the cohort rule must not glob an output directory"
    # ...and there is exactly one cohort output per arm, not one per sample
    assert 'f"{OUT}/{{arm}}/cohort.cassette_score.tsv"' in smk
    assert "{sample}" not in smk[smk.index("rule mhcmatch_cassette_score:"):]


def test_both_engines_refuse_a_cassette_unit_they_cannot_source():
    """A unit is the LONG (~27 aa) window, and neither engine may fall back to the minimal epitope.

    `mhcmatch cassette`'s own `--unit-column` fallback is `peptide`, which on a scored table is the
    minimal epitope -- and a 9-mer loads onto any cell without costimulation, which is the
    tolerising configuration (PMID 17911588). It is not a smaller version of the right thing, so
    with neither a `--context` FASTA nor a named column both engines STOP and name the sample.

    Asserted on both, because this is the one refusal that is silent when it goes missing: the run
    completes, publishes a construct, and exits 0.
    """
    nf = (INTEGRATIONS / "nextflow" / "mhcmatch" / "main.nf").read_text()
    smk = (INTEGRATIONS / "snakemake" / "mhcmatch" / "Snakefile").read_text()
    assert "if (!ctx && !ucol)" in nf and "no source for the cassette unit" in nf
    assert "raise WorkflowError(" in smk and "no `windows` FASTA to pass as `--context`" in smk
    for text, engine in ((nf, "main.nf"), (smk, "Snakefile")):
        assert "17911588" in text, f"{engine} no longer records WHY the fallback is refused"
