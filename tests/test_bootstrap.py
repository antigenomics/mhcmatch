"""pmhc HF bootstrap: ``Store.from_pmhc(None)`` with no ``$MHCMATCH_PMHC`` fetches
``pmhc/pmhc_<tier>.tsv.gz`` from the public HF dataset. The routing is tested offline (monkeypatched
``fetch_pmhc`` returning a tiny synthetic table); the real network fetch runs opt-in (``RUN_HF_FETCH=1``).
"""
import gzip
import os

import pytest

from mhcmatch import Store, mimicry
from mhcmatch import store as store_mod

_HEADER = "mhc_class\tmhc_species\tepitope\tmhc_a\tmhc_b\tweight\n"
_ROWS = ["MHCI\tHomoSapiens\tNLVPMVATV\tHLA-A*02:01\t\t1\n",
         "MHCI\tHomoSapiens\tGILGFVFTL\tHLA-A*02:01\t\t1\n"]


def test_from_pmhc_routes_to_fetch_when_no_env(monkeypatch, tmp_path):
    """No path + no MHCMATCH_PMHC -> from_pmhc must call fetch_pmhc(tier) and load its result."""
    monkeypatch.delenv("MHCMATCH_PMHC", raising=False)
    tbl = tmp_path / "pmhc_shortlist.tsv.gz"
    with gzip.open(tbl, "wt") as fh:
        fh.write(_HEADER)
        fh.writelines(_ROWS)
    seen = {}

    def fake_fetch(tier="full"):
        seen["tier"] = tier
        return str(tbl)

    monkeypatch.setattr(store_mod, "fetch_pmhc", fake_fetch)
    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    assert seen["tier"] == "shortlist"                       # routed to the HF bootstrap
    assert "HLA-A*02:01" in st.alleles("mhc1")


@pytest.mark.skipif(not os.getenv("RUN_HF_FETCH"), reason="set RUN_HF_FETCH=1 for the real HF download")
def test_fetch_pmhc_real_download():
    path = store_mod.fetch_pmhc("shortlist")
    assert path.endswith("pmhc/pmhc_shortlist.tsv.gz") and os.path.exists(path)


def test_fetch_proteome_resolves_names(monkeypatch):
    """Name -> proteome/<file> resolution (no download): human/mouse alias, pathogen stem passthrough."""
    import huggingface_hub
    monkeypatch.setattr(huggingface_hub, "hf_hub_download",
                        lambda repo_id, repo_type, filename: "/tmp/" + filename)
    assert store_mod.fetch_proteome("human").endswith("proteome/human.fasta.gz")
    assert store_mod.fetch_proteome("mouse").endswith("proteome/mouse.fasta.gz")
    assert store_mod.fetch_proteome("ecoli_K12_UP000000625").endswith(
        "proteome/ecoli_K12_UP000000625.fasta.gz")


def test_proteome_from_hf_routes_to_fetch(monkeypatch, tmp_path):
    """Proteome.from_hf(name) fetches then loads the FASTA."""
    import gzip
    from mhcmatch import Proteome
    from mhcmatch import store as sm
    fa = tmp_path / "human.fasta.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(">P1 test\nNLVPMVATVKQ\n")
    monkeypatch.setattr(sm, "fetch_proteome", lambda name="human": str(fa))
    pm = Proteome.from_hf("human")
    assert "P1" in pm.seqs


# --- release consistency -------------------------------------------------------------------
# The __init__ fallback and the nextflow container pins are hand-maintained copies of
# pyproject's version. Both drifted behind a release twice (0.14 -> 0.15, 0.15 -> 0.16), and a
# stale one silently mislabels every `versions.yml` a pipeline run emits.

def _declared_version():
    import re, pathlib
    root = pathlib.Path(__file__).resolve().parents[1]
    return re.search(r'^version = "([^"]+)"', (root / "pyproject.toml").read_text(),
                     re.M).group(1), root


def test_fallback_version_matches_pyproject():
    import re
    want, root = _declared_version()
    src = (root / "src" / "mhcmatch" / "__init__.py").read_text()
    got = re.search(r'__version__ = "([^"]+)"', src).group(1)
    assert got == want, f"__init__.py fallback {got} != pyproject {want}"


def test_installed_metadata_matches_pyproject():
    # mhcmatch.__version__ reads the *install* metadata, not pyproject -- so a stale editable install
    # makes every version-keyed test compare a stale value against itself and pass. That is how a
    # 0.25.0 -> 0.26.0 bump shipped vendored models stamped 0.25.0 with a green local suite and two
    # red CI runs: CI installs fresh, so only CI saw the mismatch. Fail here instead, where the
    # message says what to do.
    from importlib.metadata import PackageNotFoundError, version
    want, _ = _declared_version()
    try:
        got = version("mhcmatch")
    except PackageNotFoundError:              # source tree with no install -- the fallback governs
        return
    assert got == want, (
        f"installed mhcmatch metadata is {got} but pyproject declares {want}; this environment is "
        f"stale and every version-keyed test is comparing it against itself. Run: pip install -e .")


def test_no_shipped_artifact_is_stale():
    # `mhcmatch build --check` as a unit test, and the reason it is *here* rather than beside the
    # artifacts it checks: the two existing guards both need the HF deposit, carry @pytest.mark.hfdata,
    # and conftest skips them whenever it is not already staged -- which is always, in CI. This one
    # reads only shipped files, so it is the guard CI actually runs.
    #
    # Both older guards were correct and both were defeated at 0.26.0 by the same thing: a stale
    # editable install made mhcmatch.__version__ report 0.25.0, so each compared a stale artifact
    # against a stale expectation and passed. test_installed_metadata_matches_pyproject closes that;
    # this closes the coverage half.
    from mhcmatch import _build
    stale = _build.check()
    assert not stale, (
        "shipped artifacts are behind __version__: "
        + "; ".join(f"{t}/{f} is {got}, want {want}" for t, f, got, want in stale)
        + ". Run: mhcmatch build")


def test_every_build_target_owns_files_that_exist():
    # A target whose file list drifts from what is actually shipped makes --check silently vacuous.
    import os
    from mhcmatch import _build
    for name, (_label, _fn, files) in _build.TARGETS.items():
        assert files, name
        for f in files:
            assert os.path.exists(os.path.join(_build.DATA, f)), f"{name}: {f} is not shipped"


def test_nextflow_pins_match_pyproject():
    """The container pins must name the version this checkout builds -- unless it is a dev version.

    A ``.devN`` suffix means no wheel has been published and no image has been pushed, so there is
    nothing for `mhcmatch==<version>` to resolve to; pinning it would produce a module that cannot
    build. On a dev version the pins are therefore allowed to lag by exactly one patch-level bump,
    which is the release they were last valid for. They are checked again at release, when the
    suffix is dropped.
    """
    import re
    want, root = _declared_version()
    if ".dev" in want:
        return
    nf = root / "integrations" / "nextflow" / "mhcmatch"
    if not nf.is_dir():                       # sdist/wheel checkouts do not carry integrations/
        return
    # **Match the pin, not any version-shaped string.** This scan used to be
    # ``re.findall(r"\b0\.\d+\.\d+\b", ...)``, which made it VACUOUS the day 1.0.0 shipped: every
    # pin it guards has been ``1.x.y`` since, so it found nothing and passed. It also never opened
    # the two files whose pins actually drifted -- ``nextflow.config`` is not ``*.nf``, so
    # ``params.mhcmatch_container`` went unchecked (it sat on 1.6.0 while the rest were on 1.6.1),
    # and ``templates/*.sbatch`` was never in the list at all, though ``setup.sbatch`` asserts the
    # installed version equals its own ``VERSION=`` and so installs the wrong release when stale.
    # Anchoring on the four spellings of a *mhcmatch* pin also keeps Nextflow's own ``21.10.6`` in
    # main.nf from reading as a stale pin, which a bare ``\d+\.\d+\.\d+`` would.
    PINS = re.compile(r"(?:mhcmatch==|mhcmatch:|MHCMATCH_VERSION=|^VERSION=)(\d+\.\d+\.\d+)", re.M)
    scanned, stale = [], {}
    for p in sorted(list(nf.rglob("*.nf")) + list(nf.rglob("*.config"))
                    + list(nf.glob("templates/*.sbatch"))
                    + [nf / "Dockerfile", nf / "environment.yml"]):
        if not p.is_file():
            continue
        for v in set(PINS.findall(p.read_text())):
            scanned.append(str(p.relative_to(root)))
            if v != want:
                stale.setdefault(str(p.relative_to(root)), set()).add(v)
    assert not stale, f"version pins behind pyproject {want}: {stale}"
    # A guard that matches nothing is the failure mode this test just had. Fail loudly instead.
    assert scanned, ("found no mhcmatch version pin at all under integrations/nextflow/mhcmatch -- "
                     "the pin spelling changed and this guard has gone vacuous again")


# --- what `bootstrap` stages must remain ingestible by ------------------------------------------
# `mhcmatch bootstrap --reference` stages a fixed list of files from `isalgo/pmhc_data`, and each
# one has exactly one module function that reads it. Staging and reading are tested together
# because a schema drift on the HF side passes the download and fails somewhere downstream as an
# empty channel or an imputed expression -- a wrong number, not an error.

@pytest.mark.hfdata
def test_every_bootstrapped_reference_is_ingestible_by_its_consumer():
    """Each file in ``cli.REFERENCE_FILES`` stages, and the function that consumes it reads rows.

    The pairing is the contract. A renamed column or a moved path on the HF side would otherwise
    surface as ``C_corpus_viral`` quietly going to zero, or every candidate taking imputed
    expression under the largest coefficient in the fitted model.
    """
    from mhcmatch import cli, expression, known, mimics
    from mhcmatch.store import fetch_file

    for rel in cli.REFERENCE_FILES:
        assert os.path.exists(fetch_file(rel)), rel

    # the two mimicry deposits -> the peptide loader the corpus and index paths both use
    for rel in ("thymus/thymus_immunopeptidome.tsv.gz", "ligandome/viral_foreign_iedb.tsv.gz"):
        peps = mimics.load_peptides(None, rel, "mhc1")
        assert len(peps) > 1000 and all(p.isalpha() for p in peps[:50]), rel

    # the thymic source-protein column, which `safety_profile` joins on and which is the one field
    # `_sources` reads beyond the peptide itself
    src = mimicry._sources(None, "thymus/thymus_immunopeptidome.tsv.gz")
    assert len(src) > 10_000 and any(v for v in src.values())

    # the expression reference -> the table `rank` reads `expr` off. `load()` returns
    # {(key_type, key, context): {stat: value}}, so the schema check is on the value dict.
    tbl = expression.load()
    assert len(tbl) > 1000
    stats = {"median_tpm", "q25_tpm", "q75_tpm"} & set(expression.COLUMNS)
    (kt, key, ctx), vals = next(iter(tbl.items()))
    assert kt and key and ctx and stats <= set(vals), (kt, key, ctx, sorted(vals))

    # the known-epitope sets -> every declared (file, label column, hit values) triple resolves
    sets = known.load()
    assert set(sets) == set(known.SET_NAMES)
    for name, peps in sets.items():
        assert peps, name


@pytest.mark.hfdata
def test_the_bootstrapped_proteomes_reach_the_functions_that_window_them():
    """``bootstrap --proteome human,mouse`` stages what ``self`` and ``self_mouse`` are built from.

    Both species matter: the shipped corpus tables carry a ``self`` channel for each, and a mouse
    run that silently fell back to the human proteome would be a differently-scaled feature under
    the same fitted weight.
    """
    from mhcmatch import mimics
    from mhcmatch.store import fetch_proteome

    for name in ("human", "mouse"):
        assert os.path.exists(fetch_proteome(name)), name
    for cat in ("self", "self_mouse"):
        w = mimics.proteome_window_array(cat, 9)
        assert len(w) > 1_000_000, (cat, len(w))
        assert w.dtype.itemsize == 9
