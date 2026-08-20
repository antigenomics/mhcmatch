"""The mimicry module's two guarantees: the channels partition the peptide, and the
tested-neoantigen database is reachable as an annotation without a fitted model."""
import pytest

from mhcmatch import mimicry


@pytest.mark.parametrize("length", [8, 9, 10, 11, 15])
def test_channels_partition_the_peptide(length):
    """anchor and tcr must cover every position exactly once -- a position counted twice would be
    double-weighted in the aggregate, and one missed would be silently invisible to both channels."""
    m = mimicry.masks(length)
    assert sorted(m["anchor"] + m["tcr"]) == list(range(length))
    assert not set(m["anchor"]) & set(m["tcr"])


def test_anchor_channel_is_the_complement_role_scheme():
    """The anchor channel is complement.ANCHORS, not a second opinion about where anchors are."""
    from mhcmatch.complement import ANCHORS
    assert set(mimicry.masks(9)["anchor"]) == {i % 9 for i in ANCHORS}


def test_annotate_finds_a_known_neoantigen_and_needs_no_fitted_model():
    """annotate() is prior evidence and must work whether or not mimicry_mhc1.json ships."""
    out = {r["peptide"]: r for r in mimicry.annotate(["KLVVVGACGV", "AAAWYLWEV"])}
    kras = out["KLVVVGACGV"]                     # KRAS G12C, assayed and immunogenic
    assert kras["known"] and kras["neoag_distance"] == 0
    assert kras["neoag_nearest"] == "KLVVVGACGV"
    assert not out["AAAWYLWEV"]["known"]


def test_score_says_what_is_missing_rather_than_failing_obscurely():
    """Until the aggregate is fitted and shipped, score() must name the missing artifact."""
    mimicry._CACHE.pop("mhc1", None)
    try:
        mimicry.params("mhc1")
    except FileNotFoundError as e:
        assert "annotate() and masks() do not need it" in str(e)


def test_score_refuses_a_silently_smaller_model():
    """A component missing from refs standardizes to zero, so the aggregate would quietly become a
    different model. That must raise, and the message must name what is gone."""
    with pytest.raises(ValueError, match="self"):
        mimicry.score(["GILGFVFTL"], {("viral", "anchor", 9): None}, cls="mhc1")


def test_rank_extended_appends_columns_without_moving_the_ranking(tmp_path, capsys):
    """--extended/--annotate are columns, not a re-score. If they ever move the order, mimicry has
    silently entered the gate, which is the thing that is explicitly not settled."""
    from mhcmatch.cli import main
    p = tmp_path / "c.scored.csv"
    p.write_text("epitope,best_allele,gene_name,tpm,score\n"
                 "GILGFVFTL,HLA-A*02:01,PMEL,3.0,0.9\n"
                 "KLVVVGACGV,HLA-A*02:01,KRAS,12.0,0.8\n")

    # `--score gate` throughout: since 0.20.0 the aggregate computes all nine of its features
    # before scoring, which loads the self-mimicry reference (6 min 15 s, ~7.5 GB). The property
    # under test is that --extended/--annotate append columns without moving the order, and that is
    # a property of those flags, not of the scorer.
    def run(*flags):
        main(["rank", "table", str(p), "--tumor", "SKCM", "--score", "gate", *flags])
        body = [ln.split("\t") for ln in capsys.readouterr().out.strip().splitlines()]
        return body[0], body[1:]

    head, rows = run()
    ehead, erows = run("--extended", "--annotate", "--no-self")
    assert [r[:len(head)] for r in erows] == rows          # same order, same values
    assert ehead[:len(head)] == head                       # base schema is a strict prefix
    for c, ch in (("viral", "anchor"), ("self", "tcr"), ("thymus", "anchor")):
        assert f"{c}_{ch}" in ehead and f"source_{c}_{ch}" in ehead
    assert "neoag_distance" in ehead


def test_no_self_cannot_be_combined_with_the_aggregate(tmp_path):
    """--no-self drops the self mimicry reference, which supplies `self_tcr` -- BOECRT's
    second-largest coefficient at +0.3154 of 1.3875 total absolute weight. Until 0.20.0 the
    combination ran and silently scored a five-feature model while reporting the nine-feature one.
    """
    from mhcmatch.cli import main
    p = tmp_path / "c.scored.csv"
    p.write_text("epitope,best_allele,gene_name,tpm,score\n"
                 "GILGFVFTL,HLA-A*02:01,PMEL,3.0,0.9\n")
    with pytest.raises(SystemExit, match="self_tcr"):
        main(["rank", "table", str(p), "--score", "aggregate", "--no-self"])


def test_cli_coefficients_reports_both_aurocs(capsys):
    """`mhcmatch mimicry --coefficients` must print the six shipped coefficients and BOTH AUROCs --
    the pooled/within-screen gap is the pooling artifact and hiding it is how it gets misquoted."""
    from mhcmatch.cli import main
    main(["mimicry", "--coefficients"])
    out, err = capsys.readouterr()
    p = mimicry.params("mhc1")
    for f in p["features"]:
        comp, ch = f.rsplit("_", 1)
        assert f"{comp}\t{ch}\t" in out
    assert f"{p['fit']['auroc_pooled']:.3f}" in err
    assert f"{p['fit']['auroc_within_screen_median']:.3f}" in err


def test_self_is_the_recipients_proteome_not_a_constant():
    """`self` means "the recipient's own proteins", so a mouse run must be able to say so.

    Scoring mouse tumour epitopes against the *human* proteome asks whether a mouse peptide
    resembles a human self protein, which is not a tolerance statement about that mouse. The
    category exists so the mouse deliverable can carry both channels rather than mislabel one.
    """
    from mhcmatch import mimics

    assert mimics.PROTEOME_REFS["self"] == ("human",)
    assert mimics.PROTEOME_REFS["self_mouse"] == ("mouse",)

    import inspect
    from mhcmatch import mimicry
    p = inspect.signature(mimicry.load_references).parameters
    assert "self_species" in p
    assert p["self_species"].default == "human", "the fitted coefficients are the human ones"


# ------------------------------------------------------------------ the reference cache (0.20.0)
def test_reference_cache_round_trips_and_agrees_with_a_fresh_build(tmp_path):
    """A cached reference set must score identically to the one it was built from.

    Exercised with `with_self=False` so it stays a unit test: the self proteome is the 6-minute,
    ~7.5 GB part, and what is under test here is the persistence, not the size.
    """
    peps = ["GILGFVFTL", "NLVPMVATV"]
    fresh = mimicry.load_references(cls="mhc1", with_self=False)
    a = mimicry.score(peps, fresh, cls="mhc1", allow_missing=True)

    built = mimicry.load_references(cls="mhc1", with_self=False, cache=tmp_path)
    b = mimicry.score(peps, built, cls="mhc1", allow_missing=True)
    loaded = mimicry.load_references(cls="mhc1", with_self=False, cache=tmp_path)
    c = mimicry.score(peps, loaded, cls="mhc1", allow_missing=True)

    assert list(tmp_path.iterdir()), "nothing was written to the cache directory"
    for x, y, z in zip(a, b, c):
        assert x.logodds == pytest.approx(y.logodds)
        assert x.logodds == pytest.approx(z.logodds), "cached load disagrees with a fresh build"
        for comp in mimicry.COMPONENTS:
            if comp in x.components:
                for ch in mimicry.CHANNELS:
                    assert x.components[comp][ch] == pytest.approx(z.components[comp][ch])


def test_reference_cache_key_changes_when_the_projection_does(tmp_path, monkeypatch):
    """A cache entry is keyed on what produced it. Bumping CACHE_VERSION must miss, not reuse --
    a stale projection is a silently wrong feature."""
    lengths = [9]
    before = mimicry._fingerprint(None, "mhc1", False, "human", lengths)
    monkeypatch.setattr(mimicry, "CACHE_VERSION", mimicry.CACHE_VERSION + 1)
    assert mimicry._fingerprint(None, "mhc1", False, "human", lengths) != before
    # and the same inputs are stable
    monkeypatch.undo()
    assert mimicry._fingerprint(None, "mhc1", False, "human", lengths) == before


def test_backing_reads_str_pairs_from_arrays():
    """`features` indexes the backing only for a best hit, so it stays memory-mapped; it still has
    to hand back plain strings."""
    import numpy as np
    b = mimicry._Backing(np.array([b"GILGFVFTL"], dtype="S9"), np.array([b"FLU"], dtype="S3"))
    assert len(b) == 1
    assert b[0] == ("GILGFVFTL", "FLU")


def test_corpus_R_is_the_fitted_luksza_form():
    """`corpus_R` must compute Z = sum_d n_d exp(-k(a0 - (L - d))), not a rescaling of it.

    Two mistakes are silent and both change the ranking rather than the calibration: flipping the
    sign on `d` (which upweights *distant* neighbours) and dropping the `L` term (which is a real
    per-row factor because peptide length varies). Either one saturates R -- measured against the
    fitted column over 328,276 cached peptides, the flipped variant has mean R 0.771 with 77.2% of
    peptides above 0.5, against a fitted mean of 3.29e-5 and none above 0.5, ranking differently at
    Spearman +0.705. So the formula is pinned here against its closed form, not merely exercised.
    """
    import math

    class _Hit:
        def __init__(self, score): self.score = score

    counts = {0: 1, 1: 28, 2: 585}                      # a real thymus row: AAAAAAAVL
    hits = [_Hit(d) for d, c in counts.items() for _ in range(c)]

    class _Index:
        def search(self, q, params): return hits

    pep = "AAAAAAAVL"                                    # L = 9
    refs = {("thymus", "tcr", len(pep)): (_Index(), 0, None)}
    got = mimicry.corpus_R([pep], refs, components=("thymus",))[0]

    k, a0 = mimicry.SHAPES["thymus"]
    z = sum(c * math.exp(-k * (a0 - (len(pep) - d))) for d, c in counts.items())
    assert got["thymus"] == pytest.approx(z / (1.0 + z), rel=1e-12)
    assert [got["thymus_n0"], got["thymus_n1"], got["thymus_n2"]] == [1, 28, 585]

    # the linear regime the appendix documents: Z stays far below 1, so R ~ Z
    assert z < 1.4e-3
    assert got["thymus"] == pytest.approx(z, rel=1e-3)

    # nearer neighbours must weigh more -- the sign check that the flipped variant fails
    solo = {d: {("thymus", "tcr", len(pep)): (type("I", (), {"search": lambda s, q, p, d=d: [_Hit(d)]})(), 0, None)}
            for d in (0, 1, 2)}
    r = [mimicry.corpus_R([pep], solo[d], components=("thymus",))[0]["thymus"] for d in (0, 1, 2)]
    assert r[0] > r[1] > r[2]


@pytest.mark.parametrize("peptide", ["AAAKFVAAWTLKAAA", "PKYVKQNTLKLATGM", "GELIGILNAAKVPAD"])
def test_class_ii_masks_follow_the_register_not_the_length(peptide):
    """A class-II ligand is anchored by a 9-mer core that floats, so its face is a function of the
    register. `masks` took only a length until 0.21.0 and `corpus_R` accepted `cls` and ignored it,
    which read every class-II ligand on the class-I layout -- a confident, wrong face."""
    from mhcmatch import complement

    m = mimicry.masks(len(peptide), "mhc2", peptide)
    assert sorted(m["anchor"] + m["tcr"]) == list(range(len(peptide)))
    assert m["anchor"] == sorted(complement.mhc2_anchors(peptide))
    assert m["anchor"] != mimicry.masks(len(peptide))["anchor"]      # not the class-I split

    # pinning the register moves the face; the class-I mask cannot
    shifted = mimicry.masks(len(peptide), "mhc2", peptide, register=1)
    assert shifted["anchor"] == [1, 4, 6, 9]
