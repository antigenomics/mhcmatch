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

    def run(*flags):
        main(["rank", "table", str(p), "--tumor", "SKCM", *flags])
        body = [ln.split("\t") for ln in capsys.readouterr().out.strip().splitlines()]
        return body[0], body[1:]

    head, rows = run()
    ehead, erows = run("--extended", "--annotate", "--no-self")
    assert [r[:len(head)] for r in erows] == rows          # same order, same values
    assert ehead[:len(head)] == head                       # base schema is a strict prefix
    for c, ch in (("viral", "anchor"), ("self", "tcr"), ("thymus", "anchor")):
        assert f"{c}_{ch}" in ehead and f"source_{c}_{ch}" in ehead
    assert "neoag_distance" in ehead


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
