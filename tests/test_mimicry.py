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
