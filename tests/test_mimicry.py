"""The mimicry module's two guarantees: the channels partition the peptide, and the
tested-neoantigen database is reachable as an annotation without a fitted model."""
import os

import pytest

from mhcmatch import mimics, mimicry


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


@pytest.mark.hfdata
def test_annotate_finds_a_known_neoantigen_and_needs_no_fitted_model():
    """annotate() is prior evidence and must work whether or not mimicry_mhc1.json ships.

    Needs the deposit: it builds the tested-neoantigen index from `isalgo/pmhc_data`."""
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


@pytest.mark.hfdata
def test_rank_extended_appends_columns_without_moving_the_ranking(tmp_path, capsys):
    """--extended/--annotate are columns, not a re-score. If they ever move the order, mimicry has
    silently entered the gate, which is the thing that is explicitly not settled.

    **The slowest test in the suite (13.1 s).** The `--score gate` note below buys back the 6 min
    15 s / ~7.5 GB self reference; what remains is building the viral and thymic reference indexes
    and running the whole CLI end to end, twice."""
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
        # rstrip newlines only. `.strip()` would eat a trailing TAB, which is exactly what an
        # empty last column (`variant_type`) emits, and the base and extended runs would then
        # disagree on width for a reason that is not about the flags under test.
        body = [ln.split("\t") for ln in capsys.readouterr().out.rstrip("\n").splitlines()]
        return body[0], body[1:]

    head, rows = run()
    ehead, erows = run("--extended", "--annotate", "--no-self")
    assert [r[:len(head)] for r in erows] == rows          # same order, same values
    assert ehead[:len(head)] == head                       # base schema is a strict prefix
    for c, ch in (("viral", "anchor"), ("self", "tcr"), ("thymus", "anchor")):
        assert f"{c}_{ch}" in ehead and f"source_{c}_{ch}" in ehead
    assert "neoag_distance" in ehead


def test_the_aggregate_no_longer_needs_the_host_proteome_index():
    """`--no-self` and `--score aggregate` were mutually exclusive through 0.20.0, because BOECRT
    scored on `self_tcr` -- its second-largest coefficient at +0.3154 of 1.3875 total absolute
    weight -- and that forced the host-proteome index: 6 min 15 s and ~7.5 GB, the largest single
    cost in the package.

    EPIC v3 scores `C_corpus_self` again -- but from a 64 KB k-mer count table, not a trie. The
    proteome index is still off the ranking path, and the property to pin is that rather than the
    absence of a `self` feature: the channels the CLI builds must be computable without ever
    calling `load_references`."""
    import mhcmatch.mimicry as MM
    from mhcmatch import rank
    from mhcmatch.cli import _aggregate_channels

    feats = rank.aggregate()["features"]
    assert "C_corpus_self" in feats and "C_corpus_thymus" in feats
    assert set(rank.CHANNEL_COLUMNS) <= set(feats)

    def _refuse(*a, **k):                       # the index build, if the corpus path ever calls it
        raise AssertionError("corpus channels must not build a reference index")

    real, MM.load_references = MM.load_references, _refuse
    try:
        got = _aggregate_channels("mhc1", no_self=True)(["GILGFVFTL", "KLVVVGACGV"])
    finally:
        MM.load_references = real
    assert set(got) == set(rank.CHANNEL_COLUMNS)
    assert all(len(v) == 2 and all(0.0 <= x <= 1.0 for x in v) for v in got.values())


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
@pytest.mark.hfdata
def test_backing_reads_str_pairs_from_arrays():
    """`features` indexes the backing only for a best hit, so it stays memory-mapped; it still has
    to hand back plain strings.

    The windows are a fixed-width byte array and the sources are a str list, because a source is a
    ``;``-joined accession list of mean length 7.5 and maximum 2,141 -- padding every row to the
    longest cost 56.3 MB for 26,302 thymic entries, ~99 % of it NUL."""
    import numpy as np
    b = mimicry._Backing(np.array([b"GILGFVFTL"], dtype="S9"), ["FLU"])
    assert len(b) == 1
    assert b[0] == ("GILGFVFTL", "FLU")


def _spectrum_from(peptides, kappa, k=3, cls="mhc1"):
    """A spectrum dict built by hand from a peptide list -- the same shape `corpus_spectrum` returns."""
    import numpy as np
    T = np.zeros(20 ** k)
    for p in peptides:
        idx = mimicry.face_kmers(p, cls, k)
        if idx.size:
            np.add.at(T, idx, 1.0)
    beta = np.exp(-kappa)
    K = (1 - beta) * np.eye(20) + beta * np.ones((20, 20))
    C = T.reshape((20,) * k)
    for ax in range(k):
        C = np.moveaxis(np.tensordot(K, C, axes=([1], [ax])), 0, ax)
    return {"thymus": (C.ravel(), float(T.sum()), k)}, T


def test_corpus_R_is_the_exact_all_vs_all_sum_not_a_truncation():
    """The contraction must reproduce a brute-force sum over EVERY reference k-mer.

    This is the property the radius-capped search never had: it walked a Hamming ball and stopped,
    capturing a median 0.4999 of the true sum on real 9-mers. Here the decaying exponent is the
    threshold and nothing is dropped, so the check is equality with the definition, not closeness.
    """
    import numpy as np
    rng = np.random.default_rng(0)
    AA = mimicry.AA
    ref = ["".join(rng.choice(list(AA), size=n)) for n in rng.integers(8, 12, size=400)]
    kappa, k = 2.25, 3
    spec, T = _spectrum_from(ref, kappa, k)
    beta = np.exp(-kappa)

    nz = np.nonzero(T)[0]
    dig = (nz[:, None] // (20 ** np.arange(k)[::-1])) % 20          # every occupied k-mer, unpacked
    for pep in ["GILGFVFTL", "SIINFEKL", "KLINSQINL", "GILGFVFTLAV"]:
        idx = mimicry.face_kmers(pep, "mhc1", k)
        qd = (idx[:, None] // (20 ** np.arange(k)[::-1])) % 20
        brute = sum(float((T[nz] * beta ** ((dig != qd[i]).sum(1))).sum()) for i in range(len(idx)))
        got = mimicry.corpus_R([pep], spec)[0]["thymus"] * idx.size * spec["thymus"][1]
        assert got == pytest.approx(brute, rel=1e-12), pep


def test_corpus_R_is_a_density_per_query_window():
    """`rho` divides by both the query's window count and the reference's total mass.

    Without the first it is a length detector; without the second `thymus` (26,513 peptides) and
    `self` (12 M proteome windows) are not on one scale and 'thymus makes the others redundant' is
    a statement about set size rather than about biology."""
    import numpy as np
    ref = ["GILGFVFTL", "GILGFVFTLA", "SIINFEKL"]
    spec, T = _spectrum_from(ref, 2.25)
    table, n, k = spec["thymus"]
    pep = "KLINSQINL"
    idx = mimicry.face_kmers(pep, "mhc1", k)
    assert mimicry.corpus_R([pep], spec)[0]["thymus"] == pytest.approx(
        float(table[idx].sum()) / (idx.size * n))
    # doubling the reference set leaves the density unchanged: N_k doubles with the sum
    spec2, _ = _spectrum_from(ref * 2, 2.25)
    assert mimicry.corpus_R([pep], spec2)[0]["thymus"] == pytest.approx(
        mimicry.corpus_R([pep], spec)[0]["thymus"])


def test_the_face_is_contiguous_and_k_is_the_only_width_every_length_supports():
    """k=3 is a consequence, not a preference: the face is L-5 wide and the shortest ligand is 8."""
    for L in range(8, 16):
        sel = mimicry.masks(L, "mhc1")["tcr"]
        assert sel == list(range(3, L - 2)), L                       # contiguous, width L-5
    assert mimicry.CORPUS_K == 3
    assert mimicry.face_kmers("SIINFEKL").size == 1                  # W=3, exactly one 3-mer
    assert mimicry.face_kmers("SIINFEKL", k=4).size == 0             # W=3 cannot carry a 4-mer
    assert mimicry.face_kmers("GILGFVFTLAV").size == 4               # W=6 -> 6-3+1


def test_corpus_R_averages_over_query_windows_so_length_alone_cannot_move_it():
    """The shipped fixed-face column varied 17x in mean across lengths 8-11 (Spearman -0.502 with
    length); the `m_k(q)` divisor is what removes that. The face is `L - 5` wide so it cannot be
    held constant across lengths -- but a *homopolymer* face is every window at once, so its density
    must be exactly the single-window value at every length."""
    spec, _ = _spectrum_from(["AAAWWWAAA", "CCCYYYCCC", "GGGAAAGGG"], 2.25)
    table, n, k = spec["thymus"]
    one = float(table[mimicry.face_kmers("AAAAAAAA", "mhc1", k)].sum()) / n     # W=3, a single AAA
    for L in range(8, 13):
        pep = "A" * L                                   # face is A*(L-5): L-7 windows, all AAA
        got = mimicry.corpus_R([pep], spec)[0]["thymus"]
        assert got == pytest.approx(one), (L, got, one)


def test_corpus_R_takes_any_position_additive_kernel():
    """BLOSUM62 is the same contraction with a different 20x20, so the graded Luksza form is exact
    too. Only gapped alignment fails to factorise, which is why `features` keeps its index."""
    import numpy as np
    from seqtree import SubstitutionMatrix
    AA = mimicry.AA
    S = np.array([[SubstitutionMatrix.blosum62().similarity(a, b) for b in AA] for a in AA], float)
    rng = np.random.default_rng(1)
    ref = ["".join(rng.choice(list(AA), size=9)) for _ in range(200)]
    k, kb = 3, 0.30
    T = np.zeros(20 ** k)
    for p in ref:
        np.add.at(T, mimicry.face_kmers(p, "mhc1", k), 1.0)
    C = T.reshape((20,) * k)
    for ax in range(k):
        C = np.moveaxis(np.tensordot(np.exp(kb * S), C, axes=([1], [ax])), 0, ax)
    spec = {"thymus": (C.ravel(), float(T.sum()), k)}

    nz = np.nonzero(T)[0]
    dig = (nz[:, None] // (20 ** np.arange(k)[::-1])) % 20
    pep = "GILGFVFTL"
    idx = mimicry.face_kmers(pep, "mhc1", k)
    qd = (idx[:, None] // (20 ** np.arange(k)[::-1])) % 20
    brute = sum(float((T[nz] * np.exp(kb * S[np.repeat(qd[i][None, :], len(nz), 0), dig].sum(1))).sum())
                for i in range(len(idx)))
    got = mimicry.corpus_R([pep], spec)[0]["thymus"] * idx.size * spec["thymus"][1]
    assert got == pytest.approx(brute, rel=1e-12)


def test_load_references_no_longer_takes_a_cache(monkeypatch):
    """The reference cache is gone (0.24.0). A caller who passes `cache=` gets a TypeError naming
    it, which is the deprecation; the env var is simply unread and has no failure mode."""
    import inspect
    assert "cache" not in inspect.signature(mimicry.load_references).parameters
    assert not [n for n in mimicry.__all__ if "CACHE" in n]
    assert not hasattr(mimicry, "CACHE_VERSION") and not hasattr(mimicry, "REFERENCE_CACHE_ENV")


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


def test_an_unmeasured_component_reports_nan_not_zero():
    """`autoimmune` is a safety read-out. Under `allow_missing` the self component has no reference
    index, and standardizing its absence to the training mean made it contribute exactly zero --
    which prints as `0` and reads as "no self-similarity found" when the truth is "never looked".

    Only reachable through `allow_missing=True`, but reachable by default since 0.21.0: EPIC does
    not score on `self_tcr`, so `--no-self --score aggregate` is now allowed where 0.20.0 refused
    it, and `--extended` prints these columns beside a score that is perfectly well defined.
    """
    import math

    pep = "GILGFVFTL"
    # a refs dict carrying viral and thymus but not self -- exactly load_references(with_self=False)
    refs = {(c, ch, 9): (_EmptyIndex(), 1000, None)
            for c in ("viral", "thymus") for ch in mimicry.CHANNELS}

    with pytest.raises(ValueError, match="would standardize to"):
        mimicry.score([pep], refs)

    s, = mimicry.score([pep], refs, allow_missing=True)
    assert all(math.isnan(v) for v in s.components["self"].values()), s.components["self"]
    assert math.isnan(s.autoimmune), "an unmeasured self channel must not read as zero risk"
    assert math.isnan(s.logodds), "the fitted coefficients describe the full component set"
    # the measured components are still numbers
    assert all(math.isfinite(v) for v in s.components["viral"].values())
    assert all(math.isfinite(v) for v in s.components["thymus"].values())


class _EmptyIndex:
    """A reference index that matches nothing -- enough to exercise the component bookkeeping."""

    def search(self, q, params):
        return []


# ------------------------------------------------------------------ the safety read-out
def _score(peptide, nearest, autoimmune=-0.5):
    """A `MimicryScore` built by hand. `safety` reads only `.nearest`, `.peptide` and `.autoimmune`,
    so no reference index is needed to exercise the join."""
    return mimicry.MimicryScore(peptide=peptide, components={}, logodds=0.0,
                                autoimmune=autoimmune, nearest=nearest)


def _hit(pep, source, subs=1, n=1):
    return {"subs": subs, "peptide": pep, "source": source, "n": n}


def test_safety_needs_the_symbol_map_to_resolve_an_accession(monkeypatch):
    """The thymus deposit names sources as UniProt accessions and `expression.safety_profile` is
    keyed on HGNC symbols. Without a map the accession does not resolve and `profile` comes back
    empty -- which reads as *no risk found*, the dangerous direction to be wrong in."""
    from mhcmatch import expression as EX
    monkeypatch.setattr(EX, "safety_profile",
                        lambda gene, top=10: [("Heart", 812.0)] if gene == "TTN" else [])

    s = _score("ESDPIVAQY", {"thymus": {"tcr": _hit("ESDPIVAQF", "Q8WZ42")}})

    with_map, = mimicry.safety([s], symbols={"Q8WZ42": "TTN"})
    hit, = with_map["hits"]
    assert hit["gene"] == "TTN"                       # resolved through the map
    assert hit["source"] == "Q8WZ42"                  # ... and the deposit's own id is preserved
    assert hit["profile"] == [("Heart", 812.0)]
    assert hit["component"] == "thymus" and hit["mimic"] == "ESDPIVAQF"

    without, = mimicry.safety([s])
    assert without["hits"][0]["profile"] == [], "an unresolved accession must not invent a profile"
    assert without["hits"][0]["gene"] == "Q8WZ42", "the raw source is still reported"


def test_safety_returns_the_self_mimic_even_though_it_has_no_source(monkeypatch):
    """The `self` component is built from proteome windows, which carry no source column, so it is
    not resolvable to a gene however good the map is. The mimic itself is still returned."""
    from mhcmatch import expression as EX
    monkeypatch.setattr(EX, "safety_profile", lambda gene, top=10: [("Lung", 1.0)])

    s = _score("GILGFVFTL", {"self": {"tcr": _hit("GILGFVFTV", "")}})
    out, = mimicry.safety([s])
    hit, = out["hits"]
    assert hit["component"] == "self" and hit["mimic"] == "GILGFVFTV"
    assert hit["gene"] == "" and hit["profile"] == []
    assert out["autoimmune_logodds"] == -0.5


def test_safety_reports_only_the_tolerance_side():
    """`viral` is a priming argument, not a withdrawal one, so it is not a safety hit."""
    s = _score("GILGFVFTL", {"viral": {"tcr": _hit("GILGFVFTV", "P03485")}})
    out, = mimicry.safety([s])
    assert out["hits"] == []


def test_safety_takes_no_tumour_argument():
    """`safety` accepted `tumor` at positional #2 through 1.4.0 and never read it, because
    `expression.safety_profile` conditions on no context at all. `safety(scores, "SKCM")` therefore
    returned the pooled tolerance profile while reading as if it had been narrowed to melanoma --
    on the read-out whose job is to say which tissue you cannot afford to damage, that is the
    dangerous direction. The parameter is gone; #2 is `top`."""
    import inspect

    params = list(inspect.signature(mimicry.safety).parameters)
    assert "tumor" not in params, "a tumour argument that nothing reads must not be accepted"
    assert params == ["scores", "top", "symbols"]


# ------------------------------------------------------------------ the calibrated read-out
def test_probability_is_in_the_unit_interval_and_monotone_in_the_logodds():
    """`probability` maps the aggregate log-odds against a *named* corpus. Whatever the corpus, the
    map is a probability and it must not reorder -- ranking is `logodds`' job and stays its job."""
    scores = [_score(f"P{i}", {}, autoimmune=0.0) for i in range(5)]
    for s, lo in zip(scores, (-8.0, -2.0, 0.0, 1.5, 6.0)):
        s.logodds = lo
    ps = mimicry.probability(scores, corpus="screens")
    assert all(0.0 < p < 1.0 for p in ps)
    assert all(ps[i] < ps[i + 1] for i in range(len(ps) - 1))


def test_probability_refuses_an_unnamed_corpus():
    """An absolute probability is a property of the corpus's prevalence, not of the peptide, so a
    corpus this model was never calibrated on has to raise rather than pick one."""
    with pytest.raises(ValueError, match="no calibration for corpus"):
        mimicry.probability([_score("GILGFVFTL", {})], corpus="nosuchcorpus")


# ------------------------------------------------------------------ locus weighting
def test_locus_weights_pool_overlapping_registers_of_one_hotspot():
    """The KRAS G12 family is one locus at seven residues, not a dozen independent observations.

    This is the correction that matters for any statement about a peptide set's k-mer vocabulary:
    a recurrent public hotspot enters an assayed deposit many times over as overlapping registers,
    and an unweighted table reads that as evidence about biology when it is evidence about what got
    tested.
    """
    kras = ["KLVVVGAVGV", "KLVVVGADGV", "KLVVVGACGV", "VVGAVGVGK", "VVGADGVGK",
            "GADGVGKSAL", "GADGVGKSA"]
    others = ["SIINFEKLAA", "GILGFVFTLA", "NLVPMVATVA"]
    w = mimicry.locus_weights(kras + others)
    assert all(x == pytest.approx(1.0 / len(kras)) for x in w[:len(kras)])
    assert all(x == 1.0 for x in w[len(kras):])
    assert sum(w) == pytest.approx(1.0 + len(others))     # one locus for KRAS, one each for the rest


def test_locus_weights_width_is_the_knob_and_a_short_match_is_not_a_locus():
    a, b = "AAAAAACCCCC", "CCCCCGGGGGG"                   # share exactly 5
    assert mimicry.locus_weights([a, b], w=5) == [0.5, 0.5]
    assert mimicry.locus_weights([a, b], w=6) == [1.0, 1.0]


def test_locus_weighting_shrinks_a_deposit_and_leaves_the_proteome_alone():
    """`weights="locus"` corrects a SAMPLING bias, so it applies to assayed deposits and not to a
    proteome, where every window appears exactly once by construction."""
    _t, n_plain = mimicry.corpus_counts(comp="thymus")
    _t, n_locus = mimicry.corpus_counts(comp="thymus", weights="locus")
    assert 0.0 < n_locus < n_plain                        # recurrence really was there
    _s, s_plain = mimicry.corpus_counts(comp="self")
    _s, s_locus = mimicry.corpus_counts(comp="self", weights="locus")
    assert s_locus == s_plain                             # and the proteome is untouched


def test_weighted_density_is_still_a_density():
    spec = mimicry.corpus_spectrum(components=("thymus",), weights="locus")
    v = mimicry.corpus_R(["GILGFVFTL", "SIINFEKL", "NLVPMVATV"], spec)
    assert all(0.0 <= r["thymus"] <= 1.0 for r in v)


def test_an_unknown_weighting_is_refused_rather_than_ignored():
    with pytest.raises(ValueError, match="weights must be"):
        mimicry.corpus_counts(comp="thymus", weights="expression")


# ------------------------------------------------------------ the division of labour, pinned
def test_the_kmer_tables_are_for_the_corpus_term_and_nothing_else():
    """`corpus_R` is a weighted sum and cannot name a hit. Everything that must name one still
    searches.

    This is a contract, not an implementation note. The contraction gives
    `sum_r w(q, r)` over the whole reference set, which is the right object for a scored feature and
    the wrong one for safety: `mimicry.safety` has to say *which* self protein a candidate resembles
    so `expression.safety_profile` can say what tissue presents it, and a density cannot be
    interrogated for that. `annotate` has the same requirement for the tested-neoantigen listing.
    So the indexed search stays on exactly those paths and this test fails if it is ever quietly
    swapped for a lookup.
    """
    import inspect

    import mhcmatch.mimicry as MM
    import mhcmatch.proteome as PR

    # the search paths, by the seqtree symbols they name
    for fn in (MM.features, MM.annotate):
        src = inspect.getsource(fn)
        assert "seqtree" in src or "index.search" in src, fn.__name__
    assert "Index" in inspect.getsource(PR)                    # self-origin safety scan

    # and the corpus path, which must NOT
    for fn in (MM.corpus_counts, MM.contract, MM.corpus_R, MM.corpus_spectrum):
        assert "seqtree" not in inspect.getsource(fn), fn.__name__


def test_annotate_names_the_neoantigen_it_matched():
    """The known-neoantigen listing reports an identity and a distance, which is what a search is
    for. A density would give neither."""
    got = mimicry.annotate(["KLVVVGACGV", "AAAAAAAAA"], cls="mhc1")
    assert set(got[0]) >= {"neoag_distance", "neoag_nearest", "neoag_n_within"}
    assert got[0]["neoag_nearest"] and got[0]["neoag_distance"] == 0     # it is a tested neoantigen
    assert isinstance(got[0]["neoag_n_within"], int)


def test_features_keeps_the_source_protein_that_safety_needs():
    """`safety` resolves a tissue from the mimic's SOURCE, so `features` must carry it through."""
    refs = mimicry.load_references(cls="mhc1", with_self=False)
    f = mimicry.features(["GILGFVFTL"], refs, cls="mhc1")[0]
    near = f.get("nearest_thymus_tcr")
    assert near is not None and near.get("peptide")
    assert "source" in near                    # the field safety() joins on


@pytest.mark.hfdata
def test_the_class_two_self_table_is_the_proteomes_own_kmers_and_class_one_is_untouched():
    """A proteome window has no register, so `self` at class II is not a projected face.

    ``thymus`` and ``viral`` are ligand deposits and do have registers, so they keep the per-window
    class-II face. ``self`` is a *proteome* -- nothing in it is presented -- so resolving a register
    for each of its ~192 M windows is the model applied outside its domain, and it measured >25 min
    and ~10.7 GB against 1.7 s for the thymic deposit. It is the window's own k-mers instead.

    Two things are asserted, and the second is the one that matters: **class I is bit-identical**.
    ``C_corpus_self`` is a fitted feature of the shipped model and it was fitted on class-I rows, so
    a class-II definition may change and a class-I table may not.
    """
    import numpy as np

    t2, n2 = mimicry.corpus_counts(None, "mhc2", "self")
    assert n2 > 1e8, "the class-II self table must carry multiplicity, not one count per k-mer"
    assert (t2 > 0).all(), "every 3-mer occurs somewhere in a proteome"

    # class I: the fixed positional face, unchanged. Hashed rather than eyeballed.
    for comp, want in (("thymus", 140482.0), ("viral", 136618.0), ("self", 121968158.0)):
        t1, n1 = mimicry.corpus_counts(None, "mhc1", comp)
        assert n1 == want, comp
        assert np.isfinite(t1).all()

    # and the densities stay in [0, 1] and on a comparable scale across the two classes
    spec = mimicry.corpus_spectrum(cls="mhc2", components=("thymus", "self", "viral"))
    got = mimicry.corpus_R(["AAAKFVAAWTLKAAA"], spec, cls="mhc2")[0]
    assert set(got) == {"thymus", "self", "viral"}
    assert all(0.0 < v < 1.0 for v in got.values()), got


@pytest.mark.hfdata
def test_the_vendored_corpus_tables_are_current_and_rebuild_bit_identically():
    """The shipped ``corpus_tables.npz`` is what the current code and the current script produce.

    Two independent guards, because the artifact can go stale in two different ways:

    1. **Provenance** -- every combination ``mhcmatch build corpus`` declares is present, and
       the stamped version is this one. A release that bumps ``__version__`` without rerunning the
       builder fails here, exactly as the vendored anchor models do.
    2. **Content** -- a live rebuild of the four *deposit* channels is bit-identical to the shipped
       table. Those are 0.7-1.6 s each; the two ``self`` channels are 51.4 s and 14.5 s and are
       checked on their totals instead, which is the number a wrong face or a changed proteome
       moves first.

    The point of the artifact is that it is **145.4 kB against 115.6 s** of per-process rebuild, so
    this test is where the price of shipping it gets paid, once, instead of in every consumer.
    """
    import numpy as np

    import mhcmatch
    from mhcmatch import _build as build          # the builder lives in the package; `mhcmatch
    #                                               build corpus` is the only entry point

    mimicry._VENDORED = None                        # force a load off disk, not a warm process
    assert mimicry._vendored_counts("mhc1", "thymus", mimicry.CORPUS_K, "human") is not None, \
        "no vendored corpus tables in this build; run: mhcmatch build corpus"
    meta = mimicry._VENDORED["_meta"]
    assert meta["_"]["version"] == mhcmatch.__version__, \
        "vendored corpus tables are stale for this version; rerun: mhcmatch build corpus"
    assert meta["_"]["k"] == mimicry.CORPUS_K
    for cls, comp, sp in build.CORPUS_COMBOS:               # every combination the script declares
        assert mimicry._vendored_counts(cls, comp, mimicry.CORPUS_K, sp) is not None, (cls, comp, sp)

    for cls, comp in (("mhc1", "thymus"), ("mhc1", "viral"),
                      ("mhc2", "thymus"), ("mhc2", "viral")):
        shipped, _ = mimicry.corpus_counts(None, cls, comp)
        mimicry._COUNTS.clear()
        mimicry._VENDORED = {}                       # empty, not None: skip the artifact entirely
        fresh, _ = mimicry.corpus_counts(None, cls, comp)
        mimicry._COUNTS.clear()
        mimicry._VENDORED = None
        assert np.array_equal(shipped, fresh), \
            f"vendored {cls}/{comp} differs from a live rebuild; rerun: mhcmatch build corpus"

    # the two expensive channels, on their totals -- see the docstring
    for cls, comp, sp, want in (("mhc1", "self", "human", 121_968_158.0),
                                ("mhc1", "self", "mouse", 112_565_681.0),
                                ("mhc2", "self", "human", 110_932_623.0),
                                ("mhc2", "self", "mouse", 101_989_053.0)):
        T, n = mimicry.corpus_counts(None, cls, comp, self_species=sp)
        assert n == want, (cls, comp, sp, n)
        assert (T > 0).all()


@pytest.mark.hfdata
def test_the_vendored_table_is_bypassed_for_anything_off_the_default_path():
    """A custom deposit, locus weighting or a different ``k`` each define a *different* table.

    The shipped artifact is the default-path table only; silently returning it for a caller who
    asked for something else would be a wrong answer delivered fast, which is the failure mode a
    cache layer usually ships with.
    """
    import numpy as np

    default, _ = mimicry.corpus_counts(None, "mhc1", "thymus")
    weighted, _ = mimicry.corpus_counts(None, "mhc1", "thymus", weights="locus")
    assert not np.array_equal(default, weighted), "locus weighting must not read the vendored table"
    k4, _ = mimicry.corpus_counts(None, "mhc1", "thymus", k=4)
    assert k4.size == 20 ** 4, "a non-default k must not be served the k=3 table"


@pytest.mark.hfdata
def test_load_references_builds_only_the_lengths_asked_for():
    """The index is per-length and so is its cost, so a run pays for the lengths it queries.

    One proteome pass is ~11 s for ``window_array``, and at class II a further ~1.0 min to resolve a
    register for each of 12,685,964 windows. The class admits **fifteen** lengths, so the
    unrestricted build is **~19 min** against **83.1 s** for one; class I is **65.4 s** against
    **16.5 s**. That is the difference between the annotation path being usable and not.

    ``None`` still means every admitted length, so an existing caller is unchanged, and a length
    the class never admitted is still absent rather than newly an error.
    """
    r = mimicry.load_references(None, "mhc2", with_self=False, lengths=[15])
    assert {L for _, _, L in r} == {15}
    assert len(r) == len(mimicry.COMPONENTS[:2]) * len(mimicry.CHANNELS), r.keys()

    r2 = mimicry.load_references(None, "mhc2", with_self=False, lengths=[15, 3, 99])
    assert {L for _, _, L in r2} == {15}, "a length the class does not admit is dropped, not built"

    full = mimicry.load_references(None, "mhc2", with_self=False)
    assert {L for _, _, L in full} == set(mimics._LEN["mhc2"]), "None still means every length"


# --------------------------------------------------------------- the wildcard mask and BLOSUM62

def test_the_wildcard_mask_keeps_every_window_where_the_slice_deletes_the_anchors():
    """`slice` is why `k = 3` was forced; `wildcard` is what lifts the constraint.

    The class-I face is `L - 5` wide under `slice`, so an 8-mer supplies three residues and has no
    4-mer window at all -- a structural zero that reads as a low score. Wildcarding the anchors in
    place leaves `L - k + 1` windows at every length, so `k = 4` and `k = 5` become available.
    """
    assert mimicry.face_kmers("SIINFEKL").size == 1                       # W = 3, one 3-mer
    assert mimicry.face_kmers("SIINFEKL", k=4).size == 0                  # W = 3 cannot carry it
    assert mimicry.face_kmers("SIINFEKL", k=4, mask="wildcard").size == 5     # 8 - 4 + 1
    assert mimicry.face_kmers("SIINFEKL", k=5, mask="wildcard").size == 4     # 8 - 5 + 1
    for pep in ("SIINFEKL", "GILGFVFTL", "GILGFVFTLA", "GILGFVFTLAV"):
        for k in (3, 4, 5):
            assert mimicry.face_kmers(pep, k=k, mask="wildcard").size == len(pep) - k + 1


def test_the_wildcard_is_the_neutral_element_of_the_normalised_blosum_kernel():
    """`S(X, a) = S(a, a)` is well posed only against a self-score-normalised kernel.

    Normalised, `K[u, u] = 1` for every residue, so "matches perfectly" *is* a factor of one and the
    wildcard row and column are all ones -- a masked position drops out of the product exactly.
    Raw, BLOSUM62's diagonal runs 4..11 half-bits, so the same rule makes a masked position weight
    every reference k-mer by what sits at the position the mask was supposed to remove.
    """
    import numpy as np

    K = mimicry.blosum62_kernel(0.4, normalise=True)
    assert K.shape == (21, 21)
    assert (K.diagonal() == 1.0).all()
    assert (K[mimicry.WILDCARD] == 1.0).all() and (K[:, mimicry.WILDCARD] == 1.0).all()
    assert K[:20, :20].max() == 1.0                      # nothing beats an identical residue

    raw = mimicry.blosum62_kernel(0.4, normalise=False)
    d = raw.diagonal()[:20]
    assert d.max() / d.min() == pytest.approx(np.exp(0.4 * (11 - 4)))     # the Trp/Ala spread


def test_the_wildcard_contraction_is_the_exact_all_vs_all_sum():
    """Same guarantee as the slice path, on the wider alphabet and a graded kernel.

    The factorisation `prod_p K[u_p, x_p]` holds for any position-additive ungapped score, so the
    contraction is not an approximation of the Luksza sum -- it is the sum.
    """
    import numpy as np

    rng = np.random.default_rng(0)
    k, code = 3, {a: i for i, a in enumerate(mimicry.AA)}
    refs = ["".join(rng.choice(list(mimicry.AA), 9)) for _ in range(120)]
    T = np.zeros(21 ** k)
    for r in refs:
        np.add.at(T, mimicry.face_kmers(r, "mhc1", k, mask="wildcard"), 1.0)

    def wc(p):
        c = [code[a] for a in p]
        for i in mimicry.masks(len(p), "mhc1")["anchor"]:
            c[i] = mimicry.WILDCARD
        return c

    for kappa, norm in ((0.4, True), (0.15, False)):
        K = mimicry.blosum62_kernel(kappa, normalise=norm)
        got = float(mimicry.contract(T, kappa, k, K)[
            mimicry.face_kmers("GILGFVFTL", "mhc1", k, mask="wildcard")].sum())
        qc = wc("GILGFVFTL")
        brute = sum(
            np.prod([K[qc[i + p], rc[j + p]] for p in range(k)])
            for r in refs for rc in (wc(r),)
            for i in range(len(qc) - k + 1) for j in range(len(rc) - k + 1))
        assert got == pytest.approx(brute, rel=1e-12)


def test_a_table_cannot_be_indexed_with_the_wrong_mask():
    """The two conventions pack into different alphabets, so a cross-index is a wrong answer.

    `corpus_R` reads the mask back off the spectrum tuple rather than taking it as an argument,
    which is what makes the mismatch unrepresentable instead of merely discouraged.
    """
    import numpy as np

    k = 3
    peps = ["GILGFVFTL", "NLVPMVATV", "SIINFEKLA"]
    for mask in ("slice", "wildcard"):
        A = mimicry.alphabet(mask)
        T = np.zeros(A ** k)
        for p in peps:
            np.add.at(T, mimicry.face_kmers(p, "mhc1", k, mask=mask), 1.0)
        spec = {"thymus": (mimicry.contract(T, 3.0, k), T.sum(), k, mask)}
        rho = mimicry.corpus_R(["GILGFVFTL"], spec)[0]["thymus"]
        assert 0.0 <= rho <= 1.0
    with pytest.raises(ValueError, match="kernel is"):
        mimicry.contract(np.zeros(21 ** k), 1.0, k, np.ones((20, 20)))
    with pytest.raises(ValueError, match="not A"):
        mimicry.contract(np.zeros(1234), 1.0, k)
    with pytest.raises(ValueError, match="mask must be"):
        mimicry.face_kmers("GILGFVFTL", mask="anchors")


def test_the_slice_path_is_score_identical_to_pre_0_27():
    """The mask argument is additive: the default reproduces the shipped column bit for bit."""
    import numpy as np

    k = 3
    peps = ["GILGFVFTL", "NLVPMVATV", "SIINFEKL", "GILGFVFTLAV"]
    T = np.zeros(20 ** k)
    for p in peps:
        np.add.at(T, mimicry.face_kmers(p, "mhc1", k), 1.0)
    beta = float(np.exp(-3.0))
    hand = (1.0 - beta) * np.eye(20) + beta * np.ones((20, 20))
    assert np.array_equal(mimicry.contract(T, 3.0, k), mimicry.contract(T, 3.0, k, hand))
    three = {"thymus": (mimicry.contract(T, 3.0, k), T.sum(), k)}      # a pre-0.27 3-tuple
    four = {"thymus": (mimicry.contract(T, 3.0, k), T.sum(), k, "slice")}
    assert mimicry.corpus_R(peps, three) == mimicry.corpus_R(peps, four)


def test_the_artifact_defines_the_corpus_geometry_not_a_module_default():
    """`kappa` already came from the artifact; `k`, the mask and the kernel now do too.

    Three halves of one definition. A `kappa` fitted against a graded kernel and then scored under
    the Hamming one is a *different feature*, not a smaller effect, so a refit that changes the
    kernel has to move the scored column the same way a refit that changes `kappa` does.
    """
    import numpy as np

    g = mimicry.corpus_geometry()                    # the shipped artifact
    assert g["k"] == mimicry.CORPUS_K
    assert g["mask"] == "slice"
    assert g["family"] == "blosum62_normalised"
    K = g["kernel"](1.65)
    assert K.shape == (20, 20) and (K.diagonal() == 1.0).all()

    wild = mimicry.corpus_geometry({"corpus_k": 4, "corpus_mask": "wildcard",
                                    "corpus_kernel": "blosum62_raw"})
    assert wild["k"] == 4 and wild["mask"] == "wildcard"
    assert wild["kernel"](0.1).shape == (21, 21)

    # **No default and no fallback.** An artifact that names no kernel, or an unknown one, or an
    # unknown face, raises rather than being scored under whatever the module last preferred --
    # which is exactly how a graded `kappa` would end up contracted against a Hamming kernel.
    with pytest.raises(ValueError, match="unknown corpus kernel"):
        mimicry.corpus_geometry({"corpus_kernel": "levenshtein", "corpus_mask": "slice"})
    with pytest.raises(ValueError, match="unknown corpus kernel"):
        mimicry.corpus_geometry({"corpus_mask": "slice"})
    with pytest.raises(ValueError, match="unknown corpus face mask"):
        mimicry.corpus_geometry({"corpus_kernel": "blosum62_normalised", "corpus_mask": "tcr5"})


def test_a_v4_shaped_artifact_scores_through_the_library_unchanged():
    """The candidate artifact must be loadable, not merely well-formed.

    `aggregate_score` reads only the names the artifact asks for and `rank._finish` supplies every
    one of them, so this is the check that the two actually meet: a ten-term feature list -- the
    nine that ship plus `d_occupancy`, which is emitted and not fitted -- scores every row without
    a code change and without touching the shipped file.
    """
    import numpy as np

    from mhcmatch import rank as R

    v4 = {"model": "EPIC", "version": 4,
          "features": ["pres", "occupancy", "d_occupancy", "expr", "expr_missing",
                       "C_phys_buried", "C_phys_charge",
                       "C_corpus_thymus", "C_corpus_self", "C_corpus_viral"],
          "coef": [0.2382, 0.1270, -0.0184, 0.3303, 0.1064, 0.1593, 0.0559,
                   0.1867, -0.2625, 0.1293],
          "mu": [1.4, 0.02, 0.0, 2.2, 0.017, 0.72, 0.095, 1.2e-3, 2.9e-4, 2.0e-4],
          "sigma": [0.53, 0.063, 0.05, 1.54, 0.13, 0.055, 0.49, 4.7e-4, 1.7e-4, 1.4e-4]}
    rows = [R.Ranked(peptide=p, allele="HLA-A*02:01", presentation=2.3, binder=2.1,
                     occupancy=0.77, d_occupancy=0.12, wt_absent=0.0, expression=3.0)
            for p in ("GILGFVFTL", "NLVPMVATV", "SIINFEKLA")]
    for r in rows:
        r.components.update({c: 1e-3 for c in R.CHANNEL_COLUMNS})

    # `_AGG` is a cache keyed on `(cls, species, mode)` since the mode split, so the injection
    # point is a slot rather than the whole global. What the test asserts is unchanged: one
    # library, two artifact generations, no branch.
    saved = R._AGG
    try:
        R._AGG = {("mhc1", "human", "neoantigen"): v4}
        done = R._finish(rows, None, score="aggregate")
        assert all(np.isfinite(r.score) for r in done)
        assert all(r.imputed == "" for r in done), "a v4 term fell back to its training mean"
        assert {r.components["model"] for r in done} == {"EPIC"}
    finally:
        R._AGG = saved
