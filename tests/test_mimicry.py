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
        body = [ln.split("\t") for ln in capsys.readouterr().out.strip().splitlines()]
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

    GRAND takes its corpus term from the thymic channel alone (26,513 peptides), so the proteome
    index is off the ranking path. This pins that: the model must declare no `self` feature, and
    the channels the CLI builds must be buildable with `with_self=False`."""
    from mhcmatch import rank

    feats = rank.aggregate()["features"]
    assert not any("self" in f for f in feats), feats
    assert "C_corpus_thymus" in feats
    assert set(rank.CHANNEL_COLUMNS) <= set(feats)


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
    to hand back plain strings."""
    import numpy as np
    b = mimicry._Backing(np.array([b"GILGFVFTL"], dtype="S9"), np.array([b"FLU"], dtype="S3"))
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

    Only reachable through `allow_missing=True`, but reachable by default since 0.21.0: GRAND does
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
