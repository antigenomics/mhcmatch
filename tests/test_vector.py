"""Cassette assembly: unit construction, the per-allotype stopping rule, and junction layout."""
import math

import pytest

from mhcmatch import vector


# --------------------------------------------------------------------------- units

def test_unit_centres_the_mutation():
    ctx = "".join("ABCDEFGHIJ"[i % 10] for i in range(60))
    u = vector.unit(ctx, 30, length=27, gene="G", allele="HLA-A*02:01", p=0.2)
    assert len(u.peptide) == 27
    assert u.mutation_index == 13                      # position 14, 1-based
    assert u.peptide[u.mutation_index] == ctx[30]


def test_unit_carries_every_register_containing_the_mutation():
    """The reason for 27/position-14: all 8..14-mers spanning the mutation fit inside the window."""
    ctx = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(80))
    u = vector.unit(ctx, 40, length=27, gene="G", allele="A", p=0.1)
    for L in range(8, 15):
        # every register of length L containing the mutation must lie wholly inside the unit
        for off in range(L):
            start = u.mutation_index - off
            assert 0 <= start and start + L <= len(u.peptide)


def test_unit_clamps_at_the_terminus_and_records_where_the_mutation_landed():
    ctx = "MKVLAAGIVGL"                                  # 11 residues, shorter than the window
    u = vector.unit(ctx, 1, length=27, gene="G", allele="A", p=0.1)
    assert u.peptide == ctx                              # clamped, cannot pad
    assert u.mutation_index == 1
    assert u.peptide[u.mutation_index] == ctx[1]

    ctx2 = "".join("ACDEFGHIKL"[i % 10] for i in range(40))
    u2 = vector.unit(ctx2, 38, length=27, gene="G", allele="A", p=0.1)
    assert len(u2.peptide) == 27
    assert u2.peptide[u2.mutation_index] == ctx2[38]     # off-centre but still correct


def test_unit_rejects_an_offset_outside_the_context():
    with pytest.raises(ValueError):
        vector.unit("MKVL", 9, gene="G", allele="A", p=0.1)


def test_unit_validates_its_probability():
    with pytest.raises(ValueError):
        vector.Unit("MKVLAAGIV", 3, "G", "A", p=1.4)


# ------------------------------------------------------------------- junction windows

def test_junction_windows_span_the_boundary_only():
    left, right = "AAAAAAAAAA", "CCCCCCCCCC"
    wins = vector.junction_windows(left, right, None, lengths=(9,))
    peps = [w for w, _ in wins]
    assert peps, "a 10+10 join must produce 9-mer windows"
    for p in peps:
        assert "A" in p and "C" in p, f"{p} does not span the boundary"
    # a 9-mer spanning a 10|10 join starts at 2..9 -> 8 windows
    assert len(peps) == 8


def test_junction_windows_count_spacer_only_windows():
    """A window of spacer plus one side existed in no genome, so it is a junction window."""
    wins = vector.junction_windows("AAAAAAAAAA", "CCCCCCCCCC", "GPGPG", lengths=(9,))
    peps = [w for w, _ in wins]
    assert any("A" in p and "C" not in p and "G" in p for p in peps)
    assert all(("G" in p) or ("A" in p and "C" in p) for p in peps)


def test_junction_windows_offsets_point_at_the_right_residues():
    left, right = "AAAAAAAAAA", "CCCCCCCCCC"
    joined = left + right
    for pep, off in vector.junction_windows(left, right, None, lengths=(9, 10)):
        assert joined[off:off + len(pep)] == pep


# --------------------------------------------------------- the per-allotype stopping rule

def _u(gene, allele, p):
    return vector.Unit("A" * 27, 13, gene, allele, p)


def test_stopping_rule_matches_the_hand_computation():
    """n0=4, p=[.5,.3,.2,.1]: thresholds 0, .100, .133, .143 -> the .1 candidate is dropped."""
    cands = [_u("g1", "A", 0.5), _u("g2", "A", 0.3), _u("g3", "A", 0.2), _u("g4", "A", 0.1)]
    sel = vector.select(cands, n0=4)
    assert [u.gene for u in sel.units] == ["g1", "g2", "g3"]
    assert [u.gene for u in sel.dropped] == ["g4"]
    thr = [round(t["threshold"], 4) for t in sel.trace]
    assert thr == [0.0, 0.1, 0.1333, 0.1429]


def test_the_rule_diversifies_across_allotypes():
    """The claim the module exists to make: a *weaker* candidate on an empty allotype beats a
    stronger one on a crowded allotype, with no quota imposed anywhere."""
    cands = [_u("a1", "A", 0.50), _u("a2", "A", 0.30), _u("a3", "A", 0.20),
             _u("a4", "A", 0.12),                       # stronger, dropped
             _u("b1", "B", 0.11)]                       # weaker, kept
    sel = vector.select(cands, n0=4)
    kept = {u.gene for u in sel.units}
    assert "b1" in kept, "first unit on an empty allotype must always be taken"
    assert "a4" not in kept, "crowded allotype must saturate"
    dropped_p = sel.dropped[0].p
    kept_p = min(u.p for u in sel.units)
    assert kept_p < dropped_p, "the kept candidate really is the weaker one"


def test_first_unit_on_an_allotype_is_always_taken():
    sel = vector.select([_u("g", "A", 1e-9)], n0=0.5)
    assert len(sel.units) == 1


def test_larger_n0_admits_more_units():
    cands = [_u(f"g{i}", "A", p) for i, p in enumerate([0.5, 0.3, 0.2, 0.12, 0.08])]
    small = vector.select(cands, n0=2)
    large = vector.select(cands, n0=20)
    assert len(large.units) >= len(small.units)


def test_expected_yield_and_per_allele_agree():
    cands = [_u("a1", "A", 0.4), _u("a2", "A", 0.3), _u("b1", "B", 0.5)]
    sel = vector.select(cands, n0=6)
    assert sel.expected_yield == pytest.approx(sum(v[2] for v in sel.per_allele().values()))
    assert sel.expected_yield < sum(u.p for u in sel.units), "saturation must discount the sum"


def test_select_requires_a_stated_capacity():
    with pytest.raises(ValueError):
        vector.select([_u("g", "A", 0.5)], n0=0)


def test_select_filters_by_class():
    c1 = vector.Unit("A" * 27, 13, "g1", "A", 0.5, cls="mhc1")
    c2 = vector.Unit("C" * 27, 13, "g2", "A", 0.5, cls="mhc2")
    assert [u.gene for u in vector.select([c1, c2], n0=4, cls="mhc2").units] == ["g2"]


# ------------------------------------------------------------------------- layout

def _clean(peps, alleles=None):
    return [0.0] * len(peps)


def test_order_prefers_no_spacer_when_junctions_are_clean():
    units = [_u("g1", "A", 0.3), vector.Unit("C" * 27, 13, "g2", "A", 0.2)]
    cas = vector.order(units, binder=_clean)
    assert cas.spacer is None
    assert cas.sequence == units[0].peptide + units[1].peptide
    assert cas.boundaries == [(0, 27), (27, 54)]


def test_order_reconstructs_its_own_sequence_from_boundaries():
    units = [vector.Unit(c * 27, 13, f"g{i}", "A", 0.3 - 0.05 * i)
             for i, c in enumerate("ACDE")]
    cas = vector.order(units, binder=_clean, spacers=("GPGPG",))
    for u, (a, b) in zip(cas.units, cas.boundaries):
        assert cas.sequence[a:b] == u.peptide
    assert cas.sequence.count("GPGPG") == len(units) - 1


def test_order_picks_the_spacer_that_kills_the_junctional_binder():
    """A binder that only fires on the 'AY' seam must drive the layout away from AAY."""
    def binder(peps, alleles=None):
        return [1.0 if "AY" in p else 0.0 for p in peps]

    units = [vector.Unit("K" * 27, 13, "g1", "A", 0.3),
             vector.Unit("R" * 27, 13, "g2", "A", 0.2)]
    cas = vector.order(units, binder=binder, spacers=("AAY", "GPGPG"))
    assert cas.spacer == "GPGPG"
    assert cas.worst_junction == 0.0


def test_order_lowers_junction_cost_against_a_deliberately_bad_ordering():
    """Two units that clash only with each other must be separated by the ordering."""
    def binder(peps, alleles=None):
        return [1.0 if ("W" in p and "Y" in p) else 0.0 for p in peps]

    units = [vector.Unit("W" * 27, 13, "w", "A", 0.4),
             vector.Unit("Y" * 27, 13, "y", "A", 0.3),
             vector.Unit("K" * 27, 13, "k", "A", 0.2)]
    cas = vector.order(units, binder=binder, spacers=(None,))
    naive = vector.scan_junctions(units, binder, None)
    assert cas.cost <= sum(j["score"] for j in naive)
    assert cas.worst_junction == 0.0, "W and Y must not end up adjacent"
    genes = [u.gene for u in cas.units]
    assert abs(genes.index("w") - genes.index("y")) == 2


def test_order_threshold_short_circuits_on_the_first_acceptable_spacer():
    calls = []

    def binder(peps, alleles=None):
        calls.append(len(peps))
        return [0.0] * len(peps)

    units = [vector.Unit(c * 27, 13, f"g{i}", "A", 0.3) for i, c in enumerate("AC")]
    vector.order(units, binder=binder, spacers=(None, "GPGPG", "AAY"), threshold=0.5)
    n_none_only = len(calls)
    calls.clear()
    vector.order(units, binder=binder, spacers=(None, "GPGPG", "AAY"))
    assert n_none_only < len(calls), "threshold must stop after the first acceptable spacer"


def test_order_handles_degenerate_payloads():
    assert vector.order([], binder=_clean).sequence == ""
    one = [_u("g", "A", 0.3)]
    cas = vector.order(one, binder=_clean)
    assert cas.sequence == one[0].peptide and cas.junctions == []


def test_scan_junctions_reports_the_strongest_window_not_the_mean():
    def binder(peps, alleles=None):
        return [5.0 if p.startswith("AAC") else 0.0 for p in peps]

    units = [vector.Unit("A" * 10, 5, "g1", "A", 0.3), vector.Unit("C" * 10, 5, "g2", "A", 0.2)]
    js = vector.scan_junctions(units, binder, None, lengths=(9,))
    assert len(js) == 1
    assert js[0]["score"] == 5.0
    assert js[0]["peptide"].startswith("AAC")
    assert units[0].peptide + units[1].peptide != ""     # sanity


def test_scan_junctions_rejects_a_binder_of_the_wrong_arity():
    units = [vector.Unit("A" * 10, 5, "g1", "A", 0.3), vector.Unit("C" * 10, 5, "g2", "A", 0.2)]
    with pytest.raises(ValueError):
        vector.scan_junctions(units, lambda peps, alleles=None: [0.0], None)


def test_from_sequence_round_trips_a_single_spacer_grammar():
    units = [vector.Unit(c * 12, 6, f"g{i}", "A", 0.3) for i, c in enumerate("ACD")]
    cas = vector.order(units, binder=_clean, spacers=("GPGPG",))
    back = vector.from_sequence(cas.sequence, "GPGPG")
    assert [u.peptide for u in back] == [u.peptide for u in cas.units]


def test_rebuild_keeps_the_payload_and_relays_it():
    def binder(peps, alleles=None):
        return [1.0 if "AY" in p else 0.0 for p in peps]

    units = [vector.Unit("K" * 27, 13, "g1", "A", 0.3), vector.Unit("R" * 27, 13, "g2", "A", 0.2)]
    bad = vector.order(units, binder=binder, spacers=("AAY",))
    good = vector.rebuild(bad, binder=binder, spacers=("GPGPG",))
    assert {u.gene for u in good.units} == {u.gene for u in bad.units}
    assert good.worst_junction < bad.worst_junction


# ------------------------------------------------- m1-pseudouridine frameshifting

def _translate(cds):
    """Minimal codon table, enough to assert deslip is synonymous."""
    import itertools
    bases = "TCAG"
    aas = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
    tab = {"".join(c): a for c, a in zip(itertools.product(bases, repeat=3), aas)}
    return "".join(tab[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))


def test_slippery_sites_finds_the_published_motif_and_nothing_else():
    """TTT followed by a codon starting T or C (Mulroney 2024, PMID 38057663)."""
    #        0:AAA 1:TTT 2:TGG 3:TTT 4:CAA 5:TTT 6:AAA 7:GGG
    cds = "AAA" "TTT" "TGG" "TTT" "CAA" "TTT" "AAA" "GGG"
    hits = vector.slippery_sites(cds)
    assert [h["codon_index"] for h in hits] == [1, 3], "TTT+A must not count"
    assert [h["next_codon"][0] for h in hits] == ["T", "C"]
    assert hits[0]["nt_offset"] == 3


def test_slippery_sites_accepts_rna_and_is_case_insensitive():
    assert vector.slippery_sites("aaauuuugg") == vector.slippery_sites("AAATTTTGG")


def test_slippery_sites_ignores_a_ttt_run_that_is_not_codon_aligned():
    """Only the codon-aligned motif was characterised, so only it is reported."""
    assert vector.slippery_sites("ATT" "TTG" "GGG") == []


def test_deslip_removes_every_site_synonymously_and_is_idempotent():
    cds = "AAA" "TTT" "TGG" "TTT" "CAA" "TTT" "AAA"
    fixed, n = vector.deslip(cds)
    assert n == 2
    assert vector.slippery_sites(fixed) == []
    assert _translate(fixed) == _translate(cds), "the fix must not change the protein"
    again, n2 = vector.deslip(fixed)
    assert n2 == 0 and again == fixed, "idempotent"


def test_deslip_leaves_a_clean_sequence_untouched():
    cds = "AAAGGGCCC"
    assert vector.deslip(cds) == (cds, 0)


def test_a_gly_ser_linker_can_manufacture_a_slippery_site():
    """The concrete reason this belongs in a cassette module: linker codons put U-runs at seams."""
    # ...Phe | Gly-Ser linker starting TCN -> TTT followed by TCT
    cds = "GCA" "TTT" "TCT" "GGA" "AGC"
    assert len(vector.slippery_sites(cds)) == 1
    fixed, n = vector.deslip(cds)
    assert n == 1 and _translate(fixed) == _translate(cds)


def test_store_binder_reads_a_field_restriction_actually_has():
    """`store_binder` reaches into `Restriction`, so a rename there must fail here and not in a
    four-minute analysis run. (It did: the first version read `percent_rank`, which does not exist.)"""
    import dataclasses

    from mhcmatch.store import Restriction
    fields = {f.name for f in dataclasses.fields(Restriction)}
    assert "rank" in fields, f"store_binder reads .rank; Restriction has {sorted(fields)}"


def test_store_binder_scales_rank_so_stronger_binders_score_higher():
    class _R:
        def __init__(self, rank): self.rank = rank

    class _Store:
        def restriction(self, p, cls=None, alleles=None, calibrated=None):
            return {"S": [_R(0.1)], "W": [_R(5.0)], "N": [], "U": [_R(None)]}[p[0]]

    b = vector.store_binder(_Store(), ["HLA-A*02:01"])
    strong, weak, none_, uncal = b(["Sxx", "Wxx", "Nxx", "Uxx"])
    assert strong > weak > none_, "higher must mean a stronger predicted binder"
    assert none_ == uncal == pytest.approx(-math.log10(100.0)), "no usable rank -> non-binder"


def test_cassette_carries_no_backbone():
    """The cassette is epitopes only -- no ATG, no stop -- so it can be cloned into a backbone
    that already supplies them."""
    units = [vector.unit("M" + "ACDEFGHIKL" * 6, 20, gene="g", allele="A", p=0.3)]
    cas = vector.order(units, binder=_clean)
    assert not cas.sequence.startswith("M") or len(cas.sequence) == len(units[0].peptide)
    assert "*" not in cas.sequence

def test_order_objectives_can_disagree_and_rate_is_length_neutral():
    """The reason `objective` is explicit: summing per-junction maxima penalises a longer spacer for
    having more registers, so it drifts toward no spacer even when the *rate* of binders is worse."""
    # No spacer: every junction register is all-Z, a *moderate* binder (0.6) -- low maximum, but a
    # 100% rate. With the spacer: only the registers covering its "W" bind, and they bind *strongly*
    # (1.0) -- higher maximum, but a much lower rate. So sum-of-maxima prefers no spacer while the
    # rate prefers the spacer, on the same payload.
    def binder(peps, alleles=None):
        return [1.0 if "W" in p else (0.6 if set(p) == {"Z"} else 0.0) for p in peps]

    units = [vector.Unit("Z" * 12, 6, "g1", "A", 0.3), vector.Unit("Z" * 12, 6, "g2", "A", 0.2)]
    by_sum = vector.order(units, binder=binder, spacers=(None, "WQQQQ"), objective="sum")
    by_rate = vector.order(units, binder=binder, spacers=(None, "WQQQQ"),
                           objective="rate", binder_threshold=0.5)
    assert by_sum.spacer is None, "sum-of-maxima takes the lower peak"
    assert by_rate.spacer == "WQQQQ", "rate takes the lower share of binding registers"
    assert by_sum.spacer != by_rate.spacer, "the objectives genuinely disagree here"


def test_order_rate_objective_requires_its_threshold():
    units = [vector.Unit(c * 12, 6, f"g{i}", "A", 0.3) for i, c in enumerate("AC")]
    with pytest.raises(ValueError, match="binder_threshold"):
        vector.order(units, binder=_clean, objective="rate")
    with pytest.raises(ValueError, match="objective"):
        vector.order(units, binder=_clean, objective="median")


def test_scan_junctions_counts_binders_only_when_asked():
    def binder(peps, alleles=None):
        return [1.0 if p.startswith("AAC") else 0.0 for p in peps]

    units = [vector.Unit("A" * 10, 5, "g1", "A", 0.3), vector.Unit("C" * 10, 5, "g2", "A", 0.2)]
    assert vector.scan_junctions(units, binder, None, lengths=(9,))[0]["n_over"] is None
    j = vector.scan_junctions(units, binder, None, lengths=(9,), binder_threshold=0.5)[0]
    assert 0 < j["n_over"] <= j["n_windows"]


# ------------------------------------------------- essential-tissue exclusion
#: The MAGE-A3 epitope an affinity-enhanced HLA-A*01:01 TCR was raised against, and the titin peptide
#: it actually killed two patients through (PMID 23770775, PMID 23926201).
MAGEA3 = "EVDPIGHLY"
TITIN = "ESDPIVAQY"


def test_the_titin_pair_is_anchor_matched_and_tcr_divergent():
    """Why no distance threshold reaches titin from MAGE-A3, and the screen does not claim to."""
    same = [i for i, (a, b) in enumerate(zip(MAGEA3, TITIN), 1) if a == b]
    assert same == [1, 3, 4, 5, 9], "P3 and P9 are the A*01:01 anchors, and both are conserved"
    assert len(MAGEA3) - len(same) == 4, "four TCR-facing substitutions -- no radius reaches it"


def _titin_risk(unit, registers):
    """A `risk` callable standing in for `self_origin_risk`: only the titin register is dangerous."""
    return [{"clause": "unrelated self origin", "register": p, "protein": "sp|Q8WZ42|TITIN_HUMAN",
             "gene": "TTN", "subs": 0, "position": 24336,
             "tissue": "Heart - Left Ventricle", "tpm": 64.41}
            for p in registers if p == TITIN]


def test_screen_drops_a_unit_whose_buried_register_matches_an_essential_tissue_protein():
    """The dangerous register need not be the one carrying the mutation."""
    danger = vector.unit("K" * 9 + TITIN + "K" * 9, 4, gene="MAGEA3", allele="HLA-A*01:01", p=0.9)
    safe = vector.unit("W" * 27, 13, gene="SAFE", allele="HLA-A*01:01", p=0.4)
    assert TITIN in danger.peptide and danger.peptide[danger.mutation_index] == "K"

    kept, rejected = vector.screen([danger, safe], _titin_risk)
    assert [u.gene for u in kept] == ["SAFE"]
    assert len(rejected) == 1
    unit_, register, reason = rejected[0]
    assert unit_.gene == "MAGEA3" and register == TITIN
    assert reason["gene"] == "TTN" and reason["tissue"].startswith("Heart")
    assert reason["protein"].endswith("TITIN_HUMAN"), "the record has to name what withdrew it"


def test_screening_first_frees_capacity_for_the_next_safe_candidate():
    """Exclusion before selection, not after: the point of the ordering."""
    danger = vector.unit("K" * 9 + TITIN + "K" * 9, 4, gene="MAGEA3", allele="A", p=0.9)
    others = [vector.unit(c * 27, 13, gene=f"G{i}", allele="A", p=0.5 - 0.1 * i)
              for i, c in enumerate("WMF")]

    unscreened = vector.select([danger] + others, n0=2.0)
    kept, rejected = vector.screen([danger] + others, _titin_risk)
    screened = vector.select(kept, n0=2.0)

    assert "MAGEA3" in [u.gene for u in unscreened.units], "unscreened, the lethal one ranks first"
    assert "MAGEA3" not in [u.gene for u in screened.units]
    assert len(screened.units) > len(unscreened.units) - 1, "its slot goes to the next safe unit"
    assert len(rejected) == 1


def test_a_unit_level_reason_records_no_register():
    """The MAGE-A12 shape: the target gene itself is the hazard, no register search needed."""
    u = vector.unit("W" * 27, 13, gene="MAGEA12", allele="A", p=0.8)
    kept, rejected = vector.screen(
        [u], lambda unit, regs: [{"clause": "target gene", "gene": unit.gene,
                                  "tissue": "Brain - Caudate (basal ganglia)", "tpm": 0.33}])
    assert kept == []
    assert rejected[0][1] is None, "unit-level reasons name no register"
    assert rejected[0][2]["clause"] == "target gene"


def test_self_origin_risk_refuses_to_run_without_the_accession_map():
    """No map means nothing resolves, which reads as `no risk` -- refuse instead of reporting safe."""
    with pytest.raises(ValueError, match="symbols"):
        vector.self_origin_risk(proteome=None, symbols=None)
    with pytest.raises(ValueError, match="symbols"):
        vector.self_origin_risk(proteome=None, symbols={})


class _Hit:
    def __init__(self, protein, n_subs=0, position=24336):
        self.protein, self.n_subs, self.position = protein, n_subs, position


class _Proteome:
    """Stands in for `Proteome.find_sources` -- the batch form, which is what `risk` must call."""

    def __init__(self, per_peptide):
        self.per_peptide, self.calls = per_peptide, 0

    def find_sources(self, peptides, max_subs=1, **kw):
        self.calls += 1
        return {p: [h for h in self.per_peptide(p) if h.n_subs <= max_subs] for p in peptides}


def _stub_expression(monkeypatch):
    """TTN's real GTEx profile; everything else silent."""
    from mhcmatch import expression as EX
    monkeypatch.setattr(EX, "safety_profile",
                        lambda g, **kw: [("Muscle - Skeletal", 351.35),
                                         ("Heart - Left Ventricle", 64.41),
                                         ("Testis", 1.51)] if g == "TTN" else [])


def test_self_origin_risk_joins_an_unrelated_source_hit_to_its_tissue(monkeypatch):
    """The whole join, with the proteome and the expression table stubbed at their real values."""
    _stub_expression(monkeypatch)
    p = _Proteome(lambda pep: [_Hit("sp|Q8WZ42|TITIN_HUMAN")] if pep == TITIN else [])

    risk = vector.self_origin_risk(p, {"Q8WZ42": "TTN"})
    u = vector.Unit("K" * 9 + TITIN + "K" * 9, 4, "MAGEA3", "HLA-A*01:01", 0.9)
    hits = risk(u, [MAGEA3, TITIN])
    assert p.calls == 1, "one batch query for the whole register list, not one per register"
    assert {h["tissue"] for h in hits} == {"Muscle - Skeletal", "Heart - Left Ventricle"}
    assert "Testis" not in {h["tissue"] for h in hits}, "1.51 TPM in testis is not the hazard"
    assert all(h["clause"] == "unrelated self origin" and h["gene"] == "TTN" for h in hits)
    assert all(h["register"] == TITIN for h in hits), "MAGE-A3 resolves to nothing here"


def test_a_units_own_parent_protein_is_native_context_not_a_hazard(monkeypatch):
    """Without this exclusion the screen rejects every unit of every cassette.

    A 27-mer is native sequence by design: the flanks *are* self peptides from the parent, and the
    mutated register is one substitution from the parent's wild type. Both resolve to the unit's own
    gene, and tolerance already covers them.
    """
    _stub_expression(monkeypatch)
    # every register traces back to one protein, the unit's own parent
    parent = _Proteome(lambda pep: [_Hit("sp|P00000|KRAS_HUMAN", n_subs=0 if pep != TITIN else 1)])

    risk = vector.self_origin_risk(parent, {"P00000": "KRAS", "Q8WZ42": "TTN"})
    own = vector.Unit("K" * 9 + TITIN + "K" * 9, 4, "KRAS", "HLA-A*01:01", 0.9)
    assert risk(own, [MAGEA3, TITIN]) == [], "its own parent is native context, at 0 subs and at 1"

    # the same registers under a unit from a different gene: now the match is unrelated, and KRAS
    # being silent is what keeps this test about the parent exclusion rather than about tissue.
    titin = _Proteome(lambda pep: [_Hit("sp|Q8WZ42|TITIN_HUMAN")] if pep == TITIN else [])
    unrelated = vector.self_origin_risk(titin, {"Q8WZ42": "TTN"})
    other = vector.Unit("K" * 9 + TITIN + "K" * 9, 4, "MAGEA3", "HLA-A*01:01", 0.9)
    assert unrelated(other, [MAGEA3, TITIN]), "a MAGE-A3 unit matching titin still fires"
    assert vector.screen([own], risk)[0] == [own]


def test_max_subs_is_pushed_into_the_search_not_filtered_after(monkeypatch):
    """Searching at the library default and discarding is the same answer for twice the work."""
    _stub_expression(monkeypatch)
    seen = {}

    class _P(_Proteome):
        def find_sources(self, peptides, max_subs=1, **kw):
            seen["max_subs"] = max_subs
            return super().find_sources(peptides, max_subs, **kw)

    p = _P(lambda pep: [_Hit("sp|Q8WZ42|TITIN_HUMAN", n_subs=1)])
    risk = vector.self_origin_risk(p, {"Q8WZ42": "TTN"}, max_subs=0)
    hits = risk(vector.Unit("W" * 27, 13, "G", "A", 0.5), [TITIN])
    assert seen["max_subs"] == 0
    assert hits == [], "a 1-substitution hit must not survive max_subs=0"


def test_the_target_gene_clause_fires_without_any_register_match(monkeypatch):
    """MAGE-A12's shape: nothing in the cassette resembles anything, the gene itself is the hazard."""
    from mhcmatch import expression as EX
    monkeypatch.setattr(EX, "safety_profile",
                        lambda g, **kw: [("Brain - Caudate (basal ganglia)", 0.33)]
                        if g == "MAGEA12" else [])

    risk = vector.self_origin_risk(_Proteome(lambda pep: []), {"x": "y"})
    hits = risk(vector.Unit("W" * 27, 13, "MAGEA12", "A", 0.8), ["WWWWWWWWW"])
    assert [h["clause"] for h in hits] == ["target gene"]
    assert hits[0]["tpm"] == 0.33 and "register" not in hits[0]
    assert risk(vector.Unit("W" * 27, 13, "MAGEA3", "A", 0.8), ["WWWWWWWWW"]) == []


def test_the_default_radius_is_exact_coincidence():
    """Per-unit, not per-register: a 27-mer carries ~70 registers and any one withdraws it.

    Measured (`bench/results/vector_screen_radius.md`): at `max_subs=1` over 8-11mers, 3 of 6
    random 27-mers are withdrawn by chance -- an 8-mer plus its 152 neighbours is ~153 of 20^8
    against 68 M proteome windows. Radius 0 is clean at every length and still catches titin.
    """
    import inspect

    r = inspect.signature(vector.self_origin_risk).parameters["max_subs"].default
    assert r == 0, "raising this needs 8-mers dropped from `lengths` in the same change"


def test_the_default_floor_catches_the_lower_of_the_two_fatal_precedents():
    """A screen tuned to the cardiac case alone would have passed the neurological one.

    GTEx medians, measured: titin 64.41 TPM in heart left ventricle, MAGE-A12 **0.33** in brain
    caudate. The conventional 5-TPM 'is it expressed' cut sits between them, so the floor is set
    under the lower one. MAGE-A3's own 0.00 elsewhere is what keeps 0.25 from meaning 'everything'.
    """
    import inspect

    floor = inspect.signature(vector.self_origin_risk).parameters["min_tpm"].default
    assert floor <= 0.31, "must catch MAGE-A12 in brain putamen (PMID 23377668)"
    assert floor > 0.0, "0.00 is MAGE-A3's own median outside testis -- silent must stay silent"


def test_essential_tissues_are_real_gtex_names_and_cover_both_fatal_precedents():
    """Prefixes must match the GTEx SMTSD vocabulary, or the screen silently matches nothing."""
    gtex = ["Heart - Left Ventricle", "Heart - Atrial Appendage", "Brain - Cortex", "Lung",
            "Liver", "Kidney - Cortex", "Muscle - Skeletal", "Testis", "Whole Blood", "Ovary"]
    hit = [t for t in gtex if t.startswith(vector.ESSENTIAL_TISSUES)]
    assert "Heart - Left Ventricle" in hit, "the titin/cardiac precedent (PMID 23770775)"
    assert "Brain - Cortex" in hit, "the MAGE-A12/brain precedent (PMID 23377668)"
    assert "Testis" not in hit, "cancer-testis antigens are the target class, not the hazard"


# ------------------------------------------------------------------- back-translation

AA20 = "ACDEFGHIKLMNPQRSTVWY"


def test_back_translation_is_synonymous():
    """The whole point: the CDS must encode the peptide it was built from, for every residue."""
    import random

    rng = random.Random(0)
    for _ in range(200):
        pep = "".join(rng.choice(AA20) for _ in range(rng.randint(8, 60)))
        assert vector.translate(vector.back_translate(pep)) == pep


def test_back_translation_leaves_no_slippery_site():
    """A concatemer that frameshifts translates a second, unscreened antigen payload (PMID 38057663).

    The default table cannot emit ``TTT`` at all -- ``TTC`` is the more used phenylalanine codon --
    so this holds by construction and the deslip pass is there for a caller's own table. Asserted
    rather than assumed, because a future table edit would break it silently.
    """
    import random

    rng = random.Random(1)
    for _ in range(200):
        pep = "".join(rng.choice(AA20) for _ in range(rng.randint(8, 60)))
        assert vector.slippery_sites(vector.back_translate(pep)) == []


def test_back_translation_shortens_homopolymers_it_cannot_always_remove():
    """`max_run` is a target, not a bound, and the test says which is which.

    Poly-proline is the case that pins it: all four proline codons begin ``CC``, so no synonymous
    choice brings consecutive prolines below a 5-run. What the backoff must do is beat the
    most-frequent-codon baseline, which reaches 13 on the same peptides.
    """
    import random

    syn = vector._synonyms(vector.CODON_USAGE_HUMAN)
    rng = random.Random(2)
    peps = ["".join(rng.choice(AA20) for _ in range(40)) for _ in range(300)]
    ours = max(vector._longest_run(vector.back_translate(p)) for p in peps)
    naive = max(vector._longest_run("".join(syn[a][0] for a in p)) for p in peps)
    assert ours < naive, "the backoff must buy something"
    assert ours <= 6
    assert vector._longest_run(vector.back_translate("P" * 12)) == 5


def test_codon_table_is_complete_and_is_the_standard_code():
    """A missing or mistyped codon is a silently wrong protein, so the table is checked as a code."""
    assert len(vector.CODON_USAGE_HUMAN) == 64
    aas = {aa for aa, _ in vector.CODON_USAGE_HUMAN.values()}
    assert aas == set(AA20) | {"*"}
    # the three stops and the two single-codon residues, which pin the reading of the table
    assert {c for c, (aa, _) in vector.CODON_USAGE_HUMAN.items() if aa == "*"} == {"TAA", "TAG", "TGA"}
    assert [c for c, (aa, _) in vector.CODON_USAGE_HUMAN.items() if aa == "M"] == ["ATG"]
    assert [c for c, (aa, _) in vector.CODON_USAGE_HUMAN.items() if aa == "W"] == ["TGG"]


def test_back_translation_rejects_a_residue_it_has_no_codon_for():
    with pytest.raises(ValueError):
        vector.back_translate("SIINFEKLX")


# ------------------------------------------------------------- the rank -> unit join

def _record(gene, ctx, mut_at, tpm="1.0"):
    marked = ctx[:mut_at] + "(" + ctx[mut_at] + ")" + ctx[mut_at + 1:]
    header = (f"Somatic:chr1:100:C:T:missense_variant:{marked}:{marked}:{tpm}:"
              f"ENSG0:ENST0:{gene}:P0:0.9:1:1")
    return header, ctx[mut_at - 4:mut_at + 5]


def test_units_from_context_recentres_minimal_epitopes_on_their_mutation():
    ctx = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(55))
    header, epitope = _record("GENE1", ctx, 27)
    units = vector.units_from_context(
        [{"peptide": epitope, "gene": "GENE1", "allele": "HLA-A*02:01", "p": "0.7"}],
        [(header, epitope)])
    assert len(units) == 1
    u = units[0]
    assert len(u.peptide) == 27
    assert u.peptide[u.mutation_index] == ctx[27]
    assert u.gene == "GENE1" and u.allele == "HLA-A*02:01" and u.p == 0.7


def test_units_from_context_collapses_registers_of_one_variant_to_one_unit():
    """Twenty registers of one mutation are one thing to put in a cassette, and `select` spends
    capacity per unit -- so the group must collapse, keeping its best-scoring row's allotype."""
    ctx = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(55))
    header, _ = _record("GENE1", ctx, 27)
    rows = [{"peptide": ctx[27 - k:27 + 9 - k], "gene": "GENE1",
             "allele": f"HLA-A*02:0{k + 1}", "p": str(0.1 * (k + 1))} for k in range(4)]
    units = vector.units_from_context(rows, [(header, "")])
    assert len(units) == 1
    assert units[0].p == pytest.approx(0.4)
    assert units[0].allele == "HLA-A*02:04"


def test_units_from_context_skips_windows_with_no_marked_mutation():
    """Fusion and CNV headers use different delimiters and have no single mutated residue."""
    units = vector.units_from_context(
        [{"peptide": "SIINFEKL", "gene": "X", "allele": "H2-Kb", "p": "0.9"}],
        [("Fusion:A--B:INFRAME:SIINFEKLSIINFEKL|X:E1--E2:G1--G2:--:0.3:1:0", "SIINFEKL")])
    assert units == []


# --------------------------------------------------------------------------- the cassette map

# The SIM2 27-mer from PMID 24690990: its HLA-A*02:01 9-mer and a class-II epitope overlap by
# design, which is what let the long peptide replace an exogenous HBVcore helper outright.
SIM2 = "NMFMFRASLDLKLIFLDSRVTEVTGYE"
ERG = "AAYQIVGLVAVQEHVLKAMKQLGLSKD"
SIM2_I, SIM2_II = "KLIFLDSRV", "LKLIFLDSRVTEVTG"


def _map_cassette(spacer_first=None):
    units = [vector.Unit(peptide=SIM2, mutation_index=13, gene="SIM2", allele="HLA-A*02:01", p=0.8),
             vector.Unit(peptide=ERG, mutation_index=13, gene="ERG", allele="HLA-B*07:02", p=0.6)]
    spacers = (spacer_first,) if spacer_first is not None else vector.SPACERS
    return vector.order(units, binder=lambda peps, alleles=None: [0.0] * len(peps), spacers=spacers)


def _r1(peps):
    return [[("HLA-A*02:01", 0.31), ("HLA-B*07:02", 1.8)] if p == SIM2_I else [] for p in peps]


def _r2(peps):
    return [[("HLA-DRB1*01:01", 0.67)] if p == SIM2_II else [] for p in peps]


def test_map_duplicates_an_epitope_per_presenting_allele():
    """A heterozygote presents the same peptide on two molecules: two rows, not one."""
    feats = vector.epitope_map(_map_cassette(), _r1, _r2)
    hits = [f for f in feats if f.kind == "epitope" and f.seq == SIM2_I]
    assert len(hits) == 2
    assert sorted(f.allele for f in hits) == ["HLA-A*02:01", "HLA-B*07:02"]
    assert {f.start for f in hits} == {(f.start for f in hits).__next__()}, "same span, two alleles"


def test_map_finds_the_overlapping_class_i_and_class_ii_pair_both_ways():
    feats = vector.epitope_map(_map_cassette(), _r1, _r2)
    e1 = [f for f in feats if f.cls == "mhc1"]
    e2 = [f for f in feats if f.cls == "mhc2"]
    assert e1 and e2
    for a in e1:
        assert {f.id for f in e2} == set(a.overlaps)
    for b in e2:
        assert {f.id for f in e1} == set(b.overlaps)
    # ...and the class-II core is resolved into cassette coordinates, 9 residues wide.
    assert all(b.core_end - b.core_start == 8 for b in e2)
    assert all(b.start <= b.core_start and b.core_end <= b.end for b in e2)


def test_map_reports_which_units_carry_their_own_class_ii_help():
    cas = _map_cassette()
    s = vector.map_summary(cas, vector.epitope_map(cas, _r1, _r2))
    assert s["n_units_with_self_help"] == 1
    assert s["units"][0]["gene"] == "SIM2" and s["units"][0]["self_help"] is True
    assert s["units"][1]["gene"] == "ERG" and s["units"][1]["self_help"] is False


def test_map_carries_the_binding_core_for_both_classes():
    """`core_start`/`core_end` were cassette coordinates with no residues beside them, and stay
    class-II-only: a class-I core drops the bulge, so it is not a contiguous span."""
    from mhcmatch.store import binding_core
    rows = vector.map_rows(vector.epitope_map(_map_cassette(), _r1, _r2))
    eps = [r for r in rows if r["kind"] == "epitope"]
    assert eps, "the fixture must produce at least one epitope for this to mean anything"
    assert "core" in vector.MAP_COLUMNS
    for r in eps:
        assert r["core"] == binding_core(r["seq"], r["cls"])[0]
        assert len(r["core"]) == 9
        if r["cls"] == "mhc2":                       # span given, and it locates the core in seq
            off = r["core_start"] - r["start"]
            assert r["seq"][off:off + 9] == r["core"]
        else:
            assert r["core_start"] is None and r["core_end"] is None


def test_map_without_a_class_ii_ranker_reports_no_help_rather_than_guessing():
    cas = _map_cassette()
    s = vector.map_summary(cas, vector.epitope_map(cas, _r1, None))
    assert s["n_mhc2"] == 0 and s["n_units_with_self_help"] == 0


def test_units_and_linkers_tile_the_cassette_exactly():
    for spacer in (None, "GPGPG"):
        cas = _map_cassette(spacer)
        feats = vector.epitope_map(cas, None, None)
        tiles = sorted((f.start, f.end, f.seq) for f in feats if f.kind in ("unit", "linker"))
        assert "".join(t[2] for t in tiles) == cas.sequence
        at = 1
        for start, end, seq in tiles:
            assert start == at and end == at + len(seq) - 1
            at = end + 1
        assert at == len(cas.sequence) + 1
        if spacer:
            assert [f.seq for f in feats if f.kind == "linker"] == [spacer]


def test_map_coordinates_are_one_based_and_slice_back_to_the_peptide():
    cas = _map_cassette("GPGPG")
    for f in vector.epitope_map(cas, _r1, _r2):
        assert cas.sequence[f.start - 1:f.end] == f.seq


def test_an_epitope_spanning_a_junction_is_flagged_as_unit_zero():
    cas = _map_cassette(None)
    lo = cas.boundaries[0][1] - 4                     # straddles the unit-1/unit-2 boundary
    junctional = cas.sequence[lo:lo + 9]
    feats = vector.epitope_map(cas, lambda peps: [[("HLA-A*02:01", 0.4)] if p == junctional else []
                                             for p in peps], None)
    hits = [f for f in feats if f.kind == "epitope"]
    assert hits and all(f.unit == 0 and f.gene == "" for f in hits)
    assert vector.map_summary(cas, feats)["n_junction_spanning"] == len(hits)


def test_write_map_round_trips_through_both_formats(tmp_path):
    import csv
    import json
    cas = _map_cassette("GPGPG")
    feats = vector.epitope_map(cas, _r1, _r2)
    tsv, js = tmp_path / "m.tsv", tmp_path / "m.json"
    summary = vector.write_map(cas, feats, str(tsv), str(js))
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert len(rows) == len(feats)
    assert list(rows[0]) == list(vector.MAP_COLUMNS)
    assert all(len(r) == len(vector.MAP_COLUMNS) for r in rows)          # one value per cell
    blob = json.loads(js.read_text())
    assert blob["sequence"] == cas.sequence and blob["summary"] == summary
    assert len(blob["features"]) == len(feats)
    ov = [r for r in rows if r["overlaps"]]
    assert ov and all("," not in r["id"] for r in rows)


# ------------------------------------------------------- `vector --quota` writes what it composes
def _quota_case(tmp_path, n):
    """`n` candidate windows through the CLI's own readers: a context FASTA and a unit table."""
    from mhcmatch.predict import parse_variant_header, variant_product
    fa, tsv = tmp_path / "ctx.fasta", tmp_path / "cand.tsv"
    subs = ("missense_variant", "frameshift_variant")
    rows = ["peptide\tgene\tallele\tp\tcls\tvariant_type"]
    with open(fa, "w") as fh:
        for i in range(n):
            left = "ACDEFGHIKLMNPQRSTVWYACDEFGH"
            right = "YWVTSRQPNMLKIHGFEDCAYWVTSRQ"
            ctx = f"{left}({'ACDEFGHIKL'[i % 10]}){right}"
            h = (f"Somatic:chr1:{100 + i}:G:C:{subs[i % 2]}:{ctx}:{ctx}:9.9:"
                 f"ENSG{i}:ENST{i}:GENE{i}:U{i}:0.9:5:9")
            seq = ctx.replace("(", "").replace(")", "")
            fh.write(f">{h}\n{seq}\n")
            rows.append(f"{seq[24:33]}\tGENE{i}\tHLA-A*02:01\t{0.30 - 0.02 * i:.4f}\tmhc1\t"
                        f"{variant_product(parse_variant_header(h))}")
    tsv.write_text("\n".join(rows) + "\n")
    return str(tsv), str(fa)


def _run_vector(tmp_path, n, extra):
    from mhcmatch import cli
    cand, ctx = _quota_case(tmp_path, n)
    faa, fna = tmp_path / "c.faa", tmp_path / "c.fna"
    cli.main(["vector", "--candidates", cand, "--context", ctx, "--n0", "8",
              "--fasta", str(faa), "--fasta-nt", str(fna),
              "--out", str(tmp_path / "r.tsv")] + extra)
    heads = [ln[1:].split()[0] for ln in faa.read_text().splitlines() if ln.startswith(">")]
    nt = [ln[1:].split()[0] for ln in fna.read_text().splitlines() if ln.startswith(">")]
    return heads, nt, (tmp_path / "r.tsv").read_text()


def test_quota_emits_the_composed_cassette_and_the_score_only_one(tmp_path):
    """Through 0.24.0 ``--quota`` composed a set and then built the sequence from ``select`` anyway.

    It reported and did not act. Both cassettes now reach the FASTA, so the caller can compare the
    layout the portfolio chose against the one a ranked list gives on the *same* slot budget.
    """
    heads, nt, report = _run_vector(tmp_path, 4, ["--quota", "mhc1=1:1,nonconventional=1:1"])
    assert heads == ["cassette_composed", "cassette_topk"]
    assert nt == ["cassette_composed_cds", "cassette_topk_cds"]
    # both arms filled, so the non-conventional quota actually bit
    assert "\tnonconventional\t" in report and "\tmhc1\t" in report


def test_without_quota_the_output_is_the_single_cassette_it_always_was(tmp_path):
    heads, nt, _ = _run_vector(tmp_path, 1, [])
    assert heads == ["cassette"] and nt == ["cassette_cds"]
