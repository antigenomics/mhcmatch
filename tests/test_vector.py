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
