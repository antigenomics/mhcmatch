"""Runnable checks for mhcmatch v0. Mirror seqtree's positive controls where applicable."""
import random

import pytest

import mhcmatch
from mhcmatch import Pseudoseq, Store, search
from mhcmatch.pseudoseq import learn_anchor_weights

# Dolton et al. (Cell 2023) HLA-A*02:01 cross-reactive trio -- TCR-facing positive control.
TRIPLE = ["EAAGIGILTV", "LLLGIGILVL", "NLSALGIFST"]

# Synthetic class-I panel: each allele fixes its 5-residue presentation signature
# (P1,P2,P3,PΩ-1,PΩ = peptide indices 0,1,2,7,8); the TCR-facing middle is random.
_SIG = {"HLA-A*02:01": "ALACV", "HLA-B*07:02": "GPRWL", "HLA-A*01:01": "STDEY"}
_AA = "ACDEFGHIKLMNPQRSTVWY"


def _make_store(n=40, seed=0):
    rng = random.Random(seed)
    recs = []
    for allele, sig in _SIG.items():
        for _ in range(n):
            mid = "".join(rng.choice(_AA) for _ in range(4))
            pep = sig[0] + sig[1] + sig[2] + mid + sig[3] + sig[4]  # length 9
            recs.append({"epitope": pep, "mhc_a": allele, "mhc_class": "MHCI"})
    return Store.from_records(recs)


# -- restriction / presentation (forward + reverse problem) -------------------
def test_restriction_recovers_allele():
    store = _make_store()
    query = "ALA" + "EEEE" + "CV"  # A*02:01 signature, novel middle
    res = store.restriction(query, cls="mhc1")
    assert res[0].allele == "HLA-A*02:01"
    assert res[0].binder
    assert store.is_binder(query, "HLA-A*02:01")


def test_nonbinder_rejected():
    store = _make_store()
    # signature WWWWW matches no allele -> no presentation-signature neighbours
    assert not store.is_presented("WWW" + "EEEE" + "WW", cls="mhc1")


def test_allele_subset_filter():
    store = _make_store()
    query = "ALA" + "KKKK" + "CV"
    res = store.restriction(query, cls="mhc1", alleles=["HLA-B*07:02"])
    assert [r.allele for r in res] == ["HLA-B*07:02"]
    assert not res[0].binder  # B*07:02 does not present an A*02:01-signature peptide


# -- anchor / TCR-facing decomposition (X masks) -----------------------------
def test_decompose_mhc1():
    d = mhcmatch.Store().decompose("KLEEEEEEV", cls="mhc1")
    assert d.anchors == (1, 8)                  # P2 + PΩ for a 9-mer
    assert d.tcr_facing == "K" + "X" + "EEEEEE" + "X"
    assert d.presentation == "X" + "L" + "XXXXXX" + "V"


def test_decompose_mhc2():
    d = mhcmatch.Store().decompose("AAAYAAKAAVAAAAA", cls="mhc2")
    assert len(d.anchors) == 4                   # P1/P4/P6/P9 of the anchored core
    assert d.tcr_facing.count("X") == 4
    assert d.presentation.count("X") == len(d.peptide) - 4


# -- large-scale similarity search (TCR-facing) ------------------------------
def test_dolton_trio_mutual_homologs():
    for q in TRIPLE:
        hits = {m.peptide for m in search.search(q, TRIPLE, mode="tcr", cls="mhc1",
                                                 k=4, max_subs=2, min_shared=1)}
        assert hits == set(TRIPLE) - {q}, f"{q}: {hits}"


def test_find_mimics_evalue():
    res = search.find_mimics("EAAGIGILTV", ["LLLGIGILVL", "KLGGALQAK", "GILGFVFTL"],
                             bacterial_sets={"ecoli": ["NLSALGIFST", "MMMMMMMMM"]},
                             max_subs=2, min_shared=1)
    assert "LLLGIGILVL" in {h.epitope for h in res["self"]["hits"]}
    assert any(h.epitope == "NLSALGIFST" for h in res["ecoli"]["hits"])


# -- motif logos --------------------------------------------------------------
def test_logo_pwm_and_lengths():
    store = _make_store()
    m = mhcmatch.logo.motif(store, "HLA-A*02:01", "mhc1")
    assert m["width"] == 9
    for col in m["pwm"]:
        assert abs(sum(col.values()) - 1.0) < 1e-9
    assert all(0.0 <= b <= __import__("math").log2(20) + 1e-9 for b in m["bits"])
    assert m["length_hist"] == {9: 40}


# -- pseudosequence kernel + diffusion ---------------------------------------
def test_kernel_symmetry_and_neighbours():
    ps = Pseudoseq("mhc1", h=2.0)
    assert ps.kernel("HLA-A*02:01", "HLA-A*02:06") == ps.kernel("HLA-A*02:06", "HLA-A*02:01")
    # a close relative scores higher than a distant allele
    assert ps.kernel("HLA-A*02:01", "HLA-A*02:06") > ps.kernel("HLA-A*02:01", "HLA-B*07:02")


def test_shrink_limits():
    prefs = {"A": __import__("collections").Counter({"L": 8, "I": 2}),
             "B": __import__("collections").Counter({"V": 10})}
    # craft two pseudosequences differing in 1 position so the kernel is well-defined
    ps = Pseudoseq("mhc1")
    ps.seqs = {"A": "A" * 34, "B": "C" + "A" * 33}
    near = Pseudoseq("mhc1", h=1e-9)
    near.seqs = ps.seqs
    far = Pseudoseq("mhc1", h=1e9)
    far.seqs = ps.seqs
    own = near.shrink(prefs, "A")              # h->0: no borrowing
    assert abs(own["L"] - 0.8) < 1e-6 and "V" not in own
    pooled = far.shrink(prefs, "A")            # h->inf: borrows B fully
    assert "V" in pooled and pooled["V"] == pytest.approx(10 / 20)


def test_learn_anchor_weights_mi():
    # position 0 perfectly predicts the anchor residue; the rest are constant (no info)
    seqs = {f"al{i}": ("L" if i % 2 else "P") + "A" * 33 for i in range(10)}
    anchor = {f"al{i}": ("V" if i % 2 else "K") for i in range(10)}
    w = learn_anchor_weights(seqs, anchor)
    assert w[0] == max(w) and w[0] > 0
    assert all(w[p] == 0 for p in range(1, 34))


# -- diffusion: rare-allele rescue -------------------------------------------
def _build(rng, allele, p2, pomega, n, recs):
    for _ in range(n):
        mid = [rng.choice(_AA) for _ in range(7)]
        recs.append({"epitope": mid[0] + p2 + "".join(mid[1:]) + pomega,  # 9-mer, P2=idx1, PΩ=idx8
                     "mhc_a": allele, "mhc_class": "MHCI"})


def test_diffusion_rescues_rare_allele():
    rng = random.Random(1)
    recs = []
    _build(rng, "HLA-A*02:01", "L", "V", 40, recs)   # frequent, canonical A*02 anchors
    _build(rng, "HLA-B*07:02", "P", "L", 40, recs)   # frequent, distant groove
    _build(rng, "HLA-A*02:06", "L", "L", 1, recs)    # RARE: never shows PΩ=V on its own
    store = Store.from_records(recs)
    am = store.anchor_model("mhc1", h=2.0, anchors=(2, -1))  # isolate the rescue on the primary anchors
    q = "AL" + "EEEEEE" + "V"                          # P2=L, PΩ=V -- a classic A*02 peptide
    raw = am.score(q, "HLA-A*02:06", raw=True)         # off its 1 peptide: PΩ=V unseen
    diffused = am.score(q, "HLA-A*02:06")              # borrows PΩ=V from groove-neighbour A*02:01
    assert diffused > raw
    assert diffused > am.score(q, "HLA-B*07:02")       # A*02-like query prefers the A*02 groove


# -- near-exact source lookup -------------------------------------------------
def test_proteome_source_lookup():
    from mhcmatch import Proteome
    prot = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    pm = Proteome({"P1": prot})
    wild = prot[5:14]                  # a 9-mer self peptide
    mutant = wild[:4] + ("W" if wild[4] != "W" else "Y") + wild[5:]
    hits = pm.find_source(mutant, max_subs=1)
    top = next(h for h in hits if h.protein == "P1" and h.position == 5)
    assert top.n_subs == 1
    assert top.mutations[0][0] == 4   # the mutated position within the peptide
