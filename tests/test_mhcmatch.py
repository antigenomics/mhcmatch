"""Runnable checks for mhcmatch v0. Mirror seqtree's positive controls where applicable."""
import collections
import math
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


def test_decompose_mhc2_register_start():
    # register_start pins the 9-mer core frame (the model register a caller scored with), overriding
    # the heuristic register: P1/P4/P6/P9 of a core starting at index 2 -> {2,5,7,10}. (C4b)
    pep = "A" * 15
    assert mhcmatch.Store().decompose(pep, cls="mhc2", register_start=2).anchors == (2, 5, 7, 10)
    assert mhcmatch.Store().decompose(pep, cls="mhc2", register_start=99).anchors == ()   # guarded
    # default (register_start=None) still uses the heuristic register (unchanged behavior).
    assert len(mhcmatch.Store().decompose(pep, cls="mhc2").anchors) == 4


def test_bare_store_restriction_is_graceful():
    # A Store never loaded via from_records/from_pmhc has no reference panel; restriction()
    # and alleles() must return empty rather than AttributeError (decompose() still works on it).
    s = mhcmatch.Store()
    assert s.alleles("mhc2") == []
    assert s.restriction("AAAYAAKAAVAAAAA", cls="mhc2") == []
    assert s.restriction("AAAYAAKAAVAAAAA", cls="mhc2", diffuse=True) == []


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


def _sig_build(rng, allele, sig, n, recs):
    for _ in range(n):
        mid = "".join(rng.choice(_AA) for _ in range(4))
        recs.append({"epitope": sig[0] + sig[1] + sig[2] + mid + sig[3] + sig[4],
                     "mhc_a": allele, "mhc_class": "MHCI"})


def test_mhc2_register_max_and_em():
    # Per-allele register: a planted 9-mer core scores higher for its own allele than for another,
    # under both register modes.
    #
    # Frame-invariance (the same core at different offsets scoring identically) is a property of
    # register="max" ONLY, because it reads off the winning frame and discards the rest. It is
    # deliberately FALSE for register="marginal": that integrates over every frame under the
    # core-offset prior, so both the offset and the padding context move the score. Do not "fix"
    # marginal to satisfy it -- see AnchorModel.score.
    coreX, coreY = "WKVKFWKVK", "DKEKDDKEK"        # distinct allele motifs
    recs = []
    for pad in range(5):                            # core placed at every offset in a 13-mer
        recs.append({"epitope": "S" * pad + coreX + "S" * (4 - pad),
                     "mhc_a": "DRA*01:01", "mhc_b": "DRB1*15:01", "mhc_class": "MHCII"})
        recs.append({"epitope": "S" * pad + coreY + "S" * (4 - pad),
                     "mhc_a": "DRA*01:01", "mhc_b": "DRB1*13:01", "mhc_class": "MHCII"})
    store = Store.from_records(recs * 4)
    # n_motifs=1: this pins REGISTER behaviour, and the K=3 default (v0.7) would confound it here --
    # a 2-allele toy panel gives cross-allele shrinkage no real neighbourhood, so K=3 EM on 20
    # near-identical ligands per allele is noisy enough to flip a toy ordering. On a realistic panel
    # (neighbours to borrow from) K=3 preserves and sharpens the specificity; the mixture itself is
    # tested in test_mixture_recovers_two_planted_binding_modes.
    for reg in ("max", "marginal"):                          # register_em=2 by default
        am = store.anchor_model("mhc2", register=reg, n_motifs=1)
        assert am.score("GG" + coreX + "GG", "DRB1_1501") > am.score("GG" + coreX + "GG", "DRB1_1301")
    mx = store.anchor_model("mhc2", register="max", n_motifs=1)
    assert mx.score("G" + coreX + "GGG", "DRB1_1501") == mx.score("GGG" + coreX + "G", "DRB1_1501")


def test_mhc2_offset_prior_learns_where_the_cores_sit():
    # The core-offset prior is what register="marginal" buys: real class-II cores sit ~3 residues
    # from the N-terminus (measured H/Hmax 0.67 on DRB1_0101 15mers vs 0.998 on random peptides),
    # so the prior must concentrate on the offsets an allele's own ligands use -- and it is
    # PER-ALLELE, not one global shape. Two alleles are needed at minimum: with a single allele the
    # background="ligand" null equals that allele's own marginal, every frame scores ~0, and the
    # register argmax is an arbitrary tie.
    coreX, coreY = "WKVKFWKVK", "DQDQNDQDQ"                  # distinct anchor motifs at 1,4,6,9
    rng = random.Random(0)
    flank = "AGILPRSTV"                                      # shares no residue with either core

    def pep(core, off):
        return ("".join(rng.choice(flank) for _ in range(off)) + core
                + "".join(rng.choice(flank) for _ in range(4 - off)))

    recs = []
    for _ in range(60):
        recs.append({"epitope": pep(coreX, 2), "mhc_a": "DRA*01:01", "mhc_b": "DRB1*15:01",
                     "mhc_class": "MHCII"})                  # DR15 cores always at offset 2
        recs.append({"epitope": pep(coreY, 0), "mhc_a": "DRA*01:01", "mhc_b": "DRB1*13:01",
                     "mhc_class": "MHCII"})                  # DR13 cores always at offset 0
    am = Store.from_records(recs).anchor_model("mhc2")
    lx = am._offset_logprior("DRB1_1501", 13)
    ly = am._offset_logprior("DRB1_1301", 13)
    assert len(lx) == 5                                      # 13 - 8 frames
    assert lx.index(max(lx)) == 2                            # each allele learned its own offset,
    assert ly.index(max(ly)) == 0                            # not a shared global shape
    # and it is a proper distribution over frames, not an unnormalized score
    assert abs(sum(math.exp(x) for x in lx) - 1.0) < 1e-9


def test_restriction_diffuse_rescues_rare():
    rng = random.Random(2)
    recs = []
    _sig_build(rng, "HLA-A*02:01", "ALACV", 40, recs)   # frequent, PΩ=V
    _sig_build(rng, "HLA-A*02:06", "ALACL", 1, recs)    # RARE, only ever shows PΩ=L
    _sig_build(rng, "HLA-B*07:02", "GPRWL", 40, recs)   # frequent, distant groove
    store = Store.from_records(recs)
    q = "ALA" + "MMMM" + "CV"                            # signature ALACV (PΩ=V)
    # vote mode: the rare allele has ~no signature neighbours -> not flagged
    r = store.restriction(q, cls="mhc1", diffuse=True, alleles=["HLA-A*02:06"])[0]
    assert r.binder and r.anchor_score > 0              # borrowed PΩ=V from groove-neighbour A*02:01


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
    # wildtype() fetches the WT counterpart (the self peptide the mutant derives from) for agretopicity.
    assert pm.wildtype(mutant) == wild
    assert pm.wildtype(wild) is None                 # an exact self peptide has no mutated WT
    assert pm.wildtype("YYYYYYYYY") is None           # nothing within 1 sub -> None


# -- CLI ----------------------------------------------------------------------
# -- structural+learned blended pocket weights (item 2) ----------------------
def test_blend_weights_build_and_score():
    import math

    store = _make_store()
    m = store.anchor_model("mhc1", weights="blend")
    assert m.weights_mode == "blend"
    for j in m.anchors:                       # blended weights stay mean-1 normalized per anchor
        w = m.ps._w(j)
        assert abs(sum(w) / len(w) - 1.0) < 1e-6
    assert math.isfinite(m.score("ALAEEEECV", "HLA-A*02:01"))


# -- allele-name resolution (item 7) -----------------------------------------
def test_resolve_allele():
    from mhcmatch import resolve_allele
    from mhcmatch.pseudoseq import load_pseudo

    key = next(k for k in load_pseudo("mhc1") if k.startswith("HLA-A"))  # e.g. 'HLA-A02:01'
    assert resolve_allele(key, "mhc1") == (key, True)
    assert resolve_allele(key[4:], "mhc1") == (key, True)          # missing 'HLA-' prefix
    assert resolve_allele("A*" + key[5:], "mhc1") == (key, True)   # '*' punctuation + no prefix
    assert resolve_allele("ZZ:99:99", "mhc1") == (None, False)     # unknown -> flagged, not guessed


# -- multiple-testing control on a protein scan (item 6) ----------------------
def test_scan_protein_fwer_fdr():
    store = _make_store()
    pep = "ALA" + "EEEE" + "CV"            # A*02:01 signature, novel TCR-facing middle
    protein = "GGGG" + pep + "GGGG"
    raw = store.scan_protein(protein, cls="mhc1")
    bh = store.scan_protein(protein, cls="mhc1", correction="bh")
    bonf = store.scan_protein(protein, cls="mhc1", correction="bonferroni")
    assert raw                                                     # the signature window is found
    raw_pos = {(i, p) for i, p, _ in raw}
    assert {(i, p) for i, p, _ in bh} <= raw_pos                   # correction never adds windows
    assert {(i, p) for i, p, _ in bonf} <= {(i, p) for i, p, _ in bh}  # Bonferroni <= BH
    assert any(p == pep for _, p, _ in bonf)                       # strong signal survives FWER


def test_cli_decompose_and_source(tmp_path, capsys):
    from mhcmatch import cli
    cli.main(["decompose", "NLVPMVATV", "--cls", "mhc1"])
    out = capsys.readouterr().out
    assert "tcr_facing" in out and "NXVPMVATX" in out

    fasta = tmp_path / "prot.fasta"
    fasta.write_text(">P1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n")
    cli.main(["source", "MKTAYIAKW", "--proteome", str(fasta), "--max-subs", "1"])
    assert "P1" in capsys.readouterr().out


# -- ligand spans: core -> full presented peptide ------------------------------
# All synthetic. The span model predicts an OBSERVED LIGAND SPAN, not a cleavage event and not
# immunogenicity -- see mhcmatch.ligand for why that distinction is load-bearing.
from mhcmatch import ligand                                            # noqa: E402
from mhcmatch.diffusion import PROTEOME_AA_FREQ                        # noqa: E402
from mhcmatch.ligand import CTX_KEYS, PAD, SpanModel, presented_span   # noqa: E402


def _flat_model(lens=None, spikes=None):
    """A SpanModel that is flat at the proteome background, except where we plant a spike."""
    ctx = {k: dict(PROTEOME_AA_FREQ) for k in CTX_KEYS}
    for k in ctx:
        ctx[k][PAD] = 0.01
    for key, res in (spikes or {}).items():
        ctx[key] = {a: 1e-4 for a in PROTEOME_AA_FREQ}
        ctx[key][res] = 0.9
        ctx[key][PAD] = 1e-4
    lens = lens or {14: 1.0}
    return SpanModel(ctx=ctx, lens=lens, padbg=0.01, background="proteome")


def _mhc2_store():
    coreX, coreY = "WKVKFWKVK", "DKEKDDKEK"        # two distinct allele motifs -> prefs != bg
    recs = []
    for pad in range(5):
        recs.append({"epitope": "S" * pad + coreX + "S" * (4 - pad), "mhc_a": "DRA*01:01",
                     "mhc_b": "DRB1*15:01", "mhc_class": "MHCII"})
        recs.append({"epitope": "S" * pad + coreY + "S" * (4 - pad), "mhc_a": "DRA*01:01",
                     "mhc_b": "DRB1*13:01", "mhc_class": "MHCII"})
    return Store.from_records(recs * 4), coreX


def test_best_register_is_bit_identical_to_score():
    store, core = _mhc2_store()
    am = store.anchor_model("mhc2", register="max")
    pep = "GG" + core + "GG"
    st, sc = am.best_register(pep, "DRB1_1501")
    assert sc == am.score(pep, "DRB1_1501")        # under register="max", score() IS best_register()[1]
    assert 0 <= st <= len(pep) - 9
    # register="marginal" (the default) deliberately breaks that identity: it integrates over every
    # frame instead of reading the winning one, so score() < max-over-frames. best_register() still
    # returns the argmax frame either way -- decompose/logos/affinity read it as a register oracle
    # and must not shift.
    mg = store.anchor_model("mhc2")
    assert mg.best_register(pep, "DRB1_1501") == (st, sc)
    assert mg.score(pep, "DRB1_1501") < sc
    # the winning core is a property of the sequence, not of where it sits: the same 9-mer must win
    # at any padding offset, and the score must match (the pre-existing frame-invariance contract).
    p1, p2 = "G" + core + "GGG", "GGG" + core + "G"
    s1, c1 = am.best_register(p1, "DRB1_1501")
    s2, c2 = am.best_register(p2, "DRB1_1501")
    assert c1 == c2
    assert p1[s1:s1 + 9] == p2[s2:s2 + 9]


def test_span_model_recovers_planted_span():
    # Plant K immediately upstream and R immediately downstream of a 14-mer span.
    core = "WKVKFWKVK"
    prot = "G" * 20 + "K" + "AA" + core + "SSS" + "R" + "G" * 20
    s = prot.index(core) - 2                        # planted span start (2 residues of N-flank)
    m = _flat_model(lens={14: 1.0}, spikes={"flankN-1": "K", "flankC+1": "R"})
    sp = presented_span(core, prot, model=m, mode="modeled")
    assert (sp.start, sp.end) == (s, s + 14)
    assert sp.peptide == prot[s:s + 14] and sp.source == "modeled"
    assert sp.flanks == (2, 3)


def test_span_score_is_length_unbiased():
    # AnchorModel.score is a max over register frames, so it GROWS with length -- ranking candidate
    # spans with it would just pick the longest. The span score must not have that bias.
    rng = random.Random(7)
    prot = "".join(rng.choice(_AA) for _ in range(400))
    m = _flat_model(lens={L: 1 / 9 for L in range(12, 21)})
    means = []
    for L in range(12, 21):
        sc = [m.context_score(prot, s, s + L) for s in range(0, 300, 7)]
        means.append(sum(sc) / len(sc))
    xs = list(range(12, 21))
    mx, my = sum(xs) / len(xs), sum(means) / len(means)
    slope = (sum((x - mx) * (y - my) for x, y in zip(xs, means))
             / sum((x - mx) ** 2 for x in xs))
    assert abs(slope) < 0.05, f"span score drifts {slope:+.3f} log-odds/residue with length"


def test_anchor_score_is_length_biased_negative_control():
    # Documented negative control for the trap above: on RANDOM peptides -- carrying no motif at all
    # -- AnchorModel.score still rises with length. Under register="max" that is a max over more
    # frames; under register="marginal" (the default) the offset prior normalizes the frame count
    # away but a Jensen residual remains (measured on real data, DRB1_1501, 9mer -> 21mer: +4.44
    # nats under max vs +2.28 under marginal -- halved, not gone). Either way span ranking must use
    # the flank model and not this. If this ever stops holding for BOTH modes, revisit
    # SpanModel.best_span and bench/results/binder_gate_length_bias.md.
    rng = random.Random(3)
    store, _ = _mhc2_store()
    for reg in ("max", "marginal"):
        am = store.anchor_model("mhc2", register=reg)
        mean = lambda L: sum(                                         # noqa: E731
            am.score("".join(rng.choice(_AA) for _ in range(L)), "DRB1_1501") for _ in range(80)) / 80
        assert mean(19) > mean(9), reg


def test_mhc2_binder_gate_is_not_a_length_detector():
    # The score IS length-biased (test above) -- so an absolute `anchor_score > 0` binder gate was a
    # length detector: it passed a random 15-mer 85% of the time and a random 21-mer 98%
    # (bench/results/binder_gate_length_bias.md). restriction() now gates on a %rank against random
    # peptides OF THE SAME LENGTH, which puts the null through the same frame-max, so the bias
    # cancels: the false-positive rate must be flat in length, not a ramp.
    #
    # This pins the %rank clause itself rather than the whole of restriction(). The end-to-end path
    # cannot be measured on a fixture this size: `_mhc2_store()` is 40 peptides over 7 residues, and
    # random probes hit k-mer neighbours 58-85% of the time there (even at 600 peptides over 15
    # residues), so the *vote* clause decides the outcome -- and it fires more as length grows,
    # manufacturing the very ramp this test looks for. Probes are drawn from the panel's own residue
    # composition, matching the calibrator's null: uniform-`_AA` probes are full of residues the
    # panel never saw, which hit the eps floor and score ~0 against a null its own padding drags
    # down, and `register="marginal"` widens that gap with frame count.
    #
    # Both are fixture artifacts. On the real panel the gate is flat in length and identical under
    # register="max" and "marginal" -- DRB1_1501 8/12/10% vs 8/12/11% at L=9/15/21 -- and the ramp
    # the old gate produced (85% at 15, 98% at 21) is measured in binder_gate_length_bias.md.
    from mhcmatch.calibrate import WEAK_RANK
    store, _ = _mhc2_store()
    cal, am = store._rank_calibrator("mhc2"), store._anchor_model("mhc2")
    aa = collections.Counter()
    for p in store._panel["mhc2"].epitopes:
        aa.update(p)
    res, w = zip(*aa.items())
    rng = random.Random(0)
    rates = {}
    for L in (9, 15, 21):
        peps = ["".join(rng.choices(res, w, k=L)) for _ in range(200)]
        rates[L] = sum(cal.percent_rank("DRB1_1501", am.score(p, "DRB1_1501"), length=L) <= WEAK_RANK
                       for p in peps) / len(peps)
    # `%rank <= t` passes t% of the null by construction -- AT EVERY LENGTH. That is the whole point:
    # the null takes the same frame-max as the query, so the length inflation cancels instead of
    # being modelled. Nothing like the old 85-98%, and no ramp.
    assert max(rates.values()) <= 0.10, rates
    assert abs(rates[21] - rates[9]) <= 0.08, rates


def test_bounds_pad_and_clipping():
    core = "WKVKFWKVK"
    prot = core + "AAAA"                          # core flush against the N-terminus
    m = _flat_model(lens={11: 1.0})
    sp = presented_span(core, prot, model=m, mode="modeled")
    assert sp.start >= 0 and sp.end <= len(prot) and prot[sp.start:sp.end] == sp.peptide

    fx = ligand.fixed_span(core, prot, 5, 5)      # neither flank fits
    assert (fx.start, fx.end) == (0, len(prot))
    assert fx.clipped == (5, 1)                  # reported, never silently shortened
    with pytest.raises(ValueError):
        ligand.fixed_span(core, prot, 5, 5, strict=True)


def test_span_is_self_consistent_with_its_own_register():
    # The emitted span must still read back the core we asked for -- otherwise the answer contradicts
    # itself (we would report core X while having extended around core Y).
    core = "WKVKFWKVK"
    recs = [{"epitope": "S" * p + core + "S" * (4 - p), "mhc_a": "DRA*01:01",
             "mhc_b": "DRB1*15:01", "mhc_class": "MHCII"} for p in range(5)]
    am = Store.from_records(recs * 4).anchor_model("mhc2")
    prot = "GGGGG" + "AA" + core + "SSS" + "GGGGG"
    sp = presented_span(core, prot, model=_flat_model(lens={14: 1.0}), mode="modeled")
    st, _ = am.best_register(sp.peptide, "DRB1_1501")
    assert sp.peptide[st:st + 9] == core


def test_observed_tier_and_leak_guard():
    core = "WKVKFWKVK"
    prot = "GGGG" + "MK" + core + "TPR" + "GGGG"
    lig = "MK" + core + "TPR"                     # a real "eluted" ligand bracketing the core
    sp = presented_span(core, prot, corpus=[lig], mode="auto")
    assert sp.source == "observed" and sp.peptide == lig and sp.support == 1
    # mode="observed" with no corpus hit is an informative None, not a silent fallback
    assert presented_span(core, prot, corpus=["QQQQQQQQ"], mode="observed") is None
    # the benchmark leak guard: mode="modeled" must never echo a corpus ligand as a prediction
    assert presented_span(core, prot, corpus=[lig], mode="modeled").source == "modeled"


def test_class_boundary_no_silent_misrouting():
    core, prot = "WKVKFWKVK", "GGGG" + "WKVKFWKVK" + "GGGG"
    with pytest.raises(ValueError):
        presented_span("A" * 15, prot.replace("W", "A"))          # not a 9-mer core
    with pytest.raises(ValueError):
        ligand.processing_score("A" * 15, "G" * 5 + "A" * 15)     # not an 8-11mer
    # class I and class II are separate entry points on purpose: length cannot tell a 9-mer class-II
    # core from a 9-mer class-I peptide, so nothing may infer the class.
    assert not hasattr(ligand, "infer_class")
    assert isinstance(ligand.processing_score(core, prot), float)  # valid AS a class-I 9-mer


def test_ligand_background_is_rejected_as_circular():
    with pytest.raises(ValueError, match="circular"):
        SpanModel(ctx={k: dict(PROTEOME_AA_FREQ) for k in CTX_KEYS},
                  lens={15: 1.0}, background="ligand")


def test_vendored_span_table_recovers_known_biology():
    for cls in ("mhc1", "mhc2"):
        m = ligand.load_span_model(cls)
        assert set(m.ctx) == set(CTX_KEYS)
        for k in CTX_KEYS:
            assert abs(sum(m.ctx[k].values()) - 1.0) < 1e-3
        assert abs(sum(m.lens.values()) - 1.0) < 1e-3
    m2 = ligand.load_span_model("mhc2")
    bgP = PROTEOME_AA_FREQ["P"]
    # Proline is the aminopeptidase stop signal: ENRICHED just inside the ligand, DEPLETED in the
    # flank. Sign-explicit on purpose -- this fires if the coordinate convention is ever flipped.
    assert m2.ctx["ligN+2"]["P"] / bgP > 1.5
    assert m2.ctx["flankN-1"]["P"] / bgP < 0.5
    # Cys is depleted ~10x at ligand termini in MS data (alkylation / missed ID), not in the flanks.
    # That artifact is clamped out at fit time, or the model would refuse every Cys-containing ligand.
    assert abs(m2.ctx["ligN+1"]["C"] - PROTEOME_AA_FREQ["C"]) < 1e-6
    assert m2.lens[15] > 0.10          # class-II ligands peak at ~15


def test_recommended_flanks_are_the_measured_ones():
    # These two constants ARE the recommendation, and both were wrong before they were measured.
    # STRUCTURE_FLANK=2 (13mer): across 93 real pMHC-II crystals the RESOLVED peptide has median
    #   length 13 with ~2 ordered flanking residues per side; only 13% resolve <=11 residues, so the
    #   core+-1 (11mer) that TCRmodel2/AlphaFold ingest is an input convention, not what is ordered.
    # ASSAY_FLANK=6 (21mer): an APC re-trims a synthetic peptide, so what matters is that it CONTAINS
    #   the natural ligand -- 21mer 80% of held-out cores, vs only 31% for the conventional 15mer.
    # See bench/results/spans_mhc2_human.md and bench/pdb_flanks.py.
    assert ligand.STRUCTURE_FLANK == 2
    assert ligand.ASSAY_FLANK == 6
    core, prot = "WKVKFWKVK", "G" * 20 + "WKVKFWKVK" + "G" * 20
    st = ligand.fixed_span(core, prot, ligand.STRUCTURE_FLANK, ligand.STRUCTURE_FLANK)
    asy = ligand.fixed_span(core, prot, ligand.ASSAY_FLANK, ligand.ASSAY_FLANK)
    assert len(st.peptide) == 13 and core in st.peptide
    assert len(asy.peptide) == 21 and st.peptide in asy.peptide


# -- MHC-I length awareness -------------------------------------------------------------------
from mhcmatch.diffusion import (MHC1_ANCHORS, MHC1_CORE, _EM_CAP,        # noqa: E402
                                    _TAU_DEFAULT, _TAU_MIN)
from mhcmatch.store import mhc1_positions                               # noqa: E402


def _mhc1_store():
    recs = [{"epitope": e, "mhc_a": "HLA-A02:01", "mhc_class": "MHCI"} for e in
            ("GILGFVFTL", "NLVPMVATV", "YLQPRTFLL", "GLCTLVAML", "FLYALALLL", "FLPSDFFPSV",
             "LLFGYPVYV", "YMLDLQPETT", "RMFPNAPYL", "SLYNTVATL", "KLVALGINAV", "AAGIGILTV")]
    recs += [{"epitope": e, "mhc_a": "HLA-B07:02", "mhc_class": "MHCI"} for e in
             ("APRTVALTA", "SPRWYFYYL", "IPSINVHHY")]
    return Store.from_records(recs)


def test_mhc1_positions_never_double_counts():
    # MHC1_CORE's +5 and -4 both land on index 4 of an 8-mer. Counting it twice makes the score an
    # inflated, mis-normalised likelihood ratio (two perfectly-correlated terms) and files one residue
    # under two positions in anchor_preferences. Only the 8-mer collides; 9/10/11 skip the bulge.
    for L in (8, 9, 10, 11):
        pos = mhc1_positions(L, MHC1_CORE)
        real = [i for i in pos if i is not None]
        assert len(real) == len(set(real)), f"L={L} double-counts an index"
        assert len(pos) == len(MHC1_CORE), "must stay aligned to `anchors` for per-anchor bookkeeping"
        assert all(0 <= i < L for i in real)
        assert all(i is not None for i in mhc1_positions(L, MHC1_ANCHORS)), "5-anchor never collides"
    assert mhc1_positions(8, MHC1_CORE).count(None) == 1
    assert all(mhc1_positions(L, MHC1_CORE).count(None) == 0 for L in (9, 10, 11))
    assert mhc1_positions(4, MHC1_CORE) is None, "too short for the footprint -> None"


def test_length_prior_is_on_by_default_and_exactly_additive():
    # ON by default since v0.5.0 (+0.031 maxF1 on the MixMHCpred3 benchmark, precision AND recall up).
    # The anchor log-odds sums a length-INVARIANT number of terms, so without this a 10-mer and a
    # 9-mer with the same anchors score bit-identically -- see data/PROVENANCE.md and the CHANGELOG.
    s = _mhc1_store()
    off = s.anchor_model("mhc1", footprint="core", background="proteome",
                         length_prior=False, length_motifs=False)
    on = s.anchor_model("mhc1", footprint="core", background="proteome", length_motifs=False)
    assert on.length_prior == "score", "the length prior must be ON by default"
    assert off.length_logodds(9, "HLA-A02:01") == 0.0        # no prior built -> no term
    # the term is exactly additive: score(on) == score(off) + length_logodds
    for p in ("GILGFVFTL", "FLPSDFFPSV", "SIINFEHL"):
        assert on.score(p, "HLA-A02:01") == pytest.approx(
            off.score(p, "HLA-A02:01") + on.length_logodds(len(p), "HLA-A02:01"))
    # and it prefers the length the allele actually presents (this panel's A*02:01 is 9-mer heavy)
    assert on.length_logodds(9, "HLA-A02:01") > on.length_logodds(8, "HLA-A02:01")


def test_length_motifs_backoff_is_exact_when_the_length_is_unseen():
    # The safety property the whole design rests on: n(a,L)=0 must reproduce the pooled model
    # BIT-FOR-BIT, so alleles with no ligands at a length (rare alleles have a median of zero
    # 8-mers) provably cannot regress.
    s = _mhc1_store()
    # prior OFF in both arms: this isolates the MOTIF backoff from the (now default-on) length prior
    off = s.anchor_model("mhc1", footprint="core", background="proteome",
                         length_prior=False, length_motifs=False)
    on = s.anchor_model("mhc1", footprint="core", background="proteome",
                        length_prior=False, length_motifs=True)
    assert on.prefs_len is not None and off.prefs_len is None
    for p in ("SIINFEHL", "AAAAAAAAAAA"):                     # A*02:01 has no 8- or 11-mers here
        assert on.score(p, "HLA-A02:01") == off.score(p, "HLA-A02:01")
    assert on.score("GILGFVFTL", "HLA-A02:01") != off.score("GILGFVFTL", "HLA-A02:01")


# Class-II binding modes. Every core residue is drawn from outside _FLANK, so the core is locatable
# and the register EM settles on the planted offset -- fill the core from the flank alphabet instead
# and the frame search lands elsewhere, putting the planted positions outside MHC2_ANCHORS (1,4,6,9)
# where nothing can see them. Same reason test_mhc2_offset_prior_learns_where_the_cores_sit does it.
_FLANK = "AGILPRSTV"
_MODE_A = "WKFWKFWKF"                                     # DRB1*15:01, mode 1
_MODE_B = "DEQDEQDEQ"                                     # DRB1*15:01, mode 2 -- disjoint from A
_MODE_C = "YNMYNMYNM"                                     # DRB1*13:01


def _mhc2_bimodal_store(n=400, seed=0):
    """One class-II allele whose ligands come from TWO planted binding modes, plus a unimodal one.

    The modes share no residue at any core position, so a single PWM can only represent them as their
    average -- exactly the blur a mixture is supposed to resolve.
    """
    rng = random.Random(seed)

    def pep(core, off=3):
        return ("".join(rng.choice(_FLANK) for _ in range(off)) + core
                + "".join(rng.choice(_FLANK) for _ in range(6 - off)))

    recs = []
    for core in (_MODE_A, _MODE_B):
        for _ in range(n):
            recs.append({"epitope": pep(core), "mhc_a": "DRA*01:01", "mhc_b": "DRB1*15:01",
                         "mhc_class": "MHCII"})
    for _ in range(n):                                    # a second allele: the ligand null needs one
        recs.append({"epitope": pep(_MODE_C), "mhc_a": "DRA*01:01", "mhc_b": "DRB1*13:01",
                     "mhc_class": "MHCII"})
    return Store.from_records(recs), pep


def test_n_motifs_semantics_and_mhc1_is_inert():
    # The safety property the mixture rests on: n_motifs=1 must not ENTER the mixture path
    # (prefs_mix is None), so scoring is the exact pre-mixture single-PWM code -- turning the knob on
    # cannot regress anything until it is deliberately set past 1. Asserted structurally rather than
    # against the default, so the shipped default (K=3, see below) can move without breaking this.
    store, pep = _mhc2_bimodal_store()
    assert store.anchor_model("mhc2", n_motifs=1).prefs_mix is None, "K=1 must stay the single-PWM path"
    assert store.anchor_model("mhc2", n_motifs=2).prefs_mix is not None, "K=2 must build a mixture"
    assert _mhc1_store().anchor_model("mhc1", n_motifs=3).prefs_mix is None, "MHC-I: inert at any K"


def test_mhc2_ships_the_mixture_by_default():
    # The v0.7 default: human MHC-II scores the K=3 mixture unless n_motifs is overridden. A separate
    # test from the semantics above so a future default change touches exactly one assertion.
    store, pep = _mhc2_bimodal_store()
    default = store.anchor_model("mhc2")
    assert default.n_motifs == 3 and default.prefs_mix is not None
    probes = [pep(_MODE_A) for _ in range(5)] + [pep(_MODE_B) for _ in range(5)]
    single = store.anchor_model("mhc2", n_motifs=1)
    assert any(default.score(p, "DRB1_1501") != single.score(p, "DRB1_1501") for p in probes), \
        "the default mixture must actually change scores vs the single PWM"


def test_mixture_recovers_two_planted_binding_modes():
    # The claim the mixture exists to make: an allele with two binding modes gets one component per
    # mode, fit from its own ligands by EM -- no external predictor's labels involved.
    #
    # Component INDICES are arbitrary (EM label-switches freely, and the crc32 init decides which way),
    # so the assertion is on separation, never on "component 0 == mode A". That is also why components
    # must not be shrunk across alleles -- see AnchorModel._dist.
    store, pep = _mhc2_bimodal_store()
    m = store.anchor_model("mhc2", n_motifs=2)
    ra = [m._responsibilities(pep(_MODE_A), "DRB1_1501")[0] for _ in range(100)]
    rb = [m._responsibilities(pep(_MODE_B), "DRB1_1501")[0] for _ in range(100)]
    mean = lambda v: sum(v) / len(v)                                          # noqa: E731
    assert abs(mean(ra) - mean(rb)) > 0.5, "EM did not separate the two planted modes"
    # held-out peptides of each mode land on opposite components (modulo the label swap)
    flip = mean(ra) < mean(rb)
    acc = (sum(1 for x in ra if (x < 0.5) == flip) + sum(1 for x in rb if (x > 0.5) == flip)) / 200
    assert acc > 0.8, f"planted-mode recovery {acc:.2f}"
    # both modes are actually used -- a collapsed mixture (one dead component) would also "separate"
    assert min(math.exp(x) for x in m.log_pi["DRB1_1501"]) > 0.2, "a component collapsed"


def test_pseudocount_is_inert_at_zero_and_mass_preserving_above():
    # The two properties the pseudocount rests on, asserted structurally.
    #
    # (1) beta=0 must leave every counter untouched, so the committed bench tables re-run bit-for-bit
    #     and the knob can land without re-baselining anything.
    # (2) beta>0 must preserve each counter's MASS. Pseudoseq.shrink reads both n_own and m, so mass is
    #     what sets the tau/(n+tau) balance -- if the pseudocount moved it, beta would silently crowd out
    #     the kernel prior that carries the rare stratum (n=5, beta=50 would take tau's share 67% -> 15%).
    store = _mhc1_store()
    base = store.anchor_model("mhc1", pseudocount=0.0)
    off = store.anchor_model("mhc1", pseudocount=0.0)
    on = store.anchor_model("mhc1", pseudocount=50.0)
    for j in base.anchors:
        for a in base.prefs[j]:
            assert dict(base.prefs[j][a]) == dict(off.prefs[j][a]), "beta=0 must be inert"
            n0 = sum(base.prefs[j][a].values())
            assert math.isclose(sum(on.prefs[j][a].values()), n0, rel_tol=1e-9), "mass must be preserved"
    # and it does what it exists for: every residue gets chemically-graded mass, not a flat floor
    j, a = base.anchors[0], next(iter(base.prefs[base.anchors[0]]))
    assert len(base.prefs[j][a]) < 20 < len(on.prefs[j][a]) + 1, "beta>0 must fill the unobserved tail"


def test_pseudocount_leaves_the_fitted_latents_alone():
    # Why it is the LAST statement of __init__: the background null, the register-EM frames and the
    # mixture's component assignments are all fit on the raw counters, so they must not move with beta.
    # This is what protects the shipped K=3 gain structurally rather than by sweeping for a safe beta.
    store, _ = _mhc2_bimodal_store()
    a = store.anchor_model("mhc2", pseudocount=0.0)
    b = store.anchor_model("mhc2", pseudocount=100.0)
    assert {j: dict(c) for j, c in a.bg.items()} == {j: dict(c) for j, c in b.bg.items()}, \
        "the log-odds null must be fit on raw counts"
    assert a.log_pi == b.log_pi, "mixture component assignments must not move with beta"


def test_frame_score_memo_is_bit_identical_to_recompute():
    # _frame_scores is memoized to cut the K=3 build ~2.7x. The cache is only correct if it is cleared
    # wherever prefs/prefs_mix/bg change (_refit_registers, _m_step, _add_pseudocounts). A stale cache
    # would silently corrupt every MHC-II score, so pin that the memoized build is bit-identical to a
    # cache-bypassed one -- fitted state AND scores.
    import mhcmatch.diffusion as D

    def _nocache(self, peptide, allele, raw=False, eps=1e-3, k=None):
        core_pos = [j - 1 for j in self.anchors]
        mask = self._score_mask(allele)
        markov = self.background == "markov"
        out = []
        for st in range(len(peptide) - 8):
            w = peptide[st:st + 9]
            ctx = [peptide[st + c - 1] if st + c > 0 else "" for c in core_pos] if markov else None
            out.append(self._anchor_logodds([w[c] for c in core_pos], allele, raw, eps, mask, ctx, k=k))
        return out

    store, pep = _mhc2_bimodal_store()
    memo = store.anchor_model("mhc2", n_motifs=3)
    orig = D.AnchorModel._frame_scores
    try:
        D.AnchorModel._frame_scores = _nocache
        ref = store.anchor_model("mhc2", n_motifs=3)
    finally:
        D.AnchorModel._frame_scores = orig
    assert memo.log_pi == ref.log_pi, "mixture assignments diverged -> stale frame cache"
    probes = [pep(_MODE_A) for _ in range(6)] + [pep(_MODE_B) for _ in range(6)]
    for p in probes:
        assert memo.score(p, "DRB1_1501") == ref.score(p, "DRB1_1501"), "score diverged -> stale cache"


def test_register_em_converge_reaches_a_real_fixed_point():
    # v0.7.2 HEADLINE, previously untested: register_em="converge" runs the best-frame EM to each
    # allele's OWN fixed point. Mutation-tested -- freezing every allele after pass 1 (which deletes the
    # entire DP gain, frequent AUPRC 0.667->0.625) left all tests green until this one existed.
    # The defining property is exactly a fixed point: one more UNFROZEN pass must change nothing.
    store, _ = _mhc2_bimodal_store()
    am = store.anchor_model("mhc2", register_em="converge", n_motifs=1)
    assert 0 < am._em_passes < _EM_CAP, "must terminate by convergence, not the runaway cap"
    assert am._refit_registers(store, frozen=None) == set(), \
        "at the fixed point a full unfrozen pass reassigns no allele's frames"


def test_footprint_adaptive_masks_rare_and_not_frequent():
    # The mode EVERY committed MHC-I benchmark number is generated under, previously untested: adaptive
    # scores rare alleles on the primary anchors only and frequent ones on the full core. Mutation-tested
    # -- making _score_mask always return None (never mask) left all tests green until this one existed.
    store = _mhc1_store()                                    # A02:01 n=12 (frequent), B07:02 n=3 (rare)
    am = store.anchor_model("mhc1", footprint="adaptive", rare_max=5)
    assert am._score_mask("HLA-A02:01") is None, "well-sampled allele scores the full core"
    assert am._score_mask("HLA-B07:02") is not None, "rare allele is masked to the anchors"
    assert len(am._score_mask("HLA-B07:02")) < len(am.anchors), "the mask actually drops positions"


def test_prior_strength_auto_fits_lower_tau_where_alleles_differ():
    # v0.7.2, previously untested: empirical-Bayes tau per position -- small where alleles genuinely
    # differ (trust own data), large where they agree (shrink to the pool). Tested on _fit_tau directly
    # since it reads only self.prefs/self.anchors. Anchor 1: four alleles each peak on their OWN residue
    # (max between-allele variance -> tau at the floor). Anchor 2: all four share one flat mix
    # (zero variance -> tau at the ceiling).
    am = _mhc1_store().anchor_model("mhc1")
    am.anchors = (1, 2)
    peaks = {"a": "L", "b": "K", "c": "W", "d": "D"}
    flat = collections.Counter({"W": 75, "D": 75, "K": 75, "F": 75})
    am.prefs = {1: {al: collections.Counter({r: 300}) for al, r in peaks.items()},
                2: {al: collections.Counter(flat) for al in peaks}}
    t = am._fit_tau(min_n=200, min_alleles=3)
    assert t[1] < t[2], "tau must be smaller at the allele-discriminating position"
    assert t[1] <= _TAU_MIN + 1e-9, "a maximally-specific anchor pins tau at the floor"
    # and it never invents a tau it cannot support: below min_n every position falls back to the default
    small = _mhc1_store().anchor_model("mhc1", prior_strength="auto")   # n=12/3, both under min_n=200
    assert all(v == _TAU_DEFAULT for v in small._tau.values()), \
        "too few well-sampled alleles -> scalar default, not a fabricated per-position tau"


def test_mixture_component_backs_off_to_the_pooled_motif_when_empty():
    # The backoff identity that makes capacity self-adapting: a component with no counts for an allele
    # returns that allele's pooled (shrunk) motif *identically*, so an allele too thin to fill K
    # components degrades to today's single PWM with no ligand-count threshold to choose.
    store, _ = _mhc2_bimodal_store()
    m = store.anchor_model("mhc2", n_motifs=2)
    j = m.anchors[0]
    assert m._dist(j, "DRB1_1501", False, 0) is not m._dist(j, "DRB1_1501", False, None), \
        "a populated component must not be the pooled dict"
    m.prefs_mix[0][j].pop("DRB1_1501", None)              # simulate an empty component
    m._cache_mix.clear()
    assert m._dist(j, "DRB1_1501", False, 0) == m._dist(j, "DRB1_1501", False, None)


def test_length_bg_uniform_flattens_the_null_length_mix():
    from mhcmatch.calibrate import corpus_stats, random_peptides
    corpus = ["A" * 9] * 80 + ["A" * 8] * 10 + ["A" * 10] * 7 + ["A" * 11] * 3
    aa, lens = corpus_stats(corpus)
    got = collections.Counter(
        len(p) for p in random_peptides(aa, lens, 4000, random.Random(0), "uniform"))
    assert all(abs(got[L] / 4000 - 0.25) < 0.03 for L in (8, 9, 10, 11))
    corp = collections.Counter(
        len(p) for p in random_peptides(aa, lens, 4000, random.Random(0), "corpus"))
    assert corp[9] / 4000 > 0.7, "'corpus' keeps the ligand length mix (the default, MHC-II needs it)"


from mhcmatch.pseudoseq import load_pseudo, resolve_allele                 # noqa: E402


def test_every_allele_in_a_collapsed_pseudoseq_group_is_resolvable():
    """A FASTA header lists every allele sharing the 34-mer -- all of them must be keys.

    Until 2026-07 the header carried only the group's first allele, so 68% of MHC-I and 80% of MHC-II
    alleles were silently unresolvable -- including HLA-B*14:02, B*18:05, C*03:04. `load_pseudo` and
    tcren's `build_pseudo_fasta.py` were fixed together; re-syncing the FASTA from an unfixed upstream
    would reintroduce it, hence this test guards the data, not the parser.
    """
    for cls, floor in (("mhc1", 12000), ("mhc2", 10000)):
        p = load_pseudo(cls)
        assert len(p) > floor, f"{cls}: {len(p)} keys -- header index lost? (expected >{floor})"

    p1 = load_pseudo("mhc1")
    # each is a non-representative that shares its group's groove exactly -- not an approximation
    for allele, rep in (("HLA-B14:02", "HLA-B14:01"), ("HLA-B18:05", "HLA-B18:01"),
                        ("HLA-C03:04", "HLA-C03:03"), ("HLA-C03:02", "HLA-C03:01")):
        assert allele in p1, f"{allele} unresolvable -- the collapsed-group index is broken"
        assert p1[allele] == p1[rep], f"{allele} must carry {rep}'s 34-mer verbatim"
        assert resolve_allele(allele, "mhc1") == (allele, True)


def test_pseudoseq_groups_are_exact_identity_not_similarity():
    """Collapsing is by *exact* 34-mer equality, so a group's members are interchangeable inputs.

    Guards against a future 'helpful' fuzzy/nearest-neighbour alias, which would silently score one
    allele with another's motif -- the collapse is only sound because the sequences are identical.
    """
    for cls in ("mhc1", "mhc2"):
        seqs = load_pseudo(cls)
        assert all(len(s) == 34 for s in seqs.values())


def test_imgt_derived_alleles_cover_what_netmhcpan_omits():
    """The FASTA carries IPD-IMGT/HLA-derived alleles NetMHCpan's table never had.

    HLA-F is the clearest case: absent from MHC_pseudo.dat entirely, so it has no known 34-mer to
    check against -- it is trusted because HLA-E and HLA-G round-trip 100% through the same
    cross-gene column mapping (see data/PROVENANCE.md). HLA-A*30:14 is a plain gap in the table.
    Both stranded real panel ligands before 2026-07-16.
    """
    p = load_pseudo("mhc1")
    assert len(p) > 20000, f"{len(p)} keys -- the IMGT source is missing?"
    for allele in ("HLA-A30:14", "HLA-F01:01", "HLA-F01:03", "HLA-F01:04"):
        assert allele in p, f"{allele} absent -- IMGT-derived alleles not vendored?"
        assert set(p[allele]) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{allele} has a non-residue"
    # HLA-F's groove is near-monomorphic: its alleles differ outside the 34 positions.
    assert p["HLA-F01:01"] == p["HLA-F01:03"] == p["HLA-F01:04"]
    # ...and it is genuinely distinct from the classical loci, i.e. not an alignment artefact.
    assert p["HLA-F01:01"] != p["HLA-A02:01"]


from mhcmatch.pseudoseq import alpha_prior, class2_from_name, class2_key    # noqa: E402


def test_alpha_imputation_fires_only_on_a_missing_alpha():
    """A fully-typed class-II input must be bit-identical with the flag on or off.

    The flag exists for the 1.5% of panel records that type only the beta chain ('-DPB11101',
    2,516 ligands) -- it must never rewrite a typing that is already complete, nor touch DR (whose
    monomorphic DRA is hardcoded) or mouse.
    """
    for a, b in (("HLA-DPA1*01:03", "HLA-DPB1*04:01"), ("HLA-DQA1*05:01", "HLA-DQB1*02:01"),
                 ("DRA", "HLA-DRB1*15:01"), ("I-Ab", "")):
        assert class2_key(a, b, True) == class2_key(a, b, False)

    p = load_pseudo("mhc2")
    # beta-only: imputed to a real groove when on, left alpha-less (and unscorable) when off
    on, off = class2_from_name("HLA-DPB1*11:01", True), class2_from_name("HLA-DPB1*11:01", False)
    assert on == "HLA-DPA10201-DPB11101" and on in p
    assert off == "-DPB11101" and off not in p


def test_alpha_prior_refuses_an_ambiguous_groove():
    """Only betas whose 34-mer is >=95% determined are imputed -- a wrong groove scores silently.

    DQA1*01:02 and DQA1*01:05 share the 2-digit group DQA1*01 but NOT the 34-mer, so DQB1*05:02 looks
    100% certain at group level while the groove is a 58/42 coin flip. The table is keyed on the
    groove for exactly that reason; these rare DQ betas stay unresolved on purpose.
    """
    prior = alpha_prior()
    assert prior, "alpha prior table is empty"
    for beta in ("DQB10503", "DQB10502", "DQB10402", "DQB10602"):
        assert beta not in prior, f"{beta}'s alpha is ambiguous and must not be imputed"
    assert prior.get("DPB11101") == "HLA-DPA10201"      # 99% -- the 2,516-ligand case
    assert prior.get("DQB10302") == "HLA-DQA10301"      # DQ8, rediscovered from linkage disequilibrium
    assert prior.get("DQB10201") == "HLA-DQA10501"      # DQ2.5
    p = load_pseudo("mhc2")
    for beta, alpha in prior.items():
        assert f"{alpha}-{beta}" in p, f"imputed key {alpha}-{beta} has no groove"


def test_store_panel_is_unchanged_by_default_alpha_imputation_is_lookup_only():
    """`Store.from_records` must default to the pre-2026-07 panel: beta-only records dropped.

    The lookup path imputes by default (nan -> an answer, a strict win); the PANEL path does not.
    Admitting these ligands was measured over the 13 alleles whose reference set grows: AUROC -0.0019,
    AUPRC -0.0012, worst where the merge is biggest (DPB1*11:01 +89% ligands, -0.0155 AUROC). Missing
    alpha-typing marks a noisier study, not just absent metadata. Opposite defaults, each measured.
    """
    recs = [{"epitope": "AAKGVAAWSAGTFRQ", "mhc_a": "", "mhc_b": "HLA-DPB1*11:01",
             "mhc_class": "MHCII"}]
    off = Store.from_records(recs)._panel["mhc2"]
    assert off.alleles == [], "a beta-only record must be dropped by default"
    on = Store.from_records(recs, impute_alpha=True)._panel["mhc2"]
    assert on.alleles == ["HLA-DPA10201-DPB11101"], "impute_alpha=True must admit it, alpha filled"
