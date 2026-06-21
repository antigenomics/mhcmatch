"""MHC restriction & presentation from a reference epitope panel.

Productionizes the validated reverse-problem method (seqtree ``bench/bench_mhc_guess.py``):
index reference peptides by their anchored *presentation* signature
(:func:`seqtree.layout.presentation_features`), widen the search scope around a query until it
has enough neighbours, then rank presenting alleles by neighbour **vote fraction** and score
**confidence** by a binomial-tail enrichment over the panel background. The vote fraction is the
ranking statistic (robust to panel skew); the enrichment is the non-binder filter.

Significance theory: ``appendix/mhcmatch.tex`` §2-3 (forward per-allele E-value + reverse problem).
"""
from __future__ import annotations

import csv
import gzip
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass

from seqtree import KmerIndex, SearchParams, layout

_CLASS = {"MHCI": "mhc1", "I": "mhc1", "mhc1": "mhc1",
          "MHCII": "mhc2", "II": "mhc2", "mhc2": "mhc2"}
_SPECIES = {"human": "HomoSapiens", "mouse": "MusMusculus"}
_AA = set("ACDEFGHIKLMNPQRSTVWY")
_SCOPES = (0, 1, 2, 3)
_DEFAULT_LENGTHS = {"mhc1": (8, 9, 10, 11), "mhc2": (13, 14, 15, 16, 17, 18)}


def infer_class(peptide: str) -> str:
    """Heuristic class from length: MHC-I if <=11, else MHC-II. Pass ``cls`` to override."""
    return "mhc1" if len(peptide) <= 11 else "mhc2"


@dataclass
class Restriction:
    allele: str
    vote: float        # neighbour vote fraction P(allele | neighbours) -- ranking score
    enrichment: float  # -log10 binomial-tail p vs panel background -- confidence
    n_votes: int
    binder: bool

    def __iter__(self):
        return iter((self.allele, self.vote, self.enrichment, self.binder))


@dataclass
class Decomposition:
    peptide: str
    tcr_facing: str    # anchors masked with X  (recognition readout)
    presentation: str  # TCR-facing masked with X (anchor readout)
    anchors: tuple     # 0-based anchor indices


def _binom_sf(k, n, p):
    """P(Binomial(n, p) >= k) -- upper tail."""
    if k <= 0:
        return 1.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    return min(1.0, sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1)))


def _mhc2_core_anchors(peptide: str) -> tuple:
    """0-based P1/P4/P6/P9 indices of the register-anchored 9-mer core (one-pass register trick)."""
    if len(peptide) < 9:
        return ()
    best_s = max(range(len(peptide) - 8),
                 key=lambda s: layout._core_anchor_score(peptide[s:s + 9]))
    return tuple(best_s + j for j in (0, 3, 5, 8))


def anchor_indices(peptide: str, cls: str) -> tuple:
    """0-based anchor positions for a peptide: class-I P2/PΩ, class-II core P1/P4/P6/P9."""
    if cls == "mhc2":
        return _mhc2_core_anchors(peptide)
    return tuple(sorted(layout.spec_for(cls).resolve(len(peptide))))


def resolve_anchor_index(peptide: str, cls: str, anchor: int):
    """0-based index of a scoring ``anchor`` in ``peptide`` (or None if out of range).

    MHC-I: ``anchor`` is a 1-based peptide position (negatives count from the C-terminus).
    MHC-II: ``anchor`` is a 1-based position *within the register-anchored 9-mer core* (P1..P9).
    """
    if cls == "mhc2":
        if len(peptide) < 9:
            return None
        s = max(range(len(peptide) - 8),
                key=lambda i: layout._core_anchor_score(peptide[i:i + 9]))
        idx = s + (anchor - 1)
        return idx if s <= idx < s + 9 else None
    idx = (anchor - 1) if anchor > 0 else (len(peptide) + anchor)
    return idx if 0 <= idx < len(peptide) else None


class _Panel:
    """One MHC class: presentation-signature KmerIndex + allele bookkeeping."""

    def __init__(self, cls):
        self.cls = cls
        self.epitopes = []
        self.alleles = []
        self.weights = []

    def add(self, epitope, allele, weight=1.0):
        self.epitopes.append(epitope)
        self.alleles.append(allele)
        self.weights.append(weight)

    def build(self):
        feats = [layout.presentation_features(e, self.cls, register="anchored")
                 for e in self.epitopes]
        self.allele_to_id = {}
        ids = []
        for a in self.alleles:
            self.allele_to_id.setdefault(a, len(self.allele_to_id))
            ids.append(self.allele_to_id[a])
        self.index = KmerIndex.build(feats, alphabet="aa", allele_ids=ids) if feats else None
        counts = Counter(self.alleles)
        total = len(self.alleles) or 1
        self.panel = sorted(counts)
        self.freq = {a: counts[a] / total for a in self.panel}

    def tally(self, query, lo=10, hi=100):
        """Counter(allele -> votes) from the query's anchored-signature neighbours, scope-widened."""
        if self.index is None:
            return None
        feats = layout.presentation_features(query, self.cls, register="anchored")
        cands = []
        for sc in _SCOPES:
            p = SearchParams(max_subs=sc, engine="seqtm")
            cands = [c for c in self.index.seed_and_gather([feats], p, 1, -1, 1)[0]
                     if self.epitopes[c.peptide_id] != query]
            if len(cands) >= lo:
                break
        if not cands:
            return None
        return Counter(self.alleles[c.peptide_id] for c in cands[:hi])


class Store:
    """Searchable reference panel of presented peptides, partitioned by MHC class."""

    def __init__(self):
        self._panel = {"mhc1": _Panel("mhc1"), "mhc2": _Panel("mhc2")}

    # -- construction ---------------------------------------------------------
    @classmethod
    def from_records(cls, records):
        """records: dicts with ``epitope``, ``mhc_a`` (or ``mhc``), ``mhc_class``; optional
        ``weight`` (default 1.0) confidence-weights the peptide in anchor-preference estimation."""
        from .pseudoseq import class2_key
        store = cls()
        for r in records:
            c = _CLASS.get(str(r.get("mhc_class", "")).strip())
            ep = str(r.get("epitope", "")).strip().upper()
            allele = str(r.get("mhc_a") or r.get("mhc") or "").strip()
            if c is None or not ep or not allele or not all(x in _AA for x in ep):
                continue
            if c == "mhc2":  # key class II by the alpha-beta pair (locus-aware)
                allele = class2_key(allele, str(r.get("mhc_b") or "").strip())
            store._panel[c].add(ep, allele, float(r.get("weight", 1.0) or 1.0))
        for p in store._panel.values():
            p.build()
        return store

    @classmethod
    def from_pmhc(cls, path=None, tier="full", species=None, classes=("mhc1", "mhc2")):
        """Load the isalgo/pmhc_data TSV(.gz). ``species`` filters the *MHC* species
        (``"human"`` / ``"mouse"``). If ``path`` is None, uses ``$MHCMATCH_PMHC/pmhc_<tier>.tsv.gz``."""
        if path is None:
            base = os.environ.get("MHCMATCH_PMHC")
            if not base:
                raise ValueError("pass path= or set MHCMATCH_PMHC to the pmhc_data directory")
            path = os.path.join(base, f"pmhc_{tier}.tsv.gz")
        sp = _SPECIES.get(species) if species else None
        keep = {_CLASS[c] for c in classes}
        csv.field_size_limit(10 ** 7)
        op = gzip.open if str(path).endswith(".gz") else open
        recs = []
        with op(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                c = _CLASS.get(str(row.get("mhc_class", "")).strip())
                if c is None or c not in keep:
                    continue
                if sp and row.get("mhc_species") != sp:
                    continue
                recs.append(row)
        return cls.from_records(recs)

    def __len__(self):
        return sum(len(p.epitopes) for p in self._panel.values())

    def alleles(self, cls):
        return list(self._panel[cls].panel)

    # -- forward problem: restriction / presentation --------------------------
    def _allele_set(self, panel, alleles):
        if alleles == "all":
            return panel.panel
        if isinstance(alleles, str):
            alleles = [alleles]
        return [a for a in alleles if a in panel.freq]

    def restriction(self, peptide, cls=None, alleles="all", top=10, alpha=0.05):
        """Rank presenting alleles for ``peptide`` (vote fraction), flag binders (enrichment).

        ``alleles``: ``"all"``, a single allele, or a list. ``alpha``: per-allele significance for
        the non-binder flag (binder iff binomial-tail p <= alpha and the allele got votes).
        Returns ``[]`` when the peptide has no presentation-signature neighbours (treat as unknown).
        """
        peptide = peptide.strip().upper()
        cls = cls or infer_class(peptide)
        panel = self._panel[cls]
        tally = panel.tally(peptide)
        if tally is None:
            return []
        n = sum(tally.values())
        thr = -math.log10(alpha)
        out = []
        for a in self._allele_set(panel, alleles):
            k = tally.get(a, 0)
            enr = -math.log10(max(_binom_sf(k, n, panel.freq[a]), 1e-300)) if k else 0.0
            out.append(Restriction(a, k / n, enr, k, enr >= thr and k > 0))
        out.sort(key=lambda r: (r.vote, r.enrichment), reverse=True)
        return out[:top]

    def is_binder(self, peptide, allele, cls=None, alpha=0.05):
        res = self.restriction(peptide, cls=cls, alleles=[allele], top=1, alpha=alpha)
        return bool(res and res[0].binder)

    def is_presented(self, peptide, cls=None, alpha=0.05):
        """Overall presentation: does any panel allele present this peptide?"""
        return any(r.binder for r in self.restriction(peptide, cls=cls, alpha=alpha))

    def scan_protein(self, protein, cls="mhc1", alleles="all", lengths=None, alpha=0.05, top=3):
        """Slide all binding-length windows over ``protein`` and return presented peptides.

        Returns ``[(position, peptide, [Restriction, ...]), ...]`` for windows with >=1 binder.
        Per-window thresholding only (FWER/FDR over windows x panel is appendix §5, Phase 1).
        """
        protein = "".join(protein.split()).upper()
        lengths = lengths or _DEFAULT_LENGTHS[cls]
        out = []
        for L in lengths:
            for i in range(len(protein) - L + 1):
                pep = protein[i:i + L]
                if not all(c in _AA for c in pep):
                    continue
                binders = [r for r in self.restriction(pep, cls, alleles, top=top, alpha=alpha)
                           if r.binder]
                if binders:
                    out.append((i, pep, binders))
        return out

    # -- anchor / TCR-facing split -------------------------------------------
    def decompose(self, peptide, cls=None, allele=None):
        """Split ``peptide`` into anchor and TCR-facing parts, each masked with ``X``.

        ``tcr_facing``: anchors -> X (the recognition readout). ``presentation``: TCR-facing -> X
        (the anchor readout). ``allele`` is accepted for forward-compat (allele-specific learned
        anchors, Phase 1); v0 uses class-default anchor positions.
        """
        peptide = peptide.strip().upper()
        cls = cls or infer_class(peptide)
        anchors = set(anchor_indices(peptide, cls))
        tcr = "".join(layout.MASK if i in anchors else c for i, c in enumerate(peptide))
        present = "".join(c if i in anchors else layout.MASK for i, c in enumerate(peptide))
        return Decomposition(peptide, tcr, present, tuple(sorted(anchors)))

    # -- diffusion-powered forward scorer -------------------------------------
    def anchor_model(self, cls="mhc1", h=2.0, prior_strength=10.0, anchors=None, learn_weights=True):
        """Anchor-factored presentation model with cross-allele kernel-shrinkage diffusion.

        See :class:`mhcmatch.diffusion.AnchorModel`. The diffusion rescues rare alleles by borrowing
        anchor preferences from groove-similar frequent ones, with a bounded prior strength so a
        large neighbour cannot swamp a rare allele's own peptides.
        """
        from .diffusion import AnchorModel
        return AnchorModel(self, cls=cls, anchors=anchors, h=h,
                           prior_strength=prior_strength, learn_weights=learn_weights)

    # -- per-allele anchor preferences (feeds pseudoseq diffusion) ------------
    def anchor_preferences(self, cls, anchor):
        """{allele: Counter(residue)} at a 1-based ``anchor`` position (negative from C-term)."""
        panel = self._panel[cls]
        prefs = defaultdict(Counter)
        for ep, a, w in zip(panel.epitopes, panel.alleles, panel.weights):
            idx = resolve_anchor_index(ep, cls, anchor)
            if idx is not None:
                prefs[a][ep[idx]] += w
        return prefs
