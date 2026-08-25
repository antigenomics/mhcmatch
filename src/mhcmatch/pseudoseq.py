"""MHC pseudosequence allele-similarity & cross-allele diffusion.

Each allele is a 34-residue groove **pseudosequence** (NetMHCpan-style; vendored in
``data/{mhci,mhcii}_pseudo.fa``). Allele similarity is an **anchor-factored kernel** over these
positions: ``K_j(a,b) = exp(-d_j(a,b)/h)`` where ``d_j`` is a position-weighted Hamming distance and
the per-anchor weights ``w_j`` say which groove residues govern peptide anchor ``j`` (e.g. MHC-I P2
vs PΩ). :func:`learn_anchor_weights` learns ``w_j`` from data (mutual information between a groove
position and the allele's anchor-residue choice) -- the "feature importance" of each pocket.

Kernel-weighted **shrinkage** (:meth:`Pseudoseq.shrink`) borrows presented-peptide statistics from
similar alleles to rescue rare ones, lifting the seqtree limitation "distinct alleles are distinct
nulls". See the theory appendix §4.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from functools import lru_cache
from importlib import resources

_FA = {"mhc1": "mhci_pseudo.fa", "mhc2": "mhcii_pseudo.fa"}
_LEN = 34


#: A class-I HLA name written without separators, as the deposited screens write it: ``A0201``,
#: ``Cw0401``, ``A*02:01``, ``HLA A0201``. Anchored and locus-restricted so a class-II name
#: (``DRB1*01:01``) and a non-human genus (``BoLA-1:00101``, ``DLA-88*501:01``) cannot match.
_BARE_I = re.compile(r"^(?:HLA[- ])?([ABC])w?\*?(\d{2,3}):?(\d{2,3})[A-Z]?$")


def normalize_allele(a: str) -> str:
    """pmhc allele name -> pseudosequence-FASTA key.

    Drops the ``*`` (``'HLA-A*02:01'`` -> ``'HLA-A02:01'``) and repairs the mouse H-2 dash
    (pmhc ``'H-2Kb'`` -> FASTA ``'H-2-Kb'``).

    Takes one allele name. A cell naming several (``'B0801,C0701'``) is a genotype, not an allele;
    :func:`mhcmatch.rank.split_alleles` is what splits it, and normalising the cell whole is the
    defect that produced ``'HLA-B08:010701'`` -- two names run together into a spelling no table
    has, which then resolves to nothing.
    """
    a = a.replace("*", "")
    if a.startswith("H-2") and len(a) > 3 and a[3] != "-":  # mouse: 'H-2Kb' -> 'H-2-Kb'
        a = "H-2-" + a[3:]
    return a


def hla_spellings(name: str) -> list:
    """Both spellings of a class-I HLA name -- with the field colon and without it.

    The bundled pseudosequence table carries **both**: ``HLA-A02:01`` and ``HLA-A0115`` are each
    keys, because the source tables it was built from disagreed. So does every deposited screen,
    which writes ``A0201``, ``Cw0401``, ``HLA A0201`` or ``A*02:01`` for the same molecule.
    Offering both spellings is what lets :func:`resolve_allele` accept all of them without a
    caller normalising first -- and a benchmark that normalises in its own helper is a second
    convention nobody else can run. Returns ``[]`` for anything that is not a class-I HLA name.
    """
    m = _BARE_I.match((name or "").strip())
    if not m:
        return []
    loc, f1, f2 = m.groups()
    return [f"HLA-{loc}{f1}:{f2}", f"HLA-{loc}{f1}{f2}"]


@lru_cache(maxsize=1)
def alpha_prior() -> dict:
    """``DP/DQ beta chain -> most likely alpha chain``, for typings that omit the alpha.

    Learned from the IEDB-derived panel and vendored (``data/mhc2_alpha_prior.tsv``); a beta is
    listed only when its **34-mer groove** is >=95% determined over >=50 fully-typed ligands. See
    :func:`class2_key`.
    """
    text = resources.files("mhcmatch.data").joinpath("mhc2_alpha_prior.tsv").read_text()
    out = {}
    for line in text.splitlines():
        if line.startswith("#") or line.startswith("beta\t"):
            continue
        f = line.split("\t")
        if len(f) >= 2:
            out[f[0]] = f[1]
    return out


def class2_key(mhc_a: str, mhc_b: str = "", impute_alpha: bool = True) -> str:
    """pmhc class-II allele -> pseudosequence-FASTA key (locus-aware).

    DR (the DRA chain is monomorphic) is keyed by the beta chain alone, e.g.
    ``'HLA-DRB1*01:01' -> 'DRB1_0101'``. DP/DQ are keyed by the alpha-beta pair, e.g.
    ``('HLA-DPA1*01:03', 'HLA-DPB1*04:01') -> 'HLA-DPA10103-DPB10401'``. With no beta chain the
    input is returned unchanged (mouse H-2 and fallbacks).

    ``impute_alpha`` (default on) fills a **missing DP/DQ alpha** from :func:`alpha_prior`, so a
    beta-only typing resolves to a real groove instead of the unscorable ``'-DPB11101'``. This is the
    polymorphic-locus analogue of what DR already gets for free from monomorphic DRA. It fires only
    where the panel pins the *groove* to >=95% over >=50 ligands -- DQA1's polymorphism sits in the
    alpha1 domain the pseudosequence samples, so a name- or 2-digit-group-level rule is not a
    substitute: DQA1*01:02 and DQA1*01:05 share the group DQA1*01 but not the 34-mer, which reads as
    100% certain while the sequence is a 58/42 coin flip. Rare DQ betas are left unresolved on
    purpose -- a wrong groove scores silently, which is worse than not scoring.
    """
    b = (mhc_b or "").strip()
    if mhc_a.startswith("I-"):                        # mouse: 'I-Ab' / 'I-Ek' -> FASTA 'H-2-IAb'
        return "H-2-" + mhc_a.replace("-", "")
    if "DRB" in b:                                   # DR: beta-only, underscore form
        beta = b[4:] if b.startswith("HLA-") else b  # drop the HLA- prefix
        return beta.replace("*", "_").replace(":", "")
    if not b:
        return mhc_a
    beta = b.replace("*", "").replace(":", "").lstrip("-")   # '-DPB11101' is an alpha-less key, not a beta name
    if beta.startswith("HLA-"):
        beta = beta[4:]
    alpha = mhc_a.replace("*", "").replace(":", "")   # NB the HLA- prefix stays: keys are 'HLA-DPA10103-DPB10401'
    if not alpha and impute_alpha:
        alpha = alpha_prior().get(beta, "")
    return f"{alpha}-{beta}"


def class2_from_name(name: str, impute_alpha: bool = True) -> str:
    """Class-II allele *name* (user- or IEDB-typed) -> mhc2 pseudoseq key, locus-aware.

    Handles DR (beta-only ``'HLA-DRB1*15:01' -> 'DRB1_1501'``), the DP/DQ alpha-beta pair given as
    ``'HLA-DQA1*05:01/DQB1*03:01'``, a **DP/DQ beta given alone** (``'HLA-DPB1*11:01' ->
    'HLA-DPA10201-DPB11101'``, the alpha imputed via :func:`alpha_prior` -- see :func:`class2_key`),
    and mouse (``'H2-IAb'`` / ``'I-Ab'`` -> ``'H-2-IAb'``). Falls back to :func:`normalize_allele`
    for anything already in key form.
    """
    a = name.strip()
    au = a.upper()
    if au.startswith("H2-"):
        return "H-2-" + a[3:]
    if au.startswith("H-2"):
        return normalize_allele(a)
    if au.startswith("I-"):
        return class2_key(a)
    if "/" in a:
        x, y = a.split("/", 1)
        return class2_key(x.strip(), y.strip(), impute_alpha)
    if "DRB" in au:
        return class2_key("DRA", a)
    # a DP/DQ beta with no alpha alongside it -- the alpha is missing, not merely unwritten
    if re.search(r"D[PQ]B1", au) and not re.search(r"D[PQ]A1", au):
        return class2_key("", a, impute_alpha)
    return normalize_allele(a)


#: Class-II reporting granularities, finest first. See :func:`class2_report`.
REPORT_MODES = ("pair", "beta", "isotype")

_CHAIN = re.compile(r"^(D[PQR][AB]\d)(\d{4,})$")


def _imgt(chain: str) -> str:
    """``'DQB10301'`` -> ``'DQB1*03:01'``; anything else unchanged. Fields are two digits each, so a
    three-field key round-trips too (``'DRB1010101'`` -> ``'DRB1*01:01:01'``)."""
    m = _CHAIN.match(chain.replace("*", "").replace(":", "").replace("_", ""))
    if not m:
        return chain
    gene, d = m.groups()
    return gene + "*" + ":".join(d[i:i + 2] for i in range(0, len(d), 2))


def class2_report(key: str, mode: str = "pair") -> str:
    """Reduce a class-II key to a reporting granularity.

    - ``"pair"`` -- the key unchanged: ``'DRB1_0101'``, ``'HLA-DQA10501-DQB10301'``. This is
      NetMHCIIpan's own naming and what :func:`class2_key` produces, so it is the default and the
      only mode in which two tools' outputs are directly comparable as strings.
    - ``"beta"`` -- the beta chain alone, in IMGT form: ``'DRB1*01:01'``, ``'DQB1*03:01'``.
    - ``"isotype"`` -- ``'DR'`` / ``'DP'`` / ``'DQ'`` (mouse: ``'H-2'``).

    **Why the coarser modes exist.** A class-II key does not lead with the same chain at every
    isotype: DRA is monomorphic, so DR is keyed by its *beta*, while DP and DQ keys lead with the
    *alpha*. Any comparison that reads the leading gene out of a key is therefore matching DR's beta
    against DP/DQ's alpha -- two different genes, and the alpha is the less polymorphic half. It also
    splits DR against itself, because ``DRB1`` and ``DRB3`` are different leading genes at the same
    isotype. ``"beta"`` and ``"isotype"`` both compare like with like; ``"isotype"`` is the right
    granularity for the question "did the two callers even pick the same molecule family".

    Measured on the class-II arm of the Gamaleya ISP concordance (10,402 rows where both callers
    named an allele): leading-gene agreement 0.401, true isotype agreement **0.527**. The gap is
    1,318 DR-vs-DR pairs differing only in DRB gene.
    """
    if mode not in REPORT_MODES:
        raise ValueError(f"mode must be one of {REPORT_MODES}, got {mode!r}")
    k = (key or "").strip()
    if mode == "pair" or not k:
        return k
    up = k.upper()
    if up.startswith("H-2") or up.startswith("H2-"):        # mouse: one isotype, no alpha/beta key
        return "H-2" if mode == "isotype" else k
    if mode == "isotype":
        return next((iso for iso in ("DR", "DP", "DQ") if iso in up), k)
    tail = k.split("-")[-1]             # 'HLA-DQA10501-DQB10301' and '-DPB11101' both end in the beta
    beta = _imgt(tail)
    return beta if beta != tail else k  # not a class-II chain: hand the key back rather than a stub


@lru_cache(maxsize=8192)
def resolve_allele(name: str, cls: str):
    """Resolve a user-typed allele name to a pseudosequence key for ``cls``.

    Returns ``(key, exact)``. ``exact=True`` when ``name`` (after :func:`normalize_allele`, or the
    locus-aware :func:`class2_from_name` for ``cls=="mhc2"``) is a known key; otherwise the closest key
    by name---a missing ``HLA-`` prefix is repaired and a too-short (e.g. two-field ``'HLA-A02:01'``)
    name is completed by prefix to its first matching key---with ``exact=False``; ``(None, False)`` if
    nothing matches. Serotype names (``'HLA-A2'``) are not expanded. Lets callers accept messy input
    (``'A*02:01'``, ``'HLA-A0201'``) and report when a requested allele is unknown rather than
    silently dropping it.

    Memoised: a miss walks every key and sorts the prefix hits, which is 673 us against
    0.3 us for a hit. Callers reach it once per background peptide inside a calibration
    build, so an unresolvable name used to cost ~6.7 s per allele instead of one lookup.
    """
    seqs = load_pseudo(cls)
    cand = normalize_allele(name.strip())
    # class-I HLA spellings come first, colon form leading, so one molecule always resolves to one
    # key. The table carries both (17,472 keys with the field colon, 1,471 without) on the SAME
    # 34-mer, so without a fixed order `A0201` and `A*02:01` return different names for the same
    # allele -- two calibrators, and two groups that never merge.
    variants = ([class2_from_name(name)] if cls == "mhc2" else hla_spellings(name)) \
        + [cand] + ([] if cand.upper().startswith(("HLA-", "H-2")) else ["HLA-" + cand])
    for v in variants:
        if v in seqs:
            return v, True
    for v in variants:  # prefix completion (two-field -> first four-field key)
        hits = sorted(k for k in seqs if k.startswith(v))
        if hits:
            return hits[0], False
    return None, False


@lru_cache(maxsize=2)
def load_pseudo(cls: str) -> dict:
    """``allele-id -> 34-mer`` for the bundled pseudosequence FASTA of a class.

    Alleles sharing a 34-mer are collapsed to one FASTA record whose header lists **every** such
    allele (``>A B C|n=3``), so all of them are keys here. Listing only the first would silently make
    the rest unscorable -- they are not rare variants: 8,854 of the source table's 12,997 alleles
    (68%) are non-representatives, among them HLA-B*14:02, B*18:05 and C*03:04.
    """
    text = resources.files("mhcmatch.data").joinpath(_FA[cls]).read_text()
    out, names = {}, ()
    for line in text.splitlines():
        if line.startswith(">"):
            names = tuple(line[1:].split("|")[0].split())
        elif names:
            seq = line.strip()
            for n in names:
                out[n] = seq
    return out


def _weighted_hamming(s: str, t: str, w) -> float:
    """Sum of weights at mismatching, non-ambiguous positions (identity metric)."""
    return sum(w[i] for i in range(_LEN)
               if s[i] != t[i] and s[i] != "X" and t[i] != "X")


_AAU = "ACDEFGHIKLMNPQRSTVWY"


@lru_cache(maxsize=1)
def _blosum():
    """seqtree's BLOSUM62 matrix and the mean Gram penalty over distinct AA pairs.

    Lazy (not at import) so docs autodoc can mock ``seqtree``. The mean normalizes the penalty
    so an *average* substitution costs ~1 -- comparable to the identity (Hamming) metric, keeping
    the bandwidth ``h`` and edge thresholds on the same scale across metrics.
    """
    import seqtree

    m = seqtree.SubstitutionMatrix.blosum62()
    n = len(_AAU)
    mean = sum(m.penalty(a, b) for a in _AAU for b in _AAU if a != b) / (n * (n - 1))
    return m, mean


@lru_cache(maxsize=None)
def _pen(a: str, b: str) -> float:
    """Normalized BLOSUM62 Gram-distance penalty between two residues (0 on identity, X skipped)."""
    if a == b or a == "X" or b == "X":
        return 0.0
    m, mean = _blosum()
    return m.penalty(a, b) / mean


def _weighted_blosum(s: str, t: str, w) -> float:
    """Weighted sum of per-position BLOSUM Gram penalties (conservative subs cost less)."""
    return sum(w[i] * _pen(s[i], t[i]) for i in range(_LEN)
               if s[i] != "X" and t[i] != "X")


#: BLOSUM62's own background -- the Blocks pair marginals ``p(i,*)`` of Henikoff & Henikoff's
#: ``blosum62.qij`` (PMID 8743679). The matrix's lambda and this background are jointly determined:
#: ``s_ab = nint(2·log2(q_ab / (p_a·p_b)))`` holds only with *these* frequencies. Deliberately **not**
#: :data:`mhcmatch.diffusion.PROTEOME_AA_FREQ`, which answers a different question (the scoring null).
BLOSUM62_BG = {
    "A": .0742, "R": .0516, "N": .0446, "D": .0536, "C": .0247, "Q": .0343, "E": .0543,
    "G": .0741, "H": .0262, "I": .0679, "L": .0989, "K": .0582, "M": .0250, "F": .0474,
    "P": .0385, "S": .0572, "T": .0509, "W": .0130, "Y": .0323, "V": .0729,
}


def _scale(m, hi: float = 5.0, iters: int = 200):
    """The Karlin-Altschul scale L of a log-odds matrix: the unique L > 0 with
    ``sum_ab p_a p_b exp(L * s_ab) = 1``. ``None`` when no positive root exists.

    **This is why a matrix cannot simply be swapped in.** A log-odds matrix is
    ``s_ab = (1/L) ln(q_ab / (p_a p_b))``, and L is a property of the published table, not a
    constant: measured against ``BLOSUM62_BG`` the matrices seqtree carries come out at BLOSUM62
    0.321 and PAM100 0.332 (half-bit units, ``ln2/2 = 0.347``) but BLOSUM80 0.231, BLOSUM45 0.231
    and PAM250 0.219 (third-bit). Hardcoding ``2^(s/2)`` therefore recovers BLOSUM62 correctly and
    overstates the others' exponent by ~1.4x -- enough to make BLOSUM45 look *more* conservative
    than BLOSUM62, which inverts the very ordering a matrix sweep is asking about.

    ``structural`` has no positive root (its entries are all >= 0, mean +6.6), so it is a
    similarity score and not a log-odds matrix; the identity does not apply to it at all.
    """
    p = BLOSUM62_BG

    def f(L):
        return sum(p[a] * p[b] * math.exp(L * m.similarity(a, b))
                   for a in _AAU for b in _AAU) - 1.0

    if f(hi) < 0:                        # mean score is non-negative: no crossing, not log-odds
        return None
    lo = 1e-9
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) if 0.5 * (lo + hi) > 1e-4 else None


#: Substitution matrices seqtree carries, and therefore the whole menu for the pseudocount blend.
#: Note what is NOT here. There is no BLOSUM90 or BLOSUM100 to try, and seqtree's ``structural``
#: is excluded on purpose: its entries are all non-negative (mean +6.6), so it has no positive
#: Karlin-Altschul scale and is a similarity score rather than a log-odds matrix.
#:
#: The two families run in OPPOSITE directions. BLOSUM-n is built from blocks clustered at >= n%
#: identity, so a HIGHER number means closer relatives and a more conservative substitution model;
#: PAM-n is n accepted point mutations per 100 residues, so a HIGHER number means MORE divergence.
#: Roughly BLOSUM80 ~ PAM120, BLOSUM62 ~ PAM160-200, BLOSUM45 ~ PAM250.
PSEUDO_MATRICES: tuple = ("blosum62", "blosum45", "blosum80", "pam250", "pam100")


@lru_cache(maxsize=len(PSEUDO_MATRICES))
def substitution_conditional(matrix: str = "blosum62") -> dict:
    """``{observed: {r: P(r | observed)}}`` -- a substitution conditional from any seqtree matrix.

    The ``q(a|b)`` of Nielsen et al. 2004 (PMID 14962912), used to spread an anchor's observed residue
    counts onto chemically similar residues (see :meth:`mhcmatch.diffusion.AnchorModel._add_pseudocounts`).

    No ``q_ij`` table and no new dependency are needed. BLOSUM half-bits are
    ``s_ab = 2·log2(q_ab / (p_a·p_b))``, so ``q_ab = p_a·p_b·2^(s_ab/2)`` and

        ``P(a|b) = q_ab / p_b = p_a · 2^(s_ab/2)``   (normalized over ``a``)

    -- only the 20 background frequencies survive. Reads ``.similarity()`` (the raw signed half-bits);
    ``.penalty()`` is the Gram form ``s_aa + s_bb - 2·s_ab``, which forces the diagonal to zero and so
    cannot recover the log-odds.

    ``BLOSUM62_BG`` is used as ``p_a`` for every matrix. That is exact for BLOSUM62 and an
    approximation for the others, whose own background differs; the identity is
    ``P(a|b) ∝ p_a·2^(s_ab/2)``, so a wrong ``p`` tilts the conditional without changing which
    residues it calls similar. Stated because it bounds what a matrix sweep can conclude.
    """
    import seqtree

    if matrix not in PSEUDO_MATRICES:
        raise ValueError(f"unknown substitution matrix {matrix!r}; "
                         f"seqtree carries {', '.join(PSEUDO_MATRICES)}")
    m = getattr(seqtree.SubstitutionMatrix, matrix)()
    lam = _scale(m)
    if lam is None:
        raise ValueError(f"{matrix!r} has no positive Karlin-Altschul scale, so it is not a "
                         "log-odds matrix and P(a|b) cannot be recovered from it")
    out = {}
    for b in _AAU:
        col = {a: BLOSUM62_BG[a] * math.exp(lam * m.similarity(a, b)) for a in _AAU}
        z = sum(col.values())
        out[b] = {a: v / z for a, v in col.items()}
    return out


def mutual_information(xs, ys) -> float:
    """MI(X;Y) in bits for two aligned categorical sequences."""
    n = len(xs)
    if n == 0:
        return 0.0
    px, py, pxy = Counter(xs), Counter(ys), Counter(zip(xs, ys))
    mi = 0.0
    for (x, y), c in pxy.items():
        pj = c / n
        mi += pj * math.log2(pj / ((px[x] / n) * (py[y] / n)))
    return max(mi, 0.0)


def learn_anchor_weights(pseudo_seqs: dict, anchor_residue: dict, prune_dpi: bool = False,
                         tol: float = 0.0) -> list:
    """Per-position relevance ``w[p]`` = MI(groove position ``p`` residue ; anchor residue) across
    alleles, normalized to mean 1. ``anchor_residue``: ``{allele: residue}`` (e.g. the modal residue
    at one peptide anchor for that allele). Positions that discriminate the anchor get more weight.

    Raw MI is inflated by linkage between groove positions (they co-vary across alleles), so many
    positions look relevant and the per-pocket profile is smeared. With ``prune_dpi=True`` an ARACNE
    data-processing-inequality prune removes indirect links: position p's edge to the pocket is
    dropped if some other position q is more informative about the pocket and about p
    (I(p;pocket) <= min(I(q;pocket), I(p;q))), leaving the direct pocket positions sparse and distinct.
    """
    alleles = [a for a in anchor_residue if a in pseudo_seqs and len(pseudo_seqs[a]) == _LEN]
    if not alleles:
        return [1.0] * _LEN
    ys = [anchor_residue[a] for a in alleles]
    cols = [[pseudo_seqs[a][p] for a in alleles] for p in range(_LEN)]
    mi = [mutual_information(cols[p], ys) for p in range(_LEN)]
    w = list(mi)
    if prune_dpi:
        for p in range(_LEN):
            if mi[p] <= 0:
                continue
            for q in range(_LEN):  # q mediates p's link to the pocket -> p is indirect
                if q == p or mi[q] <= mi[p]:
                    continue
                if mi[p] <= mutual_information(cols[p], cols[q]) - tol:
                    w[p] = 0.0
                    break
    mean = sum(w) / _LEN
    return [x / mean for x in w] if mean > 0 else [1.0] * _LEN


class Pseudoseq:
    """Allele-similarity kernel and diffusion over groove pseudosequences for one MHC class."""

    def __init__(self, cls, h=2.0, weights=None, metric="blosum"):
        """``h``: kernel bandwidth. ``weights``: per-position list (one kernel) or
        ``{anchor: [34 weights]}`` (anchor-factored, from :func:`learn_anchor_weights`).
        ``metric``: ``"blosum"`` (default) scores each position by the BLOSUM62 Gram distance
        (conservative substitutions cost less); ``"identity"`` counts plain mismatches."""
        self.cls = cls
        self.seqs = load_pseudo(cls)
        self.h = h
        self.weights = weights
        self.metric = metric

    def _w(self, anchor=None):
        if isinstance(self.weights, dict):
            return self.weights.get(anchor, [1.0] * _LEN)
        return self.weights or [1.0] * _LEN

    def _lookup(self, a):
        s = self.seqs.get(a) or self.seqs.get(normalize_allele(a))
        return s if s and len(s) == _LEN else None

    def kernel(self, a, b, anchor=None) -> float:
        sa, sb = self._lookup(a), self._lookup(b)
        if sa is None or sb is None:
            return 0.0
        dist = _weighted_blosum if self.metric == "blosum" else _weighted_hamming
        return math.exp(-dist(sa, sb, self._w(anchor)) / self.h)

    def neighbors(self, allele, candidates=None, anchor=None, top=10, min_k=0.0):
        """``[(allele, kernel), ...]`` most groove-similar to ``allele`` (self excluded)."""
        cands = candidates if candidates is not None else self.seqs.keys()
        na = normalize_allele(allele)
        scored = [(b, self.kernel(allele, b, anchor)) for b in cands
                  if normalize_allele(b) != na]
        scored = [x for x in scored if x[1] > min_k]
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[:top]

    def cluster(self, alleles, anchor=None, threshold=0.5):
        """Single-linkage clusters: merge alleles with ``kernel >= threshold``. O(n^2); use on a
        panel (~hundreds of alleles), not the full 4k-allele set."""
        al = list(alleles)
        parent = {a: a for a in al}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(len(al)):
            for j in range(i + 1, len(al)):
                if self.kernel(al[i], al[j], anchor) >= threshold:
                    parent[find(al[i])] = find(al[j])
        groups = defaultdict(list)
        for a in al:
            groups[find(a)].append(a)
        return list(groups.values())

    def shrink(self, prefs, allele, anchor=None, candidates=None, prior_strength=None) -> dict:
        """Kernel-weighted empirical-Bayes pooling of a per-anchor residue distribution.

        ``prefs``: ``{allele: Counter(residue -> count)}`` for one anchor. Returns the shrunk
        probability dict for ``allele``.

        With ``prior_strength=None`` (default) this is the counts-weighted form
        ``(n_a π_a + Σ_b K_ab n_b π_b) / (n_a + Σ_b K_ab n_b)`` with limits ``h -> 0`` (raw
        per-allele) and ``h -> ∞`` (global pool). With ``prior_strength=τ`` it uses the
        fixed-concentration form ``(n_a π_a + τ m_a) / (n_a + τ)`` where ``m_a`` is the
        kernel-weighted neighbour mean -- a bounded prior that prevents one large neighbour from
        swamping a rare allele's own peptides and self-adapts to ``n_a`` (appendix §4, Prop. on
        bias--variance). The latter is the recommended default for the forward scorer.
        """
        na = normalize_allele(allele)
        own = Counter(prefs.get(allele, Counter()))
        nbr = Counter()
        cands = candidates if candidates is not None else prefs.keys()
        for b in cands:
            if normalize_allele(b) == na:
                continue
            k = self.kernel(allele, b, anchor)
            if k <= 0:
                continue
            for res, c in prefs.get(b, Counter()).items():
                nbr[res] += k * c

        if prior_strength is None:
            pooled = own + nbr
            total = sum(pooled.values())
            return {res: c / total for res, c in pooled.items()} if total > 0 else {}

        n_own, m = sum(own.values()), sum(nbr.values())
        total = n_own + (prior_strength if m > 0 else 0.0)
        if total <= 0:
            return {}
        pooled = {res: c for res, c in own.items()}
        if m > 0:
            for res, c in nbr.items():
                pooled[res] = pooled.get(res, 0.0) + prior_strength * (c / m)
        return {res: c / total for res, c in pooled.items()}


def blosum62_conditional() -> dict:
    """The BLOSUM62 conditional --- :func:`substitution_conditional` at its default.

    Kept as a name because it is what the anchor model has always called; the matrix is now a
    parameter rather than a constant.
    """
    return substitution_conditional("blosum62")
