"""Whether a unit is actually liberated from the construct it was put in.

Everything in :mod:`mhcmatch.rank` scores a peptide that is *assumed* to exist. A unit inside a
cassette is not assumed to exist: it has to be cut out of a longer translated sequence, and what
cuts it depends on the residues around it --- which, in a designed construct, are its neighbours.
That makes liberation an order-dependent property of the design, and it is the one part of the
antigen-processing pathway that a designer actually controls.

**The asymmetry between the two ends is a fact about the proteasome, not a modelling convenience.**
The proteasome generates the *exact* C-terminus of a class-I epitope and no further trimming
happens there, while the N-terminus is generated with an extension that ER aminopeptidases trim
afterwards --- more than a quarter of epitopes carry an extension of three residues or more, and
close to half carry five or more (Nielsen et al., *Immunogenetics* 2005;57(1-2):33-41). So a
C-terminal cut is *required at exactly one position* and an N-terminal cut is required *somewhere
in a window*. The two enter :func:`liberation` differently for that reason.

**Cleavage is stochastic and the published rescaling says how much.** A NetChop network emits a
score in ``[0, 1]``; Nielsen et al. found that a site predicted above threshold was actually used
in roughly one digest in two, and that multiplying the network output by :data:`SCALE` recovers the
observed per-digest frequency. That factor is what turns a score into a probability, and it is the
only reason the three factors below can be multiplied at all.

Three quantities, all in ``[0, 1]``:

``q_cterm``
    the C-terminal cut is made --- one position, so one probability.
``q_dest``
    **no** internal cut destroys the unit, ``prod(1 - p_i)`` over the positions strictly inside it.
    Nielsen's own benchmark counts a prediction as a false positive when the strongest internal
    site beats the C-terminal one, so internal cleavage is not a nuisance term; it is half of what
    the network was evaluated on.
``q_ntrim``
    a cut is made *somewhere* in the trimmable window upstream, ``1 - prod(1 - p_j)``.

Their product is the liberation probability, and it multiplies a unit's response probability rather
than being added to it --- a unit that is never cut out is not rescued by being well recognised.

The same shape appears in JessEV (Dorigatti and Schubert, *Bioinformatics* 2020;36:i643-i650) as
``R_e = C_Nt * C_Ct * prod(1 - C_interior)``, reached by Monte-Carlo simulation over a cleavage
process. What is different here is that the factors are calibrated probabilities via :data:`SCALE`
rather than raw scores, and that they multiply a *fitted* response probability rather than a
predictor score.

**Nothing in this module runs NetChop.** It consumes a parsed table, so it is testable with no
binary, no download and no cluster --- which matters, because NetChop is academic-agreement
software that cannot be vendored. :func:`parse` reads what the ``netchop`` wrapper prints.
"""
from __future__ import annotations

import math
import re

__all__ = ["SCALE", "TRIM_WINDOW", "parse", "probability", "liberation", "unit_geometry"]

#: Network output -> probability that the site is used in a given digest, from the rescaling
#: Nielsen et al. 2005 fitted to in-vitro digests: a predicted site was observed in one digest in
#: two for the MHC-ligand-trained network and three in five for the 20S-trained one. Keyed by the
#: ``netchop -v`` network: ``0`` is ``Cterm.3.0``, ``1`` is ``20S.3.0``.
SCALE = {"Cterm": 0.5, "20S": 0.6}

#: How far upstream of a unit's first residue a cut still liberates it, because ER aminopeptidases
#: trim the N-terminal extension afterwards. Nielsen et al. measured extensions of three or more
#: residues on more than a quarter of epitopes and five or more on close to half, so a window of
#: seven covers the great majority without pretending the tail does not exist. It is a parameter,
#: not a constant, because the distribution is a property of the cell and not of the peptide.
TRIM_WINDOW = 7

_ROW = re.compile(r"^\s*(\d+)\s+([A-Z])\s+([S.])\s+([0-9.eE+-]+)\s+(\S+)\s*$")


def parse(text: str) -> dict:
    """``netchop`` output -> ``{identifier: [score per residue]}``, in sequence order.

    Reads the long format (the default). Comment lines, the banner and the rules are skipped by
    not matching, rather than by counting lines --- the header is nine lines on one host and a
    different number on another, and a reader that counts is a reader that breaks quietly.

    The score at position *i* is the probability of a cut **after** residue *i*, which is why
    :func:`liberation` reads a unit's C-terminal cut at its *last* residue and not past it.

    >>> parse("   1   M  S   0.760600 gi|333\\n   2   A  .   0.483380 gi|333\\n")
    {'gi|333': [0.7606, 0.48338]}
    """
    out: dict = {}
    for line in text.splitlines():
        m = _ROW.match(line)
        if not m:
            continue
        pos, _aa, _call, score, ident = m.groups()
        seq = out.setdefault(ident, [])
        if int(pos) != len(seq) + 1:
            raise ValueError(
                f"processing.parse: {ident!r} jumps to position {pos} with {len(seq)} scores read. "
                "NetChop numbers each record from 1, so a gap means two records share an "
                "identifier and their scores are being concatenated into one sequence.")
        seq.append(float(score))
    return out


def probability(score: float, network: str = "Cterm") -> float:
    """A network score as a per-digest cleavage probability: ``SCALE[network] * score``.

    >>> round(probability(1.0), 4)
    0.5
    >>> round(probability(0.8, "20S"), 4)
    0.48
    """
    if network not in SCALE:
        raise ValueError(f"network must be one of {sorted(SCALE)}, got {network!r}")
    s = float(score)
    if not 0.0 <= s <= 1.0:
        raise ValueError(
            f"processing.probability: a NetChop score is in [0, 1], got {s:.6g}. A value outside "
            "that range is a parsing error upstream, and rescaling it here would hide it.")
    return SCALE[network] * s


def liberation(scores, start: int, length: int, *, network: str = "Cterm",
               trim_window: int = TRIM_WINDOW) -> dict:
    """Liberation of the unit at ``scores[start : start + length]`` from its construct.

    ``scores`` is one NetChop score per residue of the **whole construct**, so the same unit in two
    different neighbourhoods gets two different answers --- which is the entire point.

    ``start`` is 0-based. The unit's C-terminal cut is the score at its last residue; the internal
    positions are the ones strictly before that; the trimmable window is the ``trim_window``
    residues before ``start``, truncated at the construct's own N-terminus, where no cut is needed
    because translation already provides one.

    Returns ``q_cterm``, ``q_dest``, ``q_ntrim``, their product ``q``, and ``max_internal`` --- the
    strongest internal site, which is the quantity Nielsen's benchmark calls a false positive when
    it beats the C-terminal one and is worth reporting beside the product.

    >>> s = [0.1] * 10
    >>> s[6] = 0.9                                   # a clean cut at the unit's C-terminus
    >>> r = liberation(s, start=2, length=5)
    >>> round(r["q_cterm"], 4)
    0.45
    """
    n = len(scores)
    if length <= 0:
        raise ValueError(f"liberation: length must be positive, got {length}")
    if start < 0 or start + length > n:
        raise ValueError(
            f"liberation: unit [{start}, {start + length}) does not fit a construct of {n} "
            "residues. The scores must be for the whole construct, not for the unit alone.")
    p = [probability(x, network) for x in scores]
    end = start + length - 1                                   # index of the last residue

    q_cterm = p[end]
    internal = p[start:end]                                    # strictly inside, C-terminus excluded
    q_dest = math.prod(1.0 - x for x in internal) if internal else 1.0

    lo = max(0, start - int(trim_window))
    window = p[lo:start]
    # At the construct's own N-terminus no cut is needed: translation starts there. Treating that
    # as "no cleavage site found" would score the first unit of every cassette as unliberatable.
    q_ntrim = 1.0 if start == 0 else 1.0 - math.prod(1.0 - x for x in window)

    return {"q_cterm": q_cterm, "q_dest": q_dest, "q_ntrim": q_ntrim,
            "q": q_cterm * q_dest * q_ntrim,
            "max_internal": max(internal) if internal else 0.0}


def unit_geometry(unit_length: int, epitope_start: int, epitope_length: int) -> dict:
    """Where the minimal epitope sits inside its unit: the length terms, in residues.

    ``epitope_start`` is 0-based within the unit. A designer chooses the unit, not the epitope, so
    these are design variables: how long the unit is, and how much sequence sits either side of the
    epitope it was built around.

    >>> unit_geometry(27, 9, 9)
    {'len_unit': 27, 'len_epitope': 9, 'len_flankN': 9, 'len_flankC': 9}
    """
    if epitope_start < 0 or epitope_length <= 0 or epitope_start + epitope_length > unit_length:
        raise ValueError(
            f"unit_geometry: epitope [{epitope_start}, {epitope_start + epitope_length}) does not "
            f"fit a unit of {unit_length} residues.")
    return {"len_unit": int(unit_length), "len_epitope": int(epitope_length),
            "len_flankN": int(epitope_start),
            "len_flankC": int(unit_length - epitope_start - epitope_length)}
