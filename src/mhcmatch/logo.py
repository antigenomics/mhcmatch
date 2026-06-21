"""Per-allele motif logos (information content) + length distributions.

:func:`motif` returns the numeric logo (PWM, per-position bits, length histogram) -- pure-Python,
always available. :func:`render` draws it with ``logomaker`` (optional ``[logo]`` extra). MHC-I
logos use peptides of a fixed length (default the modal length); MHC-II uses register-anchored
9-mer cores. See ``appendix/mhcmatch.tex`` §6.
"""
from __future__ import annotations

import math
from collections import Counter

from seqtree import layout

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _mhc2_core(peptide):
    """The register-anchored 9-mer core window (one-pass register trick), or None."""
    if len(peptide) < 9:
        return None
    return max((peptide[s:s + 9] for s in range(len(peptide) - 8)),
               key=layout._core_anchor_score)


def motif(store, allele, cls, length=None):
    """Logo data for ``allele``'s presented peptides.

    Returns ``{allele, cls, width, n, pwm, bits, length_hist}`` where ``pwm[i]`` is a residue->freq
    dict (sums to 1), ``bits[i]`` the information content (``log2(20) - entropy``) in [0, log2(20)],
    and ``length_hist`` a length->count dict over all the allele's peptides.
    """
    panel = store._panel[cls]
    peps = [e for e, a in zip(panel.epitopes, panel.alleles) if a == allele]
    if not peps:
        raise ValueError(f"no peptides for allele {allele!r} in class {cls}")
    length_hist = Counter(len(p) for p in peps)
    if cls == "mhc2":
        frames = [c for c in map(_mhc2_core, peps) if c]
    else:
        L = length or length_hist.most_common(1)[0][0]
        frames = [p for p in peps if len(p) == L]
    if not frames:
        raise ValueError("no fixed-length frames to build a logo")
    width = len(frames[0])
    pwm, bits = [], []
    for i in range(width):
        col = Counter(f[i] for f in frames if f[i] in _AA)
        n = sum(col.values()) or 1
        freqs = {aa: col.get(aa, 0) / n for aa in _AA}
        entropy = -sum(p * math.log2(p) for p in freqs.values() if p > 0)
        pwm.append(freqs)
        bits.append(math.log2(20) - entropy)
    return {"allele": allele, "cls": cls, "width": width, "n": len(frames),
            "pwm": pwm, "bits": bits, "length_hist": dict(length_hist)}


def render(m, ax=None):
    """Render a :func:`motif` result as an information-content sequence logo (needs ``[logo]``)."""
    import logomaker
    import pandas as pd

    df = pd.DataFrame([{aa: m["pwm"][i][aa] * m["bits"][i] for aa in _AA}
                       for i in range(m["width"])])
    return logomaker.Logo(df, ax=ax)
