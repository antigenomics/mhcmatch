"""Near-exact source-peptide lookup against a reference proteome.

Given a query peptide (e.g. a neoantigen), find the nearly-exact self peptide it derives from and
its parent protein / position via **full-sequence** (unmasked) ``<= max_subs`` search over all
windows of the proteome of the query's length -- using the seqtree Hamming fast path. This is a
*distinct* mode from the anchor-masked TCR-facing homology and the presentation-signature searches.
See ``appendix/mhcmatch.tex`` §5 (near-exact source identification).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass

from seqtree import Index, SearchParams

_AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
_AA = set(_AA_ORDER)


def read_fasta(path):
    """``{name: sequence}`` from a (optionally gzipped) FASTA; name = first whitespace token."""
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


@dataclass
class SourceHit:
    protein: str
    position: int       # 0-based start in the protein
    ref_peptide: str
    n_subs: int
    mutations: tuple    # ((pos_in_peptide, query_aa, ref_aa), ...)


class Proteome:
    """A reference proteome with lazily-built per-length window indices."""

    def __init__(self, seqs):
        self.seqs = seqs
        self._cache = {}   # length -> (Index | None, [(protein, pos, window), ...])

    @classmethod
    def from_fasta(cls, path):
        return cls(read_fasta(path))

    @classmethod
    def from_hf(cls, name="human"):
        """Load a reference proteome by name, auto-fetched from the public HF dataset (no manual
        download). ``name`` = ``"human"`` / ``"mouse"`` (UP000005640 / UP000000589) or a pathogen
        stem; see :func:`mhcmatch.store.fetch_proteome`."""
        from .store import fetch_proteome
        return cls.from_fasta(fetch_proteome(name))

    def _index(self, L):
        if L not in self._cache:
            windows, meta = [], []
            for name, seq in self.seqs.items():
                s = seq.upper()
                for i in range(len(s) - L + 1):
                    w = s[i:i + L]
                    if all(c in _AA for c in w):
                        windows.append(w)
                        meta.append((name, i, w))
            self._cache[L] = (Index.build(windows, alphabet="aa") if windows else None, meta)
        return self._cache[L]

    def find_source(self, peptide, max_subs=1, exclude_exact=False):
        """Self peptides within ``max_subs`` substitutions of ``peptide``, nearest first.

        Returns ``[SourceHit, ...]``. ``exclude_exact=True`` drops perfect (0-mismatch) matches --
        useful to find the wild-type a mutated neoantigen derives from when the query is itself self.
        """
        q = peptide.strip().upper()
        idx, meta = self._index(len(q))
        if idx is None:
            return []
        p = SearchParams(max_subs=max_subs, engine="seqtm")
        out = []
        for hit in idx.search(q, p):
            name, pos, w = meta[hit.ref_id]
            muts = tuple((i, q[i], w[i]) for i in range(len(q)) if q[i] != w[i])
            if exclude_exact and not muts:
                continue
            out.append(SourceHit(name, pos, w, len(muts), muts))
        out.sort(key=lambda h: h.n_subs)
        return out

    def _window_set(self, L):
        """Set of all length-``L`` standard-AA proteome windows (lazy). ~1 GB/length as a Python set --
        much lighter than the seqtree index, and O(1) membership for the 1-sub wildtype fast path."""
        key = ("set", L)
        if key not in self._cache:
            s = set()
            for seq in self.seqs.values():
                seq = seq.upper()
                for i in range(len(seq) - L + 1):
                    w = seq[i:i + L]
                    if all(c in _AA for c in w):
                        s.add(w)
            self._cache[key] = s
        return self._cache[key]

    def wildtype(self, peptide, max_subs=1):
        """The wild-type self peptide a mutated ``peptide`` derives from, or ``None``.

        A self peptide exactly one substitution away (its point-mutation origin) -- the position-aligned
        WT counterpart needed for agretopicity / DAI when the caller has no WT window (e.g. a bare
        neoantigen list like TESLA). ``None`` when nothing is one sub away (indel / spliced / non-self,
        or the peptide is itself an exact self peptide with no mutated origin). Ties resolve to the
        first variant found (position, then residue order).

        For ``max_subs=1`` this uses a hash-set fast path (generate the L*19 single-sub variants and
        test proteome membership -- microseconds/peptide, so it scales to large corpora); larger
        ``max_subs`` falls back to the general :meth:`find_source` fuzzy search.
        """
        q = peptide.strip().upper()
        if max_subs == 1 and all(c in _AA for c in q):
            ws = self._window_set(len(q))
            for i in range(len(q)):
                pre, post = q[:i], q[i + 1:]
                for a in _AA_ORDER:
                    if a != q[i]:
                        v = pre + a + post
                        if v in ws:
                            return v
            return None
        hits = self.find_source(peptide, max_subs=max_subs, exclude_exact=True)
        return hits[0].ref_peptide if hits else None
