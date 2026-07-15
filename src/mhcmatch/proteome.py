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

_AA = set("ACDEFGHIKLMNPQRSTVWY")


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
