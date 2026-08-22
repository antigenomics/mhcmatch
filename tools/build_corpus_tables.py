#!/usr/bin/env python
"""Shim. The builder moved into the package so it is reachable from any install:

    mhcmatch build corpus

See ``mhcmatch._build.corpus_tables`` for what it does and ``src/mhcmatch/data/PROVENANCE.md``
for why these are vendored. Kept because PROVENANCE and older notes name this path.
"""
from mhcmatch._build import corpus_tables

if __name__ == "__main__":
    corpus_tables()
