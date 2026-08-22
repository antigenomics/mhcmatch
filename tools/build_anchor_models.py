#!/usr/bin/env python
"""Shim. The builder moved into the package so it is reachable from any install:

    mhcmatch build anchor

See ``mhcmatch._build.anchor_models`` for what it does and ``src/mhcmatch/data/PROVENANCE.md``
for why these are vendored. Kept because PROVENANCE and older notes name this path.
"""
from mhcmatch._build import anchor_models

if __name__ == "__main__":
    anchor_models()
