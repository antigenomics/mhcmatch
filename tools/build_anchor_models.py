#!/usr/bin/env python
"""Regenerate the vendored pre-fit MHC-II AnchorModels (release-time task).

The MHC-II register + K=3 motif EM is slow to fit on the full corpus (~1-5 min), and a
``mhcmatch predict`` run triggers it twice (the presentation scorer + the affinity register oracle),
so both configs are shipped pre-fit under ``mhcmatch.data`` and loaded read-only by
``Store.anchor_model`` -- no runtime writes, so concurrent pipeline tasks never race on a cache.

Rerun this whenever the vendored models would go stale -- i.e. on a **version bump** (the load guard
keys on ``mhcmatch.__version__``) or when the pmhc panel (``isalgo/pmhc_data`` full tier) changes::

    python tools/build_anchor_models.py

Then commit the regenerated ``src/mhcmatch/data/anchor_model_mhc2_*.pkl.gz`` alongside the bump.
"""
import os

import mhcmatch
from mhcmatch.diffusion import _VENDORED_MODELS, save_vendored_anchor_model

DATA = os.path.join(os.path.dirname(mhcmatch.__file__), "data")


def main():
    store = mhcmatch.Store.from_pmhc(tier="full", species="human", classes=("mhc2",))
    for (cls, footprint, background), name in _VENDORED_MODELS.items():
        path = os.path.join(DATA, name)
        save_vendored_anchor_model(store, cls, path, footprint=footprint, background=background)
        print(f"wrote {path}  ({os.path.getsize(path) / 1e6:.2f} MB)  [{cls} {footprint}/{background}]")


if __name__ == "__main__":
    main()
