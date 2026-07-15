#!/usr/bin/env python3
# 2026-07-15 — fit the FULL Potts affinity model on all measured nM and vendor its weights.
"""Fits the shipped Potts affinity model (one-hot, no held-out split) on every measured '=' nM row of
a class and writes the weight vector to ``src/mhcmatch/data/affinity_potts_<cls>.npz`` for the runtime
predictor in :mod:`mhcmatch.affinity`. Reuses the benchmarked feature construction from
``train_potts`` verbatim, so the vendored model is exactly the one the benchmark evaluated.

    conda run -n tcren-nb python bench/affinity/fit_potts.py --cls mhc1
    conda run -n tcren-nb python bench/affinity/fit_potts.py --cls mhc2
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from sklearn.linear_model import Ridge

sys.path.insert(0, os.path.dirname(__file__))
import train_potts as T                                       # noqa: E402
from mhcmatch.affinity import ic50_to_y                       # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "..", "..", "src", "mhcmatch", "data")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", required=True, choices=("mhc1", "mhc2"))
    ap.add_argument("--alpha", type=float, default=40.0)
    args = ap.parse_args()
    T.configure(args.cls)
    T.set_soft(1.0, 1)                                        # one-hot = the shipped encoding
    if args.cls == "mhc2":
        from mhcmatch import Store
        st = Store.from_pmhc(T.PMHC, tier="full", species=("human", "mouse"), classes=("mhc2",))
        T._AM = st.anchor_model("mhc2", background="proteome", footprint="core")

    points = T.load_points(args.cls)
    cols, y = [], []
    for (pep, key), nm in points.items():
        c = T._cols(T._core(args.cls, pep, key), key)
        if c is not None:
            cols.append(c)
            y.append(ic50_to_y(nm))
    X = T._matrix(cols)
    print(f"# fitting {args.cls} Potts on {X.shape[0]} pts x {T.NFEAT} params, alpha={args.alpha}",
          flush=True)
    model = Ridge(alpha=args.alpha, solver="lsqr", max_iter=1000).fit(X, np.asarray(y))

    w = model.coef_.astype(np.float32)
    nz = int((np.abs(w) > 1e-6).sum())
    path = os.path.join(DATA, f"affinity_potts_{args.cls}.npz")
    np.savez_compressed(path, w=w, b=np.float32(model.intercept_),
                        meta=np.array([T.PEPP, T.PSP, T.Q, int(args.alpha)], dtype=np.int32))
    print(f"# wrote {path}  ({nz}/{T.NFEAT} nonzero weights, {os.path.getsize(path)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
