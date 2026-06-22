#!/usr/bin/env python3
"""MHC allele similarity / promiscuity network (Graphviz) with detected communities.

Nodes are the alleles present in ``pmhc_data`` for one (class, species). Edges are
pseudosequence-kernel similarity -- the same diffusion kernel of :mod:`mhcmatch.pseudoseq` -- drawn
thin/grey above a SOFT threshold and bold above a HARD threshold. **Communities** (functional
supertypes, the structural basis of cross-allele promiscuity) are detected by greedy modularity
(networkx) and drawn as outlined Graphviz clusters; node size grows with the allele's presented-set
size. One network per class x species (MHC-I/II, human/mouse).

    python bench/promiscuity_graph.py --pmhc-dir /path/to/pmhc_data --out appendix

Needs: ``networkx`` (community detection) and Graphviz (``dot``/``fdp`` on PATH).
    pip install networkx          # the `graphviz` system package provides dot/fdp
"""
from __future__ import annotations

import argparse
import itertools
import math
import os
import subprocess
from collections import Counter, defaultdict

import networkx as nx

from mhcmatch import Pseudoseq, Store
from mhcmatch.diffusion import MHC1_ANCHORS, MHC2_ANCHORS
from mhcmatch.pseudoseq import learn_anchor_weights, load_pseudo, normalize_allele

# Colour-blind-friendly community palette (Graphviz fill = light tint, outline = full colour).
PALETTE = ["#2563eb", "#f59e0b", "#10b981", "#ef4444", "#8b5cf6", "#ec4899",
           "#0891b2", "#65a30d", "#9333ea", "#dc2626", "#0d9488", "#ca8a04"]
TINT = {c: c + "22" for c in PALETTE}   # Graphviz supports #RRGGBBAA fills


def _mi_weights(store, cls):
    """Per-position max mutual-information weight over the core anchors (max-MI positions)."""
    seqs = load_pseudo(cls)
    anchors = MHC1_ANCHORS if cls == "mhc1" else MHC2_ANCHORS
    vecs = []
    for j in anchors:
        prefs = store.anchor_preferences(cls, j)
        modal = {normalize_allele(a): c.most_common(1)[0][0] for a, c in prefs.items() if c}
        vecs.append(learn_anchor_weights(seqs, modal))
    combined = [max(v[p] for v in vecs) for p in range(len(vecs[0]))]
    mean = sum(combined) / len(combined)
    return [x / mean for x in combined] if mean > 0 else None


def build_graph(store, cls, soft, hard, h, metric="kernel"):
    """Allele graph. ``metric="kernel"``: predicted groove similarity (MI-weighted BLOSUM kernel).
    ``metric="shared"``: observed co-presentation -- overlap coefficient of presented-peptide sets,
    ``|P_a ∩ P_b| / min(|P_a|,|P_b|)`` (edge also carries the raw shared-epitope count)."""
    npep = Counter(store._panel[cls].alleles)
    G = nx.Graph()
    if metric == "shared":
        sets = defaultdict(set)
        for ep, a in zip(store._panel[cls].epitopes, store._panel[cls].alleles):
            sets[a].add(ep)
        nodes = [a for a in store.alleles(cls) if sets[a]]
        for a in nodes:
            G.add_node(a, n=npep[a])
        for a, b in itertools.combinations(nodes, 2):
            inter = len(sets[a] & sets[b])
            if inter == 0:
                continue
            ov = inter / min(len(sets[a]), len(sets[b]))
            if ov >= soft:
                G.add_edge(a, b, weight=round(ov, 3), hard=ov >= hard, shared=inter)
        return G
    # predicted groove similarity
    ps = Pseudoseq(cls, h=h, weights=_mi_weights(store, cls), metric="blosum")
    nodes = [a for a in store.alleles(cls) if ps._lookup(a) is not None]
    for a in nodes:
        G.add_node(a, n=npep[a])
    for a, b in itertools.combinations(nodes, 2):
        k = ps.kernel(a, b)
        if k >= soft:
            G.add_edge(a, b, weight=round(k, 3), hard=k >= hard)
    return G


def detect_communities(G):
    if G.number_of_edges() == 0:
        return [{n} for n in G.nodes]
    return list(nx.community.greedy_modularity_communities(G, weight="weight"))


def to_dot(G, comms, title):
    comm_of = {n: i for i, c in enumerate(comms) for n in c}
    out = ["graph G {",
           '  graph [layout=fdp, overlap=false, splines=true, fontname="Helvetica", fontsize=18];',
           f'  labelloc="t"; label="{title}";',
           '  node [shape=circle, style=filled, fontname="Helvetica", fontsize=11, penwidth=1.4];',
           '  edge [color="#d1d5db"];']
    # one Graphviz cluster per (multi-node) community -> dashed outline around the supertype
    for i, c in enumerate(comms):
        col = PALETTE[i % len(PALETTE)]
        if len(c) >= 2:
            # label="" so the cluster does not inherit (and repeat) the graph title.
            out.append(f'  subgraph cluster_{i} {{ label=""; style="rounded,dashed"; color="{col}";')
        for n in sorted(c):
            # log-scaled node size (presented-set size spans 1 to ~1e5): readable, capped.
            sz = round(min(1.5, 0.30 + 0.22 * math.log10(1 + G.nodes[n]["n"])), 2)
            label = n.replace("HLA-", "").replace("H-2-", "")
            out.append(f'    "{n}" [label="{label}", fillcolor="{TINT[col]}", color="{col}", '
                       f'width={sz}];')
        if len(c) >= 2:
            out.append("  }")
    for u, v, d in G.edges(data=True):
        if d["hard"]:
            out.append(f'  "{u}" -- "{v}" [penwidth=2.0, color="#6b7280"];')
        else:
            out.append(f'  "{u}" -- "{v}" [penwidth=0.4];')
    out.append("}")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--tier", default="full", choices=("full", "shortlist"))
    ap.add_argument("--out", default="appendix")
    ap.add_argument("--soft", type=float, default=0.25)
    ap.add_argument("--hard", type=float, default=0.5)
    ap.add_argument("--h", type=float, default=2.0)
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz")
    descs = {"kernel": "groove kernel", "shared": "shared-epitope overlap"}

    for cls, clabel in (("mhc1", "MHC-I"), ("mhc2", "MHC-II")):
        for species in ("human", "mouse"):
            store = Store.from_pmhc(path, tier=args.tier, species=species, classes=(cls,))
            graphs = {}
            for metric in ("kernel", "shared"):
                G = build_graph(store, cls, args.soft, args.hard, args.h, metric)
                graphs[metric] = G
                if G.number_of_nodes() < 2:
                    print(f"# {clabel} {species} ({metric}): <2 matched alleles, skipped")
                    continue
                comms = detect_communities(G)
                ncc = sum(1 for c in comms if len(c) >= 2)
                title = (f"{clabel} {species}: {G.number_of_nodes()} alleles, {ncc} communities "
                         f"({descs[metric]} soft {args.soft} / hard {args.hard})")
                stem = f"promiscuity_{'shared_' if metric == 'shared' else ''}{cls}_{species}"
                (open(os.path.join(out, stem + ".dot"), "w")).write(to_dot(G, comms, title))
                subprocess.run(["dot", "-Kfdp", "-Tpdf", "-o", os.path.join(out, stem + ".pdf"),
                                os.path.join(out, stem + ".dot")], check=True)
                print(f"# {clabel} {species} ({metric}): {G.number_of_nodes()} alleles, "
                      f"{G.number_of_edges()} edges, {ncc} communities -> {stem}.pdf")
            # Jaccard agreement of edge sets between predicted (kernel) and observed (shared),
            # restricted to alleles present in both graphs.
            common = set(graphs["kernel"].nodes()) & set(graphs["shared"].nodes())
            ek = {frozenset(e) for e in graphs["kernel"].edges() if set(e) <= common}
            es = {frozenset(e) for e in graphs["shared"].edges() if set(e) <= common}
            if ek or es:
                j = len(ek & es) / len(ek | es)
                print(f"#   Jaccard(kernel, shared edges | {len(common)} common alleles) = "
                      f"{j:.3f}  ({len(ek & es)}/{len(ek | es)})")


if __name__ == "__main__":
    main()
