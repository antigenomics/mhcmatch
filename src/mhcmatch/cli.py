"""Command-line interface for mhcmatch: ``mhcmatch <command> ...``.

Commands: ``decompose`` (no data needed), ``restriction``, ``scan``, ``logo`` (need a pmhc_data
table via ``--pmhc`` or ``$MHCMATCH_PMHC``), ``source`` (needs a proteome FASTA), and ``span``
(core -> full presented ligand; the panel is optional, and only supplies the observed-ligand tier).
"""
from __future__ import annotations

import argparse
import os

from . import Proteome, Store


def _add_store_opts(p):
    p.add_argument("--pmhc", help="pmhc_data TSV(.gz); else $MHCMATCH_PMHC/pmhc_<tier>.tsv.gz, "
                                  "else auto-fetched from the public HF dataset isalgo/pmhc_data")
    p.add_argument("--tier", default="full", choices=("full", "shortlist"))
    p.add_argument("--species", default="human", choices=("human", "mouse"))


def _store(a):
    return Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species)


def _read_seq(arg):
    """A raw sequence, or the concatenated sequences of a FASTA file path."""
    if os.path.exists(arg):
        from .proteome import read_fasta
        seqs = read_fasta(arg)
        if seqs:
            return "".join(seqs.values())
    return arg.strip()


def cmd_decompose(a):
    d = Store().decompose(a.peptide, cls=a.cls)
    print(f"peptide       {d.peptide}")
    print(f"anchors       {','.join(str(i + 1) for i in d.anchors)}")
    print(f"tcr_facing    {d.tcr_facing}")
    print(f"presentation  {d.presentation}")


def _resolve_panel_allele(store, name, cls):
    """Map a user-typed allele to a panel allele (exact, else prefix on punctuation-stripped name)."""
    pool = (cls,) if cls else ("mhc1", "mhc2")
    panel = {al for c in pool for al in store.alleles(c)}
    if name in panel:
        return name
    key = name.replace("*", "").replace(":", "")
    hits = sorted(a for a in panel if a.replace("*", "").replace(":", "").startswith(key))
    if hits:
        print(f"# resolved '{name}' -> '{hits[0]}'")
        return hits[0]
    print(f"# allele '{name}' not found in panel")
    return name


def cmd_restriction(a):
    store = _store(a)
    allele = _resolve_panel_allele(store, a.allele, a.cls) if a.allele else None
    res = store.restriction(a.peptide, cls=a.cls, alleles=[allele] if allele else "all",
                            top=a.top, diffuse=a.diffuse, calibrated=a.calibrated)
    if not res:
        print("no presenting allele (no presentation-signature neighbours)")
        return
    diffuse = a.diffuse or a.calibrated
    hdr = f"{'allele':<18}{'vote':>7}{'enr':>7}" + ("{:>8}".format("score") if diffuse else "")
    if a.calibrated:
        hdr += f"{'%rank':>8}{'P':>7}{'band':>12}"
    print(hdr + f"{'binder':>8}")
    for r in res:
        line = f"{r.allele:<18}{r.vote:>7.2f}{r.enrichment:>7.1f}"
        if diffuse:
            line += f"{(r.anchor_score or 0.0):>8.2f}"
        if a.calibrated:
            line += f"{r.rank:>8.2f}{r.p_present:>7.2f}{r.band:>12}"
        print(line + f"{'yes' if r.binder else 'no':>8}")


def cmd_affinity(a):
    store = _store(a)
    allele = _resolve_panel_allele(store, a.allele, a.cls)
    am = store.affinity_model(a.cls)
    nm = am.predict_ic50(a.peptide, allele)
    print(f"{a.peptide}  {allele}  predicted IC50 ~ {nm:,.0f} nM")
    if a.wt:
        nm_wt = am.predict_ic50(a.wt, allele)
        print(f"  WT {a.wt}: IC50 ~ {nm_wt:,.0f} nM   amplitude A=Kd_WT/Kd_MT = "
              f"{am.amplitude(a.wt, a.peptide, allele):.2f}   DAI = {am.dai(a.wt, a.peptide, allele):+.2f}")
    if a.structure:
        try:
            from .structure import StructureScorer
            sc = StructureScorer(pseudoseq=store.anchor_model(a.cls).ps)
            mj = sc.mj_energy(a.peptide, allele)
            if mj == mj:
                extra = f"   ΔΔG(WT→MT) = {sc.ddg(a.wt, a.peptide, allele):+.2f}" if a.wt else ""
                print(f"  structural MJ energy = {mj:.2f}{extra}")
            else:
                print("  (no structural template for this allele/length)")
        except ImportError as e:
            print(f"  (structure scoring unavailable: {e})")


def cmd_scan(a):
    hits = _store(a).scan_protein(_read_seq(a.protein), cls=a.cls or "mhc1",
                                  alleles=[a.allele] if a.allele else "all", top=a.top,
                                  correction=a.correction)
    label = f" ({a.correction} FWER/FDR)" if a.correction else ""
    print(f"# {len(hits)} presented window(s){label}")
    for pos, pep, binders in hits:
        print(f"{pos:>5}  {pep:<14}  {','.join(b.allele for b in binders)}")


def cmd_source(a):
    hits = Proteome.from_fasta(a.proteome).find_source(a.peptide, max_subs=a.max_subs)
    if not hits:
        print("# no source within max_subs")
        return
    for h in hits:
        muts = ",".join(f"{q}{i + 1}{r}" for i, q, r in h.mutations) or "exact"
        print(f"{h.protein}\tpos {h.position}\tsubs {h.n_subs}\t{h.ref_peptide}\t{muts}")


def cmd_span(a):
    from . import ligand
    prot = _read_seq(a.protein)
    corpus = None
    if a.pmhc or os.environ.get("MHCMATCH_PMHC"):
        corpus = _store(a)._panel["mhc2"].epitopes
    sp = ligand.presented_span(a.core.strip().upper(), prot, corpus=corpus, mode=a.mode,
                               flanks=tuple(int(x) for x in a.flanks.split(",")))
    if sp is None:
        print("# no reference ligand contains this core (mode=observed)")
        return
    nl, nr = sp.flanks
    print(f"tier      {sp.source}")
    print(f"core      {sp.core} @ {sp.core_start}")
    print(f"peptide   {sp.peptide}  [{sp.start}:{sp.end}]  len {len(sp.peptide)}")
    print(f"flanks    {nl} / {nr}" + (f"   clipped {sp.clipped}" if any(sp.clipped) else ""))
    print(f"score     {sp.score:+.2f}")
    print(f"alts      {sp.n_alternatives}" + (f"   support {sp.support}" if sp.support else ""))


def cmd_predict(a):
    from . import predict as P
    store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
    alleles = [x.strip() for x in a.alleles.split(",") if x.strip()]
    preds = P.predict_fasta(store, a.cls, a.fasta, alleles, rank_threshold=a.rank_threshold,
                            top=a.top, background=a.background, footprint=a.footprint, seed=a.seed)
    if a.native:
        P.write_native(preds, a.native)
        print(f"# wrote {a.native}")
    if a.scored_csv:
        P.write_scored_csv(preds, a.scored_csv)
        print(f"# wrote {a.scored_csv}")
    if not a.native and not a.scored_csv:
        print(f"# {len(preds)} predicted binder(s) (%rank <= {a.rank_threshold}) over "
              f"{len(alleles)} allele(s)")
        for p in preds[:(a.top or 20)]:
            print(f"{p.peptide:<15} {p.allele:<18} %rank={p.percent_rank:<6} {p.band:<11} "
                  f"{p.var.get('gene_name', '')}")


def cmd_logo(a):
    from . import logo
    m = logo.motif(_store(a), a.allele, a.cls or "mhc1")
    print(f"# {a.allele}  width={m['width']}  n={m['n']}  lengths={dict(sorted(m['length_hist'].items()))}")
    for i, (bits, col) in enumerate(zip(m["bits"], m["pwm"]), 1):
        top = sorted(col.items(), key=lambda x: -x[1])[:3]
        print(f"  pos {i:>2}  {bits:4.2f} bits  " + " ".join(f"{aa}:{p:.2f}" for aa, p in top))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mhcmatch", description="peptide-MHC presentation tools")
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("decompose", help="split a peptide into anchor / TCR-facing parts (X masks)")
    d.add_argument("peptide")
    d.add_argument("--cls", choices=("mhc1", "mhc2"))
    d.set_defaults(fn=cmd_decompose)

    r = sub.add_parser("restriction", help="rank presenting alleles for a peptide")
    r.add_argument("peptide")
    r.add_argument("--allele", help="restrict to a single allele")
    r.add_argument("--cls", choices=("mhc1", "mhc2"))
    r.add_argument("--diffuse", action="store_true", help="rare-allele-aware (diffusion-shrunk anchors)")
    r.add_argument("--calibrated", action="store_true",
                   help="add per-allele %%rank, P(present), and binding band (implies --diffuse)")
    r.add_argument("--top", type=int, default=10)
    _add_store_opts(r)
    r.set_defaults(fn=cmd_restriction)

    af = sub.add_parser("affinity", help="predict IC50 (nM) + neoantigen amplitude/DAI for a peptide")
    af.add_argument("peptide")
    af.add_argument("--allele", required=True)
    af.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    af.add_argument("--wt", help="wild-type peptide -> also report amplitude A=Kd_WT/Kd_MT and DAI")
    af.add_argument("--structure", action="store_true",
                    help="also compute the tcren MJ contact energy / ΔΔG (needs the [structure] extra)")
    _add_store_opts(af)
    af.set_defaults(fn=cmd_affinity)

    s = sub.add_parser("scan", help="find presented peptides in a protein (sequence or FASTA path)")
    s.add_argument("protein")
    s.add_argument("--allele")
    s.add_argument("--cls", choices=("mhc1", "mhc2"))
    s.add_argument("--top", type=int, default=3)
    s.add_argument("--correction", choices=("bonferroni", "bh"),
                   help="multiple-testing control over windows x alleles (FWER / BH-FDR)")
    _add_store_opts(s)
    s.set_defaults(fn=cmd_scan)

    so = sub.add_parser("source", help="find the self peptide a neoantigen derives from")
    so.add_argument("peptide")
    so.add_argument("--proteome", required=True, help="reference proteome FASTA(.gz)")
    so.add_argument("--max-subs", type=int, default=1)
    so.set_defaults(fn=cmd_source)

    lg = sub.add_parser("logo", help="motif logo (information content) + length distribution")
    lg.add_argument("allele")
    lg.add_argument("--cls", choices=("mhc1", "mhc2"))
    _add_store_opts(lg)
    lg.set_defaults(fn=cmd_logo)

    sp = sub.add_parser("span", help="extend an MHC-II binding core to the full presented ligand")
    sp.add_argument("core", help="the 9-mer binding core")
    sp.add_argument("--protein", required=True, help="source protein sequence, or a FASTA path")
    sp.add_argument("--mode", default="auto", choices=("auto", "observed", "modeled", "fixed"))
    sp.add_argument("--flanks", default="3,3", help="left,right sizes for --mode fixed")
    _add_store_opts(sp)                 # only used to supply the observed-ligand corpus
    sp.set_defaults(fn=cmd_span)

    pr = sub.add_parser("predict", help="score a variant peptide-window FASTA -> native + .scored.csv")
    pr.add_argument("fasta", help="a .peptide.fasta (pipeline schema)")
    pr.add_argument("--alleles", required=True, help="comma-separated HLA alleles (pipeline form)")
    pr.add_argument("--cls", required=True, choices=("mhc1", "mhc2"))
    pr.add_argument("--native", help="write the native TSV here")
    pr.add_argument("--scored-csv", dest="scored_csv", help="write the pipeline .scored.csv here")
    pr.add_argument("--rank-threshold", type=float, default=2.0, help="keep binders with %%rank <= this")
    pr.add_argument("--top", type=int, help="cap binders kept per window (strongest first)")
    pr.add_argument("--background", default="proteome", choices=("ligand", "proteome", "markov"))
    pr.add_argument("--footprint", default="adaptive", choices=("anchor", "core", "adaptive"))
    pr.add_argument("--seed", type=int, default=0)
    _add_store_opts(pr)
    pr.set_defaults(fn=cmd_predict)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
