"""Command-line interface for mhcmatch: ``mhcmatch <command> ...``.

Commands: ``decompose`` (no data needed), ``restriction``, ``scan``, ``logo`` (need a pmhc_data
table via ``--pmhc`` or ``$MHCMATCH_PMHC``), and ``source`` (needs a proteome FASTA).
"""
from __future__ import annotations

import argparse
import os

from . import Proteome, Store


def _add_store_opts(p):
    p.add_argument("--pmhc", help="pmhc_data TSV(.gz); else $MHCMATCH_PMHC/pmhc_<tier>.tsv.gz")
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


def cmd_restriction(a):
    res = _store(a).restriction(a.peptide, cls=a.cls, alleles=[a.allele] if a.allele else "all",
                                top=a.top, diffuse=a.diffuse)
    if not res:
        print("no presenting allele (no presentation-signature neighbours)")
        return
    print(f"{'allele':<18}{'vote':>7}{'enr':>7}" + ("{:>8}".format("score") if a.diffuse else "")
          + f"{'binder':>8}")
    for r in res:
        line = f"{r.allele:<18}{r.vote:>7.2f}{r.enrichment:>7.1f}"
        if a.diffuse:
            line += f"{(r.anchor_score or 0.0):>8.2f}"
        print(line + f"{'yes' if r.binder else 'no':>8}")


def cmd_scan(a):
    hits = _store(a).scan_protein(_read_seq(a.protein), cls=a.cls or "mhc1",
                                  alleles=[a.allele] if a.allele else "all", top=a.top)
    print(f"# {len(hits)} presented window(s)")
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
    r.add_argument("--top", type=int, default=10)
    _add_store_opts(r)
    r.set_defaults(fn=cmd_restriction)

    s = sub.add_parser("scan", help="find presented peptides in a protein (sequence or FASTA path)")
    s.add_argument("protein")
    s.add_argument("--allele")
    s.add_argument("--cls", choices=("mhc1", "mhc2"))
    s.add_argument("--top", type=int, default=3)
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

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
