#!/usr/bin/env python3
"""Refuse to publish anything that identifies a collaborator's patient.

**This repository is public.** It carries integrations written against a private clinical pipeline,
and the boundary is one-directional: pipeline code and file contracts may cross, a patient may not.

The collaborating repository has its own `privacy_check.fish`, and it is **not** sufficient here --
it was written for a private-to-private boundary. It derives the surname list from local run
directories rather than holding one, greps only `.md .tex .py .html .json` (so a `.tsv`, `.csv` or
`.fasta` full of variants passes), and explicitly whitelists one donor by name. A tree can pass it
and still be unsafe to publish.

What this checks, over every tracked file:

1. **Donor tokens and internal hostnames**, matched case-sensitively on a word boundary --- read
   from a file **outside this repository**, named by ``$MHCMATCH_PRIVACY_TOKENS``. A guard that
   carries its own blocklist in a public tree is the single largest concentration of exactly the
   thing it exists to keep out, which is what this file once was. Unset, these two screens
   are skipped and ``main`` says so; screens 2 and 4 do not need a name to work.
2. **Cyrillic.** No tracked file in the collaboration repository contains any; a real name is the
   only reason it would appear here.
3. **Their internal hosts.** Naming an unreachable private host in public documentation is useless
   to a reader and discloses infrastructure.
4. **Variant-shaped data files.** A tracked `.tsv`/`.csv`/`.vcf` carrying `chrom`+`pos`+`ref`/`alt`
   -- or the VEP/pVACseq spelling of the same thing, `Chromosome`+`Start`+`Reference`/`Variant` --
   is a variant table. Synthetic fixtures declare themselves with a `chrN` contig; anything naming a
   real contig has to be justified rather than assumed.

Run: `python tools/privacy_check.py` (exit 1 on any finding). Also a test, so it runs in CI.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Where the screened tokens live, **outside this repository**. One token per line; ``#`` comments
#: allowed; two optional prefixes:
#:
#: * ``host:`` --- an internal hostname, matched as a plain substring.
#: * ``ctx:`` --- a token that is an ordinary word elsewhere (a surname that is also a word, a city
#:   that is also an identifier stem), so a bare match would be noise. It fires only in a
#:   donor-shaped context: beside ``donor``/``patient``/``sample``, or as the stem of a
#:   ``*.mhc*`` / ``*.alleles*`` / ``*.vaccine*`` filename.
#:
#: Everything else is a donor token, matched case-sensitively on a word boundary.
TOKENS_ENV = "MHCMATCH_PRIVACY_TOKENS"
CYRILLIC = re.compile(r"[Ѐ-ӿ]")
#: The opt-out marker, spelled out so it is greppable across the whole repo.
ALLOW_CYRILLIC = "privacy-check: cyrillic-ok"
#: A ``ctx:`` token fires only in a donor-shaped context, so a surname that is also an ordinary
#: word does not make the check unusable. ``{}`` is the token.
_CTX = (r"\b(?:donor|patient|sample)[_ =:/-]*{0}\b|\b{0}[._-](?:mhc|alleles|vaccine)")


def _tokens():
    """``(donors, contextual, hosts)`` from ``$MHCMATCH_PRIVACY_TOKENS``, or three empty lists.

    Empty is not "clean" -- it means the donor and host screens did not run, which
    :func:`main` reports rather than letting an unset variable read as a pass.
    """
    donors, ctx, hosts = [], [], []
    path = os.environ.get(TOKENS_ENV, "").strip()
    if not path or not Path(path).is_file():
        return donors, ctx, hosts
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("host:"):
            hosts.append(line[5:].strip())
        elif line.startswith("ctx:"):
            ctx.append(line[4:].strip())
        else:
            donors.append(line)
    return donors, [c for c in ctx if c], hosts


SKIP_SUFFIX = {".png", ".jpg", ".svg", ".gz", ".pkl", ".npz", ".idx", ".pdf", ".ico", ".whl"}


def tracked() -> list[Path]:
    out = subprocess.run(["git", "-C", str(ROOT), "ls-files"], capture_output=True, text=True)
    return [ROOT / p for p in out.stdout.splitlines() if p]


def findings() -> list[str]:
    donors, ctx, hosts = _tokens()
    donor_re = re.compile(r"\b(" + "|".join(map(re.escape, donors)) + r")\b") if donors else None
    ctx_res = {c: re.compile(_CTX.format(re.escape(c))) for c in ctx}
    bad: list[str] = []
    for path in tracked():
        if path.suffix.lower() in SKIP_SUFFIX or not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        # The name itself, not only the contents.
        if donor_re is not None and donor_re.search(path.name):
            bad.append(f"{rel}: donor token in the FILENAME")
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue

        # No self-exemption any more, and that is the point of the move: this file used to be
        # skipped because it carried the blocklist, which made the guard clean by construction
        # while being the largest concentration of the thing it screens for.
        if donor_re is not None:
            for m in set(donor_re.findall(text)):
                bad.append(f"{rel}: donor token {m!r}")
        for name, pat in ctx_res.items():
            if pat.search(text):
                bad.append(f"{rel}: donor-shaped reference to {name!r}")
        # A file may opt out with an explicit, greppable marker naming why. Cyrillic is not
        # forbidden as such -- a non-ASCII COLUMN NAME is a real thing a pipeline table carries and
        # is worth a test. A surname is what must not appear, and an exception that has to be typed
        # out and justified is visible in review; a silently relaxed rule is not.
        if CYRILLIC.search(text) and ALLOW_CYRILLIC not in text:
            bad.append(f"{rel}: Cyrillic text (add `{ALLOW_CYRILLIC} <reason>` if deliberate)")
        for h in hosts:
            if h in text:
                bad.append(f"{rel}: internal host {h!r}")

        if path.suffix.lower() in (".tsv", ".csv", ".vcf"):
            head = text.split("\n", 1)[0].lower()
            cols = {c.strip() for c in re.split(r"[,\t]", head)}
            # Two vocabularies, because a variant table has two: ours (`chrom`/`pos`/`ref`/`alt`)
            # and the VEP/pVACseq one a standard upstream writes
            # (`Chromosome`/`Start`/`Reference`/`Variant`). Screening only the first would have
            # passed every table the generic path actually produces.
            ours = {"chrom", "pos"} <= cols and bool({"ref", "alt"} & cols)
            theirs = {"chromosome", "start"} <= cols and bool({"reference", "variant"} & cols)
            if ours or theirs:
                # Skip the header row: it literally contains the word `chrom`, which the contig
                # pattern would otherwise match as a contig named "chrom".
                body = text.split("\n", 1)[1] if "\n" in text else ""
                contigs = set(re.findall(r"(?m)(?:^|[,\t])(chr[0-9XYMxym][\w]*)(?=[,\t]|$)", body))
                real = {c for c in contigs if c.lower() != "chrn"}
                if real:
                    bad.append(f"{rel}: variant table naming real contigs {sorted(real)[:4]}")
    return sorted(set(bad))


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    # **CI cannot run the name screens**: the token list lives outside the repo by design, so
    # there is nothing for a public runner to read. Without this flag the run would print its
    # NOTE to stderr and exit 0, and a green check would report a pass it never performed.
    # The operator passes it before a push; CI runs the structural screens and says which.
    ap.add_argument("--require-tokens", action="store_true",
                    help=f"exit 2 unless ${TOKENS_ENV} names a readable token file, so an "
                         f"incomplete screen cannot be mistaken for a clean one")
    a = ap.parse_args()

    donors, ctx, hosts = _tokens()
    if not (donors or ctx or hosts):
        note = (f"${TOKENS_ENV} is unset or missing, so the donor and host screens did NOT "
                f"run. The Cyrillic and variant-table screens did.")
        if a.require_tokens:
            print(f"PRIVACY CHECK INCOMPLETE -- {note}", file=sys.stderr)
            return 2
        print(f"NOTE: {note}", file=sys.stderr)
    bad = findings()
    if bad:
        print("PRIVACY CHECK FAILED -- do not push:", file=sys.stderr)
        for b in bad:
            print(f"  {b}", file=sys.stderr)
        return 1
    which = "every screen" if (donors or ctx or hosts) else "the structural screens only"
    print(f"privacy check clean over {len(tracked())} tracked file(s) -- {which}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
