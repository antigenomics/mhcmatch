// One cohort-level task that computes nothing and exists only to refuse to start.
//
// **Every failure this catches is one that otherwise surfaces hours in, inside a fan-out task log.**
// Three of them, in the order they bite:
//
//   1. mhcmatch is not installed on the executor -- `command not found`, discovered at the FIRST
//      task of the run rather than at the twentieth.
//   2. the installed release is not the one this overlay was written against. The processes call
//      CLI flags BY NAME, so an older release exits 2 on an unknown flag deep in a task log and
//      `mhcmatch --version` cannot warn you because both print a string. It has happened: PyPI
//      served the previous release while the module had gained `--map-binder`, and both cassette
//      tasks died with `unrecognized arguments`.
//   3. the candidate table names no peptide or no allele column -- which, discovered later, reads
//      as "this candidate named no allele we know", a real and DIFFERENT state.
//
// In `off` mode it warns instead of failing, so an overlay that is merged but not yet enabled can
// never break a run.

process MHCMATCH_PREFLIGHT {
    tag 'mhcmatch'
    label 'process_single'

    conda "${moduleDir}/../mhcmatch/environment.yml"

    input:
    path table          // one candidate table, purely to read its header
    val  mode

    output:
    path 'mhcmatch_preflight.txt', emit: report

    script:
    def want = params.mhcmatch_overlay_require_version ?: ''
    def hard = mode == 'off' ? '0' : '1'
    """
    set +e
    STRICT=${hard}
    : > mhcmatch_preflight.txt

    fail() {
        echo "MHCMATCH PREFLIGHT: \$1" | tee -a mhcmatch_preflight.txt >&2
        if [ "\$STRICT" = "1" ]; then exit 3; fi
        echo "  (mhcmatch_overlay_mode = off, so this is a warning)" >&2
    }

    if ! command -v mhcmatch > /dev/null 2>&1; then
        fail "mhcmatch is not on PATH for this executor. Install it, or set params.mhcmatch_container."
    else
        # **The LAST whitespace-separated field, not the whole line.** `mhcmatch --version` prints
        # `mhcmatch 1.20.0`, so stripping whitespace gives `mhcmatch1.20.0` and the comparison
        # below rejects every correct install. Caught by running this for real; a stub would not
        # have, because a stub runs no command.
        GOT=\$(mhcmatch --version 2>/dev/null | tail -1 | awk '{print \$NF}')
        echo "mhcmatch \$GOT" >> mhcmatch_preflight.txt
        if [ -n "${want}" ] && [ "\$GOT" != "${want}" ]; then
            fail "installed mhcmatch is \$GOT, this overlay expects ${want}. The processes call CLI flags by name, so a mismatch fails inside a task log with 'unrecognized arguments'."
        fi
    fi

    # The two columns `rank pairs` genuinely requires, under any of the spellings it accepts. The
    # header is read with the delimiter it actually has -- a candidate table may be comma- or
    # tab-separated, and a tab-only split turns a wide CSV schema into one column.
    python - "${table}" >> mhcmatch_preflight.txt <<'PY'
import sys
p = sys.argv[1]
head = open(p, encoding="utf-8", errors="replace").readline().rstrip("\\n")
sep = "," if head.count(",") > head.count("\\t") else "\\t"
cols = [c.strip() for c in head.split(sep)]
pep = [c for c in ("peptide", "epitope") if c in cols]
alle = [c for c in ("allele", "best_allele") if c in cols]
print(f"table {p}: {len(cols)} column(s), delimiter {sep!r}")
print(f"  peptide column: {pep or 'MISSING'}")
print(f"  allele  column: {alle or 'MISSING'}")
if not pep or not alle:
    sys.exit(4)
PY
    RC=\$?
    if [ "\$RC" != "0" ]; then
        fail "the candidate table names no peptide (peptide/epitope) and/or no allele (allele/best_allele) column. Those two are the only columns required; everything else is carried through untouched."
    fi

    echo "MHCMATCH PREFLIGHT OK (mode=${mode})" >> mhcmatch_preflight.txt
    cat mhcmatch_preflight.txt
    exit 0
    """

    stub:
    """
    echo "MHCMATCH PREFLIGHT OK (stub, mode=${mode})" > mhcmatch_preflight.txt
    """
}
