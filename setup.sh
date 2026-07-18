#!/usr/bin/env bash
# Repo-local venv + editable install. POSIX syntax -- runs under bash or zsh (or sh).
# Flags:
#   --tests  also install the pytest extra
#   --logo   also install logo rendering extras (logomaker, matplotlib)
set -eu
repo="$(cd "$(dirname "$0")" && pwd)"
cd "$repo"

if [ ! -d .venv ]; then
    python3 -m venv .venv
fi
.venv/bin/python -m pip install -U pip

# Co-develop against a sibling seqtree checkout if present; else PyPI resolves it.
if [ -d ../seqtree ]; then
    .venv/bin/pip install -e ../seqtree
fi

# --tests takes precedence over --logo (matches the original else-if); quoting keeps zsh from globbing "[...]".
extras=""
case " $* " in
    *" --tests "*) extras="[test]" ;;
    *" --logo "*)  extras="[test,logo]" ;;
esac

.venv/bin/pip install -e ".$extras"

echo ""
echo "Done. Tests: pytest -s tests/"
echo "      Docs:  ROADMAP.md and appendix/mhcmatch.pdf (env -C appendix make)"
