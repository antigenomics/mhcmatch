#!/usr/bin/env fish
# Repo-local venv + editable install. Flags:
#   --tests  also install the pytest extra
#   --logo   also install logo rendering extras (logomaker, matplotlib)
set repo (dirname (status --current-filename))
cd $repo

if not test -d .venv
    python3 -m venv .venv
end
source .venv/bin/activate.fish
pip install -U pip

# Co-develop against a sibling seqtree checkout if present; else PyPI resolves it.
if test -d ../seqtree
    pip install -e ../seqtree
end

set extras ""
if contains -- --tests $argv
    set extras "[test]"
else if contains -- --logo $argv
    set extras "[test,logo]"
end

pip install -e ".$extras"

echo ""
echo "Done. Tests: pytest -s tests/"
echo "      Docs:  ROADMAP.md and appendix/mhcmatch.pdf (env -C appendix make)"
