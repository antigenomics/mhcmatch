# Minimal Sphinx makefile. Build: `env -C docs make html` from the repo root, with the repo-local
# `.venv/bin` on PATH (sphinx-build lives there, not on the system path).
#
# `-W --keep-going` is what .github/workflows/docs.yml runs, so a warning fails here too. A local
# build that tolerates warnings is a local build that passes while CI does not.
SPHINXOPTS    ?= -W --keep-going
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)

clean:
	rm -rf "$(BUILDDIR)"

.PHONY: html clean
