# Compile the mhcmatch theory appendix. Requires a TeX distribution (latexmk + lualatex + bibtex);
# lualatex is needed for the OldStandard OpenType text+math fonts and is pinned in .latexmkrc.
mhcmatch.pdf: mhcmatch.tex refs.bib
	latexmk -lualatex -interaction=nonstopmode -halt-on-error mhcmatch.tex

.PHONY: clean
clean:
	latexmk -C mhcmatch.tex
