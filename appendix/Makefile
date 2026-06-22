# Compile the mhcmatch theory appendix. Requires a TeX distribution (latexmk + lualatex + bibtex);
# lualatex is needed for the OldStandard OpenType text+math fonts and is pinned in .latexmkrc.
# The figures (diffusion_auc.pdf, mhc2_pockets.pdf) are committed; regenerate them from the
# pmhc_data tables with `python ../bench/make_figures.py --pmhc-dir <dir>` (needs gnuplot).
FIGS = diffusion_auc.pdf mhc2_pockets.pdf

mhcmatch.pdf: mhcmatch.tex refs.bib $(FIGS)
	latexmk -lualatex -interaction=nonstopmode -halt-on-error mhcmatch.tex

.PHONY: clean
clean:
	latexmk -C mhcmatch.tex
