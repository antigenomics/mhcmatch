# Compile the mhcmatch theory appendix. Requires a TeX distribution (latexmk + lualatex + bibtex);
# lualatex is needed for the OldStandard OpenType text+math fonts and is pinned in .latexmkrc.
# The figures (diffusion_auc.pdf, mhc2_pockets.pdf) are committed; regenerate them from the
# pmhc_data tables with `python ../bench/make_figures.py --pmhc-dir <dir>` (needs gnuplot).
FIGS = diffusion_auc.pdf \
       pockets_mhc1_human.pdf pockets_mhc1_mouse.pdf pockets_mhc2_human.pdf pockets_mhc2_mouse.pdf \
       promiscuity_mhc1_human.pdf \
       promiscuity_shared_mhc1_human.pdf promiscuity_shared_mhc1_mouse.pdf \
       promiscuity_shared_mhc2_human.pdf promiscuity_shared_mhc2_mouse.pdf promiscuity_mhc1_human.pdf

mhcmatch.pdf: mhcmatch.tex refs.bib $(FIGS)
	latexmk -lualatex -interaction=nonstopmode -halt-on-error mhcmatch.tex

.PHONY: clean
clean:
	latexmk -C mhcmatch.tex
