# Vector assembly: what the literature actually establishes

Retrieved from **PubMed** on 2026-08-18 for the `mhcmatch.vector` roadmap. Every PMID and DOI below
was retrieved in that scan, never recalled. Each claim carries a tier:

| tier | meaning |
|---|---|
| **[exp]** | measured in a wet-lab experiment, with the readout named |
| **[obs]** | observed in a clinical or cohort dataset, not an intervention |
| **[silico]** | asserted by a computational design paper with **no experimental validation** |
| **[open]** | no evidence retrieved either way |

The **[silico]** tier exists because it is the single biggest hazard in this literature. The
multi-epitope vaccine design field has a large body of papers that all use the same linker
conventions — `AAY` between CTL epitopes, `GPGPG` between helper epitopes, `KK` between B-cell
epitopes, `EAAAK` to fuse an adjuvant — and cite each other for them. Almost none of them tests a
linker against an alternative. Convention repeated is not evidence.

---

## 1. Spacers — the one head-to-head experiment goes against Gly/Ser

**[exp] Alanine-based spacers beat `GGGS` for MHC-I presentation, and epitope position did not
matter.** Aguilar-Gurrieri et al. built neoantigen concatemers differing only in the linker and
assayed H-2Kb presentation of SIINFEKL from the processed polypeptide. Spacers had "a large impact on
the efficiency of neoantigen processing and presentation by MHC-I molecules; **in contrast, the
peptide position and the flanking regions have a minimal impact**", and "linkers based on alanine
residues promote a more efficient peptide presentation than the commonly used GGGS linker."
*Cancer Immunol Immunother* 2023, PMID 36820900,
[DOI](https://doi.org/10.1007/s00262-023-03409-3).

This is the most directly relevant experiment retrieved, and it cuts two ways for us:

1. It **supports** alanine spacers (`AAA`, `AAY`) for class-I units, against the Gly/Pro-permissive
   argument the module docstring currently leans on — which came from a bioinformatic analysis of
   ligand-flanking composition (Martín-Galiano & López, PMID 30645615), not from a presentation assay.
2. It **undercuts** ordering as a lever for class I. Position and flanking had minimal impact on
   presentation efficiency in that system.

**[exp] `GPGPG` rescues class-II junctional suppression.** Unspaced concatenation of four HLA-DR
epitopes created a high-affinity junctional epitope that suppressed the response to all four;
inserting `GPGPG` restored all four. Livingston et al., *J Immunol* 2002, PMID 12023344.
Replicated for Aβ: a tandem `Aβ1-15` repeat without spacer created junctional epitopes and lost
anti-Aβ titre, and "the disruption of junctional epitopes through the introduction of a GPGPG spacer
restored the immunogenicity against all the epitopes". Guan et al., *Neuroreport* 2012, PMID 22922658,
[DOI](https://doi.org/10.1097/WNR.0b013e328358a044); same group, tetravalent construct,
*Hum Vaccin Immunother* 2013, PMID 23732905, [DOI](https://doi.org/10.4161/hv.24830).

**So the two results are not in conflict — they are about different classes.** Every retrieved
`GPGPG` result is class II or antibody; the only class-I head-to-head favours alanine. That gives a
class-conditional default rather than one global spacer list.

**[silico] Everything else.** `AAY` for CTL + `GPGPG` for HTL + `EAAAK` for adjuvant fusion is
asserted without comparison by Nezafat et al., *J Theor Biol* 2014, PMID 24512916,
[DOI](https://doi.org/10.1016/j.jtbi.2014.01.018); a furin-sensitive `RVRR` linker plus `GPGPG` and
`A(EAAAK)2A` by Safavi et al., *Vaccine* 2020, PMID 33082015,
[DOI](https://doi.org/10.1016/j.vaccine.2020.10.016); `AAY`/`GPGPG`/`GDGDG` variants by Roohparvar
Basmenj et al., *Sci Rep* 2023, PMID 37940672, [DOI](https://doi.org/10.1038/s41598-023-46408-1);
`GPGPG` + `KK` by Li et al., *Intervirology* 2016, PMID 27096202,
[DOI](https://doi.org/10.1159/000445059). None of these compares its linker against an alternative in
an assay. **Do not cite them as evidence for a linker; cite them only as evidence of convention.**

---

## 2. Ordering — a large free set, and little sign the order matters

**[exp] Position had minimal impact on MHC-I presentation efficiency** — Aguilar-Gurrieri 2023, above
(PMID 36820900, [DOI](https://doi.org/10.1007/s00262-023-03409-3)).

**[exp] All six orderings of three Fel d I regions produced no detectable junctional responses.**
Rogers et al., *Mol Immunol* 1994, PMID 7521933. The only systematic permutation experiment
retrieved, and it found no position effect.

**[exp] A polyepitope's antibody response contained no fraction against any interepitope junction**,
across several mouse strains and rabbits, while retaining an intrinsic B-cell immunodominance
hierarchy driven by the epitopes themselves. Kumar et al., *Vaccine* 1994, PMID 7513115,
[DOI](https://doi.org/10.1016/0264-410x(94)90203-8).

**[silico] The junction-free design space is astronomically large.** Lee et al. enumerated
*all* junction-free orderings for a given epitope set under class-I and class-II motif constraints and
found "the number of such variants of any given polyepitope can be astronomically high".
*Biomed Microdevices* 2010, PMID 20033850, [DOI](https://doi.org/10.1007/s10544-009-9376-7).

⇒ **Ordering is a constraint-satisfaction problem, not an optimisation problem.** Junction-clean
layouts are abundant; once inside that set, no retrieved experiment shows which one to prefer. A
solver that returns *one* cheap path is leaving a large feasible set unexplored, and the right use of
that freedom is a **second** objective, not a better first one.

**[obs] Prior art for the whole layout step: PolyCTLDesigner.** Given a set of epitopes it selects
N-terminal flanking sequences per epitope to optimise TAP binding, joins them so proteasomal /
immunoproteasomal processing liberates them, and minimises non-target junctional epitopes — using
TAP-binding and cleavage-specificity patterns with a genetic algorithm and graph theory. Antonets &
Bazhan, *BMC Res Notes* 2013, PMID 24107711, [DOI](https://doi.org/10.1186/1756-0500-6-407).
`mhcmatch.vector.order` currently implements only the third of those three. pVACvector is the other
prior art and does the same third piece by simulated annealing: Hundal et al.,
*Cancer Immunol Res* 2020, PMID 31907209.

---

## 3. Flanking residues — large, measured, and currently unmodelled

**[exp] TAP prefers N-terminally extended precursors over minimal epitopes.** Several melanoma
epitopes (tyrosinase `YMNGTMSQV`, MAGE-1 `EADPTGHSY`) are *poor* TAP substrates as minimals, because
they carry N-terminal residues deleterious for TAP binding. Synthesising the same epitopes with
N-terminal extensions of up to four residues showed "the longer peptides are indeed transported into
the ER at a significantly higher level than the original epitopes". Wang, Guttoh & Androlewicz,
*Melanoma Res* 1998, PMID 9764810, [DOI](https://doi.org/10.1097/00008390-199808000-00008).

**[exp] Flanking sequence can gate presentation entirely.** An influenza NP epitope was not presented
at all from one minigene construct; adding a **single** C-terminal methionine restored presentation.
Extending N-terminally instead required **more than 55 residues**, with an abrupt transition (91-155
gave little or none, 90-155 gave full presentation). Yellen-Shaw & Eisenlohr, *J Immunol* 1997,
PMID 9029109.

**[exp] Flanking prolines improve processing.** An HIV-2 gag epitope was "most efficiently processed
from precursors that contain two flanking proline residues", and the proline motif tracked lower
viral load in patients. Jallow et al., *Eur J Immunol* 2015, PMID 26018465,
[DOI](https://doi.org/10.1002/eji.201545451).

**[exp] A mutation *outside* an epitope can abolish it.** HCV genotype 1a escaped an HLA-B\*51 epitope
by a substitution **five residues upstream** of it, which "impaired recognition of target cells
presenting the endogenously processed epitope" — the epitope sequence itself was untouched. Walker et
al., *J Virol* 2015, PMID 26446603, [DOI](https://doi.org/10.1128/JVI.01993-15).

**[exp] Short natural flanks improved presentation from a VLP scaffold** relative to the epitope
inserted with no flanks. Rueda et al., *J Gen Virol* 2004, PMID 14993639,
[DOI](https://doi.org/10.1099/vir.0.19525-0).

⇒ A 27-mer centred on the mutation gets native flanks *by construction*, which is why our unit
default is defensible. But the **spacer** replaces the native flank on one side of every junction, and
no part of the module currently asks whether the replacement is a good one.

Note the tension with §1: Aguilar-Gurrieri found flanking regions had minimal impact where the
*spacer* dominated. The reconciliation is that these are different comparisons — native-flank vs
no-flank (large effect) against one native flank vs another (small effect).

---

## 4. Duplication and nesting

**[exp] A duplicated epitope was processed only when the copies were separated.** Of a set of VLP
constructs, "only PPV-VLPs carrying two copies of the OVA epitope **linked by two glycines** were able
to be properly processed, suggesting that the introduction of flexible residues between the two
consecutive OVA epitopes may be necessary for the correct presentation of these dimers".
Rueda et al. 2004, PMID 14993639, [DOI](https://doi.org/10.1099/vir.0.19525-0).

**[exp] Copy number is not monotonic.** Comparing 3 to 20 displayed copies of a malaria repeat,
"low copy number can reduce the abundance of low-affinity mAb epitopes while retaining high-affinity
mAb epitopes"; five copies was optimal and beat a near-full-length protein. Langowski et al.,
*PNAS* 2020, PMID 31988134, [DOI](https://doi.org/10.1073/pnas.1911792117). B-cell readout, so it
transfers to T-cell cassettes only as a shape — more copies is not more response.

**[exp] Nesting a helper epitope over a CTL epitope works, but so does putting it far away, and
position inside the nest does not matter.** For the immunodominant HLA-A\*02:01 pp65(495-503) epitope,
"homologous CD4 T cell help, located within an overlapping (nested) pp65(487-503) domain, facilitated
induction"; but "the position of the e6 epitope within this nested domain is not critical", and
"**distant CD4 T cell epitope(s) can thus provide efficient help**". Reiser et al., *J Immunol* 2011,
PMID 21810614, [DOI](https://doi.org/10.4049/jimmunol.1002512).

⇒ **Duplication buys nothing a centred long unit does not already give** (all registers spanning the
mutation are present), and where duplication is unavoidable it needs a flexible separator. **Nesting
is not required**; linkage to help is (§5).

---

## 5. CD8 + CD4 — linkage is required, help is per-epitope

**[exp] Physical linkage was necessary; co-delivery failed.** A GFP–TRP-2 fusion DNA vaccine
suppressed B16 melanoma growth while TRP-2 alone did not, and the antitumour immunity was abolished
in proteasome-activator-deficient mice. Decisively: "**genetic immunization with pGFP plus
pTRP-2(181-188) failed to exert the antitumour immunity**" — the same two components delivered as two
separate plasmids gave nothing. Zhang et al., *Immunology* 2004, PMID 15270727,
[DOI](https://doi.org/10.1111/j.1365-2567.2004.01916.x).

**[exp] Help-dependence is epitope-specific, not global.** Of three HLA-A\*02:01-restricted pp65
epitopes, only the immunodominant one "critically depends on CD4 T cell help"; the other two were
induced by monospecific vaccines with no help, and adding heterologous helper epitopes enhanced only
the help-dependent one. Reiser et al. 2011, PMID 21810614,
[DOI](https://doi.org/10.4049/jimmunol.1002512).

**[exp] Deleting the dominant epitope redistributed the response.** Abrogating the immunodominant
epitope "efficiently enhanced e3- and e8-specific T cell responses". Same paper.

**[exp] PADRE, a synthetic pan-DR helper epitope, raises CD8 output when co-formulated.** Liposomal
HER2-derived CTL peptide + PADRE gave "superior induction of CD4 and CD8 T cells responses and
significantly enhanced production of IFN-γ" over the CTL peptide alone, with tumour control.
Zamani et al., *Eur J Cancer* 2020, PMID 32145473,
[DOI](https://doi.org/10.1016/j.ejca.2020.01.010); same group with MPL,
*J Control Release* 2019, PMID 30999007, [DOI](https://doi.org/10.1016/j.jconrel.2019.04.019).
Carried into a clinical-stage vector as a genetic fusion, giving PADRE-specific CD4 alongside
target-specific CD8 and B responses without autoimmunity: Snook et al., *Hum Gene Ther Methods* 2016,
PMID 27903079, [DOI](https://doi.org/10.1089/hgtb.2016.114).

⇒ **CD4 and CD8 payloads belong in the same molecule.** That is an experimental result, not a
formulation preference, and it settles the open fork noted in the earlier design memo between
"separate formulations per class" and "link them".

---

## 6. Trafficking tags — real, but they do not steer the arm they are named for

**[exp] Sig/LAMP-1 routing raised both arms.** Linking HPV-16 E7 to the LAMP-1 sorting signal
"dramatically increased in vitro activation and in vivo expansion of E7-specific CD4 **and** CD8 T
cells" and improved tumour control. Kang et al., *Immunol Lett* 2006, PMID 16844231,
[DOI](https://doi.org/10.1016/j.imlet.2006.05.004).

**[exp] But the tag does not select the arm.** A rabies DNA vaccine compared tPA, ubiquitin and
LAMP-1 tags, each chosen to bias a different arm. All improved on the untagged vaccine, but "the
response elicited **did not pertain to the type of target sequence and the directed arm of
immunity**… the DNA vaccines that had been designed to generate different type of immune responses
yielded in effect similar response", and the authors conclude the effect is antigen-dependent.
Kaur, Rai & Bhatnagar, *Vaccine* 2009, PMID 19356616,
[DOI](https://doi.org/10.1016/j.vaccine.2009.01.128).

**[exp] Ubiquitin fusion works through the degradation pathway.** Antitumour immunity from the
GFP–TRP-2 fusion was "completely cancelled in mice deficient in proteasome activator PA28α/β".
Zhang et al. 2004, PMID 15270727, [DOI](https://doi.org/10.1111/j.1365-2567.2004.01916.x).

**[exp] Secretion drives the antibody arm.** A secretion-competent HBV envelope construct "elicited
strong and sustained immunity" where the retained form gave only a weak response, and blocking
secretion abolished the gain. Prange & Werr, *Vaccine* 1999, PMID 10067665,
[DOI](https://doi.org/10.1016/s0264-410x(98)00243-6). Relevant to a signal-peptide choice, but the
readout is humoral.

**[silico] MITD and tPA in current mRNA design papers.** Both appear as standard cassette furniture —
5' m7G cap, 5'UTR, Kozak, tPA signal peptide, adjuvant, MITD, stop, 3'UTR, poly(A)120 — in
Ahmed et al., *Front Immunol* 2025, PMID 40433366, [DOI](https://doi.org/10.3389/fimmu.2025.1480025)
and Ali & Luqman Ali, *Viral Immunol* 2025, PMID 40125606,
[DOI](https://doi.org/10.1089/vim.2025.0004). Both are docking-and-simulation papers with no
experimental readout. The experimental basis for MITD remains Kreiter et al., *J Immunol* 2008,
PMID 18097032, which is not part of this scan.

⇒ **A trafficking tag is worth having and is not worth trusting to select an arm.** It belongs in the
backbone as a recorded, swappable choice, not as a mechanism the selection layer reasons about.

---

## 7. The mRNA layer

**[exp] m1Ψ causes +1 ribosomal frameshifting, and synonymous editing fixes it.** N1-methyl-
pseudouridine incorporation produces +1 frameshifted product at slippery sequences, and cellular
immunity to the frameshifted products of BNT162b2 was detected in vaccinated mice and humans;
"synonymous targeting of such slippery sequences provides an effective strategy to reduce the
production of frameshifted products". Mulroney et al., *Nature* 2023, PMID 38057663,
[DOI](https://doi.org/10.1038/s41586-023-06800-3).

This one is **already implemented** — `vector.slippery_sites` / `vector.deslip`.

---

## 8. What this scan did *not* find

Recorded so nobody re-runs it expecting a different answer:

- **No published objective function for epitope count or per-allele allocation in an individual
  patient.** [open]
- **No trial delivering N versus 2N epitopes at matched total dose**, so the within-cassette dilution
  coefficient is unmeasured. [open]
- **No experimental comparison of `AAY` against `AAA`** — the alanine result (PMID 36820900) is
  alanine-based-vs-`GGGS`, not tyrosine-vs-alanine. [open]
- **No T-cell-readout study varying epitope copy number inside one cassette.** The copy-number data
  are B-cell display. [open]
- **No study of class-I and class-II junction interaction** in a mixed cassette. [open]
