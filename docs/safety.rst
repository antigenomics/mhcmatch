Safety, prior evidence, and what goes in the cassette
=====================================================

Ranking says which mutations are worth targeting. This page is the step after: **what to withdraw
on safety grounds, what prior evidence a candidate already carries, and how many units to spend.**
They are three different questions and none of them is answered by the candidate's score.

Everything here is computed, not looked up. That matters for the one caveat this page repeats:
on a corpus assembled *from* the reference sets, the prior-evidence columns are self-fulfilling and
tell you nothing; on a fresh patient's variants they are the most informative columns in the table.

.. contents::
   :local:
   :depth: 1

.. important::

   **Everything on this page is class I — CD8 epitopes — and deliberately so.** Both the safety
   screen and the near-exact known-antigen lookup are built on a CD8 mechanism: a register that *is*
   an essential-tissue self peptide killing the cell that presents it, and a minimal epitope close
   enough to a confirmed neoantigen that the same clonotype could see both.

   Neither transfers to class II by rewriting the length range. A class-II ligand's hazard is not
   direct cytotoxicity — CD4 self-reactivity runs through help, hypersensitivity and allergy, which
   is a different question with a different literature and different thresholds, and none of it has
   been measured here. Reusing these numbers on class-II ligands would assert an equivalence nobody
   has established.

   So: run ``vector``, ``neoag`` and the ``rank`` annotation on ``--cls mhc1``. ``mhcmatch vector``
   accepts ``--cls mhc2`` for the *register vocabulary* used in junction scanning and claims nothing
   beyond that; the shipped Nextflow subworkflow filters all three to class I.


The safety screen: two clauses, and it excludes rather than down-ranks
----------------------------------------------------------------------

:func:`mhcmatch.vector.screen` withdraws a unit outright. The second-best cassette is cheap and
myocarditis is not, and capacity spent on a unit that has to be withdrawn is capacity not spent on
a safe one.

The precedent is not hypothetical and it was not a binding-prediction failure:

* An affinity-enhanced TCR against the HLA-A\*01:01 MAGE-A3 epitope ``EVDPIGHLY`` killed the first
  two patients infused, by cardiogenic shock within days. Autopsy found myocardial damage with **no
  MAGE-A3 expressed in heart at all**; the off-target was ``ESDPIVAQY``, from titin (Linette *et
  al.*, *Blood* 2013;122(6):863-71, PMID 23770775; Cameron *et al.*, *Sci Transl Med*
  2013;5(197):197ra103, PMID 23926201).
* A TCR recognising MAGE-A3/A9/A12 caused necrotising leukoencephalopathy and two deaths, because
  MAGE-A12 turned out to be transcribed in human brain (Morgan *et al.*, *J Immunother*
  2013;36(2):133-51, PMID 23377668).

Both were invisible to binding prediction and visible in **expression**. So the check joins a
peptide to a *protein* to a *tissue*, and it asks two questions:

``target gene``
   Is the unit's own gene transcribed in a tissue that must not be attacked? The MAGE-A12 case. No
   register search is needed to see it — the floor is 0.25 TPM rather than the conventional 5,
   because MAGE-A12 in brain caudate is **0.33**.

   **Asked only of a product the normal proteome carries.** MAGE-A12 is a cancer-testis antigen — a
   shared, *unmutated* self protein — so its brain transcription is the hazard precisely because the
   construct encodes a sequence brain tissue also presents. A somatic neoantigen is a different
   object: a missense, a frameshift, an inframe indel or a fusion junction encodes a sequence absent
   from normal tissue **by construction**, so its parent gene's expression is not that hazard, and
   what is one is the second clause. The gate is :data:`mhcmatch.predict.NOVEL_PRODUCTS` read off
   :attr:`~mhcmatch.vector.Unit.kind`; an ``isoform``, a wild-type or overexpressed target keeps the
   clause, and **an unknown or missing kind keeps it too** — the screen does not exempt a unit
   because nobody annotated it.

   Ungated, the clause withdrew a candidate for the fact that its parent gene exists: on a 37-donor
   cohort, 10 of 37 donors lost every unit they had, and one lost 1,098 of 1,618 to this clause
   alone.

``unrelated self origin``
   Does any register of the unit coincide with a self peptide from a **different** gene, and is that
   gene transcribed in an essential tissue? The titin-shaped case.

   **Asked only of the registers that carry novel sequence**, where there is novel sequence to carry.
   A 27-mer unit is thirteen-fourteenths wild type by construction, and the unrestricted clause read
   that design as the hazard: on **178 experimentally immunogenic somatic neoantigens** rebuilt as
   cassette units, **178 of 178 trip it**, median **36** self registers each — and 36 is exactly
   ``12 + 10 + 8 + 6``, the number of 8/9/10/11-mer windows of a 27-mer that cannot contain a centred
   mutation. The measured self fraction matches that geometry at every length (60.02 / 52.6 / 44.4 /
   35.2 % against 60.0 / 52.6 / 44.4 / 35.3 predicted), reaching **99.1 % of the geometric ceiling**,
   while at the minimal-epitope level the clause is clean: 0 of 178 mutant epitopes are in the
   proteome and 178 of 178 wild types are. There were no genuine coincidences to find.

   A window with no novel residue in it is therefore **structurally exempt** — it is wild-type
   sequence, it was always going to be in the proteome, and no cassette avoids it short of not using
   long units. Which windows those are depends on the product:
   :data:`~mhcmatch.predict.TRACT_PRODUCTS` (``frameshift``, ``fusion``) are novel from the variant
   offset **to the end of the unit**; the rest of :data:`~mhcmatch.predict.NOVEL_PRODUCTS` at that
   one index. Every clause-2 reason carries ``n_registers_spanning`` and ``n_hit_spanning`` so the
   exemption is auditable. For an ``isoform``, a ``cnv`` or an unannotated unit **every** register is
   judged, as before — the same gate as clause 1, so the two rules cannot disagree.

Two floors, not one
~~~~~~~~~~~~~~~~~~~

``min_tpm = 0.25`` is a **reporting** floor: below it a finding is not recorded at all, and it sits
under MAGE-A12's 0.33 so the fatal case is always visible. It is not an exclusion line — 0.25 TPM is
"detectable somewhere", which nearly every human gene is. ``veto_tpm = 5.0``, the conventional
"is it expressed" cut, is the exclusion line, and it applies under
``self_origin_risk(..., graded=True)`` / ``mhcmatch vector --screen-mode graded``: a finding below it
is kept as a per-unit **off-target fingerprint** rather than a refusal, reported in the cassette's
``fingerprint`` rows and priced into composition by ``--weight-offtarget``
(:func:`mhcmatch.vector.offtarget_cost` → :func:`mhcmatch.portfolio.compose`). The default stays
``veto``.

**The unit's own gene is excluded from the second clause, and it has to be.** A vaccine unit is a
long window of native context: its flanking registers *are* self peptides from its own parent
protein, and its mutated register sits one substitution from that protein's wild type. Screened
naively, every unit of every cassette fires. Those matches are also the ones tolerance already
covers — the flanks are presented in normal tissue daily.


Near-exact identity, not similarity
-----------------------------------

The obvious alternative — score each register with :mod:`mhcmatch.mimicry` and flag whatever
resembles a tolerance-side reference — was built and measured against this one, on 1,000 viral
epitopes (which cannot be self, so every firing is a false positive) and 1,000 thymic peptides from
essential-tissue genes:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20

   * - route
     - false pos.
     - true pos.
   * - mimicry, anchor-masked
     - 0.693
     - 0.944
   * - self origin (this one)
     - **0.020**
     - 0.940

.. note::

   ``bench/results/...`` paths on this page resolve in the benchmark repository,
   ``2026-mhcmatch-code`` (private, released with the manuscript), not in the
   library repo.

**Equal sensitivity, 35× the false positives.** Anchor-channel similarity to a presented reference
is *presentation*, not recognition, so a masked match fires for every peptide sharing the allele's
motif — the influenza epitope ``GILGFVFTL`` draws 14 essential-tissue hits. A screen that flags
two-thirds of a candidate list excludes nothing in practice, because nobody withdraws two-thirds of
a cassette. ``bench/results/vector_safety_screen.md``.


Why the search radius is zero
-----------------------------

``max_subs=0`` — exact coincidence — because **the decision is per unit while the search is per
register, and that multiplies.** A 27-mer carries ~70 class-I registers and is withdrawn if any one
fires, so a per-register false-positive rate that reads as small is not the rate a cassette
experiences. Measured on six hazard-free 27-mers plus one burying the real titin epitope
(``bench/results/vector_screen_radius.md``) — **units falsely withdrawn, out of six**:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - ``max_subs``
     - 9-mers
     - 9-10-11
     - 8-9-10-11
   * - 0
     - 0
     - 0
     - 0
   * - 1
     - 1
     - 1
     - **4**

Radius 0 is clean at every length set; radius 1 is clean at none. It collapses once 8-mers enter: an
8-mer plus its 152 one-substitution neighbours is ~153 of 20\ :sup:`8` sequences against the
proteome's ~68 M windows, so a chance hit per register is expected and the ~20 8-mer registers in a
unit make it near-certain. **Every setting still catches the titin unit**, so radius 1 buys nothing
and costs most of the cassette.

.. warning::

   **What this does not catch**, stated because a safety screen that oversells itself is worse than
   none. It would not have caught the titin event *as it happened*: the construct contained MAGE-A3,
   whose profile is clean (13.4 TPM testis, 0.00 elsewhere), the cross-reactive titin peptide was
   never in it, and four TCR-facing substitutions separate the two — no distance threshold reaches
   one from the other. The affinity-enhanced TCR was the cause. What this catches is the adjacent
   and commoner failure.


The third clause: report near-identity, never withdraw on it
-------------------------------------------------------------

Both deaths were *near*-identity rather than identity, so an exact-only screen cannot be the whole
answer — and the section above is why a ``d = 1`` **veto** cannot be it either. ``report_subs=1``
(``mhcmatch vector --report-subs 1``) resolves that by separating the two decisions: a ``d = 1``
coincidence is **reported and the unit is kept**, arriving in ``screen``'s ``notes`` with
``"veto": False`` regardless of ``graded``, so it is a safety consideration attached to a candidate
rather than a refusal of it.

Left raw, that annotation is useless — two thirds of every cassette. Four filters, each answering a
different way of not being a hazard, take it to one unit in twelve. Measured end to end through the
shipped path on 178 experimentally immunogenic somatic neoantigens, rebuilt as 27-mer units
(``bench/results/vector_report_tier.md``):

.. list-table::
   :header-rows: 1
   :widths: 52 12 12 24

   * - layer
     - units
     - of 174
     - action
   * - **clause 2** — exact, different gene, mutation-spanning
     - 2
     - **1.1%**
     - withdrawn
   * - **clause 3** — d=1, different + expressed + non-homologous
     - 116
     - 66.7%
     - reported
   * - … and 9-11mers only
     - 27
     - 15.5%
     - reported
   * - … and the variant is itself presented
     - 14
     - **8.0%**
     - reported

**8-mers are the whole difference between 66.7% and 15.5%**, and for the same reason radius 1 is
refused above: an 8-mer's 152-neighbour ball against 68,398,087 proteome windows in 20\ :sup:`8`
expects **0.41** chance hits per register, where a 9-mer's 171 neighbours in 20\ :sup:`9` expect
**0.023** — 18× fewer. On this arm 8-mers report 101 units and 9-11mers report 25, and 76 units are
reported on the strength of an 8-mer *alone*. Exact matching keeps its 8-mers untouched: at ``d = 0``
an 8-mer expects 0.0027 hits, which is why ``max_subs=0`` can scan a length ``report_subs=1`` must
not.

The homology filter (:func:`~mhcmatch.vector.flank_identity`, cut at 0.5) is the second: a gene that
shares the unit's *flanks* as well as its register is related by descent, and a T cell that sees it
is one tolerance already had to deal with. The unit's 27-mer bounds the comparison at ±9-10
residues, so this separates loci rather than superfamilies — NRAS → KRAS survives it, correctly, and
is reported.

The last is presentation. :func:`~mhcmatch.vector.presented` asks whether the **off-target's own
sequence** is predicted presented on the allotype the unit was selected for; a variant no allotype
shows is a sequence coincidence, not a hazard. The cut is read off the positives rather than
borrowed — on this scorer the 176 assayed immunogenic peptides have a median of 0.69% rank, and
**30% rank keeps 97.2% of them** where the conventional 2% would discard three in ten. On a safety
read-out that is the expensive error, so the default is permissive by construction and it still
halves the tier, 27 units to 14.

.. note::

   **This tier annotates; it does not gate.** Nothing in it withdraws a unit, and with
   ``--weight-offtarget 0`` (the default) nothing in it changes composition either. What it changes
   is that a kept unit can now say *which* essential-tissue gene it sits one substitution from, at
   what expression, in what tissue, and whether that gene's own version is presented — which is the
   part of a safety argument that a clinician overrides or accepts.


Prior evidence: near-exact matches to known antigens
----------------------------------------------------

A candidate that is already a confirmed neoantigen, or sits one substitution from one, carries
direct evidence no model output can match. :mod:`mhcmatch.known` reports **exact** membership of
five reference sets and :mod:`mhcmatch.rank` floats those candidates into a tier of their own, with
the model score still shown beside them — burying "this is a confirmed NCI neoantigen" inside a
weighted sum lets a mediocre score dilute the one piece of direct evidence in the row.

``mhcmatch neoag`` and ``mhcmatch rank --annotate`` widen that to **near-exact**, reporting
``neoag_distance`` (substitutions to the nearest tested neoantigen), ``neoag_nearest`` (which one)
and ``neoag_n_within``. Fuzzy beats exact: held out, matching at ≤2 substitutions roughly doubles to
triples the recall of a fresh cohort's true positives over exact lookup, which is why
``--max-subs`` defaults to 2 and 1 is the tighter reading.

.. warning::

   **These columns are meaningless on a benchmark corpus and informative on a patient.** The
   reference sets *are* the deposits our own evaluation arms are built from, so a candidate drawn
   from those arms is at distance 0 to itself and the column is all-positive by construction. It
   only starts carrying information when the candidates come from somewhere the reference has never
   seen — which is exactly the delivery case. Read a filled ``neoag_distance`` column on internal
   data as a schema check, never as a result.

The same holds for the ``self`` and ``thymus`` mimicry channels, for the opposite reason: a hit
argues *reduced* immunogenicity (reactive T cells were plausibly deleted during negative selection)
**and** cross-reactivity risk, and those two conclusions pull in different directions. They are
reported, never folded into the score.


Coverage: how many units, and on which allotypes
-------------------------------------------------

:func:`mhcmatch.vector.select` spends capacity with a **per-allotype stopping rule**, not a global
budget. The clinical numbers (20 for autogene cevumeran, 34 for mRNA-4157, 20 in four pools for
NeoVax) derive from no published objective function, but the *shape* of the competition is
established: it is for the antigen-presenting cell rather than for MHC, it is strongest within an
allotype, and it can be net-positive across allotypes.

So the expected yield is a sum of independently saturating per-allotype terms, and the rule is one
line — **keep adding to an allotype while the next candidate's probability beats that allotype's
current expected yield per slot**:

.. math::

   E = \sum_a \frac{n_0\,S_a(n_a)}{n_0 + n_a},
   \qquad
   \text{add to } a \iff p_{a,n+1} > \frac{S_a(n_a)}{n_0 + n_a}

Because each allotype saturates on its own, a crowded allotype's next unit falls below an empty
allotype's first one, so **diversification across the patient's allotypes falls out of the
arithmetic** instead of being imposed as a quota. ``n0`` is per-allotype capacity, it is the one
free parameter, and it has **no default** — nothing in the public record fits it, so the value is
the caller's to defend and :attr:`~mhcmatch.vector.Selection.n0` carries it into the result so a
cassette can always name its own assumption.

Screening runs **before** selection, for the same reason: capacity spent on a unit that will be
withdrawn is capacity not spent on a safe one.

Once a set is selected, what that set is *worth* — coverage, redundancy, :math:`\Pr(\ge k)` — is
:doc:`portfolio`.


The cassette map: coordinates, linkers, and who helps whom
-----------------------------------------------------------

:func:`mhcmatch.vector.epitope_map` annotates the assembled cassette and
:func:`~mhcmatch.vector.write_map` writes it as a flat TSV and a JSON a viewer can draw from without
recomputing anything. One row per **unit**, **linker** and predicted **epitope**, in 1-based
inclusive coordinates over :attr:`~mhcmatch.vector.Cassette.sequence`:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - column
     - meaning
   * - ``id`` ``kind``
     - ``u1``/``l1``/``e1``; ``unit``, ``linker`` or ``epitope``
   * - ``start`` ``end`` ``length``
     - 1-based inclusive span. Units and linkers **tile the cassette exactly**
   * - ``seq``
     - the residues; ``sequence[start-1:end]`` reproduces it
   * - ``cls`` ``allele`` ``rank``
     - ``mhc1``/``mhc2``, the presenting allotype, its %rank
   * - ``unit`` ``gene``
     - which unit contains it; **0 means it spans a junction**
   * - ``core_start`` ``core_end``
     - class II only: the 9-mer register core, in cassette coordinates
   * - ``overlaps`` ``n_overlaps``
     - ids of the **other class's** epitopes overlapping this one

**A peptide presented by two of the recipient's alleles gets two rows.** At a heterozygous locus
those are two independent presentation events, they are what :func:`~mhcmatch.vector.select` spends
per-allotype capacity on, and collapsing them would under-count the coverage of exactly the patient
the cassette was personalised for.

**The overlap column is the point, not a decoration.** A cassette that carries a CD8 epitope and
borrows its CD4 help from an unrelated universal helper (PADRE, HBVcore) raises no class-II response
against the tumour antigen at all. Kissick *et al.* built one 27-mer around the HLA-A\*02:01
SIM2\ :sub:`237-245` epitope so a class-II epitope from the **same protein** overlapped it, and the
long peptide alone then replaced the exogenous HBVcore helper outright — a CD8 IFN-γ recall response
equal to the 9-mer-plus-helper, *and* a CD4 IL-2 response to SIM2\ :sub:`240-254`, with 137 class-II
binders predicted across DR/DP/DQ from that single 27-mer (*PLoS One* 2014;9(4):e93231, PMID
24690990, `doi:10.1371/journal.pone.0093231 <https://doi.org/10.1371/journal.pone.0093231>`_). So
:func:`~mhcmatch.vector.map_summary` reports per unit whether its class-I epitopes have overlapping
class-II epitopes — ``self_help`` — and a unit without it is the configuration that needed the
borrowed helper.

.. code-block:: python

   from mhcmatch import vector as V

   r1 = V.store_ranker(store_mhc1, ["HLA-A*02:01", "HLA-B*07:02"], cls="mhc1")
   r2 = V.store_ranker(store_mhc2, ["HLA-DRB1*01:01"],             cls="mhc2")
   feats   = V.epitope_map(cassette, r1, r2, threshold=2.0)
   summary = V.write_map(cassette, feats, "cassette.map.tsv", "cassette.map.json")
   summary["n_units_with_self_help"], summary["n_junction_spanning"]

``store_ranker`` is the per-allele adapter and is deliberately not
:func:`~mhcmatch.vector.store_binder`, which collapses to the best allele because a *layout* cost
only needs to know that some binder forms. Both rankers are injected, so the whole map is testable
with no panel and no download.

.. note::

   Map coordinates are over the **epitope cassette only** — no start codon, no stop, no leader, no
   trafficking domain, because those belong to the vector backbone. The next section adds them, and
   :func:`~mhcmatch.vector.mrna` returns its own parts map in nucleotides for exactly that reason.


Choosing a linker
-----------------

:data:`~mhcmatch.vector.LINKERS` is a table of named presets, so a construct format is named rather
than typed. ``mhcmatch cassette linkers`` prints it; ``--linker NAME`` pins one; any function taking
a spacer takes a name, explicit residues, or ``none``.

.. code-block:: bash

   mhcmatch cassette linkers
   mhcmatch cassette build --candidates units.tsv --n0 8 --linker GS10 --mrna cassette.mrna.fa

.. list-table::
   :header-rows: 1
   :widths: 20 16 14 50

   * - family
     - example
     - intended class
     - what it is trying to buy
   * - minimal
     - ``none``, ``G``
     - any
     - the shortest construct and the fewest junctional residues able to form a binder
   * - GS-rich flexible
     - ``GS10``, ``G4S``
     - any
     - flexible separation, no secondary structure of its own, low intrinsic immunogenicity
   * - class-I favouring
     - ``AAY``, ``AAA``
     - MHC-I
     - a preferred cleavage site at the seam, so the flanking units are liberated with the exact
       C-terminus class I needs
   * - class-II oriented
     - ``GPGPG``
     - MHC-II
     - independent processing of class-II epitopes
   * - protease-cleavable
     - ``furin``
     - any
     - cleavage placed rather than predicted
   * - rigid
     - ``EAAAK``
     - any
     - enforced spatial separation

**The table does not rank itself, and that is deliberate.** ``GPGPG`` is the only linker with causal
evidence behind it and the evidence is class II: unspaced concatenation of four HLA-DR epitopes
created a junctional epitope that suppressed the response to all four, and inserting ``GPGPG``
restored all four (Livingston *et al.*, *J Immunol* 2002, PMID 12023344). Against that, all six
orderings of three Fel d I regions produced no detectable junctional response at all (Rogers *et
al.*, *Mol Immunol* 1994, PMID 7521933). Both results are real.

For class I the two mechanisms that would settle a ranking act **at different positions**. Glycine
and proline are abundant in the C-terminal regions from which class-I ligands are cleaved
(Martin-Galiano & Lopez, *PLoS One* 2019, PMID 30645615), which is a processing argument for
them. The same residues immediately flanking a class-I epitope inhibit recognition of
the epitope on their *amino-terminal* side, and flanking context moved the ratio of two responses
from one construct by up to fifty-fold (Bergmann *et al.*, *J Immunol* 1996, PMID 8871618), which
argues against. A preset chosen on either citation alone is a linker chosen by reputation.

So the ``cls`` column records the class each linker is **intended** for, which is provenance, and
nothing here is a measurement of which one wins. :func:`~mhcmatch.vector.order` measures each
candidate against the recipient's own allotypes and that measurement is what selects — pinning with
``--linker`` is for the case where the construct format is already decided and only the ordering is
open.

.. warning::

   ``GS10`` is a **reconstruction**. The manufactured pentatope format is described in the
   methodological literature only as joining units with "non-immunogenic 10-mer glycine/serine
   linkers"; ``GGSGGGGSGG`` is the explicit 10-residue sequence reported for closely related RNA
   constructs in the patent literature, and it is not read off a published construct residue by
   residue. The compositional variants ``GS10b`` and ``GS10c`` are in the table because the format
   is specified by description rather than by residue. Anything conditional on the exact linker is
   conditional on that reconstruction.


Building the mRNA
-----------------

:func:`~mhcmatch.vector.mrna` assembles the molecule and returns what it assembled:

.. code-block:: python

   from mhcmatch import vector as V

   cas = V.order(units, binder=V.store_binder(store, alleles), linker="GS10")
   m   = V.mrna(cas, leader=SIGNAL, trailer=TRAFFICKING, utr5=U5, utr3=U3, poly_a=100)

   m.checks["translates"]                        # the frame survived every element
   m.checks["gc"], m.checks["longest_homopolymer"], m.checks["slippery_sites"]
   [p["name"] for p in m.of_kind("unit")]        # what is in it, in the order it is in
   m.as_rna()                                    # the same molecule in RNA letters

The construct is nine kinds of part, in this order, every one of them optional except the units:

.. list-table::
   :header-rows: 1
   :widths: 14 16 60

   * - part
     - alphabet
     - supplied by
   * - ``utr5``
     - nucleotides
     - the caller
   * - ``start``
     - —
     - the library, unless the reading frame already begins with M
   * - ``leader``
     - amino acids
     - the caller — a secretory signal peptide
   * - ``unit``
     - amino acids
     - the cassette
   * - ``linker``
     - amino acids
     - the chosen preset, between consecutive units
   * - ``trailer``
     - amino acids
     - the caller — an MHC class-I trafficking domain
   * - ``stop``
     - nucleotides
     - the library, one codon, default ``TGA``
   * - ``utr3``
     - nucleotides
     - the caller
   * - ``poly_a``
     - nucleotides
     - the caller, as a length

**What the library supplies, and what it refuses to.** The coding sequence is generated here, and
generated for the whole open reading frame **in one pass** — so homopolymer avoidance and the m1Ψ
+1-frameshift repair (:func:`~mhcmatch.vector.back_translate`, :func:`~mhcmatch.vector.deslip`) act
across the seams the designer created rather than inside each unit separately. That is not a detail:
a concatemer has far more codon-boundary junctions per kilobase than a natural ORF, the residues at
those junctions are chosen by the designer, and ``AAA`` back-translated is where a synthesis-hostile
run comes from. The backbone is **not** supplied and defaults to nothing, because a UTR belongs to a
particular vector and a plausible invented one is worse than none.

``parts`` **tiles the molecule exactly** — consecutive, non-overlapping, 0-based half-open, and the
concatenation of every part's slice is the whole sequence. An element that silently went missing
appears as a gap and one placed twice as an overlap, neither of which a length check would catch.

``checks`` reports numbers rather than a verdict, and ``translates`` is the one that must hold: the
coding sequence, read back in the frame the assembled construct sets it in, gives exactly
``protein``. It is what catches a frame broken by a supplied 5' element. The composition figures —
GC, longest homopolymer, slippery sites — are over the **coding sequence only**, since a poly(A)
tail is a homopolymer by construction and would swamp both.

.. note::

   This is not a codon optimiser and does not claim to be one. It fixes the two things that make a
   *polyepitope* construct fail where a natural ORF would not, and leaves GC content, secondary
   structure, splice sites and CpG alone. Where a manufacturer's own optimiser is available it
   should be preferred; this exists so a cassette ships with a usable coding sequence rather than
   none.

Running it
----------

.. code-block:: bash

   mhcmatch cassette build \
       --candidates ranked.tsv \
       --context windows.fasta \
       --n0 6 --alleles A*02:01,B*07:02 --cls mhc1 \
       --screen \
       --out cassette.tsv --fasta cassette.faa --fasta-nt cassette.fna \
       --map cassette.map.tsv --map-json cassette.map.json \
       --map-alleles-mhc2 DRB1*15:01,DRB1*07:01,DQB1*03:01

``--candidates`` and ``--context`` are both required and neither is redundant: ``rank`` emits
**minimal epitopes** and a vaccine unit is the long window around the mutation. Injecting a minimal
epitope is not a smaller version of the right thing — a 9-mer loads onto any cell without
costimulation and is the tolerising configuration (PMID 17911588) — so the reader refuses a table it
cannot tell apart rather than guessing.

``cassette.tsv`` carries a ``withdrawn`` section, one row per (unit, register, source gene) with the
clause that fired, the gene, the substitutions, the tissue and its TPM. **A withdrawn candidate has
to say what withdrew it**: "the screen dropped 3 of 40" is not a safety argument, and the reason is
what a clinician overrides or accepts.

Under ``--screen-mode graded`` it also carries a ``fingerprint`` section in the same shape, for the
units that were **kept**: every essential-tissue finding below ``--veto-tpm``, with the unit's total
off-target count in the index column. Why a unit was kept but discounted was previously answerable
only by re-running the screen.

.. note::

   ``--screen`` builds one whole-proteome window index **per register length**, and every unit's
   registers then resolve in a single query — so screen the whole candidate list in one process,
   never one unit per invocation. Measured on the human proteome (144,182 proteins, 69,486,637
   residues) with 3,000 units: **172.5 s and 11.1 GB peak** for all four class-I register lengths
   together. **Without the flag no safety check runs at all** and the cassette carries whatever it
   was handed. On a cluster this is the process that needs 48 GB; see
   ``integrations/nextflow/mhcmatch/README.md``.

A site with its own toxicity list substitutes the policy wholesale: ``screen`` takes
``risk(unit, registers) -> [reason, ...]`` as an argument, so the shipped
:func:`~mhcmatch.vector.self_origin_risk` is one implementation and not the interface.
