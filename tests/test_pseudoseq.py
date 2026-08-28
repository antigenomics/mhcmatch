"""Allele-name resolution: one molecule must resolve to exactly one key.

The pseudosequence tables carry several spellings of the same molecule on purpose -- deposited
screens disagree, so ``HLA-A02:01`` and ``HLA-A0201`` are both keys on the same 34-mer. That is
useful for *accepting* input and dangerous for *keying* anything, which is what these tests guard.
"""
# -- one molecule, one key (the mouse class-I allele collision) ---------------

def test_all_three_mouse_class1_spellings_fold_to_one_key():
    """``H-2Kb`` (pmhc), ``H2-Kb`` (deposits) and ``H-2-Kb`` (FASTA) are one molecule.

    ``mhci_pseudo.fa`` carries the last two as *separate* keys on a byte-identical 34-mer, so before
    1.4.1 ``resolve_allele('H2-Kb', 'mhc1')`` returned ``('H2-Kb', True)`` -- exact, and backed by
    zero panel ligands. The panel is keyed on the pmhc spelling, so the presentation head silently
    fell back to kernel shrinkage: SIINFEKL scored presentation %rank 20.19 under ``'H-2-Kb'``
    against 0.0040 under ``'H-2Kb'``, a 5,000x move driven by nothing but how the caller typed it.

    Class II has had this rule since :func:`class2_from_name`; human class I has it as
    :func:`hla_spellings`. This is the class-I mouse corner that had neither.
    """
    from mhcmatch.pseudoseq import normalize_allele, resolve_allele
    for stem in ("Kb", "Db", "Kd", "Ld", "Dd", "Kk", "Dq"):
        keys = {normalize_allele(f"H-2{stem}"), normalize_allele(f"H2-{stem}"),
                normalize_allele(f"H-2-{stem}")}
        assert keys == {f"H-2-{stem}"}, f"H-2{stem} spellings did not fold: {keys}"
        resolved = {resolve_allele(f"H-2{stem}", "mhc1"), resolve_allele(f"H2-{stem}", "mhc1"),
                    resolve_allele(f"H-2-{stem}", "mhc1")}
        assert resolved == {(f"H-2-{stem}", True)}, f"H-2{stem} did not resolve alike: {resolved}"


def test_the_two_mouse_fasta_spellings_carry_the_same_groove():
    """The premise of the fold: ``H2-Kb`` and ``H-2-Kb`` are not two molecules.

    If a future ``mhci_pseudo.fa`` ever gave them different 34-mers, folding them would be wrong and
    this test says so before the fold silently merges two grooves.
    """
    from mhcmatch.pseudoseq import load_pseudo
    seqs = load_pseudo("mhc1")
    for stem in ("Kb", "Db", "Kd", "Ld"):
        a, b = f"H-2-{stem}", f"H2-{stem}"
        if a in seqs and b in seqs:
            assert seqs[a] == seqs[b], f"{a} and {b} carry different pseudosequences"
