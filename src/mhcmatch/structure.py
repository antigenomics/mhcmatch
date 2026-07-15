"""Structure-based pMHC Miyazawa–Jernigan contact energy and WT/MT ΔΔG (**optional** ``tcren`` dep).

Threads a peptide onto a template pMHC crystal groove and sums the MJ contact potential (all shipped
by ``tcren``). For a query allele with no own template we borrow the **groove-closest** template
(the same pseudosequence kernel the diffusion uses). For a single-mutation WT/MT pair the
``ΔΔG = MJ(mut) − MJ(wt)`` on one backbone is a physics-based differential affinity estimator
(adaptive double threading, Jojic et al. 2006) -- the structural complement to the sequence-based
:mod:`mhcmatch.affinity`.

On measured HLA-A*02:01 the MJ energy tracks log-IC50 at Spearman ≈ 0.55 (see
``bench/affinity/bench_structure.py``), at ~0.02 ms/peptide after a one-time template build.

Needs the ``[structure]`` extra::

    pip install 'mhcmatch[structure]'      # pulls tcren

Template structures are **not vendored** (they live in ``tcren``'s Canonical2026 set); point
``structure_dir`` at them (default: ``$MHCMATCH_STRUCTURES`` or the sibling ``tcren-ms`` checkout).
"""
from __future__ import annotations

import json
import os
from importlib import resources

# normalized allele -> template. ``chains`` maps PDB chain id -> role (canonical Canonical2026:
# C=peptide, D=MHCα, E=β2m). Curated + extensible; see data/structure_templates.json.
_DEFAULT_DIR = os.environ.get(
    "MHCMATCH_STRUCTURES",
    os.path.expanduser("~/vcs/code/tcren-ms/data/Canonical2026"))


def _require_tcren():
    try:
        import tcren  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "mhcmatch.structure needs tcren -- install the optional extra: "
            "pip install 'mhcmatch[structure]'") from exc


def _load_templates():
    src = resources.files("mhcmatch.data").joinpath("structure_templates.json")
    return json.loads(src.read_text()) if src.is_file() else {}


class StructureScorer:
    """MJ contact-energy scorer over template pMHC structures.

    ``pseudoseq`` (a :class:`mhcmatch.Pseudoseq`) enables borrowing the groove-closest template for
    alleles without their own; omit it to restrict to exact-allele templates.
    """

    def __init__(self, structure_dir=None, templates=None, pseudoseq=None, cutoff=5.0):
        _require_tcren()
        from .pseudoseq import normalize_allele
        self._dir = structure_dir or _DEFAULT_DIR
        self._tpl = {normalize_allele(a): v for a, v in (templates or _load_templates()).items()}
        self._ps = pseudoseq
        self._cutoff = cutoff
        self._cm = {}          # pdb -> ContactMap (built once)
        self._norm = normalize_allele

    def template_for(self, allele, length):
        """``(pdb, chains)`` of the best template for ``allele`` at peptide ``length``: exact allele
        if present, else the groove-closest templated allele (needs ``pseudoseq``). ``None`` if none."""
        cand = [(a, t) for a, t in self._tpl.items() if t.get("length") == length]
        if not cand:
            return None
        a0 = self._norm(allele)
        for a, t in cand:
            if a == a0:
                return t["pdb"], t.get("chains", {})
        if self._ps is None:
            return None
        best = max(cand, key=lambda at: self._ps.kernel(a0, at[0]))   # groove-closest
        return best[1]["pdb"], best[1].get("chains", {})

    def _contactmap(self, pdb, chains):
        if pdb not in self._cm:
            from tcren.structure.io import import_structure
            from tcren.contactmap import ContactMap
            from tcren.annotation.chains import _tag_peptide
            s = import_structure(os.path.join(self._dir, f"{pdb}.pdb.gz"))
            roles = chains or {"C": "peptide", "D": "MHCa", "E": "B2M"}
            for c in s.chains:
                role = roles.get(c.chain_id)
                if role == "peptide":
                    _tag_peptide(c)
                elif role:
                    c.chain_type = role
            self._cm[pdb] = ContactMap.from_structure(s, cutoff=self._cutoff)
        return self._cm[pdb]

    def mj_energies(self, peptides, allele):
        """``{peptide: MJ energy}`` for equal-length ``peptides`` on ``allele``'s template (lower =
        stronger binding). Empty dict if no matching-length template. One batch swap."""
        peptides = [p.strip().upper() for p in peptides]
        if not peptides:
            return {}
        tpl = self.template_for(allele, len(peptides[0]))
        if tpl is None:
            return {}
        from tcren.scoring import score_peptides
        from tcren.potential import mj
        cm = self._contactmap(*tpl)
        df = score_peptides(cm, peptides, mj(), interface="peptide_mhc", substituted_side="from")
        return dict(zip(df["peptide"].to_list(), df["score"].to_list()))

    def mj_energy(self, peptide, allele):
        """MJ contact energy of ``peptide`` on ``allele``'s template (``nan`` if untemplated)."""
        return self.mj_energies([peptide], allele).get(peptide.strip().upper(), float("nan"))

    def ddg(self, wt, mut, allele):
        """``MJ(mut) − MJ(wt)`` on the shared backbone (the WT/MT differential; <0 = mutant binds
        better). ``wt`` and ``mut`` must be the same length. ``nan`` if untemplated."""
        e = self.mj_energies([wt, mut], allele)
        ew, em = e.get(wt.strip().upper()), e.get(mut.strip().upper())
        return (em - ew) if (ew is not None and em is not None) else float("nan")


if __name__ == "__main__":  # pragma: no cover - needs the [structure] extra + a template on disk
    try:
        _require_tcren()
    except ImportError as e:
        print(f"skip: {e}")
        raise SystemExit(0)
    sc = StructureScorer()
    if sc.template_for("HLA-A*02:01", 9) is None:
        print("skip: no HLA-A*02:01 template on disk (set MHCMATCH_STRUCTURES)")
        raise SystemExit(0)
    e = sc.mj_energies(["GILGFVFTL", "GILGFVFTK", "AAAAAAAAA"], "HLA-A*02:01")
    assert e["GILGFVFTL"] < e["GILGFVFTK"] < e["AAAAAAAAA"], e   # native < bad-anchor < poly-Ala
    assert sc.ddg("GILGFVFTL", "GILGFVFTK", "HLA-A*02:01") > 0   # K anchor destabilizes -> ΔΔG>0
    print(f"structure.py self-check OK  MJ(GILGFVFTL)={e['GILGFVFTL']:.2f} < "
          f"MJ(...K)={e['GILGFVFTK']:.2f} < MJ(polyA)={e['AAAAAAAAA']:.2f}")
