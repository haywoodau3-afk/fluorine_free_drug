"""Recalculate the reported 1X70 symmetry-aware aligned heavy-atom RMSDs."""
from pathlib import Path
import json

from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT = Path(__file__).resolve().parents[1]


def main():
    raw_reference = Chem.SDMolSupplier(str(ROOT / "inputs/reference-redocking/1x70_bound_ligand.sdf"), removeHs=False)[0]
    reference = Chem.RemoveHs(raw_reference)
    runs = json.loads((ROOT / "data/docking/sitagliptin_redocking.json").read_text())
    for run in runs:
        parsed = PDBQTMolecule.from_file(str(ROOT / run["docking_output"]))
        docked = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(parsed)[0])
        rmsds = [rdMolAlign.GetBestRMS(docked, reference, prbId=i, refId=0, numThreads=1)
                 for i in range(docked.GetNumConformers())]
        archived = run["symmetry_aware_heavy_atom_rmsd_angstrom"]
        assert len(rmsds) == len(archived)
        assert all(abs(a-b) < 1e-6 for a, b in zip(rmsds, archived))
        print(f"1X70 seed {run['seed']}: first-pose aligned heavy-atom RMSD {rmsds[0]:.3f} Å")


if __name__ == "__main__":
    main()
