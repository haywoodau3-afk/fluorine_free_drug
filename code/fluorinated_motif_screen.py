"""OECD 2021 structural PFAS motif screening utilities.

The primary filter looks for a saturated carbon bearing two or three fluorines,
with no attached hydrogen or chlorine, bromine, or iodine. It scans every atom
and every disconnected component in the submitted molecular structure.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from rdkit import Chem


@dataclass(frozen=True)
class MotifMatch:
    atom_index: int
    motif: str
    fluorine_neighbors: tuple[int, ...]


@dataclass(frozen=True)
class ScreeningResult:
    valid_smiles: bool
    contains_pfas_motif: bool | None
    fully_fluorine_free: bool | None
    matches: tuple[MotifMatch, ...]
    error: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "valid_smiles": self.valid_smiles,
            "contains_pfas_motif": self.contains_pfas_motif,
            "fully_fluorine_free": self.fully_fluorine_free,
            "matches": [asdict(match) for match in self.matches],
            "error": self.error,
        }


_EXCLUDED_NEIGHBOR_ATOMIC_NUMBERS = {1, 17, 35, 53}


def screen_smiles(smiles: str) -> ScreeningResult:
    """Screen an entire SMILES, including all salt/disconnected components.

    Aryl C-F bonds do not qualify because aromatic and unsaturated carbons are
    excluded from the fully fluorinated saturated methyl/methylene motif.
    Invalid input returns an explicit indeterminate result.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ScreeningResult(False, None, None, (), "RDKit could not parse SMILES")

    matches: list[MotifMatch] = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            continue
        bonds = atom.GetBonds()
        if any(bond.GetBondType() != Chem.BondType.SINGLE for bond in bonds):
            continue

        fluorine_neighbors = tuple(
            neighbor.GetIdx()
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 9
        )
        fluorine_count = len(fluorine_neighbors)
        if fluorine_count not in (2, 3) or atom.GetTotalNumHs() != 0:
            continue

        other_neighbors = [
            neighbor
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 9
        ]
        if len(other_neighbors) != 4 - fluorine_count:
            continue
        if any(
            neighbor.GetAtomicNum() in _EXCLUDED_NEIGHBOR_ATOMIC_NUMBERS
            for neighbor in other_neighbors
        ):
            continue

        matches.append(
            MotifMatch(
                atom_index=atom.GetIdx(),
                motif="CF3" if fluorine_count == 3 else "CF2",
                fluorine_neighbors=fluorine_neighbors,
            )
        )

    return ScreeningResult(
        valid_smiles=True,
        contains_pfas_motif=bool(matches),
        fully_fluorine_free=not any(atom.GetAtomicNum() == 9 for atom in molecule.GetAtoms()),
        matches=tuple(matches),
    )
