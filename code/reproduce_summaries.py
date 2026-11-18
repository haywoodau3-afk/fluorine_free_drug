"""Recompute article counts and a first-pose sensitivity from the packaged records."""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import json
import math

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
from sklearn.metrics import average_precision_score

from fluorinated_motif_screen import screen_smiles

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "reproduced"


def read(relative):
    with (ROOT / relative).open(newline="") as stream:
        return list(csv.DictReader(stream))


def put_csv(relative, rows):
    with (OUTPUT / relative).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pdb_atoms(path, first_model=False):
    for line in path.read_text().splitlines():
        if first_model and line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")):
            residue = f"{line[21:22].strip()}:{line[17:20].strip()}:{line[22:26].strip()}"
            xyz = [float(line[start:start+8]) for start in (30, 38, 46)]
            yield residue, xyz


def docking_sensitivity(zero):
    docking_rows = {row["design_id"]: row for row in read("data/docking/results.csv")}
    by_drug = defaultdict(list)
    for row in zero:
        by_drug[row["drug_id"]].append(row)

    sensitivity = []
    for drug, designs in by_drug.items():
        receptor_atoms = list(pdb_atoms(ROOT / f"inputs/docking/receptors/{drug}.pdb"))
        residues, xyz = zip(*receptor_atoms)
        tree = cKDTree(xyz)
        comparisons = []
        for row in designs:
            design = row["design_id"]
            result = docking_rows[design]
            parent = set(result["parent_contact_residues"].split(";"))
            path = ROOT / f"inputs/docking/poses/{drug}/{design}.pdbqt"
            overlaps = []
            for first in (False, True):
                coords = [point for _, point in pdb_atoms(path, first)]
                indices = tree.query_ball_point(coords, 4.0)
                contacts = {residues[index] for group in indices for index in group}
                overlaps.append(len(parent & contacts) / len(parent | contacts))
            assert abs(overlaps[0] - float(row["contact_jaccard"])) <= 0.000051, design
            comparisons.append(dict(drug_id=drug, design_id=design,
                reported_union_pose_rank=int(row["zero_fluorine_rank"]),
                union_pose_jaccard=round(overlaps[0], 4),
                first_pose_jaccard=round(overlaps[1], 4),
                vina_kcal_mol=float(row["vina_kcal_mol"])))
        comparisons.sort(key=lambda row: (-row["first_pose_jaccard"], row["vina_kcal_mol"], row["design_id"]))
        for rank, row in enumerate(comparisons, 1):
            row["first_pose_rank"] = rank
        sensitivity.extend(comparisons)
    return sensitivity


def tetrahedral_centres(molecule, identities):
    copy = Chem.Mol(molecule)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return {identities[centre.centeredOn] for centre in Chem.FindPotentialStereo(copy)
            if centre.type == Chem.StereoType.Atom_Tetrahedral}


def classify_topology(strict):
    structures = Chem.SDMolSupplier(str(ROOT / "data/atlas/mapped_atlas_and_parents.sdf"))
    parents = {m.GetProp("drug_id"): m for m in structures if m.GetProp("_Name").startswith("PARENT-")}
    burden = {r["state_id"]: r for r in read("data/atlas/edit_burden.csv")}
    rows = []
    for graph in read("data/atlas/mapped_graph_edits.csv"):
        state = graph["state_id"]
        if state not in strict:
            continue
        drug = graph["drug_id"]
        parent = parents[drug]
        product = Chem.MolFromSmiles(graph["product_smiles_mapped"])
        parent_ids = {a.GetIdx(): a.GetIdx() + 1 for a in parent.GetAtoms()}
        product_ids = {a.GetIdx(): a.GetAtomMapNum() or -(a.GetIdx() + 1) for a in product.GetAtoms()}
        retained_maps = {identity for identity in product_ids.values() if identity > 0}
        parent_ring_atoms = {parent_ids[a.GetIdx()] for a in parent.GetAtoms() if a.IsInRing()}
        product_ring_atoms = {product_ids[a.GetIdx()] for a in product.GetAtoms() if a.IsInRing()}
        ring_change = parent_ring_atoms != product_ring_atoms
        stereo_change = tetrahedral_centres(parent, parent_ids) != tetrahedral_centres(product, product_ids)
        deleted_framework = int(burden[state]["deleted_other_nonfluorine_heavy_atom_count"])
        site_count = int(burden[state]["parent_cf3_site_count"])
        anchors = {n.GetAtomMapNum() for a in product.GetAtoms() if a.GetAtomMapNum() == 0
                   for n in a.GetNeighbors() if n.GetAtomMapNum() > 0}
        anchors.update(i + 1 for bond in json.loads(graph["broken_parent_bonds_json"])
                       for i in bond["parent_atom_indices"] if i + 1 in retained_maps)
        anchors.update(i + 1 for i in json.loads(burden[state]["mapped_cf3_carbon_element_change_indices_json"])
                       if i + 1 in retained_maps)
        if ring_change or stereo_change or deleted_framework > 4:
            topology = "ring or framework"
        elif len(anchors) >= 2 or site_count >= 2:
            topology = "multi-anchor or multi-site"
        else:
            topology = "one-anchor"
        rows.append(dict(state_id=state, drug_id=drug, topology=topology, ring_atom_change=ring_change,
                         potential_tetrahedral_centre_change=stereo_change, deleted_framework_atoms=deleted_framework,
                         retained_anchor_count=len(anchors), original_cf3_sites=site_count,
                         edit_burden=burden[state]["total_deleted_plus_added_heavy_atoms"]))
    return rows


def summarize():
    membership = read("data/atlas/atlas_membership.csv")
    strict = {r["state_id"]: r for r in membership if r["primary_cohort"] == "strict_194_state"}
    assert len(membership) == 215 and len(strict) == 194
    assert len({r["connectivity_family_id"] for r in strict.values()}) == 186
    topology = classify_topology(strict)
    expected = {"one-anchor": 48, "multi-anchor or multi-site": 100, "ring or framework": 46}
    assert Counter(r["topology"] for r in topology) == expected
    assert len(read("data/atlas/protected_mapping_checks.csv")) == 220
    assert len(read("data/atlas/stereochemical_exclusions.csv")) == 21
    proposals = read("data/atlas/expansion_proposals.csv")
    assert len(proposals) == 781
    assert sum(r["chemical_disposition"] == "accepted_chemistry_only_design" for r in proposals) == 646
    fragments = [r for r in read("data/atlas/cf3_site_fragments.csv") if r["state_id"] in strict]
    sites = Counter(r["fragment_status"] for r in fragments)
    assert sorted(sites.values()) == [7, 22, 262]
    normalized = set()
    for row in fragments:
        if row["fragment_status"] != "resolved_attachment_labelled_fragment":
            continue
        fragment = Chem.MolFromSmiles(row["replacement_fragment_smiles_attachment_labelled"])
        for atom in fragment.GetAtoms():
            atom.SetAtomMapNum(0)
            if atom.GetAtomicNum() == 0:
                atom.SetIsotope(0)
        normalized.add(Chem.MolToSmiles(fragment, isomericSmiles=True))
    assert len(normalized) == 34

    route = read("data/retrosynthesis/design_outcomes.csv")
    zero = read("data/retrosynthesis/zero_fluorine_designs.csv")
    assert len(route) == 926 and len(zero) == 42
    screens = {r["design_id"]: screen_smiles(r["canonical_isomeric_smiles"]) for r in route}
    assert all(r.valid_smiles and r.contains_pfas_motif is False for r in screens.values())
    assert {key for key, result in screens.items() if result.fully_fluorine_free} == {r["design_id"] for r in zero}
    assert sum(r["stock_terminating_route"] == "True" for r in route) == 134
    assert sum(r["stock_terminating_route"] == "True" for r in zero) == 24
    assert all(r["docked"] == "True" for r in zero)
    assert Counter(r["exhaustiveness"] for r in read("data/docking/executed_runs.csv")) == {"2": 287, "8": 100}

    matched = [r for r in read("data/generator-comparisons/template_task_comparison.csv") if r["matched_allocation"]=="True"]
    checks = dict(enumeration_slots=312,random_control_slots=312,enumeration_accepted=184,
                  random_control_accepted=204,enumeration_unique_families=183,random_control_unique_families=179)
    assert len(matched)==11
    assert all(sum(int(r[k]) for r in matched) == n for k,n in checks.items())

    predictions = read("data/dpp4/held_out_predictions.csv")
    folds = defaultdict(list)
    for row in predictions:
        folds[row["fold"]].append(row)
    ap_rf, ap_nn, prevalence = [], [], []
    hits = selection = 0
    for fold in folds.values():
        labels = [int(r["observed_active_at_100nM"]) for r in fold]
        ap_rf.append(average_precision_score(labels, [float(r["random_forest_active_probability"]) for r in fold]))
        ap_nn.append(average_precision_score(labels, [float(r["nearest_neighbor_activity_score"]) for r in fold]))
        prevalence.append(np.mean(labels))
        n = math.ceil(.05*len(fold))
        shortlist = sorted(fold,key=lambda r:(-float(r["random_forest_active_probability"]),r["parent_molecule_chembl_id"]))[:n]
        hits += sum(int(r["observed_active_at_100nM"]) for r in shortlist)
        selection += n
    assert (len(predictions),hits,selection)==(3098,117,156)
    properties=Counter(r["candidate_parent_comparison_status"] for r in read("data/lipophilicity/candidate_parent_summary.csv"))
    assert properties == {"abstain_parent_reference_missing_or_ambiguous":148,"fail_model_resolution_comparison":46}
    scope=Counter(r["drug_id"] for r in strict.values())
    return dict(strict_atlas_by_parent=scope, topology_counts=expected, site_observations=sites, distinct_nonempty_fragments=len(normalized),
                target_aware_designs=len(route), zero_fluorine_designs=len(zero), stock_terminating_full=134,
                stock_terminating_zero_fluorine=24, held_out_dpp4=dict(n=len(predictions),top5_actives=hits,
                shortlist_size=selection,random_forest_ap=float(np.mean(ap_rf)),nearest_neighbour_ap=float(np.mean(ap_nn)),
                prevalence=float(np.mean(prevalence))), lipophilicity=properties), docking_sensitivity(zero), topology


if __name__ == "__main__":
    OUTPUT.mkdir(exist_ok=True)
    summary, sensitivity, topology = summarize()
    (OUTPUT / "article_checks.json").write_text(json.dumps(summary, indent=2) + "\n")
    put_csv("first_pose_sensitivity.csv", sensitivity)
    put_csv("strict_atlas_topology.csv", topology)
    print(json.dumps(summary, indent=2))
