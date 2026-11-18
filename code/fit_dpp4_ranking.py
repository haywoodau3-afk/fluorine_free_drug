#!/usr/bin/env python3
"""Run the frozen publication- and scaffold-held-out target-wide DPP-4 benchmark."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import rdkit
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import hypergeom
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = ROOT / "methods/dpp4_ranking.json"
FOLD_MANIFEST_PATH = ROOT / "methods/dpp4_folds.json"
ACTIVITY_PATH = ROOT / "inputs/chembl/dpp4_benchmark_input.csv"
OUTPUT = ROOT / "reproduced/dpp4"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_records(protocol: dict[str, object]) -> list[dict[str, object]]:
    with ACTIVITY_PATH.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    endpoint = protocol["eligible_records"]
    excluded_docs = set(protocol["excluded_prior_test_documents"])
    threshold = endpoint["positive_threshold_inclusive_nM"]
    eligible_rows = []
    for row in source_rows:
        if (
            row["target_chembl_id"] != endpoint["target_chembl_id"]
            or row["target_organism"] != endpoint["target_organism"]
            or row["standard_type"] != "IC50"
            or row["standard_relation"] != "="
            or row["standard_units"] != "nM"
            or row["parent_pfas_filter"] != "False"
            or not row["parent_canonical_smiles"]
            or row["document_chembl_id"] in excluded_docs
        ):
            continue
        molecule = Chem.MolFromSmiles(row["parent_canonical_smiles"])
        if molecule is None:
            continue
        eligible_rows.append({
            **row,
            "_mol": molecule,
            "_label": int(float(row["standard_value"]) <= threshold),
            "_scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=False),
        })

    by_parent = defaultdict(list)
    for row in eligible_rows:
        by_parent[row["parent_molecule_chembl_id"]].append(row)
    duplicates = sorted(parent_id for parent_id, rows in by_parent.items() if len(rows) > 1)
    unique = sorted(
        (rows[0] for parent_id, rows in by_parent.items() if len(rows) == 1),
        key=lambda row: row["parent_molecule_chembl_id"],
    )
    return unique, duplicates


def score_metrics(rows, scores: np.ndarray, ranking_scores: np.ndarray, top_k: int) -> dict[str, object]:
    labels = np.asarray([row["_label"] for row in rows], dtype=np.uint8)
    order = sorted(
        range(len(rows)),
        key=lambda index: (-float(ranking_scores[index]), rows[index]["parent_molecule_chembl_id"]),
    )
    top_indices = order[:top_k]
    positives = int(labels.sum())
    recovered = int(labels[top_indices].sum())
    return {
        "parent_count": len(rows),
        "positive_count": positives,
        "negative_count": len(rows) - positives,
        "prevalence": positives / len(rows),
        "average_precision": float(average_precision_score(labels, scores)),
        "roc_auc": float(roc_auc_score(labels, scores)),
        "shortlist_size": top_k,
        "recovery_at_shortlist": recovered,
        "precision_at_shortlist": recovered / top_k,
        "recall_at_shortlist": recovered / positives,
        "random_expected_recovery_at_shortlist": top_k * positives / len(rows),
        "random_probability_at_least_observed_recovery": float(
            hypergeom.sf(recovered - 1, len(rows), positives, top_k)
        ),
        "shortlist_parent_ids": [rows[index]["parent_molecule_chembl_id"] for index in top_indices],
    }


def main() -> None:
    protocol = json.loads(PROTOCOL_PATH.read_text())
    fold_manifest = json.loads(FOLD_MANIFEST_PATH.read_text())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    records, duplicate_ids = load_records(protocol)
    if len(records) != fold_manifest["deduplicated_unique_parent_count"]:
        raise ValueError("Deduplicated target-wide parent count differs from the frozen fold manifest")
    if duplicate_ids != fold_manifest["removed_duplicate_parent_ids"]:
        raise ValueError("Duplicate parent exclusion set differs from the frozen fold manifest")

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=True
    )
    fingerprints = [generator.GetFingerprint(row["_mol"]) for row in records]
    by_document = defaultdict(list)
    for index, row in enumerate(records):
        by_document[row["document_chembl_id"]].append(index)
    labels = np.asarray([row["_label"] for row in records], dtype=np.uint8)
    x = np.asarray([[int(bit) for bit in fingerprint.ToBitString()] for fingerprint in fingerprints], dtype=np.uint8)
    model_parameters = {
        key: value for key, value in protocol["primary_model"].items()
        if key not in {"estimator", "hyperparameter_tuning"}
    }

    folds_report = {}
    prediction_rows = []
    for fold in fold_manifest["folds"]:
        fold_number = fold["fold"]
        test_documents = set(fold["test_document_chembl_ids"])
        test_indices = sorted(index for doc_id in sorted(test_documents) for index in by_document[doc_id])
        if set(records[index]["parent_molecule_chembl_id"] for index in test_indices) != set(fold["test_parent_ids"]):
            raise ValueError(f"Fold {fold_number} test parent membership changed")
        test_parents = {records[index]["parent_molecule_chembl_id"] for index in test_indices}
        test_scaffolds = {records[index]["_scaffold"] for index in test_indices}
        test_fingerprints = [fingerprints[index] for index in test_indices]

        training_candidates = [
            index for index, row in enumerate(records)
            if row["document_chembl_id"] not in test_documents
            and row["parent_molecule_chembl_id"] not in test_parents
            and row["_scaffold"] not in test_scaffolds
        ]
        max_similarity = {}
        for index in training_candidates:
            max_similarity[index] = max(
                DataStructs.BulkTanimotoSimilarity(fingerprints[index], test_fingerprints)
            )
        train_indices = [
            index for index in training_candidates
            if max_similarity[index] <= protocol["split"]["maximum_train_test_morgan_tanimoto"]
        ]
        train_scaffolds = {records[index]["_scaffold"] for index in train_indices}
        if test_scaffolds & train_scaffolds:
            raise ValueError(f"Fold {fold_number} has Murcko scaffold overlap")
        if test_parents & {records[index]["parent_molecule_chembl_id"] for index in train_indices}:
            raise ValueError(f"Fold {fold_number} has parent overlap")
        if not train_indices:
            raise ValueError(f"Fold {fold_number} has no training parents")
        max_retained_similarity = max(max_similarity[index] for index in train_indices)
        if max_retained_similarity > protocol["split"]["maximum_train_test_morgan_tanimoto"]:
            raise ValueError(f"Fold {fold_number} exceeds the frozen similarity cap")

        estimator = RandomForestClassifier(**model_parameters)
        estimator.fit(x[train_indices], labels[train_indices])
        rf_scores = estimator.predict_proba(x[test_indices])[:, list(estimator.classes_).index(1)]
        nn_scores = []
        nn_similarities = []
        for test_fp in test_fingerprints:
            similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(test_fp, [fingerprints[i] for i in train_indices]), dtype=float)
            maximum = float(similarities.max())
            ties = [int(labels[train_indices[j]]) for j, value in enumerate(similarities) if value == maximum]
            nn_scores.append(float(np.mean(ties)))
            nn_similarities.append(maximum)
        nn_scores = np.asarray(nn_scores)
        test_rows = [records[index] for index in test_indices]
        top_k = max(1, math.ceil(0.05 * len(test_rows)))
        rf_metrics = score_metrics(test_rows, rf_scores, rf_scores, top_k)
        nn_metrics = score_metrics(test_rows, nn_scores, nn_scores, top_k)
        nn_order = sorted(
            range(len(test_rows)),
            key=lambda local_index: (-float(nn_scores[local_index]), -float(nn_similarities[local_index]), test_rows[local_index]["parent_molecule_chembl_id"]),
        )
        nn_top = nn_order[:top_k]
        nn_recovered = sum(test_rows[index]["_label"] for index in nn_top)
        nn_metrics["recovery_at_shortlist"] = int(nn_recovered)
        nn_metrics["precision_at_shortlist"] = nn_recovered / top_k
        nn_metrics["recall_at_shortlist"] = nn_recovered / nn_metrics["positive_count"]
        nn_metrics["shortlist_parent_ids"] = [test_rows[index]["parent_molecule_chembl_id"] for index in nn_top]
        folds_report[str(fold_number)] = {
            "held_out_document_count": len(test_documents),
            "training_parent_count": len(train_indices),
            "training_positive_count": int(labels[train_indices].sum()),
            "training_negative_count": int(len(train_indices) - labels[train_indices].sum()),
            "training_scaffold_count": len(train_scaffolds),
            "test_scaffold_count": len(test_scaffolds),
            "maximum_train_test_tanimoto_after_pruning": max_retained_similarity,
            "test_records_with_validity_comments": sum(bool(row["data_validity_comment"]) for row in test_rows),
            "random_forest": rf_metrics,
            "nearest_neighbor": nn_metrics,
        }

        for local_index, row in enumerate(test_rows):
            prediction_rows.append({
                "fold": fold_number,
                "assay_chembl_id": row["assay_chembl_id"],
                "document_chembl_id": row["document_chembl_id"],
                "parent_molecule_chembl_id": row["parent_molecule_chembl_id"],
                "parent_canonical_smiles": row["parent_canonical_smiles"],
                "observed_ic50_nM": row["standard_value"],
                "observed_active_at_100nM": row["_label"],
                "data_validity_comment": row["data_validity_comment"],
                "murcko_scaffold": row["_scaffold"],
                "random_forest_active_probability": f"{rf_scores[local_index]:.8f}",
                "nearest_neighbor_activity_score": f"{nn_scores[local_index]:.8f}",
                "nearest_neighbor_max_tanimoto": f"{nn_similarities[local_index]:.8f}",
            })

    macro_ap = float(np.mean([fold["random_forest"]["average_precision"] for fold in folds_report.values()]))
    macro_nn_ap = float(np.mean([fold["nearest_neighbor"]["average_precision"] for fold in folds_report.values()]))
    macro_prevalence = float(np.mean([fold["random_forest"]["prevalence"] for fold in folds_report.values()]))
    folds_above_prevalence = sum(
        fold["random_forest"]["average_precision"] > fold["random_forest"]["prevalence"]
        for fold in folds_report.values()
    )
    pooled_recovery = sum(fold["random_forest"]["recovery_at_shortlist"] for fold in folds_report.values())
    pooled_random_expected = sum(fold["random_forest"]["random_expected_recovery_at_shortlist"] for fold in folds_report.values())
    qualifies = (
        macro_ap > macro_nn_ap
        and macro_ap > macro_prevalence
        and folds_above_prevalence >= 4
        and pooled_recovery >= 1.25 * pooled_random_expected
    )
    report = {
        "protocol": str(PROTOCOL_PATH.relative_to(ROOT)),
        "fold_manifest": str(FOLD_MANIFEST_PATH.relative_to(ROOT)),
        "status": "completed",
        "model_versions": {
            "python": sys.version.split()[0],
            "rdkit": rdkit.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "dataset_counts": {
            "unique_parent_count": len(records),
            "positive_count": int(labels.sum()),
            "negative_count": int(len(labels) - labels.sum()),
            "source_document_count": len(set(row["document_chembl_id"] for row in records)),
            "source_assay_count": len(set(row["assay_chembl_id"] for row in records)),
            "scaffold_count": len(set(row["_scaffold"] for row in records)),
            "validity_comment_parent_count": sum(bool(row["data_validity_comment"]) for row in records),
            "excluded_duplicate_parent_ids": len(duplicate_ids),
        },
        "fold_results": folds_report,
        "aggregate_metrics": {
            "random_forest_macro_average_precision": macro_ap,
            "nearest_neighbor_macro_average_precision": macro_nn_ap,
            "mean_fold_prevalence": macro_prevalence,
            "folds_where_random_forest_ap_exceeds_prevalence": folds_above_prevalence,
            "random_forest_pooled_top5pct_recovery": pooled_recovery,
            "random_expected_pooled_top5pct_recovery": pooled_random_expected,
            "top5pct_enrichment_over_random": pooled_recovery / pooled_random_expected,
        },
        "qualified_for_targetwide_activity_support": qualifies,
        "qualification_rule": protocol["evaluation"]["qualification_rule"],
        "input_sha256": {
            "protocol": sha256(PROTOCOL_PATH),
            "fold_manifest": sha256(FOLD_MANIFEST_PATH),
            "activity_csv": sha256(ACTIVITY_PATH),
        },
        "interpretation": "Retrospective target-wide classification across heterogeneous human DPP-4 IC50 assays. No candidate potency, GLP-1 response, clinical exposure, safety, or therapeutic-equivalence claim is supported.",
        "predictions_csv": str((OUTPUT / "predictions.csv").relative_to(ROOT)),
    }
    with (OUTPUT / "predictions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(prediction_rows[0]))
        writer.writeheader()
        writer.writerows(prediction_rows)
    (OUTPUT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    write_markdown(report)
    print(json.dumps({
        "qualified": qualifies,
        "aggregate_metrics": report["aggregate_metrics"],
        "report": str((OUTPUT / "report.md").relative_to(ROOT)),
    }, indent=2))


def write_markdown(report: dict[str, object]) -> None:
    lines = [
        "# Target-wide DPP-4 activity benchmark",
        "",
        "**Protocol:** `methods/dpp4_ranking.json`; fold membership is frozen in `methods/dpp4_folds.json`.",
        "**Scope:** Five publication-held-out folds across heterogeneous human DPP-4 IC50 assays. All training compounds are scaffold-disjoint and have Morgan Tanimoto <=0.50 to the corresponding held-out fold.",
        "",
        "## Fold results",
        "",
        "| Fold | Test N | Active | Method | AP | ROC-AUC | Top-5% recovered | Top-5% size | Random expected |",
        "| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for fold_id, fold in report["fold_results"].items():
        for method, title in (("random_forest", "Morgan random forest"), ("nearest_neighbor", "Nearest-neighbor")):
            result = fold[method]
            lines.append(
                f"| {fold_id} | {result['parent_count']} | {result['positive_count']} | {title} | "
                f"{result['average_precision']:.3f} | {result['roc_auc']:.3f} | "
                f"{result['recovery_at_shortlist']}/{result['shortlist_size']} | {result['shortlist_size']} | "
                f"{result['random_expected_recovery_at_shortlist']:.2f} |"
            )
    aggregate = report["aggregate_metrics"]
    lines.extend([
        "",
        f"Macro average precision: random forest {aggregate['random_forest_macro_average_precision']:.3f}; "
        f"nearest neighbor {aggregate['nearest_neighbor_macro_average_precision']:.3f}; "
        f"mean prevalence {aggregate['mean_fold_prevalence']:.3f}.",
        f"Random forest AP exceeds fold prevalence in {aggregate['folds_where_random_forest_ap_exceeds_prevalence']}/5 folds. "
        f"Pooled top-5% recovery: {aggregate['random_forest_pooled_top5pct_recovery']}/" 
        f"{sum(fold['random_forest']['shortlist_size'] for fold in report['fold_results'].values())}; "
        f"random expected {aggregate['random_expected_pooled_top5pct_recovery']:.2f} "
        f"({aggregate['top5pct_enrichment_over_random']:.2f}x).",
        "",
        f"**Qualification:** {'Pass' if report['qualified_for_targetwide_activity_support'] else 'Fail'} under the frozen rule: {report['qualification_rule']}",
        "",
        "This benchmark does not predict clinical response, exposure, safety, or therapeutic equivalence. It is not a prospective model validation.",
        "",
        "All held-out compound predictions and data-validity comments are in `predictions.csv`; fold metrics and source checksums are in `report.json`.",
    ])
    (OUTPUT / "report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
