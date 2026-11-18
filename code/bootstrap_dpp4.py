#!/usr/bin/env python3
"""Estimate publication-clustered uncertainty for the frozen target-wide predictions."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import sklearn
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
PREDICTIONS = ROOT / "data/dpp4/held_out_predictions.csv"
OUTPUT = ROOT / "reproduced/dpp4/document_bootstrap.json"
REPORT = ROOT / "reproduced/dpp4/document_bootstrap.md"
REPLICATES = 2000
SEED = 20261003


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile_interval(values: np.ndarray) -> list[float]:
    return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]


def main() -> None:
    with PREDICTIONS.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    by_fold_doc: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_fold_doc[row["fold"]][row["document_chembl_id"]].append(row)

    rng = np.random.default_rng(SEED)
    rf_macro, nn_macro, prevalence_macro, ap_difference, top5_enrichment = [], [], [], [], []
    for _ in range(REPLICATES):
        rf_fold_aps, nn_fold_aps, prevalences = [], [], []
        recovered_total = 0
        random_expected_total = 0.0
        for fold in sorted(by_fold_doc, key=int):
            documents = sorted(by_fold_doc[fold])
            sampled_docs = rng.choice(documents, size=len(documents), replace=True)
            sampled_rows = [row for document in sampled_docs for row in by_fold_doc[fold][document]]
            labels = np.asarray([int(row["observed_active_at_100nM"]) for row in sampled_rows], dtype=np.uint8)
            rf_scores = np.asarray([float(row["random_forest_active_probability"]) for row in sampled_rows], dtype=float)
            nn_scores = np.asarray([float(row["nearest_neighbor_activity_score"]) for row in sampled_rows], dtype=float)
            rf_fold_aps.append(float(average_precision_score(labels, rf_scores)))
            nn_fold_aps.append(float(average_precision_score(labels, nn_scores)))
            prevalence = float(labels.mean())
            prevalences.append(prevalence)
            shortlist_size = max(1, int(np.ceil(0.05 * len(sampled_rows))))
            ranking = sorted(
                range(len(sampled_rows)),
                key=lambda index: (-rf_scores[index], sampled_rows[index]["parent_molecule_chembl_id"]),
            )
            recovered_total += int(labels[ranking[:shortlist_size]].sum())
            random_expected_total += shortlist_size * prevalence
        rf_macro_value = float(np.mean(rf_fold_aps))
        nn_macro_value = float(np.mean(nn_fold_aps))
        prevalence_value = float(np.mean(prevalences))
        rf_macro.append(rf_macro_value)
        nn_macro.append(nn_macro_value)
        prevalence_macro.append(prevalence_value)
        ap_difference.append(rf_macro_value - nn_macro_value)
        top5_enrichment.append(recovered_total / random_expected_total)

    arrays = {
        "random_forest_macro_average_precision": np.asarray(rf_macro),
        "nearest_neighbor_macro_average_precision": np.asarray(nn_macro),
        "mean_fold_prevalence": np.asarray(prevalence_macro),
        "macro_ap_difference_random_forest_minus_nearest_neighbor": np.asarray(ap_difference),
        "pooled_top5pct_recovery_enrichment": np.asarray(top5_enrichment),
    }
    summary = {
        key: {
            "median": float(np.median(values)),
            "percentile_95_interval": quantile_interval(values),
        }
        for key, values in arrays.items()
    }
    result = {
        "status": "exploratory uncertainty analysis after the frozen primary benchmark; does not change qualification",
        "resampling_unit": "source document within each of the five already-frozen held-out folds",
        "replicates": REPLICATES,
        "seed": SEED,
        "conditioned_on": "The already-fitted fold models and their predictions; training data are not resampled or refit.",
        "summary": summary,
        "bootstrap_fraction_random_forest_ap_exceeds_nearest_neighbor": float(np.mean(arrays["macro_ap_difference_random_forest_minus_nearest_neighbor"] > 0)),
        "bootstrap_fraction_random_forest_ap_exceeds_mean_prevalence": float(np.mean(arrays["random_forest_macro_average_precision"] > arrays["mean_fold_prevalence"])),
        "input_sha256": {"predictions_csv": sha256(PREDICTIONS)},
        "software": {"numpy": np.__version__, "scikit_learn": sklearn.__version__},
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    write_report(result)
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "summary": summary}, indent=2))


def write_report(result: dict[str, object]) -> None:
    lines = [
        "# Publication-grouped uncertainty for the target-wide DPP-4 benchmark",
        "",
        "This exploratory analysis resamples source documents within each frozen held-out fold. It is conditional on the existing fitted models and does not refit them, alter the predeclared qualification decision, or replace prospective validation.",
        "",
        f"Bootstrap replicates: {result['replicates']}; seed: {result['seed']}; grouped resampling unit: source document.",
        "",
        "| Quantity | Median | 95% percentile interval |",
        "| --- | ---: | ---: |",
    ]
    for key, value in result["summary"].items():
        label = key.replace("_", " ")
        low, high = value["percentile_95_interval"]
        lines.append(f"| {label} | {value['median']:.3f} | [{low:.3f}, {high:.3f}] |")
    lines.extend([
        "",
        f"Fraction of resamples where RF macro AP exceeds nearest-neighbor AP: {result['bootstrap_fraction_random_forest_ap_exceeds_nearest_neighbor']:.3f}.",
        f"Fraction of resamples where RF macro AP exceeds mean fold prevalence: {result['bootstrap_fraction_random_forest_ap_exceeds_mean_prevalence']:.3f}.",
        "",
        "Interpret the intervals as variation over held-out publication composition conditional on the five trained fold models. They do not include model-training uncertainty, publication-selection uncertainty, assay measurement error, or future prospective error.",
    ])
    REPORT.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
