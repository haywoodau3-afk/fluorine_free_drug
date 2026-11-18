"""Replay a recorded docking case from its archived PDBQT inputs and Vina log."""
from argparse import ArgumentParser
from pathlib import Path
import csv
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("design_id", help="Design identifier from data/docking/executed_runs.csv")
    parser.add_argument("--vina", default="vina", help="AutoDock Vina 1.2.7 executable")
    args = parser.parse_args()
    with (ROOT / "data/docking/executed_runs.csv").open(newline="") as stream:
        row = next((r for r in csv.DictReader(stream) if r["design_id"] == args.design_id), None)
    if row is None:
        parser.error(f"unknown design identifier: {args.design_id}")

    params = json.loads((ROOT / "methods/docking.json").read_text())["vina"]
    log = (ROOT / row["log_file"]).read_text()
    receptor = f"inputs/docking/receptors/{row['drug_id']}.pdbqt"
    ligand = f"inputs/docking/ligands/{row['drug_id']}/{args.design_id}.pdbqt"
    out = ROOT / "reproduced/docking" / f"{args.design_id}.pdbqt"
    out.parent.mkdir(parents=True, exist_ok=True)
    command = [
        args.vina, "--receptor", str(ROOT / receptor), "--ligand", str(ROOT / ligand),
        "--seed", row["seed"], "--exhaustiveness", row["exhaustiveness"],
        "--num_modes", str(params["num_modes"]), "--energy_range", str(params["energy_range_kcal_per_mol"]),
        "--cpu", str(params["cpu_threads"]), "--out", str(out),
    ]
    for grid in ("center", "size"):
        values = re.search(rf"^Grid {grid}\s*:\s*X ([\d.\-]+) Y ([\d.\-]+) Z ([\d.\-]+)", log, re.MULTILINE)
        if values is None:
            raise ValueError(f"Grid {grid} missing from the log for {args.design_id}")
        for axis, value in zip("xyz", values.groups()):
            command.extend((f"--{grid}_{axis}", value))
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
