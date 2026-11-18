# Computational design of fluorine-free drug analogues informs synthetic method development

Data and code accompanying *Computational design of fluorine-free drug analogues informs synthetic method development*. The manuscript and Supporting Information accompany this repository as separate submission documents.

## Article guide

| Article result | Files | Main identifiers |
| --- | --- | --- |
| Figures 1–2 and Table 1. Atlas, stereochemical audit and graph edits | [Atlas data](data/atlas), [enumeration settings](methods/atlas_enumeration.json) | `state_id`, `mask_id`, `connectivity_family_id` |
| Fixed-budget generator experiments | [Attempt and task records](data/generator-comparisons), [raw proposals](inputs/generator-proposals), [model records](methods/generator_models.json) | Task and proposal IDs |
| Figure 3 and Table 2. Literature evidence | [Six-case review](data/synthetic-precedent/six_case_review.csv), [12-state selection](data/synthetic-precedent/twelve_state_selection.csv), [24 route perspectives](data/synthetic-precedent/route_perspectives.csv), [compounds 1–4](data/atlas/illustrated_designs.csv) | `state_id`, case and selection IDs |
| Figure 4A. Sitagliptin reference pose | [Redocking files](inputs/reference-redocking), [RMSD results](data/docking/sitagliptin_redocking.json) | PDB ID, seed, pose number |
| Figure 4B–C. DPP-4 ranking and literature comparisons | [DPP-4 data](data/dpp4), [ChEMBL input](inputs/chembl), [splits](methods/dpp4_folds.json) | ChEMBL parent, document and assay IDs |
| Figures 4D and 6. Retrosynthesis | [926 design outcomes](data/retrosynthesis/design_outcomes.csv), [42 zero-fluorine designs](data/retrosynthesis/zero_fluorine_designs.csv), [figure design crosswalk](data/target-aware/illustrated_designs.csv) | `design_id` |
| Figure 5. Docking | [Manuscript v6.4 artwork](analysis/figure-5/Figure-5.png), [panel-to-source index](analysis/figure-5/panel_index.csv), [displayed atom distances](analysis/figure-5/displayed_distances.csv), [docking results](data/docking/results.csv), [first-pose sensitivity](analysis/first_pose_sensitivity.csv) | Panel letter, `design_id`, PDB ID |
| Lipophilicity comparison | [Per-state comparisons and per-seed validation](data/lipophilicity) | `state_id`, endpoint, seed |

Supporting Information S10 gives exact paths for each manuscript figure and table, along with a second map for inputs, settings and code.

For Figure 5, select a row in the [panel index](analysis/figure-5/panel_index.csv) to locate the individual panel, crystallographic mmCIF file, prepared receptor, Vina log and saved pose. The [distance table](analysis/figure-5/displayed_distances.csv) identifies the eight displayed atom pairs using the ligand PDBQT atom serial and the receptor author chain, residue number and atom name. The article composite is the exact Figure 5 embedded in manuscript v6.4. Individual panel files and [the alternative source-view composite](analysis/figure-5/Figure-5-source-view.png) provide another layout of the same coordinate views. The illustrations and dashed distances use model 1. The tabulated Jaccard values use all returned poses.

Historical source paths retained inside CSV records are resolved by [the source-path map](provenance/original_path_map.csv), with [mapping notes](provenance/README.md) and [synthetic precedent source notes](data/synthetic-precedent/source_notes.md). Original input hashes are preserved; selected-column exports have their own hashes.

## Reproduce the analysis

Use Python 3.12 with the dependencies in [requirements.txt](requirements.txt). Run the following from the repository root. Outputs are written to `reproduced/`.

```sh
python -m pip install -r requirements.txt
python code/reproduce_summaries.py
python code/fit_dpp4_ranking.py
python code/bootstrap_dpp4.py
python code/verify_reference_pose.py
```

The summary script validates cohort counts, site-fragment identities, map-free topology, generator budgets, docking settings, route outcomes, DPP-4 metrics and lipophilicity decisions, and recomputes first-pose contact sensitivity. The fitting script fits the forest and nearest-neighbour baseline on the published analysis input and fold assignments. The bootstrap operates on the distributed held-out predictions. [Analysis outputs](analysis) contain the checked numerical summaries and sensitivity results.

For target-aware docking replay, point `code/replay_docking.py` to a Vina 1.2.7 executable. It reads the grid, seed and exhaustiveness recorded for a design and uses the accompanying prepared PDBQT files.

```sh
python code/replay_docking.py DME-2c3bc93701a42d819162 --vina vina
```

The archived generation and route-search outputs can be inspected from the proposal and per-design tables. Regenerating the original libraries and fitting the lipophilicity model are distinct from reproducing the distributed summaries. The package supplies the preserved design and validation records and settings, plus retrieval pointers for external inputs. For a fresh planner run, use AiZynthFinder 4.4.1 and the matching files identified in [SOURCES.md](SOURCES.md) and [the asset manifest](methods/retrosynthesis_assets.json). Place verified assets in `external/retrosynthesis/`, generate a one-SMILES-per-line target file from `data/retrosynthesis/design_outcomes.csv` and use the settings in [methods/retrosynthesis.json](methods/retrosynthesis.json) with the `uspto` expansion policy and filter and `zinc` stock. Run and merge unresolved targets separately as described in Supporting Information S8.

## Data dictionary

- `state_id` identifies an atom-mapped atlas record. Use it for atlas joins. The same product or connectivity can occur in multiple state records.
- `connectivity_family_id` groups atom-mapped records by drug-specific heavy-atom connectivity. `candidate_id` is a design lineage key and can repeat.
- `mask_id` identifies editable and protected parent atoms. Atom index lists are zero-based. Mapped SMILES atom-map numbers encode retained parent indices plus one. New atoms have no positive parent map.
- `design_id` identifies a target-aware structure. Canonical isomeric SMILES provide a second identity check. `fluorine_atoms` counts fluorine independently of motif eligibility. `stock_terminating_route` records the bounded planner result.
- `contact_jaccard` uses the union of contacts over all saved poses per design. `first_pose_jaccard` uses the first pose only. Residue keys encode author chain, residue name and sequence number.
- Activity fields report IC₅₀ in nM unless a different endpoint is explicit. An empty measurement field indicates an unavailable value. The lipophilicity output records an interval-union decision for every strict atlas state.
- Short labels 1–4 and figure panel letters identify the illustrated molecules in the article. Their persistent data keys are listed in the two illustrated-designs tables.

Original code is licensed under [MIT](LICENSE); original study data, figures and narrative documentation under [CC BY 4.0](DATA_LICENSE.md). Third-party terms, including CC BY-SA 3.0 for ChEMBL-derived material, remain as described in [SOURCES.md](SOURCES.md). SHA256SUMS covers the files distributed with this package.
