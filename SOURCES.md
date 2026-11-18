# Source records and redistribution

## ChEMBL data

ChEMBL data are from [EMBL-EBI ChEMBL](https://www.ebi.ac.uk/chembl/). The DPP-4 REST extract used for the ranking was found to be content-equivalent to a status-bracketed ChEMBL 37 pull. The original REST payload supplied no release field. Other public retrievals retain individual ChEMBL record and source identifiers. The [ChEMBL 37 release files](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/) and the [ChEMBL REST API](https://www.ebi.ac.uk/chembl/api/data/) provide source access.

ChEMBL database extracts and their adaptations are distributed under [CC BY-SA 3.0 Unported](https://creativecommons.org/licenses/by-sa/3.0/). This applies to ChEMBL-derived tables, including ChEMBL identifiers and structures in study tables. Cite Zdrazil et al., [doi:10.1093/nar/gkad1004](https://doi.org/10.1093/nar/gkad1004), and Mendez et al., [doi:10.1093/nar/gky1075](https://doi.org/10.1093/nar/gky1075). The distributed ranking input selects the columns the script uses and retains the source row order. The activity extract includes unique activity IDs and source assay and document identifiers. It omits abstracts, assay descriptions, ontology labels and computed commercial property objects. The values have not been reattributed to the underlying papers where primary row mapping remains unresolved.

## PDB coordinates

[PDB archive data](https://www.rcsb.org/pages/usage-policy) are available under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). The original mmCIF inputs are distributed in `inputs/pdb/` as downloaded, alongside prepared structures and docking outputs. Cite the primary PDB records [1X70](https://www.rcsb.org/structure/1X70), [3LN1](https://www.rcsb.org/structure/3LN1), [1IKW](https://www.rcsb.org/structure/1IKW), [6HLO](https://www.rcsb.org/structure/6HLO) and [5HI2](https://www.rcsb.org/structure/5HI2). Prepared PDB/PDBQT files are computational derivatives. The docking log export replaces machine-specific absolute input paths with filenames while preserving numeric settings and outputs.

## Literature, assays and property inputs

The route and DPP-4 evidence tables identify publications or patents by DOI, patent number, compound number or ChEMBL document ID. They provide study-derived comparisons and do not distribute underlying article or patent PDFs. The fixed-dose labels and identifiers are from Carpio et al., [doi:10.3390/ijms26125579](https://doi.org/10.3390/ijms26125579), Table S1. The published DPP-4 and preparative comparisons are cross-referenced by the citation numbers in the article.

AstraZeneca lipophilicity source files are available through [TDC](https://tdcommons.ai/single_pred_tasks/adme/#lipophilicity) and Harvard Dataverse. Retrieve the [task table](https://dataverse.harvard.edu/api/access/datafile/4259595) and the [fixed ADMET group archive](https://dataverse.harvard.edu/api/access/datafile/4426004). Extract `lipophilicity_astrazeneca/train_val.csv` and `lipophilicity_astrazeneca/test.csv`. The study uses PyTDC 1.1.15 internal scaffold splitting and the endpoint-specific settings in `data/lipophilicity/validation_protocol.json`. Derived metrics and predictions are distributed with this package.

## Model and stock assets

The pretrained generator identifiers, commits and source locators are in `methods/generator_models.json`. The repository distributes generated proposals and analysis records. Retrieve the pretrained weights from their cited providers.

AiZynthFinder 4.4.1 uses a USPTO expansion model, USPTO reaction filter, USPTO templates and a ZINC stock. The [exact asset checksums](methods/retrosynthesis_assets.json) distinguish the release inputs from similarly named files. Obtain `uspto_model.onnx` and `uspto_filter_model.onnx` from [Zenodo 7797465](https://doi.org/10.5281/zenodo.7797465), the matching template table from [Zenodo 7341155](https://doi.org/10.5281/zenodo.7341155), and the frozen stock file `zinc_stock_17_04_20.hdf5` from [Figshare 12334577 v1](https://doi.org/10.6084/m9.figshare.12334577.v1). Match hashes before use and rename the stock and templates to the filenames in the supplied configuration. The route-status tables are derived summaries. Model, stock and detailed route-tree files remain with their respective providers or study records. Availability of exact route trees is subject to the redistribution terms for stock-derived contents.

Original code in `code/` is released under the [MIT License](LICENSE), copyright 2026 Long Hay Au. Original study data, figures and narrative documentation are released under [CC BY 4.0](DATA_LICENSE.md). These grants apply only to the author's original contributions. ChEMBL-derived extracts and adaptations retain CC BY-SA 3.0; PDB archive coordinates retain CC0 1.0; all other third-party terms above continue to apply to their respective records and assets.

The historical pair-input exports under `inputs/source-lineage/` are ChEMBL-derived study tables covered by the ChEMBL attribution and CC BY-SA 3.0 terms above. They retain the original source row order and omit free-text assay descriptions. The source-path map distinguishes original input hashes from exported-file hashes. Synthetic precedent source notes summarise existing study reviews and cite the underlying papers and patents; no article or patent PDF is distributed.

The fixed-dose input CSV transcribes only study identifiers, structures and the experimental labels used in the analysis from Carpio et al., Table S1. Source compounds and experimental outcomes remain attributed to that publication.
