# Synthetic precedent source notes

These notes accompany manuscript v6.4 and preserve source access for historical review paths retained in the data tables. They record the existing bounded reviews; they do not add a new literature search or change an evidence tier. Citation numbers below refer to the manuscript.

## Six-case synthesis evidence

The six cases and exact structures are in [six_case_review.csv](six_case_review.csv) and [six_case_structures.csv](six_case_structures.csv). The three strict-atlas case identities are also listed in [illustrated_designs.csv](../atlas/illustrated_designs.csv).

- Celecoxib-derived acid, `STATE-70d6e1863660ad083cad`: Raji et al., [DOI 10.1016/j.bmc.2016.12.032](https://doi.org/10.1016/j.bmc.2016.12.032), experimental section 4.1.4, reports acid 5 and 96% ester hydrolysis yield. Rogez-Florent et al., [DOI 10.1016/j.bmc.2012.10.029](https://doi.org/10.1016/j.bmc.2012.10.029), reports the exact ester 2c and 78% cyclocondensation yield. The printed ester charge in the acid paper has the mass/amount discrepancy described in SI S5; the published yield is retained as reported. These are manuscript references 26 and 35. The parent process, reference 36, [DOI 10.1039/D0RE00346H](https://doi.org/10.1039/D0RE00346H), supports a related reaction class, not a matched whole-route comparison.
- Sitagliptin C3-methyl design, `DESIGN-b7150737828a9c774fe4` / `STATE-361ff71ff48d8d51b659`: Kshirsagar et al., [DOI 10.14233/ajchem.2024.32667](https://doi.org/10.14233/ajchem.2024.32667), reports the 3-methyl tetrahydrotriazolopyrazine core (reference 27). [EP1412357B1](https://patents.google.com/patent/EP1412357B1/en) supports analogous protected beta-amino-acid coupling on a different, fluorinated aryl side chain (reference 28). Access to the exact nonfluorinated stereodefined side-chain partner, the exact coupling and complete target remain unverified. This is partial precedent, not an established synthesis of the design.
- Sorafenib pyridyl ether, `STATE-780a70aa4786ecb01a2c`: [US7528255B2](https://patents.google.com/patent/US7528255B2/en) supports related pyridinecarboxamide formation and ether assembly (reference 29). Exact activated-pyridine identity, monoarylation and retention of the free phenol remain candidate-specific questions.
- Efavirenz carboxamide, celecoxib thiazole and aprepitant morpholine cases retain the primary paper and patent identifiers in their respective `six_case_review.csv` rows. Those sources support scaffold or reaction classes to the stated extent, not synthesis of the exact targets.

## Frozen-shortlist source notes

The historical frozen-shortlist review was dated 11 October 2026. Its initial run was superseded by the corrected shortlist; candidate identifiers, rather than obsolete SYN labels, identify the relevant structures. The twelve-state selection later used a separate prespecified selection rule. The current selected structures, search dates, source URLs, candidate-specific risks and review status remain in [twelve_state_selection.csv](twelve_state_selection.csv); the 24 perspective-level conclusions are in [route_perspectives.csv](route_perspectives.csv).

The methyl sitagliptin evidence is identified above. For other frozen-shortlist structures, parent or class precedents support only a route concept: the exact edited core, side-chain partner, precursor access, compatibility or stereochemical operation remains unresolved where the tables say so. A source URL or shared starting-material listing does not establish access to an exact candidate or a complete synthesis. The selected aprepitant bromo/methoxymethyl structure must be read from its canonical SMILES; an older narrative called the substituent methoxy. The packaged structure and existing candidate-specific evidence tier are retained.

## Current-shortlist supplement

The 11 October 2026 supplement covered three replacement shortlist candidates:

- Sitagliptin C3-ethynyl, `DESIGN-e9c631257052d5c5d53f`: the parent process [DOI 10.1021/op0500786](https://doi.org/10.1021/op0500786) supports the assembly platform. It does not establish the exact C3-ethynyl core or the nonfluorinated chiral side chain.
- Celecoxib C3-(2-hydroxyethyl), `DESIGN-8a0023e1a976d0bba6a6`: the parent process [DOI 10.1039/D0RE00346H](https://doi.org/10.1039/D0RE00346H) supports cyclocondensation for the CF3-bearing parent. The exact hydroxyethyl precursor, regioselective use and compatibility remain unresolved.
- Celecoxib C3-dimethylamino, `DESIGN-2661ac09ff7d1e1d388c`: the same parent process is a platform precedent. No adequate exact-substrate or defining-transformation precedent was established by that review.

## Selected-state source crosswalk

This crosswalk preserves the source URLs already recorded for the twelve-state selection. Evidence grades and risks remain in the accompanying tables. Historical paths embedded in original source fields are resolved through [the original-path map](../../provenance/original_path_map.csv).

| Selection | State | Existing source URLs |
| --- | --- | --- |
| FRS-S-DPP4-01 | STATE-014d0aba511e80056de0 | http://pubs.acs.org/oprdfk/article/9/5/634/3689298/First-Generation-Process-for-the-Preparation-of; https://patents.google.com/patent/CN105017260B/en |
| FRS-S-DPP4-02 | STATE-78c6c6b4c35c41c40598 | https://patents.google.com/patent/CN105017260B/en |
| FRS-S-DPP4-03 | STATE-2eb986ac70cfa7454f97 | http://pubs.acs.org/oprdfk/article/9/5/634/3689298/First-Generation-Process-for-the-Preparation-of |
| FRS-S-DPP4-04 | STATE-bd72b76a7c3a9e72ba78 | https://patents.google.com/patent/CN105017260B/en; https://en.cnki.com.cn/Article_en/CJFDTOTAL-ZHOU201011004.htm |
| FRS-S-COX2-01 | STATE-0439b72d38063349790a | https://pubmed.ncbi.nlm.nih.gov/9135032/; https://doi.org/10.1333/s00897132496a |
| FRS-S-COX2-02 | STATE-38ca39c00e7756b40f5e | https://pubmed.ncbi.nlm.nih.gov/9135032/; https://www.mdpi.com/2813-446X/2/4/14 |
| FRS-S-COX2-03 | STATE-51bd7c3f3fe08227a4b3 | https://pubmed.ncbi.nlm.nih.gov/9135032/ |
| FRS-S-COX2-04 | STATE-84a44bf5306445e38163 | https://pubmed.ncbi.nlm.nih.gov/9135032/; https://www.mdpi.com/2813-446X/2/4/14 |
| FRS-S-NK1R-01 | STATE-0354e630cdb88c1544f9 | https://pubs.acs.org/doi/abs/10.1021/ja027458g; https://doi.org/10.1016/j.tetlet.2004.09.185 |
| FRS-S-NK1R-02 | STATE-0b9e7128bda9e246e7eb | https://pubs.acs.org/doi/abs/10.1021/ja027458g; https://doi.org/10.1016/j.tetlet.2004.09.185 |
| FRS-S-NK1R-03 | STATE-04ab44ba365d6e8495d4 | https://pubs.acs.org/doi/abs/10.1021/ja027458g |
| FRS-S-NK1R-04 | STATE-0a9de8acd430dad44a94 | https://pubs.acs.org/doi/abs/10.1021/ja027458g; https://doi.org/10.1016/j.tetlet.2004.09.185 |
