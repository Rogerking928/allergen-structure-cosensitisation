# Predicted structural similarity and allergen co-sensitisation

Analysis code for a study asking how far predicted structural similarity between
allergens corresponds to co-sensitisation measured in patients, and how that
compares with the sequence criterion used in food safety assessment.

## Data

The analysis uses the **Allergen Chip Challenge** public release (4,271 patients,
11 French university hospitals, three microarray platforms), published under
Licence Ouverte 2.0 and available from data.gouv.fr
(https://doi.org/10.60597/j5fe-g420). The file is semicolon-separated with
comma decimal marks; `analysis/common.py` reads it. Place it in `data/`.

Allergen names are resolved against the WHO/IUIS Allergen Nomenclature database
(allergen.org), predicted structures are taken from the AlphaFold Protein
Structure Database, and Pfam / SCOP superfamily annotations from UniProt
cross-references. These are downloaded by the scripts below and cached under
`data/`.

No patient-level data are redistributed here. The pair-level estimates for all
18,848 component pairs are published with the article as Supplementary Data S1.

## Pipeline

Run from `analysis/`, in this order:

| Step | Script | Output |
| --- | --- | --- |
| Allergen name to UniProt | `map_uniprot.py` | `out/component_uniprot.csv` |
| Download AlphaFold models | `fetch_structures.py` | `data/afdb/` |
| Structural and sequence similarity | `similarity.py` | `out/pair_similarity.csv` |
| Co-sensitisation per pair | `cosens.py` | `out/cosens_pooled.csv` |
| Join and classify pairs | `structure_cosens.py` | `out/pairs_structure_cosens.csv` |
| Discrimination, thresholds | `discrim.py` | `out/discrimination.json` |
| TM-score bands, permutation nulls | `bands.py` | `out/tm_bands_ci.csv` |
| Sensitivity analyses | `robust.py`, `rho_ci.py`, `plddt.py` | `out/robust_summary.csv` |
| Joint-model coefficient difference | `coef_diff.py` | `out/coef_diff.json` |
| Split-sample analysis | `holdout.py` | `out/holdout.json` |
| Collect every reported number | `summary.py` | `out/summary.json` |
| Additional sensitivity analyses | `extra.py` | `out/extra.json` |
| Figures | `figs.py` | `out/Figure_*.png` |
| Tables, figure book, source data | `book.py` | deliverables |

`overlay.py` and `overlay_render.py` produce the superposed AlphaFold cartoons
(PyMOL required). Structural comparison uses Foldseek's TM-align mode.

Scripts in `code/` build the manuscript: `build_refs.py` resolves references
against PubMed and Crossref, `number_refs.py` numbers them in order of first
citation, `build_manuscript_doc.py` renders the Word file, and
`check_manuscript.py` verifies that every number in the text can be found in a
figure or table and that figures and tables are numbered in citation order.

Set `DELIVERABLE_DIR` to choose where the figures, tables and source-data
workbook are written.

## Statistical conventions

The unit of analysis is the allergen pair. Pairs are not independent, so P
values come from permuting allergen labels and confidence intervals from
resampling allergens rather than pairs, with the weight of a pair being the
product of the resampling counts of its two members.

## Requirements

Python 3.12 with pandas, numpy, scipy, scikit-learn, biopython, matplotlib,
networkx, python-docx and openpyxl; Foldseek for structural comparison; PyMOL
for the superposition panels.

## Citation and licence

Code released under the MIT licence. If you use it, please cite the article.
