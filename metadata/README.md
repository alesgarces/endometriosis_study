# Sample metadata

`sample_metadata.csv` is a curated, analysis-ready manifest with one row per
local donor–library archive. The primary key is `sample_id`, which matches the
archive basename and the keys of the notebook's `sample_adatas` dictionary.

## Sources

- ArrayExpress E-MTAB-14039:
  https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039
- Marečková et al. (2024):
  https://doi.org/10.1038/s41588-024-01873-w
- Supplementary Table 1: harmonized donor metadata.
- Supplementary Table 2: Cell Ranger library-level QC summaries.

The sources were reconciled against all 35 archives on 2026-08-14. The
ArrayExpress broad and fine cycle stages agreed with Supplementary Table 1 for
all 35 rows.

## Metadata levels

| Level | Key | Fields |
|---|---|---|
| Donor–library sample | `sample_id` | archive, donor, library, assay |
| Donor | `donor_id` | age, disease, cycle stage, hormone exposure, endometriosis stage, pathology, tissue location, dissociation |
| Library | `library_id` | Cell Ranger called cells, reads per cell, genes per cell |

Library-level Cell Ranger values are repeated when a multiplexed library was
split into donor-specific matrices. They describe the original library and are
not the number of cells in a local donor-specific matrix.

## Column dictionary

| Column | Level | Description |
|---|---|---|
| `sample_id` | sample | `<library_id>_<donor_id>`; unique primary key |
| `archive_filename` | sample | Local source archive |
| `library_id` | library | 10x sequencing-library identifier |
| `donor_id` | donor | Biological replicate identifier |
| `source_dataset` | donor | Dataset designation in Supplementary Table 1 |
| `organism` | study | Organism profiled |
| `tissue` | study | Tissue profiled |
| `tissue_location` | donor | `endometrium` for superficial biopsy or `whole_uterus` for full-thickness uterine wall |
| `assay` | sample | Sequencing modality; all local archives are `scRNA-seq` |
| `disease` | donor | ArrayExpress value: `normal` or `endometriosis` |
| `case_control` | donor | Analysis label derived from `disease`: `control` or `case` |
| `age_years` | donor | Age at sampling |
| `cycle_context` | donor | Derived label: `natural_cycle` or `exogenous_hormones` |
| `menstrual_cycle_stage` | donor | Broad stage from ArrayExpress: Menstrual, Proliferative, Secretory, or Hormones |
| `menstrual_cycle_stage_fine` | donor | More detailed stage where available |
| `hormonal_treatment` | donor | Exogenous treatment; source `NA` was converted to `none` because the paper defines it as no hormonal treatment for at least 3 months |
| `endometriosis_stage` | donor | rASRM stage; `0` denotes a control without endometriosis |
| `endometriosis_stage_assignment` | donor | Evidence used to assign disease stage |
| `endometrial_pathology_code` | donor | Source codes: A adenomyosis, C no pathology, E endometriosis, F fibroids, H hyperplasia, P polyp |
| `endometrial_pathology` | donor | Expanded semicolon-delimited pathology description |
| `tissue_dissociation_processing` | donor | Processing method recorded in Supplementary Table 1 |
| `cellranger_called_cells` | library | Cell Ranger cells for the unsplit library |
| `cellranger_reads_per_cell` | library | Cell Ranger reads per cell for the unsplit library |
| `cellranger_genes_per_cell` | library | Cell Ranger genes per cell for the unsplit library |
| `arrayexpress_accession` | study | Public data accession |
| `publication_doi` | study | Source publication DOI |

For statistical tests, `donor_id` is the biological replicate. Multiple
`sample_id` or `library_id` values from one donor must not be treated as
independent donors.
