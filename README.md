# Endometriosis single-cell RNA-seq analysis

This project reanalyzes human endometrial single-cell RNA-seq data to identify endometriosis-associated, cell-type-specific transcriptional programs, with a focus on the early decidualized stromal state `dStromal_early`. It combines Marečková [E-MTAB-14039](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039) with the independent Huang cohort [GSE214411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214411), building on Marečková et al., [*An integrated single-cell reference atlas of the human endometrium*](https://www.nature.com/articles/s41588-024-01873-w) (*Nature Genetics*, 2024). The workflow is in [`code/notebooks/endometriosis_scRNAseq_analysis.ipynb`](code/notebooks/endometriosis_scRNAseq_analysis.ipynb).

## Research question and hypothesis

Marečková et al. identified `dStromal_early`, `dStromal_mid`, uM1 and uM2 among the states most enriched for expression of genes near endometriosis risk variants. They described `dStromal_early` as an early-secretory decidualized stromal state.

> **Working hypothesis:** because `dStromal_early` is enriched for expression of genes near endometriosis risk variants, it should show a donor-level disease-associated transcriptional program. Decidualization, stromal-remodeling and inflammatory pathways should agree in direction across Marečková and Huang even if individual genes do not replicate. Comparison with other cell types should reveal whether this program is `dStromal_early`-localized or more broadly distributed.

The fGWAS enrichment motivates the cellular focus but does not establish differential expression, pathway direction, cell-type specificity or therapeutic relevance. Marečková is therefore the primary cohort; Huang is modeled separately as independent sensitivity evidence.

## Input data

| Cohort | Local input | Input structure | Metadata |
|---|---|---|---|
| Marečková, E-MTAB-14039 | `data/work/` | Filtered 10x matrices packaged by sequencing library | [`metadata/sample_metadata.csv`](metadata/sample_metadata.csv) |
| Huang, GSE214411 | `huang_data/GSE214411_RAW.tar` | Thirteen prefixed 10x matrices: six endometriosis and seven control donors | [`metadata/huang_sample_metadata.csv`](metadata/huang_sample_metadata.csv) |

Each library is read as a separate raw-count `AnnData`, assigned namespaced identifiers and curated `.obs` metadata, then combined over shared Ensembl genes. Integer counts remain in `layers['counts']` for count-based modeling.

### Why Huang was the only external atlas dataset added

This is a focused disease comparison, not a reconstruction of all seven HECA datasets. Huang was selected because it provides:

- compatible raw counts with donor identities for pseudobulk analysis;
- six endometriosis and seven control donors after metadata harmonization;
- selectable natural-cycle, untreated, eutopic samples with compatible menstrual-stage information; and
- `dStromal_early` cells recovered through the common QC and CellTypist workflow.

More cells cannot replace independent case and control donors in donor-level DE. Additional sources lacking a valid within-study contrast would add technical and compositional heterogeneity without identifying the disease effect. Huang is therefore a deliberately limited sensitivity cohort, not evidence of generalization to every HECA dataset.

## Notebook workflow

### 1. Data loading and metadata

The notebook streams the 10x matrices, validates identifiers, constructs one `AnnData` per library, attaches donor and clinical metadata, and combines the cohorts.

### 2. Quality control and preprocessing

The primary cell filters follow the paper at a high level:

- at least 1,000 detected genes per cell;
- no more than 20% mitochondrial counts;
- Scrublet scoring performed independently for each sequencing library;
- predicted doublets retained in the general atlas but excluded from the downstream `dStromal_early` DE cohort;
- normalization to 10,000 counts per cell followed by `log1p`;
- cell-cycle scoring without regressing cell-cycle signal from expression;
- exclusion of the CDK1-associated cycling-gene module from integration features while retaining those genes for annotation and DE.

Across both cohorts, 204,385 of 249,684 cells passed the gene-count and mitochondrial filters (81.9%).

| Disease group | Cells before QC | Cells after QC | Retained |
|---|---:|---:|---:|
| Control | 114,328 | 88,364 | 77.3% |
| Endometriosis | 135,356 | 116,021 | 85.7% |
| **Total** | **249,684** | **204,385** | **81.9%** |

#### QC before filtering

The red dashed lines mark the minimum detected-gene and maximum mitochondrial-percentage thresholds.

![Detected genes, total counts, and mitochondrial content before QC](results/qc/figures/01_pre_filter_qc_distributions.png)

#### QC after filtering

![Detected genes, total counts, and mitochondrial content after QC](results/qc/figures/03_post_filter_qc_distributions.png)

### 3. scVI integration, clustering and cell typing

The notebook selects 2,000 library-aware highly variable genes, excludes the inferred cycling module from integration features, and trains scVI on raw counts with `library_id` as the batch key. Its 30-dimensional latent representation supports the neighbor graph, UMAP and Leiden clustering; scVI-decoded expression is not used for cell typing or DE.

Leiden clustering is performed at resolution 1.0:

![Leiden clusters on the library-corrected scVI UMAP](results/cell_typing/figures/01_scvi_umap_leiden.png)

Cell types are predicted from log-CP10K expression using CellTypist's `Human_Endometrium_Atlas.pkl` model and over-cluster majority voting. scVI coordinates are used only for visualization.

![CellTypist Human Endometrium Atlas labels on the scVI UMAP](results/cell_typing/figures/02_celltypist_labels_on_scvi_umap.png)

### 4. Canonical donor-level differential expression and pathway analysis

The inferential analysis retains natural-cycle, untreated, eutopic samples; excludes Scrublet-predicted doublets; and requires at least 20 cells of a given type per donor. Raw counts are summed into one pseudobulk per donor, and each model requires at least three Normal and three Endometriosis donors per study.

Eight compartments meet these criteria in both cohorts: `dStromal_early`, `eStromal`, `preGlandular`, `SOX9_functionalis_II`, `Venous`, `ePV_2`, `Immune_Myeloid` and `Immune_Lymphoid`. All 16 cell-type-by-study models fitted successfully. Each cohort is modeled separately with PyDESeq2:

```text
~ cycle_group + disease_group
```

The cycle term is retained when estimable. The contrast is Endometriosis minus Normal, with gene-level FDR controlled at 0.10 within each cell type and dataset. Weighted preranked GSEA tests 50 Hallmark pathways using the full DESeq2 Wald-statistic ranking; P values are adjusted within each cell type and across all cell-type-by-pathway tests per dataset. Cohorts are not pooled, and cells are never treated as replicates.

![Cell-type-resolved Hallmark dysregulation](results/differential_expression/revised_analysis/cell_type_pathways/figures/05_celltype_hallmark_dysregulation.png)

![Cross-study Hallmark concordance by cell type](results/differential_expression/revised_analysis/cell_type_pathways/figures/06_celltype_hallmark_cross_study_concordance.png)

Across the eight cell types, 52 pathway-cell-type pairs have within-cell-type BH-FDR < 0.25 in both cohorts and matching enrichment direction. Nine occur in `dStromal_early`, all higher in Endometriosis:

- allograft rejection;
- androgen response;
- coagulation;
- epithelial-mesenchymal transition;
- IL2-STAT5 signaling;
- KRAS signaling up;
- mTORC1 signaling;
- myogenesis; and
- UV response down.

In `dStromal_early`, 31 Hallmarks pass the exploratory threshold in Marečková and 32 in Huang. Nineteen pass in both, but only these nine agree in direction; the other ten are not considered replicated mechanisms.

Four pathways—androgen response, KRAS signaling up, mTORC1 signaling and UV response down—are concordant only in `dStromal_early`. Coagulation and IL2-STAT5 recur in one other compartment; allograft rejection also recurs in `preGlandular` and `Immune_Lymphoid`; and epithelial-mesenchymal transition and myogenesis recur with the opposite direction in `SOX9_functionalis_II`.

#### Focal `dStromal_early` gene-level findings

The focal results are a prespecified view of the canonical models, not a second DE or GSEA run. Although `dStromal_early` is early-secretory-associated, restricting donors to `menstrual_cycle_stage_fine == 'Secretory Early'` leaves too few for inference in both cohorts.

| Dataset | Donors (Normal / Endometriosis) | Genes tested | FDR < 0.10 | Higher in endometriosis | Higher in Normal |
|---|---:|---:|---:|---:|---:|
| Marečková E-MTAB-14039 | 8 (3 / 5) | 13,803 | 12 | 11 | 1 |
| Huang GSE214411 | 13 (7 / 6) | 14,606 | 46 | 32 | 14 |

The 58 FDR-significant gene-study associations are study-specific: none reaches FDR < 0.10 in both datasets. Leading genes include `MT1G`, `SPP1`, `LINC01320`, `IFITM1` and `CD74` in Marečková, and `AL592183.1`, `ADGRB3`, `NEU1`, `FMN1` and `MYC` in Huang.

The figure shows the 40 smallest adjusted P values. Positive log2 fold changes indicate higher expression in Endometriosis; gold outlines mark Open Targets drug or clinical-candidate records, not demonstrated treatments.

![Study-specific dStromal early differential-expression effects](results/differential_expression/revised_analysis/figures/01_study_specific_significant_gene_effects.png)

Six significant genes overlap Open Targets drug or clinical-candidate records:

| Dataset | Gene | log2 fold change | FDR | Open Targets drugs or candidates | Relationship to the nine concordant `dStromal_early` pathways |
|---|---|---:|---:|---|---|
| Marečková | `SPP1` | 4.390 | 0.0012 | ASK-8007 | Leading-edge gene in its source cohort |
| Marečková | `CD74` | 2.606 | 0.0355 | Milatuzumab; repotrectinib | Leading-edge gene in its source cohort |
| Huang | `SMAD7` | 1.453 | 0.0039 | Mongersen sodium | Leading-edge gene in its source cohort |
| Huang | `TGFB1` | 0.667 | 0.0169 | Bintrafusp alfa; fresolimumab; luspatercept; LY-2382770; metelimumab | Leading-edge gene in its source cohort |
| Huang | `ADM` | 1.503 | 0.0480 | Enibarcimab | DE-only candidate relative to the nine concordant sets |
| Huang | `VWF` | -1.818 | 0.0946 | Caplacizumab; egaptivon pegol; von Willebrand factor products | Coagulation-set member, but not a source-cohort leading-edge gene |

Drug annotations do not affect statistical significance. `SPP1`, `CD74`, `TGFB1` and `SMAD7` connect gene- and pathway-level results only within their source cohorts; no candidate is a leading-edge driver in both studies. Database links may represent direct binding, RNA or protein modulation, replacement therapy, or fusion activity, so they motivate follow-up rather than treatment.

#### CD74 example

In Marečková, `CD74` was higher in Endometriosis (`log2FoldChange = 2.606`, approximately 6.1-fold; `FDR = 0.0355`). CD74 is the MHC class II invariant-chain chaperone and can act as a receptor for macrophage migration inhibitory factor, linking it to antigen presentation and inflammatory signaling ([NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/972)). Its elevation is compatible with altered stromal immune interaction, but RNA abundance does not establish mechanism or cell-surface protein abundance.

Two oncology-related records were returned for `CD74`, but they have different meanings:

- [Milatuzumab](https://www.cancer.gov/publications/dictionaries/cancer-drug/def/milatuzumab) directly binds CD74 but has been studied as an antineoplastic agent, not for endometriosis.
- [Repotrectinib](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2023/218213Orig1s000Lbl.pdf) targets ROS1/TRK in selected cancers. Its CD74 link reflects oncogenic `CD74–ROS1` fusions, not inhibition of ordinary CD74 in `dStromal_early` cells.

Thus, `CD74` is a statistically supported, pharmacologically annotated candidate in Marečková, but therapeutic relevance requires protein localization, evidence of causal CD74/MIF signaling, endometrial-model perturbation and independent donor validation.

Consolidated results are in [`results/differential_expression/revised_analysis/`](results/differential_expression/revised_analysis/), with cell-type DE, Hallmark GSEA and concordance tables under [`cell_type_pathways/`](results/differential_expression/revised_analysis/cell_type_pathways/).

## Hypothesis evaluation and principal finding

The working hypothesis is **supported at the pathway level, but not by individual-gene replication**:

- **Gene level:** Marečková identifies 12 FDR-significant `dStromal_early` genes and Huang identifies 46, mostly higher in Endometriosis. None reaches FDR < 0.10 in both cohorts, so there is no replicated gene signature.
- **Pathway level:** nine `dStromal_early` Hallmarks show matching positive direction and within-cell-type FDR < 0.25 in both cohorts. Four—mTORC1 signaling, KRAS signaling up, androgen response and UV response down—are concordant only in `dStromal_early` among the eight tested states.
- **Interpretation:** the convergent evidence supports a coordinated, cell-state-dependent remodeling, signaling and immune-interaction program. The fGWAS result motivates the cellular focus; the cross-study pathway analysis supplies the disease-associated functional evidence.
- **Experimental hypothesis:** test **partial mTORC1 attenuation as a mechanistic probe** in donor-derived stromal-cell decidualization models, asking whether it reduces the disease-associated program while preserving decidual function. This is not a treatment recommendation and requires protein or activity confirmation, perturbation, dose-response testing and independent-donor replication.
- **Limitations:** the focal Marečková model has three Normal donors; fine menstrual-stage matching is incomplete; cell labels are computationally transferred; RNA does not establish protein activity; and pathway FDR < 0.25 is exploratory. A secretory-only descriptive analysis preserves direction in all eight key-gene-by-study comparisons but has too few donors for separate inference.

## Reproducible environment

Create the Conda environment from the repository root:

```bash
conda env create --file environment/environment.yml
```

Register it as a Jupyter kernel:

```bash
conda run --name endometriosis-study \
  python -m ipykernel install --user \
  --name endometriosis-study \
  --display-name "Endometriosis study"
```

Open [`code/notebooks/endometriosis_scRNAseq_analysis.ipynb`](code/notebooks/endometriosis_scRNAseq_analysis.ipynb) in JupyterLab or an IDE and select the **Endometriosis study** kernel. The focal `dStromal_early` tables and figures are derived from the same canonical cell-type models; the notebook does not rerun focal DESeq2 or GSEA.

## Main outputs

| Output | Location |
|---|---|
| QC tables and figures | [`results/qc/`](results/qc/) |
| CellTypist annotations and figures | [`results/cell_typing/`](results/cell_typing/) |
| Consolidated focal DE, pathway, sensitivity and drug-target results | [`results/differential_expression/revised_analysis/`](results/differential_expression/revised_analysis/) |
| Cell-type-resolved DE, Hallmark GSEA and cross-study concordance | [`results/differential_expression/revised_analysis/cell_type_pathways/`](results/differential_expression/revised_analysis/cell_type_pathways/) |
| Consolidation equivalence report | [`results/differential_expression/revised_analysis/consolidation_equivalence_summary.csv`](results/differential_expression/revised_analysis/consolidation_equivalence_summary.csv) |
| Filtered and integrated `AnnData` objects | `data/processed/` |
| Cached CellTypist and scVI models | `data/models/` and `data/processed/scvi_model/` |
