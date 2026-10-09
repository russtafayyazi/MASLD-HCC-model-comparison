# MASLD-HCC model comparison

Analysis code for the manuscript comparing FAT-MASH and HOT-MASH mouse models of MASLD-associated hepatocellular carcinoma.

This repository contains the analysis code used to generate the main results and figures in the manuscript.

## Analysis workflow

The scripts are numbered in the approximate order of the analysis workflow.

| Script | Description |
| --- | --- |
| `01_phenotype_statistics.R` | Statistical analysis of histology, lesion burden, fibrosis, and serum ALT measurements. |
| `02_scRNAseq_integration_annotation.R` | scRNA-seq preprocessing, quality control, integration, clustering, and primary cell-type annotation. |
| `03_tumor_cell_identification.R` | Identification of the tumor-enriched hepatocyte cluster using an HCC-associated gene signature. |
| `04_tumor_programs_hoshida.R` | Hoshida HCC subtype scoring using pseudobulk singscore and single-cell UCell analyses. |
| `05_human_MASLD_mapping.R` | Projection of mouse hepatic lineages onto the human MASLD reference dataset GSE202379. |
| `06_human_MASLD_HCC_signature.R` | Derivation of a human MASLD-HCC tumor signature from GSE164760 and scoring in mouse tumor cells. |
| `07_create_lineage_objects.R` | Preparation of refined lineage-specific RNA objects for downstream analyses. |
| `08_model_vs_control_pseudobulk_DE_fGSEA.R` | Model-vs-chow pseudobulk differential expression and fGSEA across major hepatic lineages. |
| `09_model_vs_model_pseudobulk_DE_fGSEA.R` | HOT-MASH vs FAT-MASH pseudobulk differential expression and fGSEA, including malignant hepatocytes. |
| `10_fGSEA_theme_collapsing_and_plotting.R` | Grouping of significant fGSEA pathways into broader biological themes and plotting of theme-level summaries. |
| `11_candidate_gene_analysis.R` | Candidate-gene expression analysis across hepatocytes, endothelial cells, macrophages, and HSCs. |
| `12_cellchat_analysis_and_plotting.R` | CellChat analysis and plotting of signaling roles and model-specific communication networks. |

## Data

Public human datasets used in the analysis:

- `GSE202379` — human MASLD reference used for cross-species projection.
- `GSE164760` — human MASLD-HCC dataset used to derive the tumor-associated gene signature.

The repository also includes:

- `data/tumorCells.VSTnorm.humanOrths.txt` — VST-normalized pseudobulk tumor-cell expression mapped to human ortholog symbols.

Large Seurat objects, intermediate RDS files, and analysis outputs are not included in this repository.

## Notes

Some scripts retain file paths and intermediate object names from the original analysis environment and may require path changes before use.

Cross-species mapping uses `biomaRt` to retrieve one-to-one mouse-human orthologs from Ensembl, so returned annotations may vary over time.

For fGSEA, the original analysis used `.rnk` files prepared from the DESeq2 Wald statistic. The deposited scripts reconstruct the same ranked vectors directly from the DESeq2 results.

Cytoscape EnrichmentMap was used for the pathway-network panels in Figures 5B-C and 6B.
