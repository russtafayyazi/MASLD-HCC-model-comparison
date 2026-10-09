# Create refined lineage objects for downstream pseudobulk DE and CellChat
#
# Input:
#   master_integrated_primaryAnnotated.rds
#
# Outputs:
#   Model-vs-chow RNA-only objects
#   Model-vs-model RNA-only object
#   Lineage-specific RNA-only objects used in downstream analyses
#
# Internal condition nomenclature retained from the original analysis:
#   CCl4_model = FAT-MASH
#   DEN_model  = HOT-MASH
#   RT_chow    = chow control

library(Seurat)
library(dplyr)


# =============================================================================
# Part 1: load annotated atlas and retain RNA assay
# =============================================================================

master_integrated <- readRDS(
  "master_integrated_primaryAnnotated.rds"
)

DefaultAssay(master_integrated) <- "RNA"

if ("SCT" %in% Assays(master_integrated)) {
  master_integrated[["SCT"]] <- NULL
}

if ("integrated" %in% Assays(master_integrated)) {
  master_integrated[["integrated"]] <- NULL
}

gc()

Idents(master_integrated) <- "seurat_clusters"


# =============================================================================
# Part 2: refine biological annotations for downstream analyses
# =============================================================================

master_integrated$primary_label_refined <- as.character(
  master_integrated$primary_label
)

# Cluster 13 model cells are retained as malignant hepatocytes
master_integrated$primary_label_refined[
  master_integrated$seurat_clusters == "13" &
    master_integrated$condition %in% c(
      "CCl4_model",
      "DEN_model"
    )
] <- "Malignant hepatocytes"

# Chow cells in cluster 13 remain hepatocytes
master_integrated$primary_label_refined[
  master_integrated$seurat_clusters == "13" &
    master_integrated$condition == "RT_chow"
] <- "Hepatocytes"

# Cluster 14 remains biologically annotated as cholangiocytes
master_integrated$primary_label_refined[
  master_integrated$seurat_clusters == "14"
] <- "Cholangiocytes"

master_integrated$primary_label_refined <- factor(
  master_integrated$primary_label_refined,
  levels = c(
    "Hepatocytes",
    "Malignant hepatocytes",
    "Endothelial cells",
    "Macrophages",
    "HSCs",
    "Cholangiocytes"
  )
)


# -----------------------------------------------------------------------------
# Analysis label used for DE and CellChat
# -----------------------------------------------------------------------------

master_integrated$analysis_label <- as.character(
  master_integrated$primary_label_refined
)

# Cholangiocytes were excluded from lineage-level downstream comparisons
# because of low per-mouse abundance.
master_integrated$analysis_label[
  master_integrated$primary_label_refined == "Cholangiocytes"
] <- "Other"

master_integrated$analysis_label <- factor(
  master_integrated$analysis_label,
  levels = c(
    "Hepatocytes",
    "Malignant hepatocytes",
    "Endothelial cells",
    "Macrophages",
    "HSCs",
    "Other"
  )
)


# Basic annotation checks
table(
  master_integrated$primary_label_refined,
  master_integrated$condition
)

table(
  master_integrated$analysis_label,
  master_integrated$condition
)


# =============================================================================
# Part 3: create model-vs-chow objects
# =============================================================================

CCl4_RT_obj <- subset(
  master_integrated,
  subset = condition %in% c(
    "CCl4_model",
    "RT_chow"
  )
)

DEN_RT_obj <- subset(
  master_integrated,
  subset = condition %in% c(
    "DEN_model",
    "RT_chow"
  )
)

saveRDS(
  CCl4_RT_obj,
  "CCl4model_RTchow_refinedAnnotated_RNAonly.rds"
)

saveRDS(
  DEN_RT_obj,
  "DENmodel_RTchow_refinedAnnotated_RNAonly.rds"
)


# =============================================================================
# Part 4: create lineage-specific FAT-MASH vs chow objects
# =============================================================================

DefaultAssay(CCl4_RT_obj) <- "RNA"
Idents(CCl4_RT_obj) <- "analysis_label"

CCl4_hepatocytes_obj <- subset(
  CCl4_RT_obj,
  subset = analysis_label == "Hepatocytes"
)

CCl4_endothelial_obj <- subset(
  CCl4_RT_obj,
  subset = analysis_label == "Endothelial cells"
)

CCl4_macrophages_obj <- subset(
  CCl4_RT_obj,
  subset = analysis_label == "Macrophages"
)

CCl4_HSC_obj <- subset(
  CCl4_RT_obj,
  subset = analysis_label == "HSCs"
)

saveRDS(
  CCl4_hepatocytes_obj,
  "CCl4_RT_hepatocytes_RNAonly.rds"
)

saveRDS(
  CCl4_endothelial_obj,
  "CCl4_RT_endothelial_RNAonly.rds"
)

saveRDS(
  CCl4_macrophages_obj,
  "CCl4_RT_macrophages_RNAonly.rds"
)

saveRDS(
  CCl4_HSC_obj,
  "CCl4_RT_HSC_RNAonly.rds"
)


# =============================================================================
# Part 5: create lineage-specific HOT-MASH vs chow objects
# =============================================================================

DefaultAssay(DEN_RT_obj) <- "RNA"
Idents(DEN_RT_obj) <- "analysis_label"

DEN_hepatocytes_obj <- subset(
  DEN_RT_obj,
  subset = analysis_label == "Hepatocytes"
)

DEN_endothelial_obj <- subset(
  DEN_RT_obj,
  subset = analysis_label == "Endothelial cells"
)

DEN_macrophages_obj <- subset(
  DEN_RT_obj,
  subset = analysis_label == "Macrophages"
)

DEN_HSC_obj <- subset(
  DEN_RT_obj,
  subset = analysis_label == "HSCs"
)

saveRDS(
  DEN_hepatocytes_obj,
  "DEN_RT_hepatocytes_RNAonly.rds"
)

saveRDS(
  DEN_endothelial_obj,
  "DEN_RT_endothelial_RNAonly.rds"
)

saveRDS(
  DEN_macrophages_obj,
  "DEN_RT_macrophages_RNAonly.rds"
)

saveRDS(
  DEN_HSC_obj,
  "DEN_RT_HSC_RNAonly.rds"
)


# =============================================================================
# Part 6: create model-vs-model objects
# =============================================================================

CCl4_DEN_obj <- subset(
  master_integrated,
  subset = condition %in% c(
    "CCl4_model",
    "DEN_model"
  )
)

saveRDS(
  CCl4_DEN_obj,
  "Models_refinedAnnotated_RNAonly.rds"
)

DefaultAssay(CCl4_DEN_obj) <- "RNA"
Idents(CCl4_DEN_obj) <- "analysis_label"


# Shared non-malignant lineages
Models_hepatocytes_obj <- subset(
  CCl4_DEN_obj,
  subset = analysis_label == "Hepatocytes"
)

Models_endothelial_obj <- subset(
  CCl4_DEN_obj,
  subset = analysis_label == "Endothelial cells"
)

Models_macrophages_obj <- subset(
  CCl4_DEN_obj,
  subset = analysis_label == "Macrophages"
)

Models_HSC_obj <- subset(
  CCl4_DEN_obj,
  subset = analysis_label == "HSCs"
)

# Tumor-enriched hepatocytes from cluster 13
Models_malignantHeps_obj <- subset(
  CCl4_DEN_obj,
  subset = analysis_label == "Malignant hepatocytes"
)


saveRDS(
  Models_hepatocytes_obj,
  "Models_hepatocytes_RNAonly.rds"
)

saveRDS(
  Models_endothelial_obj,
  "Models_endothelial_RNAonly.rds"
)

saveRDS(
  Models_macrophages_obj,
  "Models_macrophages_RNAonly.rds"
)

saveRDS(
  Models_HSC_obj,
  "Models_HSC_RNAonly.rds"
)

saveRDS(
  Models_malignantHeps_obj,
  "Models_malignantHeps_RNAonly.rds"
)


# =============================================================================
# Part 7: final checks
# =============================================================================

table(
  CCl4_RT_obj$analysis_label,
  CCl4_RT_obj$condition
)

table(
  DEN_RT_obj$analysis_label,
  DEN_RT_obj$condition
)

table(
  CCl4_DEN_obj$analysis_label,
  CCl4_DEN_obj$condition
)
