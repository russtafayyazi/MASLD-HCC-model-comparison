# Cross-species projection of mouse hepatic lineages onto a human MASLD reference
# Supports Figure 4A.
#
# Human reference:
#   GSE202379
#
# Mouse lineages:
#   Hepatocytes
#   Endothelial cells
#   Macrophages
#   Hepatic stellate cells (HSCs)
#
# Internal condition nomenclature retained from the original analysis:
#   CCl4_model = FAT-MASH
#   DEN_model  = HOT-MASH
#
# File paths reflect the original analysis environment and should be adapted
# as needed.

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(tibble)
library(biomaRt)
library(sparseMatrixStats)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(grid)


# =============================================================================
# Part 1: prepare human GSE202379 reference
# =============================================================================

# Original GSE202379 Seurat object
Gribben <- readRDS(
  "/Users/russtafayyaziv2/Desktop/Medical Science MSc/Thesis/Human data papers/Gribben Nature/GSE202379_SeuratObject_AllCells.rds"
)

DefaultAssay(Gribben) <- "SCT"

# Retain only the SCT assay and remove unnecessary reductions/assays
Gribben@reductions$umap <- NULL
Gribben@reductions$harmony_t.0 <- NULL
Gribben@reductions$umap_harmony_t.0 <- NULL
Gribben[["RNA"]] <- NULL

Gribben <- DietSeurat(
  Gribben,
  assays = "SCT",
  dimreducs = "pca",
  graphs = NULL,
  misc = FALSE,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

# Restrict human reference to the four lineages used for projection
keep_types <- c(
  "Hepatocytes",
  "Endothelial",
  "Stellate",
  "Macrophages"
)

stopifnot(
  all(
    c(
      "cell.annotation",
      "Patient.ID",
      "Disease.status"
    ) %in% colnames(Gribben@meta.data)
  )
)

cells_keep <- colnames(Gribben)[
  Gribben$cell.annotation %in% keep_types
]

Gribben_4L <- subset(
  Gribben,
  cells = cells_keep
)

# Retain metadata used for mapping and interpretation
meta_keep <- c(
  "cell.annotation",
  "Patient.ID",
  "Disease.status",
  "orig.ident",
  "nCount_SCT",
  "nFeature_SCT"
)

meta_keep <- intersect(
  meta_keep,
  colnames(Gribben_4L@meta.data)
)

Gribben_4L@meta.data <- Gribben_4L@meta.data[
  ,
  meta_keep,
  drop = FALSE
]

Gribben_4L$cell.annotation <- factor(
  Gribben_4L$cell.annotation,
  levels = keep_types
)

rm(Gribben)
gc()


# =============================================================================
# Part 2: identify strict one-to-one mouse-human orthologs
# =============================================================================

mm <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

orth_raw <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "with_hsapiens_homolog",
  values = TRUE,
  mart = mm
)

orth_1to1 <- orth_raw %>%
  filter(
    hsapiens_homolog_orthology_type == "ortholog_one2one",
    external_gene_name != "",
    hsapiens_homolog_associated_gene_name != ""
  ) %>%
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  ) %>%
  distinct(
    hsapiens_homolog_ensembl_gene,
    .keep_all = TRUE
  ) %>%
  transmute(
    mouse_ensembl = ensembl_gene_id,
    mouse_symbol = external_gene_name,
    human_ensembl = hsapiens_homolog_ensembl_gene,
    human_symbol = hsapiens_homolog_associated_gene_name
  )


# =============================================================================
# Part 3: align human and mouse objects in ortholog space
# =============================================================================

master_integrated <- readRDS(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/My Paper - Thesis Plus/revision 2 - finding tumor cells/master integration and annotation/master_integrated_primaryAnnotated.rds"
)

DefaultAssay(master_integrated) <- "SCT"

# Only SCT is required for projection mapping
master_integrated[["RNA"]] <- NULL
master_integrated[["integrated"]] <- NULL

# Model cells only
CCl4_model_SCT <- subset(
  master_integrated,
  subset = condition == "CCl4_model"
)

DEN_model_SCT <- subset(
  master_integrated,
  subset = condition == "DEN_model"
)

rm(master_integrated)
gc()

human_genes <- rownames(
  Gribben_4L[["SCT"]]
)

CCl4_genes <- rownames(
  CCl4_model_SCT[["SCT"]]
)

DEN_genes <- rownames(
  DEN_model_SCT[["SCT"]]
)

mouse_union <- union(
  CCl4_genes,
  DEN_genes
)

# Retain strict 1:1 orthologs represented in the datasets
orth_use <- orth_1to1 %>%
  filter(
    mouse_symbol %in% mouse_union,
    human_symbol %in% human_genes
  )

human_keep <- orth_use$human_symbol

CCl4_keep <- orth_use$mouse_symbol[
  orth_use$mouse_symbol %in% CCl4_genes
]

DEN_keep <- orth_use$mouse_symbol[
  orth_use$mouse_symbol %in% DEN_genes
]

Gribben_orth <- subset(
  Gribben_4L,
  features = human_keep
)

CCl4_model_orth <- subset(
  CCl4_model_SCT,
  features = CCl4_keep
)

DEN_model_orth <- subset(
  DEN_model_SCT,
  features = DEN_keep
)

rm(
  Gribben_4L,
  CCl4_model_SCT,
  DEN_model_SCT
)

gc()


# -----------------------------------------------------------------------------
# Rename mouse SCT features using human ortholog symbols
# -----------------------------------------------------------------------------

orth_use$mouse_symbol <- as.character(
  orth_use$mouse_symbol
)

orth_use$human_symbol <- as.character(
  orth_use$human_symbol
)

rename_sct_features_v5 <- function(
    seu,
    map_tbl,
    assay = "SCT"
) {

  stopifnot(
    assay %in% Assays(seu)
  )

  DefaultAssay(seu) <- assay

  old <- rownames(
    seu[[assay]]
  )

  idx <- match(
    old,
    map_tbl$mouse_symbol
  )

  new <- old

  hit <- !is.na(idx)

  new[hit] <- map_tbl$human_symbol[
    idx[hit]
  ]

  new <- make.unique(new)

  counts_mat <- GetAssayData(
    seu,
    assay = assay,
    layer = "counts"
  )

  data_mat <- GetAssayData(
    seu,
    assay = assay,
    layer = "data"
  )

  rownames(counts_mat) <- new
  rownames(data_mat) <- new

  ass <- CreateAssayObject(
    counts = counts_mat
  )

  ass <- SetAssayData(
    ass,
    layer = "data",
    new.data = data_mat
  )

  seu[[assay]] <- ass

  seu
}

CCl4_model_orth <- rename_sct_features_v5(
  CCl4_model_orth,
  orth_use,
  assay = "SCT"
)

DEN_model_orth <- rename_sct_features_v5(
  DEN_model_orth,
  orth_use,
  assay = "SCT"
)


# Enforce identical feature space and feature order
genes_ref <- rownames(
  Gribben_orth[["SCT"]]
)

genes_CCl4 <- rownames(
  CCl4_model_orth[["SCT"]]
)

genes_DEN <- rownames(
  DEN_model_orth[["SCT"]]
)

genes_common <- Reduce(
  intersect,
  list(
    genes_ref,
    genes_CCl4,
    genes_DEN
  )
)

Gribben_orth <- Gribben_orth[
  genes_common,
]

CCl4_model_orth <- CCl4_model_orth[
  genes_common,
]

DEN_model_orth <- DEN_model_orth[
  genes_common,
]

stopifnot(
  identical(
    rownames(Gribben_orth[["SCT"]]),
    rownames(CCl4_model_orth[["SCT"]])
  ),
  identical(
    rownames(Gribben_orth[["SCT"]]),
    rownames(DEN_model_orth[["SCT"]])
  )
)


# =============================================================================
# Part 4: lineage-specific human reference preparation
# =============================================================================

# The original analysis rebuilt the human reference independently for each
# lineage using the 2,500 most variable SCT features followed by scaling,
# 50-component PCA, neighbors, clustering, and UMAP.

build_reference <- function(
    human_obj,
    human_lineage
) {

  Idents(human_obj) <- "cell.annotation"

  ref <- subset(
    human_obj,
    idents = human_lineage
  )

  meta_keep <- c(
    "cell.annotation",
    "Disease.status",
    "Patient.ID"
  )

  meta_keep <- intersect(
    meta_keep,
    colnames(ref@meta.data)
  )

  ref@meta.data <- ref@meta.data[
    ,
    meta_keep,
    drop = FALSE
  ]

  DefaultAssay(ref) <- "SCT"

  VariableFeatures(ref) <- character(0)
  ref@reductions <- list()
  ref@graphs <- list()

  sct_data <- GetAssayData(
    ref,
    assay = "SCT",
    layer = "data"
  )

  gene_variances <- sparseMatrixStats::rowVars(
    sct_data
  )

  gene_variances[
    !is.finite(gene_variances)
  ] <- NA

  gene_variances <- gene_variances[
    !is.na(gene_variances)
  ]

  hvgs <- names(
    sort(
      gene_variances,
      decreasing = TRUE
    )
  )[
    seq_len(
      min(
        2500,
        length(gene_variances)
      )
    )
  ]

  VariableFeatures(ref) <- hvgs

  ref <- ScaleData(
    ref,
    features = hvgs,
    verbose = TRUE
  )

  ref <- RunPCA(
    ref,
    features = hvgs,
    npcs = 50,
    verbose = TRUE
  )

  ref <- FindNeighbors(
    ref,
    dims = 1:50,
    verbose = TRUE
  )

  ref <- FindClusters(
    ref,
    resolution = 0.3,
    verbose = TRUE
  )

  ref <- RunUMAP(
    ref,
    dims = 1:50,
    return.model = TRUE,
    verbose = TRUE
  )

  list(
    reference = ref,
    hvgs = hvgs
  )
}


hep_ref <- build_reference(
  Gribben_orth,
  "Hepatocytes"
)

endo_ref <- build_reference(
  Gribben_orth,
  "Endothelial"
)

macro_ref <- build_reference(
  Gribben_orth,
  "Macrophages"
)

hsc_ref <- build_reference(
  Gribben_orth,
  "Stellate"
)


# =============================================================================
# Part 5: project mouse lineages onto human MASLD disease states
# =============================================================================

# This helper reproduces the repeated lineage-matched mapping workflow used in
# the original analysis: subset lineage, use reference HVGs, scale, run PCA,
# identify SCT transfer anchors, and transfer Disease.status labels.

map_lineage <- function(
    mouse_obj,
    mouse_lineage,
    reference_obj,
    reference_hvgs
) {

  DefaultAssay(mouse_obj) <- "SCT"
  Idents(mouse_obj) <- "primary_label"

  query <- subset(
    mouse_obj,
    idents = mouse_lineage
  )

  VariableFeatures(query) <- reference_hvgs

  query <- ScaleData(
    query,
    features = reference_hvgs,
    verbose = TRUE
  )

  query <- RunPCA(
    query,
    features = reference_hvgs,
    npcs = 50,
    verbose = TRUE
  )

  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query = query,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:50,
    verbose = TRUE
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_obj$Disease.status,
    dims = 1:50,
    verbose = TRUE
  )

  AddMetaData(
    query,
    metadata = predictions
  )
}


# -----------------------------------------------------------------------------
# Hepatocytes
# -----------------------------------------------------------------------------

hep_CCl4 <- map_lineage(
  CCl4_model_orth,
  "Hepatocytes",
  hep_ref$reference,
  hep_ref$hvgs
)

hep_DEN <- map_lineage(
  DEN_model_orth,
  "Hepatocytes",
  hep_ref$reference,
  hep_ref$hvgs
)


# -----------------------------------------------------------------------------
# Endothelial cells
# -----------------------------------------------------------------------------

endo_CCl4 <- map_lineage(
  CCl4_model_orth,
  "Endothelial cells",
  endo_ref$reference,
  endo_ref$hvgs
)

endo_DEN <- map_lineage(
  DEN_model_orth,
  "Endothelial cells",
  endo_ref$reference,
  endo_ref$hvgs
)


# -----------------------------------------------------------------------------
# Macrophages
# -----------------------------------------------------------------------------

macro_CCl4 <- map_lineage(
  CCl4_model_orth,
  "Macrophages",
  macro_ref$reference,
  macro_ref$hvgs
)

macro_DEN <- map_lineage(
  DEN_model_orth,
  "Macrophages",
  macro_ref$reference,
  macro_ref$hvgs
)


# -----------------------------------------------------------------------------
# Hepatic stellate cells
# -----------------------------------------------------------------------------

hsc_CCl4 <- map_lineage(
  CCl4_model_orth,
  "HSCs",
  hsc_ref$reference,
  hsc_ref$hvgs
)

hsc_DEN <- map_lineage(
  DEN_model_orth,
  "HSCs",
  hsc_ref$reference,
  hsc_ref$hvgs
)


# =============================================================================
# Part 6: summarize human disease-stage prediction scores
# =============================================================================

summarize_projection <- function(
    query,
    model,
    lineage
) {

  score_cols <- grep(
    "^prediction.score\\.",
    colnames(query@meta.data),
    value = TRUE
  )

  # Exclude Seurat's max prediction score, if present
  score_cols <- setdiff(
    score_cols,
    "prediction.score.max"
  )

  model_summary <- query@meta.data %>%
    summarise(
      across(
        all_of(score_cols),
        ~ mean(.x, na.rm = TRUE)
      )
    ) %>%
    mutate(
      Model = model,
      CellType = lineage,
      .before = 1
    )

  per_mouse <- query@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      across(
        all_of(score_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Model = model,
      CellType = lineage,
      .after = "orig.ident"
    )

  list(
    model = model_summary,
    per_mouse = per_mouse
  )
}


projection_summaries <- list(
  summarize_projection(
    hep_CCl4,
    "FAT-MASH",
    "Hepatocyte"
  ),
  summarize_projection(
    endo_CCl4,
    "FAT-MASH",
    "Endothelial"
  ),
  summarize_projection(
    macro_CCl4,
    "FAT-MASH",
    "Macrophage"
  ),
  summarize_projection(
    hsc_CCl4,
    "FAT-MASH",
    "HSC"
  ),
  summarize_projection(
    hep_DEN,
    "HOT-MASH",
    "Hepatocyte"
  ),
  summarize_projection(
    endo_DEN,
    "HOT-MASH",
    "Endothelial"
  ),
  summarize_projection(
    macro_DEN,
    "HOT-MASH",
    "Macrophage"
  ),
  summarize_projection(
    hsc_DEN,
    "HOT-MASH",
    "HSC"
  )
)

projection_model_means <- bind_rows(
  lapply(
    projection_summaries,
    `[[`,
    "model"
  )
)

projection_per_mouse <- bind_rows(
  lapply(
    projection_summaries,
    `[[`,
    "per_mouse"
  )
)

projection_model_means
projection_per_mouse


# =============================================================================
# Part 7: Figure 4A summary heatmap
# =============================================================================
#
# The final figure was assembled in PowerPoint. The values below are the
# rounded model-by-lineage mean prediction scores used in the original
# Figure 4A heatmap.

figure4A_values <- tribble(
  ~Model, ~CellType, ~Healthy, ~NAFLD,
  ~`NASH w/o cirrhosis`, ~`NASH w cirrhosis`, ~`End stage`,

  "FAT-MASH", "Hepatocyte",
  0.0489, 0.0104, 0.5787, 0.0048, 0.3571,

  "FAT-MASH", "Endothelial",
  0.0468, 0.0606, 0.4788, 0.4092, 0.0045,

  "FAT-MASH", "Macrophage",
  0.0356, 0.0543, 0.7961, 0.0658, 0.0482,

  "FAT-MASH", "HSC",
  0.0182, 0.1015, 0.6963, 0.1081, 0.0759,

  "HOT-MASH", "Hepatocyte",
  0.0441, 0.0053, 0.8626, 0.0011, 0.0870,

  "HOT-MASH", "Endothelial",
  0.0173, 0.0055, 0.7199, 0.2573, 0.0000,

  "HOT-MASH", "Macrophage",
  0.0103, 0.0109, 0.9441, 0.0278, 0.0070,

  "HOT-MASH", "HSC",
  0.0098, 0.0573, 0.7170, 0.2147, 0.0011
)

stage_cols <- c(
  "Healthy",
  "NAFLD",
  "NASH w/o cirrhosis",
  "NASH w cirrhosis",
  "End stage"
)

figure4A_values <- figure4A_values %>%
  select(
    Model,
    CellType,
    all_of(stage_cols)
  )

figure4A_rows <- figure4A_values %>%
  mutate(
    Row = paste(
      Model,
      CellType,
      sep = " — "
    )
  )

figure4A_mat <- figure4A_rows %>%
  column_to_rownames("Row") %>%
  select(
    -Model,
    -CellType
  ) %>%
  as.matrix()

row_split <- factor(
  figure4A_rows$Model,
  levels = c(
    "FAT-MASH",
    "HOT-MASH"
  )
)

names(row_split) <- rownames(
  figure4A_mat
)

lineage_order <- c(
  "Hepatocyte",
  "Endothelial",
  "Macrophage",
  "HSC"
)

row_order <- order(
  row_split,
  match(
    figure4A_rows$CellType,
    lineage_order
  )
)

figure4A_mat <- figure4A_mat[
  row_order,
  ,
  drop = FALSE
]

row_split <- row_split[
  row_order
]

col_fun <- colorRamp2(
  seq(
    0,
    1,
    length.out = 7
  ),
  magma(7)
)

ht <- Heatmap(
  figure4A_mat,
  name = "Mean probability",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  row_title = NULL,
  column_title = NULL,
  column_names_rot = 45,
  border = TRUE,
  heatmap_legend_param = list(
    at = seq(
      0,
      1,
      by = 0.25
    ),
    labels = seq(
      0,
      1,
      by = 0.25
    )
  )
)

draw(
  ht,
  heatmap_legend_side = "right"
)
