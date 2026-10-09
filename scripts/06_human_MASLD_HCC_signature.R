# Human MASLD-HCC tumor-signature derivation and scoring
# Supports Figure 4B-D using GSE164760.
#
# The human signature was derived from HCC tumors compared with a composite
# non-tumor background of adjacent liver, cirrhotic liver, and NASH liver.
# Genes were retained at adjusted P < 0.05 and log2 fold change > 1.
#
# The resulting signature was scored in mouse tumor-cell pseudobulk profiles
# using singscore and in individual tumor cells using UCell.
#
# Internal condition nomenclature retained from the original analysis:
#   CCl4_model = FAT-MASH
#   DEN_model  = HOT-MASH
#   RT_chow    = chow control
#
# Paths to public GEO data may need to be updated.

library(hgu219.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(tidyr)
library(limma)
library(singscore)
library(pheatmap)
library(Seurat)
library(UCell)
library(ggplot2)


# =============================================================================
# Part 1: derive the human MASLD-HCC tumor signature from GSE164760
# =============================================================================

# Expression matrix downloaded from GEO.
# First column contains Affymetrix probe IDs.
GSE164760_expr_matrix <- read.delim(
  "path/to/GSE164760_expression_matrix.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

probe_ids <- rownames(
  GSE164760_expr_matrix
)


# -----------------------------------------------------------------------------
# Map Affymetrix probe IDs to gene symbols
# -----------------------------------------------------------------------------

gene_symbols <- mapIds(
  hgu219.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

GSE164760_expr_matrix$gene <- gene_symbols

GSE164760_expr_matrix <- GSE164760_expr_matrix[
  !is.na(GSE164760_expr_matrix$gene),
]


# Collapse multiple probes mapping to the same gene by mean expression
GSE164760_expr_gene <- GSE164760_expr_matrix %>%
  group_by(gene) %>%
  summarise(
    across(
      where(is.numeric),
      mean
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("gene")


# -----------------------------------------------------------------------------
# Define sample groups from GEO column names
# -----------------------------------------------------------------------------

group <- case_when(
  grepl(
    "^HCC.tum",
    colnames(GSE164760_expr_gene)
  ) ~ "HCC_tumor",

  grepl(
    "^HCC.nontumadj",
    colnames(GSE164760_expr_gene)
  ) ~ "adjacent",

  grepl(
    "^cirrhotic",
    colnames(GSE164760_expr_gene)
  ) ~ "cirrhotic",

  grepl(
    "^nash",
    colnames(GSE164760_expr_gene)
  ) ~ "nash",

  grepl(
    "^healthy",
    colnames(GSE164760_expr_gene)
  ) ~ "healthy",

  TRUE ~ NA_character_
)

stopifnot(
  !anyNA(group)
)

group <- factor(group)

table(group)


# -----------------------------------------------------------------------------
# Log2 transform expression matrix
# -----------------------------------------------------------------------------

# The downloaded matrix was analyzed after log2 transformation in the
# original workflow.
expr <- log2(
  GSE164760_expr_gene + 1
)


# -----------------------------------------------------------------------------
# Differential expression: HCC tumor vs composite non-tumor background
# -----------------------------------------------------------------------------

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- levels(
  group
)

contrast_matrix <- makeContrasts(
  HCC_vs_nonTumor =
    HCC_tumor -
    (adjacent + cirrhotic + nash) / 3,
  levels = design
)

fit <- lmFit(
  expr,
  design
)

fit <- contrasts.fit(
  fit,
  contrast_matrix
)

fit <- eBayes(
  fit
)

res <- topTable(
  fit,
  number = Inf,
  sort.by = "P"
)


# Final human MASLD-HCC-associated tumor signature
MASLD_HCC_human <- rownames(res)[
  res$logFC > 1 &
    res$adj.P.Val < 0.05
]

MASLD_HCC_human


# Exact final signature retained from the original analysis
expected_MASLD_HCC_human <- c(
  "CAP2",
  "AKR1B15",
  "TXNRD1",
  "AKR1C3",
  "AKR1B10",
  "SPINK1",
  "FAT1",
  "GGH",
  "CSTB",
  "GLUL",
  "UGT2B11",
  "LCN2",
  "FGGY",
  "TUBA1C",
  "ANXA2P2",
  "TUBA1B",
  "FAM72B",
  "PDIA3",
  "CLVS1",
  "ARHGAP27P1",
  "C8orf76",
  "VMP1",
  "SPP1",
  "TBCA",
  "CUTA",
  "HNRNPA2B1",
  "IGFBP1"
)

stopifnot(
  setequal(
    MASLD_HCC_human,
    expected_MASLD_HCC_human
  )
)


# =============================================================================
# Part 2: pseudobulk singscore in mouse tumor cells
# =============================================================================

# VST-normalized tumor-cell pseudobulk matrix with mouse genes mapped to
# human ortholog symbols. This derived matrix is deposited with the repository.
tumorCells.VSTnorm.humanOrths <- read.delim(
  "data/tumorCells.VSTnorm.humanOrths.txt",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  all(
    c(
      "gene",
      "ortholog_name"
    ) %in% colnames(
      tumorCells.VSTnorm.humanOrths
    )
  )
)

expr_mat <- tumorCells.VSTnorm.humanOrths %>%
  select(-gene) %>%
  column_to_rownames("ortholog_name") %>%
  as.matrix()


# Number of human signature genes represented in the mouse pseudobulk matrix
length(
  intersect(
    MASLD_HCC_human,
    rownames(expr_mat)
  )
)


# Rank genes within each pseudobulk sample and calculate singscore
rankData <- rankGenes(
  expr_mat
)

MASLD_HCC_score <- simpleScore(
  rankData,
  upSet = MASLD_HCC_human
)

pseudobulk_scores <- data.frame(
  sample = colnames(expr_mat),
  MASLD_HCC_score = MASLD_HCC_score$TotalScore
)


# Recover model identity from original sample names
pseudobulk_scores$condition <- case_when(
  grepl(
    "^CCl4_model_",
    pseudobulk_scores$sample
  ) ~ "CCl4_model",

  grepl(
    "^DEN_model_",
    pseudobulk_scores$sample
  ) ~ "DEN_model",

  TRUE ~ NA_character_
)

pseudobulk_scores$condition <- factor(
  pseudobulk_scores$condition,
  levels = c(
    "CCl4_model",
    "DEN_model"
  )
)

pseudobulk_scores


# -----------------------------------------------------------------------------
# Figure 4D: pseudobulk tumor-cell MASLD-HCC signature scores
# -----------------------------------------------------------------------------

score_mat <- matrix(
  pseudobulk_scores$MASLD_HCC_score,
  ncol = 1,
  dimnames = list(
    pseudobulk_scores$sample,
    "MASLD-HCC signature"
  )
)

pheatmap(
  score_mat,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)


# =============================================================================
# Part 3: single-cell UCell scoring
# =============================================================================

master_integrated <- readRDS(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/My Paper - Thesis Plus/revision 2 - finding tumor cells/master integration and annotation/master_integrated_primaryAnnotated.rds"
)

DefaultAssay(
  master_integrated
) <- "RNA"

# Only the RNA assay is required for UCell scoring
if ("SCT" %in% Assays(master_integrated)) {
  master_integrated[["SCT"]] <- NULL
}

if ("integrated" %in% Assays(master_integrated)) {
  master_integrated[["integrated"]] <- NULL
}

gc()

Idents(
  master_integrated
) <- "seurat_clusters"


# -----------------------------------------------------------------------------
# Mouse orthologs of the final human MASLD-HCC signature
# -----------------------------------------------------------------------------

# These are the mouse orthologs used in the original single-cell analysis.
MASLD_HCC_mouse <- c(
  "Cap2",
  "Txnrd1",
  "Akr1b8",
  "Akr1b7",
  "Akr1b10",
  "Spink1",
  "Fat1",
  "Ggh",
  "Cstb",
  "Glul",
  "Lcn2",
  "Fggy",
  "Tuba1c",
  "Tuba1b",
  "Fam72a",
  "Pdia3",
  "Clvs1",
  "9130401M01Rik",
  "Vmp1",
  "Spp1",
  "Tbca",
  "Cuta",
  "Hnrnpa2b1",
  "Igfbp1"
)

obj_genes <- unique(
  rownames(
    master_integrated[["RNA"]]
  )
)

MASLD_HCC_mouse_present <- intersect(
  MASLD_HCC_mouse,
  obj_genes
)

c(
  signature_total = length(
    MASLD_HCC_mouse
  ),
  signature_present = length(
    MASLD_HCC_mouse_present
  )
)


# -----------------------------------------------------------------------------
# Calculate single-cell UCell scores
# -----------------------------------------------------------------------------

master_integrated <- AddModuleScore_UCell(
  obj = master_integrated,
  features = list(
    MASLD_HCC_signature =
      MASLD_HCC_mouse_present
  ),
  assay = "RNA",
  slot = "counts",
  name = NULL
)


# =============================================================================
# Part 4: tumor-cell single-cell and per-mouse analyses
# =============================================================================

# Cluster 13 is the tumor-enriched hepatocyte cluster identified previously.
cl13 <- subset(
  master_integrated,
  idents = "13"
)

# Only disease-model cluster-13 cells were retained as tumor cells.
tumor_cells <- subset(
  cl13,
  subset = condition %in%
    c(
      "CCl4_model",
      "DEN_model"
    )
)

tumor_cells$condition <- factor(
  tumor_cells$condition,
  levels = c(
    "CCl4_model",
    "DEN_model"
  )
)


# -----------------------------------------------------------------------------
# Figure 4B: single-cell MASLD-HCC signature scores
# -----------------------------------------------------------------------------

p_single_cell <- VlnPlot(
  tumor_cells,
  features = "MASLD_HCC_signature",
  group.by = "condition",
  pt.size = 0.15
)

p_single_cell &
  scale_fill_manual(
    values = c(
      "CCl4_model" = "#9EC5DC",
      "DEN_model" = "#EEC464"
    ),
    labels = c(
      "CCl4_model" = "FAT-MASH",
      "DEN_model" = "HOT-MASH"
    )
  ) &
  theme(
    axis.text.y = element_text(
      size = 12
    )
  )


# -----------------------------------------------------------------------------
# Figure 4C: per-mouse mean MASLD-HCC UCell scores
# -----------------------------------------------------------------------------

mouse_scores <- tumor_cells@meta.data %>%
  group_by(
    sample_id,
    condition
  ) %>%
  summarise(
    mean_UCell = mean(
      MASLD_HCC_signature,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

mouse_scores


plot_df <- mouse_scores %>%
  mutate(
    condition = factor(
      condition,
      levels = c(
        "CCl4_model",
        "DEN_model"
      ),
      labels = c(
        "FAT-MASH",
        "HOT-MASH"
      )
    )
  )

ggplot(
  plot_df,
  aes(
    x = condition,
    y = mean_UCell
  )
) +
  geom_boxplot(
    alpha = 0.3,
    outlier.shape = NA
  ) +
  geom_point(
    size = 3,
    position = position_jitter(
      width = 0.08
    )
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean UCell score"
  )
