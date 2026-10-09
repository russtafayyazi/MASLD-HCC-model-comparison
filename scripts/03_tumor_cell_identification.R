# Tumor-cell identification using an HCC-associated gene signature
# Primarily supports Figure 2E-G.
#
# Internal condition nomenclature retained from the original analysis:
#   CCl4_model = FAT-MASH
#   DEN_model  = HOT-MASH
#   RT_chow    = chow control
#
# Original metadata values are preserved because downstream scripts depend
# on them. Manuscript terminology is used only for displayed plot labels.

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(UCell)

# -----------------------------------------------------------------------------
# Load primary annotated atlas
# -----------------------------------------------------------------------------

master_integrated <- readRDS(
  "master_integrated_primaryAnnotated.rds"
)

# Keep RNA assay only to reduce object size
DefaultAssay(master_integrated) <- "RNA"

if ("SCT" %in% Assays(master_integrated)) {
  master_integrated[["SCT"]] <- NULL
}

if ("integrated" %in% Assays(master_integrated)) {
  master_integrated[["integrated"]] <- NULL
}

gc()

Idents(master_integrated) <- "seurat_clusters"


# -----------------------------------------------------------------------------
# HCC-associated UCell signature
# -----------------------------------------------------------------------------

# Literature-informed HCC-associated signature:
# oncofetal, progenitor/cancer stem-cell, and proliferation-associated genes
hcc_module <- c(
  "Afp", "Gpc3", "Igf2", "Dlk1", "Tff3",
  "Epcam", "Prom1", "Cd24a", "Sox9", "Krt19",
  "Mki67", "Top2a", "Ube2c", "Cdk1", "Ccnb1"
)

# Retain genes represented in the dataset
hcc_module <- hcc_module[
  hcc_module %in% rownames(master_integrated)
]

hcc_module
length(hcc_module)

# Calculate single-cell UCell scores
master_integrated <- AddModuleScore_UCell(
  master_integrated,
  features = list(HCC_UCell = hcc_module)
)

# Confirm generated metadata column
grep(
  "HCC_UCell",
  colnames(master_integrated@meta.data),
  value = TRUE
)


# -----------------------------------------------------------------------------
# Figure 2E: HCC-signature enrichment across hepatocyte clusters
# -----------------------------------------------------------------------------

# Restrict visualization to clusters assigned as hepatocytes in the primary
# annotation
hepatocyte_obj <- subset(
  master_integrated,
  subset = primary_label == "Hepatocytes"
)

Idents(hepatocyte_obj) <- "seurat_clusters"

VlnPlot(
  hepatocyte_obj,
  features = "HCC_UCell_UCell",
  group.by = "seurat_clusters"
)

# Cluster-level summary of HCC UCell scores
hcc_cluster_summary <- hepatocyte_obj@meta.data %>%
  mutate(
    cluster = as.character(seurat_clusters)
  ) %>%
  group_by(cluster) %>%
  summarise(
    median_ucell = median(
      HCC_UCell_UCell,
      na.rm = TRUE
    ),
    mean_ucell = mean(
      HCC_UCell_UCell,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(
    desc(median_ucell)
  )

View(hcc_cluster_summary)


# -----------------------------------------------------------------------------
# Figure 2F: HCC-signature enrichment within cluster 13 by condition
# -----------------------------------------------------------------------------

cl13 <- subset(
  hepatocyte_obj,
  idents = "13"
)

cl13$condition <- factor(
  cl13$condition,
  levels = c(
    "RT_chow",
    "CCl4_model",
    "DEN_model"
  )
)

VlnPlot(
  cl13,
  features = "HCC_UCell_UCell",
  group.by = "condition"
) +
  scale_fill_manual(
    values = c(
      "CCl4_model" = "#9EC5DC",
      "DEN_model" = "#EEC464",
      "RT_chow" = "grey70"
    ),
    labels = c(
      "CCl4_model" = "FAT-MASH",
      "DEN_model" = "HOT-MASH",
      "RT_chow" = "Chow"
    )
  )

# Condition-level summary for cluster 13
cl13_condition_summary <- cl13@meta.data %>%
  group_by(condition) %>%
  summarise(
    mean_ucell = mean(
      HCC_UCell_UCell,
      na.rm = TRUE
    ),
    median_ucell = median(
      HCC_UCell_UCell,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

View(cl13_condition_summary)


# -----------------------------------------------------------------------------
# Figure 2G: differential expression of cluster 13 vs other hepatocytes
# -----------------------------------------------------------------------------

hep_obj <- subset(
  master_integrated,
  subset = primary_label == "Hepatocytes"
)

hep_obj$group_13 <- ifelse(
  hep_obj$seurat_clusters == "13",
  "cluster13",
  "other_hepatocytes"
)

Idents(hep_obj) <- "group_13"

# Positive markers enriched in cluster 13 relative to remaining hepatocytes
tumor_markers <- FindMarkers(
  hep_obj,
  ident.1 = "cluster13",
  ident.2 = "other_hepatocytes",
  assay = "RNA",
  only.pos = TRUE
)

View(tumor_markers)

# Candidate genes displayed in Figure 2G
genes_keep <- c(
  "Tff3", "Afp", "Nid1", "Gpc3", "Tspan8",
  "Gldn", "Igdcc4", "Aepb1", "Slpi", "Pnpla5",
  "Col4a3", "Cpe", "Spink1", "Cacna1b", "Inhbb"
)

genes_keep <- genes_keep[
  genes_keep %in% rownames(hep_obj)
]

# Average log-normalized expression in cluster 13 and remaining hepatocytes
avg_exp <- AverageExpression(
  hep_obj,
  assays = "RNA",
  features = genes_keep,
  group.by = "group_13",
  layer = "data"
)$RNA

avg_exp_df <- as.data.frame(avg_exp) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "group",
    values_to = "avg_expression"
  )

# Preserve display order from the candidate-gene list
avg_exp_df$gene <- factor(
  avg_exp_df$gene,
  levels = genes_keep
)

ggplot(
  avg_exp_df,
  aes(
    x = gene,
    y = avg_expression,
    fill = group
  )
) +
  geom_col(
    position = "dodge"
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Average normalized expression",
    fill = NULL
  ) +
  theme_classic()
