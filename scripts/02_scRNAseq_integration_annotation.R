# scRNA-seq preprocessing, QC, integration, clustering, and primary annotation
# Primarily supports Figure 2 and provides the integrated object used downstream.
#
# Internal sample nomenclature retained from the original analysis:
#   CCl4_model = FAT-MASH
#   DEN_model  = HOT-MASH
#   RT_chow    = chow control
#
# Original sample IDs and metadata values are intentionally preserved because
# downstream scripts depend on them.

library(Seurat)
library(dplyr)
library(ggplot2)

# helper: wrap a matrix into a Seurat object with consistent metadata + QC fields
make_obj <- function(mat,
                     sample_id,
                     role,
                     group,
                     housing_temp,
                     min_features = 200) {
  # basic checks to catch input mistakes early
  stopifnot(role %in% c("model", "control"))
  stopifnot(group %in% c("CCl4", "DEN", "RT_chow"))
  stopifnot(housing_temp %in% c("RT", "TN"))

  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.features = min_features
  )

  # core metadata
  obj$sample_id    <- sample_id
  obj$role         <- role
  obj$group        <- group
  obj$model_label  <- group
  obj$housing_temp <- housing_temp

  # useful combined metadata for later plotting / tables
  obj$condition <- dplyr::case_when(
    group == "RT_chow" ~ "RT_chow",
    group == "CCl4"    ~ "CCl4_model",
    group == "DEN"     ~ "DEN_model"
  )

  # QC metrics (mouse conventions)
  obj[["pct.mt"]]   <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["pct.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rps|^Rpl|^RPS|^RPL")

  return(obj)
}

# ---- read matrices from disk ----
CCl4_model_45.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/CCl4_V/Feature:Barcode Matrix (sample)/CCl4_V_45_sample.h5"
)

CCl4_model_66.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/CCl4_V/Feature:Barcode Matrix (sample)/CCl4_V_66_sample.h5"
)

CCl4_model_67.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/CCl4_V/Feature:Barcode Matrix (sample)/CCl4_V_67_sample.h5"
)

DEN_model_451.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/DEN_V/Feature:Barcode Matrix (sample)/DEN_V_451_sample.h5"
)

DEN_model_487.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/DEN_V/Feature:Barcode Matrix (sample)/DEN_V_487_sample.h5"
)

DEN_model_508.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/DEN_V/Feature:Barcode Matrix (sample)/DEN_V_508_sample.h5"
)

RT_chow_1.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/Ctrl_RT/Feature:Barcode Matrix (sample)/Ctrl_RT_1_sample.h5"
)

RT_chow_2.data <- Read10X_h5(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/GMEL single cell project/Demultiplexed data files/Ctrl_RT/Feature:Barcode Matrix (sample)/Ctrl_RT_2_sample.h5"
)

# ---- create Seurat objects with original sample naming ----
CCl4_model_45 <- make_obj(
  CCl4_model_45.data,
  sample_id = "CCl4_model_45",
  role = "model",
  group = "CCl4",
  housing_temp = "RT"
)

CCl4_model_66 <- make_obj(
  CCl4_model_66.data,
  sample_id = "CCl4_model_66",
  role = "model",
  group = "CCl4",
  housing_temp = "RT"
)

CCl4_model_67 <- make_obj(
  CCl4_model_67.data,
  sample_id = "CCl4_model_67",
  role = "model",
  group = "CCl4",
  housing_temp = "RT"
)

DEN_model_451 <- make_obj(
  DEN_model_451.data,
  sample_id = "DEN_model_451",
  role = "model",
  group = "DEN",
  housing_temp = "TN"
)

DEN_model_487 <- make_obj(
  DEN_model_487.data,
  sample_id = "DEN_model_487",
  role = "model",
  group = "DEN",
  housing_temp = "TN"
)

DEN_model_508 <- make_obj(
  DEN_model_508.data,
  sample_id = "DEN_model_508",
  role = "model",
  group = "DEN",
  housing_temp = "TN"
)

RT_chow_1 <- make_obj(
  RT_chow_1.data,
  sample_id = "RT_chow_1",
  role = "control",
  group = "RT_chow",
  housing_temp = "RT"
)

RT_chow_2 <- make_obj(
  RT_chow_2.data,
  sample_id = "RT_chow_2",
  role = "control",
  group = "RT_chow",
  housing_temp = "RT"
)

sample_list <- list(
  CCl4_model_45 = CCl4_model_45,
  CCl4_model_66 = CCl4_model_66,
  CCl4_model_67 = CCl4_model_67,
  DEN_model_451 = DEN_model_451,
  DEN_model_487 = DEN_model_487,
  DEN_model_508 = DEN_model_508,
  RT_chow_1     = RT_chow_1,
  RT_chow_2     = RT_chow_2
)

names(sample_list)

rm(
  CCl4_model_45, CCl4_model_45.data,
  CCl4_model_66, CCl4_model_66.data,
  CCl4_model_67, CCl4_model_67.data,
  DEN_model_451, DEN_model_451.data,
  DEN_model_487, DEN_model_487.data,
  DEN_model_508, DEN_model_508.data,
  RT_chow_1, RT_chow_1.data,
  RT_chow_2, RT_chow_2.data
)

gc()

# ---- pre-QC summary table ----
pre_qc_tbl <- bind_rows(
  lapply(sample_list, function(obj) {
    obj@meta.data %>%
      summarize(
        sample_id         = dplyr::first(sample_id),
        role              = dplyr::first(role),
        group             = dplyr::first(group),
        condition         = dplyr::first(condition),
        housing_temp      = dplyr::first(housing_temp),
        cells_pre_qc      = n(),
        med_genes_pre_qc  = median(nFeature_RNA),
        med_umis_pre_qc   = median(nCount_RNA),
        med_mt_pre_qc     = median(pct.mt),
        med_ribo_pre_qc   = median(pct.ribo)
      )
  })
) %>%
  arrange(group, sample_id)

View(pre_qc_tbl)

# ---- QC violin plots before filtering ----
plot_qc <- function(obj, title_text) {
  VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "pct.mt"),
    pt.size = 0.1,
    ncol = 3
  ) +
    ggtitle(title_text)
}

# generate and print one plot per sample
qc_plot_list <- lapply(names(sample_list), function(nm) {
  plot_qc(sample_list[[nm]], paste0(nm, " (pre-QC)"))
})

names(qc_plot_list) <- names(sample_list)

for (nm in names(qc_plot_list)) {
  print(qc_plot_list[[nm]])
}

# ---- adaptive QC filtering across all samples ----
filter_one <- function(obj) {
  md <- obj@meta.data

  nf_lo_q <- quantile(md$nFeature_RNA, 0.01, na.rm = TRUE)
  nf_hi_q <- quantile(md$nFeature_RNA, 0.995, na.rm = TRUE)

  nc_lo_q <- quantile(md$nCount_RNA, 0.01, na.rm = TRUE)
  nc_hi_q <- quantile(md$nCount_RNA, 0.995, na.rm = TRUE)

  keep <- with(
    md,
    nFeature_RNA >= max(200, nf_lo_q) &
      nFeature_RNA <= nf_hi_q &
      nCount_RNA >= nc_lo_q &
      nCount_RNA <= nc_hi_q &
      pct.mt <= 25
  )

  subset(obj, cells = rownames(md)[keep])
}

sample_list_qc <- lapply(sample_list, filter_one)

# ---- cells post-QC summary ----
post_qc_tbl <- bind_rows(
  lapply(sample_list_qc, function(obj) {
    obj@meta.data %>%
      summarize(
        sample_id = dplyr::first(sample_id),
        cells_post_qc = n()
      )
  })
)

print(post_qc_tbl)

qc_retention <- pre_qc_tbl %>%
  select(sample_id, cells_pre_qc) %>%
  left_join(post_qc_tbl, by = "sample_id") %>%
  mutate(
    retained_pct = round(
      100 * cells_post_qc / cells_pre_qc,
      1
    )
  )

print(qc_retention)

# ---- save and load checkpoint ----
saveRDS(
  sample_list_qc,
  "all_samples_post_QC.rds"
)

rm(
  sample_list,
  sample_list_qc,
  qc_plot_list
)

gc()

sample_list_qc <- readRDS(
  "all_samples_post_QC.rds"
)

# ---- doublet detection per sample ----
library(SingleCellExperiment)
library(scDblFinder)

run_scdbl <- function(seu, sample_name, seed = 123) {

  sce <- as.SingleCellExperiment(seu)

  set.seed(seed)

  sce <- scDblFinder(sce)

  df <- as.data.frame(
    colData(sce)
  )

  seu$doublet_class <- df$scDblFinder.class
  seu$doublet_score <- df$scDblFinder.score

  singlet_cells <- colnames(seu)[
    seu$doublet_class == "singlet"
  ]

  seu_singlets <- subset(
    seu,
    cells = singlet_cells
  )

  message(
    sample_name,
    ": removed ",
    sum(seu$doublet_class != "singlet"),
    " doublets (",
    round(
      mean(seu$doublet_class != "singlet") * 100,
      1
    ),
    "%)"
  )

  return(seu_singlets)
}

sample_list_singlets <- lapply(
  names(sample_list_qc),
  function(nm) {
    run_scdbl(
      sample_list_qc[[nm]],
      sample_name = nm
    )
  }
)

names(sample_list_singlets) <- names(sample_list_qc)

# ---- post-doublet summary ----
post_doublet_tbl <- bind_rows(
  lapply(sample_list_singlets, function(obj) {
    obj@meta.data %>%
      summarize(
        sample_id = dplyr::first(sample_id),
        cells_post_doublet = n()
      )
  })
)

print(post_doublet_tbl)

doublet_retention <- qc_retention %>%
  left_join(
    post_doublet_tbl,
    by = "sample_id"
  ) %>%
  mutate(
    retained_pct_after_doublet = round(
      100 * cells_post_doublet / cells_post_qc,
      1
    )
  )

View(doublet_retention)

# ---- save and load checkpoint ----
saveRDS(
  sample_list_singlets,
  "all_samples_post_QC_post_doublet.rds"
)

rm(
  sample_list_qc,
  sample_list_singlets
)

gc()

sample_list_singlets <- readRDS(
  "all_samples_post_QC_post_doublet.rds"
)

# ---- SCTransform per sample ----
library(future)

plan("sequential")

options(
  future.globals.maxSize = 8 * 1024^3
)

sample_list_sct <- lapply(
  sample_list_singlets,
  function(obj) {

    SCTransform(
      obj,
      method = "glmGamPoi",
      vars.to.regress = "pct.mt",
      verbose = TRUE
    )
  }
)

# ---- save and load checkpoint ----
saveRDS(
  sample_list_sct,
  "all_samples_post_SCT.rds"
)

rm(
  sample_list_singlets,
  sample_list_sct
)

gc()

sample_list_sct <- readRDS(
  "all_samples_post_SCT.rds"
)

# ---- integration setup ----
library(Seurat)

obj_list <- sample_list_sct

rm(sample_list_sct)

gc()

integration_features <- SelectIntegrationFeatures(
  object.list = obj_list,
  nfeatures = 3000
)

obj_list <- PrepSCTIntegration(
  object.list = obj_list,
  anchor.features = integration_features
)

obj_list <- lapply(
  obj_list,
  function(x) {

    DefaultAssay(x) <- "SCT"

    RunPCA(
      x,
      features = integration_features,
      npcs = 30,
      verbose = FALSE
    )
  }
)

# ---- pre-integration merged object for diagnostic UMAP ----
master_merged_pre <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  project = "MASLD_HCC_master_preintegration"
)

DefaultAssay(master_merged_pre) <- "SCT"

master_merged_pre <- RunPCA(
  master_merged_pre,
  features = integration_features,
  npcs = 30,
  verbose = FALSE
)

master_merged_pre <- RunUMAP(
  master_merged_pre,
  dims = 1:30
)

DimPlot(
  master_merged_pre,
  reduction = "umap",
  group.by = "sample_id",
  pt.size = 0.1
) +
  ggtitle("Master atlas pre-integration — by sample")

DimPlot(
  master_merged_pre,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.1
) +
  ggtitle("Master atlas pre-integration — by condition")

rm(master_merged_pre)

gc()

# ---- RPCA anchor finding ----
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reduction = "rpca",
  dims = 1:30
)

rm(obj_list)

gc()

# ---- integrate all 8 samples ----
# This step was completed on the workstation for memory purposes.
master_integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30
)

# ---- post-integration dimensional reduction and clustering ----
DefaultAssay(master_integrated) <- "integrated"

master_integrated <- RunPCA(
  master_integrated,
  npcs = 30,
  verbose = FALSE
)

master_integrated <- FindNeighbors(
  master_integrated,
  dims = 1:30
)

master_integrated <- FindClusters(
  master_integrated,
  resolution = 0.4
)

master_integrated <- RunUMAP(
  master_integrated,
  dims = 1:30
)

DimPlot(
  master_integrated,
  reduction = "umap",
  group.by = "sample_id",
  pt.size = 0.1
)

DimPlot(
  master_integrated,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.1
)

DimPlot(
  master_integrated,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.1
)

# ---- save and load checkpoint ----
saveRDS(
  master_integrated,
  "master_integrated_allConditions.rds"
)

rm(anchors)

gc()

master_integrated <- readRDS(
  "master_integrated_allConditions.rds"
)

table(
  Idents(master_integrated)
)

table(
  master_integrated$condition,
  Idents(master_integrated)
)

# ---- add normalized layer to RNA assay ----
DefaultAssay(master_integrated) <- "RNA"

master_integrated <- JoinLayers(
  master_integrated,
  assay = "RNA"
)

master_integrated <- NormalizeData(
  master_integrated,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE
)

# ---- save and load checkpoint ----
saveRDS(
  master_integrated,
  "master_integrated_joined_logged.rds"
)

gc()

master_integrated <- readRDS(
  "master_integrated_joined_logged.rds"
)

# -----------------------------------------------------------------------------
# Primary annotation
# -----------------------------------------------------------------------------

DefaultAssay(master_integrated) <- "RNA"

Idents(master_integrated) <- "seurat_clusters"

# ---- broad lineage marker panel ----
marker_panel <- c(

  # Hepatocytes
  "Alb", "Ttr", "Apoa1", "Ass1",
  "Cyp2e1", "Glul", "Cyp2f2",

  # Cholangiocytes
  "Krt19", "Krt7", "Sox9", "Epcam",

  # Endothelial
  "Kdr", "Pecam1", "Cdh5", "Tek",
  "Klf2", "Dnase1l3", "Clec4g",

  # HSC / fibroblast / mesenchymal
  "Rgs5", "Pdgfrb", "Des", "Reln",
  "Col1a1", "Col3a1", "Acta2", "Tagln",

  # Kupffer / macrophage / monocyte
  "Clec4f", "Marco", "Adgre1", "Lyz2",
  "Cd68", "Ccr2", "Itgam",

  # T / NK
  "Cd3d", "Cd3e", "Cd3g", "Trbc1",
  "Trbc2", "Nkg7", "Xcl1", "Klrb1c",

  # B / plasma
  "Cd79a", "Ms4a1", "Cd74",
  "H2-Ab1", "Jchain", "Mzb1",

  # Dendritic
  "Flt3", "Itgax", "Xcr1",
  "Clec10a", "Ccr9",

  # Neutrophils
  "S100a8", "S100a9", "Ly6g",
  "Retnlg", "Ccl4",

  # HCC-associated flags
  "Afp", "Gpc3"
)

marker_panel <- marker_panel[
  marker_panel %in% rownames(master_integrated)
]

p_dot <- DotPlot(
  master_integrated,
  features = marker_panel
) +
  RotatedAxis() +
  ggtitle(
    "Broad liver lineage markers + HCC-associated markers"
  )

print(p_dot)

# ---- identify cluster-enriched markers for annotation review ----
all_markers <- FindAllMarkers(
  master_integrated,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.1
)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 30
  ) %>%
  ungroup()

str(top_markers)

write.table(
  top_markers,
  file = "top_markers.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

library(tidyr)

# ---- cluster-by-sample representation ----
cluster_sample_long <- master_integrated@meta.data %>%
  as.data.frame() %>%
  mutate(
    cluster = as.character(seurat_clusters)
  ) %>%
  count(
    cluster,
    sample_id,
    name = "n"
  ) %>%
  group_by(cluster) %>%
  mutate(
    total_in_cluster = sum(n),
    pct_in_cluster = round(
      100 * n / total_in_cluster,
      1
    )
  ) %>%
  ungroup()

sample_order <- sort(
  unique(master_integrated$sample_id)
)

cluster_sample_counts <- cluster_sample_long %>%
  select(
    cluster,
    sample_id,
    n
  ) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = n,
    values_fill = 0
  ) %>%
  select(
    cluster,
    all_of(sample_order)
  )

cluster_sample_pct <- cluster_sample_long %>%
  select(
    cluster,
    sample_id,
    pct_in_cluster
  ) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = pct_in_cluster,
    values_fill = 0,
    names_glue = "{sample_id}_pct"
  ) %>%
  select(
    cluster,
    all_of(
      paste0(sample_order, "_pct")
    )
  )

cluster_totals <- cluster_sample_long %>%
  distinct(
    cluster,
    total_in_cluster
  )

cluster_by_mouse <- cluster_totals %>%
  left_join(
    cluster_sample_counts,
    by = "cluster"
  ) %>%
  left_join(
    cluster_sample_pct,
    by = "cluster"
  ) %>%
  mutate(
    cluster_num = suppressWarnings(
      as.numeric(cluster)
    )
  ) %>%
  arrange(cluster_num) %>%
  select(-cluster_num)

View(cluster_by_mouse)

# ---- apply primary annotations ----
Idents(master_integrated) <- "seurat_clusters"

cl_map <- c(
  `0`  = "Hepatocytes",
  `1`  = "Hepatocytes",
  `2`  = "Hepatocytes",
  `3`  = "Hepatocytes",
  `4`  = "Hepatocytes",
  `5`  = "Hepatocytes",
  `6`  = "Hepatocytes",
  `7`  = "Endothelial cells",
  `8`  = "Hepatocytes",
  `9`  = "Macrophages",
  `10` = "Macrophages",
  `11` = "HSCs",
  `12` = "Hepatocytes",
  `13` = "Hepatocytes",
  `14` = "Cholangiocytes"
)

clusters_present <- sort(
  unique(
    as.character(
      master_integrated$seurat_clusters
    )
  )
)

if (
  length(
    setdiff(
      clusters_present,
      names(cl_map)
    )
  ) > 0
) {
  stop(
    "Missing clusters in map: ",
    paste(
      setdiff(
        clusters_present,
        names(cl_map)
      ),
      collapse = ", "
    )
  )
}

master_integrated$primary_label <- unname(
  cl_map[
    as.character(
      master_integrated$seurat_clusters
    )
  ]
)

master_integrated$lineage <- master_integrated$primary_label

master_integrated$lineage <- factor(
  master_integrated$lineage,
  levels = c(
    "Hepatocytes",
    "Endothelial cells",
    "Macrophages",
    "HSCs",
    "Cholangiocytes"
  )
)

# ---- sanity checks ----
table(
  master_integrated$primary_label
)

table(
  master_integrated$primary_label,
  master_integrated$condition
)

DimPlot(
  master_integrated,
  reduction = "umap",
  group.by = "primary_label",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  NoLegend()

# ---- save primary annotated object ----
saveRDS(
  master_integrated,
  "master_integrated_primaryAnnotated.rds"
)

master_integrated <- readRDS(
  "master_integrated_primaryAnnotated.rds"
)
