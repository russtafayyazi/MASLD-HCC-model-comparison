# CellChat intercellular communication analysis
#
# Supports Figure 7:
#   A - outgoing and incoming inferred interaction strength
#   B - FAT-MASH inferred communication network
#   C - HOT-MASH inferred communication network
#
# CellChat was run separately for:
#   FAT-MASH (CCl4_model) vs chow
#   HOT-MASH (DEN_model) vs chow
#
# Log-normalized RNA expression was used as input.
# CellChatDB.mouse was used with:
#   computeCommunProb(type = "truncatedMean", trim = 0.1)
#   filterCommunication(min.cells = 10)
#
# Pathway-level communication probabilities were calculated with
# computeCommunProbPathway(), networks were aggregated with aggregateNet(),
# and signaling-role centrality was calculated with
# netAnalysis_computeCentrality(slot.name = "netP").
#
# During the original workflow, the outgoing/incoming values returned by
# get_signaling_role_df() were copied into the final Figure 7A plotting script.
# Those exact values are retained below.
#
# Final figure assembly and minor formatting were performed outside R.

library(Seurat)
library(CellChat)
library(dplyr)
library(tibble)
library(ggplot2)
library(methods)


# =============================================================================
# Shared settings
# =============================================================================

keep_labels <- c(
  "Hepatocytes",
  "Endothelial cells",
  "Macrophages",
  "HSCs",
  "Malignant hepatocytes"
)

cell_colors <- c(
  "Hepatocytes" = "#009E73",
  "Endothelial cells" = "#56B4E9",
  "Macrophages" = "#E69F00",
  "HSCs" = "#F6768E",
  "Malignant hepatocytes" = "#F0E442"
)

CellChatDB <- CellChatDB.mouse


# =============================================================================
# Helper: run CellChat for one model/control object
# =============================================================================

run_cellchat_model <- function(
    object_file
) {

  obj <- readRDS(
    object_file
  )

  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "analysis_label"

  # Retain the cell populations used for CellChat analysis.
  obj <- subset(
    obj,
    subset = analysis_label %in% keep_labels
  )

  obj$analysis_label <- droplevels(
    obj$analysis_label
  )

  # Split pooled cells into chow control and disease-model groups.
  ctrl <- subset(
    obj,
    subset = role == "control"
  )

  model <- subset(
    obj,
    subset = role == "model"
  )

  ctrl$analysis_label <- droplevels(
    ctrl$analysis_label
  )

  model$analysis_label <- droplevels(
    model$analysis_label
  )

  table(ctrl$analysis_label)
  table(model$analysis_label)


  # ---------------------------------------------------------------------------
  # Create CellChat objects
  # ---------------------------------------------------------------------------

  data_ctl <- Seurat::GetAssayData(
    ctrl,
    layer = "data",
    assay = "RNA"
  )

  data_model <- Seurat::GetAssayData(
    model,
    layer = "data",
    assay = "RNA"
  )

  meta_ctl <- data.frame(
    labels = ctrl$analysis_label,
    row.names = colnames(ctrl)
  )

  meta_model <- data.frame(
    labels = model$analysis_label,
    row.names = colnames(model)
  )

  cellchat_ctl <- createCellChat(
    object = data_ctl,
    meta = meta_ctl,
    group.by = "labels"
  )

  cellchat_model <- createCellChat(
    object = data_model,
    meta = meta_model,
    group.by = "labels"
  )


  # ---------------------------------------------------------------------------
  # CellChat database and communication inference
  # ---------------------------------------------------------------------------

  cellchat_ctl@DB <- CellChatDB

  cellchat_ctl <- subsetData(
    cellchat_ctl
  )

  cellchat_model@DB <- CellChatDB

  cellchat_model <- subsetData(
    cellchat_model
  )


  cellchat_ctl <- identifyOverExpressedGenes(
    cellchat_ctl
  )

  cellchat_ctl <- identifyOverExpressedInteractions(
    cellchat_ctl
  )

  cellchat_model <- identifyOverExpressedGenes(
    cellchat_model
  )

  cellchat_model <- identifyOverExpressedInteractions(
    cellchat_model
  )


  cellchat_ctl <- computeCommunProb(
    cellchat_ctl,
    type = "truncatedMean",
    trim = 0.1
  )

  cellchat_model <- computeCommunProb(
    cellchat_model,
    type = "truncatedMean",
    trim = 0.1
  )


  cellchat_ctl <- filterCommunication(
    cellchat_ctl,
    min.cells = 10
  )

  cellchat_model <- filterCommunication(
    cellchat_model,
    min.cells = 10
  )


  # Infer communication at the signaling-pathway level.
  cellchat_ctl <- computeCommunProbPathway(
    cellchat_ctl
  )

  cellchat_model <- computeCommunProbPathway(
    cellchat_model
  )


  # Aggregate the inferred cell-cell communication networks.
  cellchat_ctl <- aggregateNet(
    cellchat_ctl
  )

  cellchat_model <- aggregateNet(
    cellchat_model
  )


  # Calculate outgoing/incoming signaling-role centrality.
  cellchat_ctl <- netAnalysis_computeCentrality(
    cellchat_ctl,
    slot.name = "netP"
  )

  cellchat_model <- netAnalysis_computeCentrality(
    cellchat_model,
    slot.name = "netP"
  )


  list(
    control = cellchat_ctl,
    model = cellchat_model
  )
}


# =============================================================================
# Run CellChat
# =============================================================================

# Inputs generated in 07_create_lineage_objects.R

CCl4_cellchat <- run_cellchat_model(
  "CCl4model_RTchow_refinedAnnotated_RNAonly.rds"
)

DEN_cellchat <- run_cellchat_model(
  "DENmodel_RTchow_refinedAnnotated_RNAonly.rds"
)


CCl4_cellchat_ctl <- CCl4_cellchat$control
CCl4_cellchat_model <- CCl4_cellchat$model

DEN_cellchat_ctl <- DEN_cellchat$control
DEN_cellchat_model <- DEN_cellchat$model


# =============================================================================
# Figure 7A: outgoing and incoming signaling roles
# =============================================================================

# This helper reproduces the values used by
# netAnalysis_signalingRole_scatter(), but returns the underlying data frame
# rather than the CellChat plot.

get_signaling_role_df <- function(
    object,
    slot.name = "netP",
    signaling = NULL,
    x.measure = "outdeg",
    y.measure = "indeg"
) {

  if (
    length(
      methods::slot(
        object,
        slot.name
      )$centr
    ) == 0
  ) {

    stop(
      "Run netAnalysis_computeCentrality() before extracting signaling roles."
    )
  }


  centr <- methods::slot(
    object,
    slot.name
  )$centr


  cell.levels <- levels(
    object@idents
  )


  outgoing <- matrix(
    0,
    nrow = length(cell.levels),
    ncol = length(centr)
  )

  incoming <- matrix(
    0,
    nrow = length(cell.levels),
    ncol = length(centr)
  )


  dimnames(outgoing) <- list(
    cell.levels,
    names(centr)
  )

  dimnames(incoming) <- dimnames(
    outgoing
  )


  for (i in seq_along(centr)) {

    outgoing[, i] <- centr[[i]][[x.measure]]

    incoming[, i] <- centr[[i]][[y.measure]]
  }


  if (!is.null(signaling)) {

    signaling <- signaling[
      signaling %in% object@netP$pathways
    ]

    if (length(signaling) == 0) {

      stop(
        "No significant communication for the requested pathways."
      )
    }


    outgoing <- outgoing[
      ,
      signaling,
      drop = FALSE
    ]

    incoming <- incoming[
      ,
      signaling,
      drop = FALSE
    ]
  }


  outgoing.cells <- rowSums(
    outgoing
  )

  incoming.cells <- rowSums(
    incoming
  )


  agg <- aggregateNet(
    object,
    signaling = signaling,
    return.object = FALSE,
    remove.isolate = FALSE
  )


  num.link.mat <- agg$count

  num.link <- rowSums(
    num.link.mat
  ) +
    colSums(
      num.link.mat
    ) -
    diag(
      num.link.mat
    )


  df <- data.frame(
    celltype = names(
      incoming.cells
    ),
    outgoing = as.numeric(
      outgoing.cells
    ),
    incoming = as.numeric(
      incoming.cells
    ),
    Count = as.numeric(
      num.link
    ),
    stringsAsFactors = FALSE
  )


  df$celltype <- factor(
    df$celltype,
    levels = names(
      incoming.cells
    )
  )


  df
}


# Values underlying the final scatterplot can be reproduced/inspected here.

CCl4_control_roles <- get_signaling_role_df(
  CCl4_cellchat_ctl
)

CCl4_model_roles <- get_signaling_role_df(
  CCl4_cellchat_model
)

DEN_control_roles <- get_signaling_role_df(
  DEN_cellchat_ctl
)

DEN_model_roles <- get_signaling_role_df(
  DEN_cellchat_model
)


CCl4_control_roles
CCl4_model_roles
DEN_control_roles
DEN_model_roles


# -----------------------------------------------------------------------------
# Final Figure 7A plotting data
# -----------------------------------------------------------------------------

# In the original workflow, the values printed above were copied into the
# plotting script. The exact values used for the manuscript plot are retained.

condition <- c(
  rep(
    "Chow",
    4
  ),
  rep(
    "FAT-MASH",
    5
  ),
  rep(
    "WD-DEN",
    5
  )
)


celltype <- c(
  "Hepatocytes",
  "Endothelial Cells",
  "Macrophages",
  "HSCs",

  "Hepatocytes",
  "Endothelial Cells",
  "Macrophages",
  "HSCs",
  "Tumor Cells",

  "Hepatocytes",
  "Endothelial Cells",
  "Macrophages",
  "HSCs",
  "Tumor Cells"
)


Outgoing <- c(
  2.029243,
  1.982788,
  1.095001,
  2.472504,

  0.6947052,
  0.4175429,
  0.0845471,
  1.4143774,
  0.8832953,

  3.017859,
  2.879671,
  1.304002,
  5.426911,
  2.792968
)


Incoming <- c(
  1.467139,
  2.201749,
  2.321097,
  1.589551,

  0.9435673,
  0.1626463,
  0.5091359,
  0.3572382,
  1.5218802,

  3.617480,
  2.746184,
  3.964792,
  2.168402,
  2.924552
)


Count <- c(
  276,
  613,
  447,
  570,

  140,
  138,
  125,
  217,
  197,

  383,
  771,
  561,
  752,
  488
)


PlotData <- data.frame(
  condition,
  celltype,
  Outgoing,
  Incoming,
  Count
)


PlotData <- PlotData %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c(
        "Chow",
        "FAT-MASH",
        "WD-DEN"
      )
    ),
    celltype = factor(
      celltype,
      levels = c(
        "Hepatocytes",
        "Endothelial Cells",
        "Macrophages",
        "HSCs",
        "Tumor Cells"
      )
    )
  )


p_signaling_roles <- ggplot2::ggplot(
  PlotData,
  ggplot2::aes(
    x = Outgoing,
    y = Incoming
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = celltype,
      size = Count,
      shape = condition
    ),
    color = "black",
    alpha = 0.9,
    stroke = 0.8
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Chow" = 21,
      "FAT-MASH" = 23,
      "WD-DEN" = 24
    )
  ) +
  ggplot2::scale_fill_manual(
    name = "Cell Type",
    values = c(
      "Hepatocytes" = "#009E73",
      "Endothelial Cells" = "#56B4E9",
      "Macrophages" = "#E69F00",
      "HSCs" = "#F6768E",
      "Tumor Cells" = "#F0E442"
    )
  ) +
  ggplot2::scale_size_continuous(
    name = "Number of Interactions",
    range = c(
      4,
      12
    )
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        size = 6,
        color = "black"
      )
    )
  ) +
  ggplot2::labs(
    x = "Outgoing Interaction Strength",
    y = "Incoming Interaction Strength",
    shape = "Condition",
    fill = "Cell Type"
  ) +
  ggplot2::theme_classic(
    base_size = 14
  ) +
  ggplot2::theme(
    axis.title =
      ggplot2::element_text(
        face = "bold"
      ),
    legend.title =
      ggplot2::element_text(
        face = "bold"
      ),
    legend.text =
      ggplot2::element_text(
        size = 11
      )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      6
    ),
    breaks = seq(
      0,
      6,
      by = 1
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      4.25
    ),
    breaks = seq(
      0,
      4.5,
      by = 1
    )
  )


p_signaling_roles


# =============================================================================
# Figure 7B: FAT-MASH communication network
# =============================================================================

# Network edge width represents inferred interaction strength.
# Edge color represents the sending/source cell population.

CCl4_model_cols <- cell_colors[
  colnames(
    CCl4_cellchat_model@net$weight
  )
]


netVisual_circle(
  CCl4_cellchat_model@net$weight,
  weight.scale = TRUE,
  label.edge = FALSE,
  edge.width.max = 12,
  title.name = "Interaction strength - model",
  color.use = CCl4_model_cols,
  vertex.label.cex = 0.01
)


# =============================================================================
# Figure 7C: HOT-MASH communication network
# =============================================================================

DEN_model_cols <- cell_colors[
  colnames(
    DEN_cellchat_model@net$weight
  )
]


netVisual_circle(
  DEN_cellchat_model@net$weight,
  weight.scale = TRUE,
  label.edge = FALSE,
  edge.width.max = 12,
  title.name = "Interaction strength - model",
  color.use = DEN_model_cols,
  vertex.label.cex = 0.01
)
