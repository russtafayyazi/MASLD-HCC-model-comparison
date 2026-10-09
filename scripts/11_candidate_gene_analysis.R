# Candidate gene expression analysis
#
# Supports Figure 6C-F.
#
# Candidate-gene expression was evaluated in hepatocytes, endothelial cells,
# macrophages, and hepatic stellate cells (HSCs).
#
# Log-normalized RNA expression was summarized per mouse from the annotated
# Seurat object. During the original analysis, these per-mouse values were
# normalized to the mean of the FAT-MASH mice for each gene outside the scripted
# workflow (e.g., by copying console output into Excel). The exact normalized
# values used for downstream statistics and plotting are retained explicitly
# below rather than reconstructing an undocumented normalization step.
#
# Statistical comparisons between models used two-sided Welch t-tests.
# Final plotted values are mean +/- SEM across mice.
#
# Minor final figure assembly/annotation was performed outside R.

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# =============================================================================
# Part 1: load annotated object and define refined labels
# =============================================================================

master_integrated <- readRDS(
  "/Users/russtafayyaziv2/Desktop/Steinberg Lab/My Paper - Thesis Plus/revision 2 - finding tumor cells/master integration and annotation/master_integrated_primaryAnnotated.rds"
)

DefaultAssay(master_integrated) <- "RNA"

if ("SCT" %in% Assays(master_integrated)) {
  master_integrated[["SCT"]] <- NULL
}

if ("integrated" %in% Assays(master_integrated)) {
  master_integrated[["integrated"]] <- NULL
}

gc()

DefaultAssay(master_integrated) <- "RNA"
Idents(master_integrated) <- "seurat_clusters"

master_integrated$seurat_clusters <- as.character(
  master_integrated$seurat_clusters
)

master_integrated$primary_label_refined <- as.character(
  master_integrated$primary_label
)

# Cluster 13 is treated as malignant hepatocytes only in the disease models.
master_integrated$primary_label_refined[
  master_integrated$seurat_clusters == "13" &
    master_integrated$condition %in% c(
      "DEN_model",
      "CCl4_model"
    )
] <- "Malignant hepatocytes"

# Chow cluster 13 remains non-malignant hepatocytes.
master_integrated$primary_label_refined[
  master_integrated$seurat_clusters == "13" &
    master_integrated$condition == "RT_chow"
] <- "Hepatocytes"

# Cluster 14 remains cholangiocytes.
master_integrated$primary_label_refined[
  master_integrated$seurat_clusters == "14"
] <- "Cholangiocytes"


# =============================================================================
# Part 2: candidate gene sets
# =============================================================================

hep_candidates <- c(
  "Acly", "Acaca", "Acacb", "Fasn", "Scd1", "Acsl4",
  "Hmgcr", "Sqle", "Soat1", "Srebf2", "Txnrd1", "Slc7a11",
  "Gpx4", "Lcn2", "Glul", "Spink1", "Ephb2"
)

endo_candidates <- c(
  "Plvap", "Kdr", "Angpt2", "Tek", "Esm1", "Clec14a",
  "Flt1", "Fabp4", "Abca1", "Vwf", "Hmox1"
)

macro_candidates <- c(
  "Trem2", "Gpnmb", "Cd9", "Cd63", "Fabp5", "Lpl",
  "Spp1", "Lgals3", "Mertk", "Ccl2", "Ccr2", "Csf1r",
  "C5ar1", "Hmox1", "Folr2", "Cd163", "Mrc1"
)

hsc_candidates <- c(
  "Col1a1", "Fn1", "Serpinh1", "Loxl2", "Postn", "Timp1",
  "Mmp2", "Mmp14", "Tgfb1", "Ccn2", "Serpine1", "Pdgfrb",
  "Itga11", "Fap", "Acta2", "Gdf15", "Ccne1"
)


# =============================================================================
# Part 3: extract per-mouse log-normalized expression
# =============================================================================

extract_mouse_expression <- function(
    object,
    lineage_label,
    candidate_genes,
    mouse_col = "orig.ident"
) {
  
  lineage_obj <- subset(
    object,
    subset = primary_label_refined == lineage_label
  )
  
  DefaultAssay(lineage_obj) <- "RNA"
  
  rna_genes <- rownames(
    lineage_obj[["RNA"]]
  )
  
  gene_check <- data.frame(
    gene = candidate_genes,
    present = candidate_genes %in% rna_genes
  )
  
  print(gene_check)
  
  genes_present <- gene_check$gene[
    gene_check$present
  ]
  
  genes_missing <- gene_check$gene[
    !gene_check$present
  ]
  
  if (length(genes_missing) > 0) {
    message(
      "Missing genes: ",
      paste(
        genes_missing,
        collapse = ", "
      )
    )
  }
  
  expr_mat <- GetAssayData(
    lineage_obj,
    assay = "RNA",
    layer = "data"
  )[
    genes_present,
    ,
    drop = FALSE
  ]
  
  expr_mat <- as.matrix(
    expr_mat
  )
  
  meta <- lineage_obj@meta.data
  
  mean_expr_mouse <- lapply(
    unique(meta[[mouse_col]]),
    function(m) {
      
      cells <- rownames(meta)[
        meta[[mouse_col]] == m
      ]
      
      rowMeans(
        expr_mat[
          ,
          cells,
          drop = FALSE
        ]
      )
    }
  )
  
  mean_expr_mouse <- do.call(
    cbind,
    mean_expr_mouse
  )
  
  colnames(mean_expr_mouse) <- unique(
    meta[[mouse_col]]
  )
  
  mean_expr_mouse
}


hep_mean_expr_mouse <- extract_mouse_expression(
  master_integrated,
  lineage_label = "Hepatocytes",
  candidate_genes = hep_candidates
)

endo_mean_expr_mouse <- extract_mouse_expression(
  master_integrated,
  lineage_label = "Endothelial cells",
  candidate_genes = endo_candidates
)

macro_mean_expr_mouse <- extract_mouse_expression(
  master_integrated,
  lineage_label = "Macrophages",
  candidate_genes = macro_candidates
)

hsc_mean_expr_mouse <- extract_mouse_expression(
  master_integrated,
  lineage_label = "HSCs",
  candidate_genes = hsc_candidates
)


# Per-mouse values inspected during the original analysis.
hep_mean_expr_mouse
endo_mean_expr_mouse
macro_mean_expr_mouse
hsc_mean_expr_mouse


# =============================================================================
# Part 4: normalized per-mouse expression used in the final analysis
# =============================================================================

# The tables below contain the exact values used in the original statistical
# analysis and Figure 6 candidate-gene plots.
#
# Each gene was normalized so that the mean expression across the three
# FAT-MASH mice is 1.


# -----------------------------------------------------------------------------
# Hepatocytes
# -----------------------------------------------------------------------------

hep_rel_expr <- tibble::tribble(
  ~gene, ~CCl4_45, ~CCl4_66, ~CCl4_67, ~DEN_451, ~DEN_487, ~DEN_508,
  "Acly",    0.986258592, 0.852825805, 1.160915603, 1.363389570, 1.084359059, 1.799101472,
  "Acaca",   0.903141505, 0.811012073, 1.285846422, 1.160793934, 1.108717515, 1.328003852,
  "Acacb",   1.065572344, 0.593239212, 1.341188444, 2.609629403, 2.082786510, 2.885381440,
  "Fasn",    0.991709995, 0.797125747, 1.211164258, 1.723302694, 1.202033991, 1.619432816,
  "Scd1",    0.932816298, 1.126363278, 0.940820423, 1.264701212, 1.499434078, 1.367988010,
  "Acsl4",   1.052492264, 0.628634643, 1.318873093, 0.552768308, 0.431063268, 0.415771746,
  "Hmgcr",   1.044052595, 0.742952190, 1.212995215, 0.953839468, 0.726376075, 1.325170405,
  "Sqle",    1.479445564, 0.322714992, 1.197839444, 1.993989171, 0.445543202, 1.399190839,
  "Soat1",   0.810046171, 1.146245556, 1.043708273, 0.417137155, 1.316212601, 0.410920617,
  "Srebf2",  0.997594679, 0.890268133, 1.112137188, 1.168167207, 1.141887820, 1.338733236,
  "Txnrd1",  0.955069897, 0.924655923, 1.120274180, 0.993212682, 1.163394156, 1.149508468,
  "Slc7a11", 0.602861769, 0.825832559, 1.571305672, 0.049258434, 0.395604816, 0.265967711,
  "Lcn2",    0.659923214, 1.071579038, 1.268497749, 0.556646689, 7.508532499, 6.705716758,
  "Glul",    1.233650219, 0.623890712, 1.142459069, 1.896140510, 2.479892339, 2.205648451,
  "Spink1",  0.114554357, 1.601613519, 1.283832124, 0.935673450, 0.366502293, 0.232731177,
  "Ephb2",   1.559802993, 0.574446524, 0.865750483, 1.239457255, 0.122402472, 0.140534905
)


# -----------------------------------------------------------------------------
# Endothelial cells
# -----------------------------------------------------------------------------

endo_rel_expr <- tibble::tribble(
  ~gene, ~CCl4_45, ~CCl4_66, ~CCl4_67, ~DEN_451, ~DEN_487, ~DEN_508,
  "Plvap",   0.406693800, 1.504592974, 1.088713226, 1.572043475, 3.701039638, 1.781693906,
  "Kdr",     1.019189739, 0.927730221, 1.053080041, 1.952095291, 0.978353601, 1.115693876,
  "Angpt2",  1.457567764, 0.158124659, 1.384307577, 3.563676021, 6.628543832, 1.564776639,
  "Tek",     0.874532429, 0.880845047, 1.244622525, 1.504853039, 0.756616401, 0.671279817,
  "Esm1",    1.041518964, 1.084473809, 0.874007228, 1.662236810, 1.564546554, 0.699256810,
  "Clec14a", 0.738873563, 0.998892909, 1.262233529, 1.481326397, 1.927794612, 1.120061639,
  "Flt1",    1.057946206, 0.968129206, 0.973924588, 1.524485888, 1.313548612, 1.156966477,
  "Fabp4",   0.981984640, 1.118763520, 0.899251840, 0.922602733, 1.562268915, 0.499079860,
  "Abca1",   0.627674292, 1.370812585, 1.001513123, 1.383219445, 2.944660348, 2.081435094,
  "Vwf",     1.257092137, 1.385520003, 0.357387860, 0.695264525, 0.798278965, 0.328750522,
  "Hmox1",   0.478667507, 1.662834510, 0.858497983, 0.591818224, 13.73808609, 2.303276688
)


# -----------------------------------------------------------------------------
# Macrophages
# -----------------------------------------------------------------------------

macro_rel_expr <- tibble::tribble(
  ~gene, ~CCl4_45, ~CCl4_66, ~CCl4_67, ~DEN_451, ~DEN_487, ~DEN_508,
  "Trem2",  1.209721138, 0.551771268, 1.238507594, 0.940287077, 0.281911343, 0.245528401,
  "Gpnmb",  1.076491779, 0.902996615, 1.020511605, 0.306959859, 0.157287662, 0.049540364,
  "Cd9",    0.690395236, 1.285331465, 1.024273299, 0.957421965, 1.840922339, 1.149499783,
  "Cd63",   0.873792138, 0.895327685, 1.230880178, 0.680545708, 0.373526696, 0.307918797,
  "Fabp5",  1.168897113, 0.837016136, 0.994086752, 0.514211970, 1.893332340, 0.757153825,
  "Spp1",   0.877119457, 1.213642712, 0.909237831, 0.358891010, 0.390623221, 0.425259448,
  "Lgals3", 0.957718510, 0.957130359, 1.085151132, 0.924524605, 0.724894968, 0.515193976,
  "Mertk",  0.959779212, 0.990455636, 1.049765152, 1.638229872, 1.322769441, 1.176154342,
  "Ccl2",   1.124587385, 1.276451149, 0.598961465, 1.052433448, 0.914131732, 0.440651049,
  "Ccr2",   0.810848826, 0.853444735, 1.335706439, 2.197532608, 0.382986388, 0.675775091,
  "Csf1r",  0.857064483, 0.933255861, 1.209679656, 1.740141297, 1.234676724, 0.917520526,
  "C5ar1",  0.818822100, 0.973659648, 1.207518252, 1.205116725, 0.407994134, 0.479432085,
  "Folr2",  0.662420093, 1.429971982, 0.907607926, 2.274084605, 6.666950134, 3.355055027,
  "Mrc1",   1.149685422, 0.961381788, 0.888932790, 1.688733406, 1.650221216, 1.431947298,
  "Cd163",  0, 0, 0, 0, 0, 0
)


# -----------------------------------------------------------------------------
# HSCs
# -----------------------------------------------------------------------------

hsc_rel_expr <- tibble::tribble(
  ~gene, ~CCl4_45, ~CCl4_66, ~CCl4_67, ~DEN_451, ~DEN_487, ~DEN_508,
  "Col1a1",   0.912557308, 1.110346400, 0.977096293, 0.585369966, 0.198788213, 0.256646235,
  "Fn1",      0.610184091, 1.213970203, 1.175845706, 1.265969472, 2.177617153, 2.129231291,
  "Serpinh1", 0.721096404, 1.266032131, 1.012871465, 0.846919977, 0.887410306, 0.843817058,
  "Loxl2",    0.746621016, 1.359305182, 0.894073801, 1.993289251, 1.327402209, 0.275887092,
  "Postn",    0.893531720, 0.728031393, 1.378436887, 2.510267531, 2.505123871, 1.143637780,
  "Timp1",    0.436945442, 0.636228929, 1.926825629, 0.615985098, 0.160649887, 0.189620242,
  "Mmp2",     1.145471957, 0.940757473, 0.913770570, 0.581490565, 0.072896875, 0.235323479,
  "Tgfb1",    0.615662175, 1.333822403, 1.050515421, 1.670081243, 1.858656246, 1.831780818,
  "Ccn2",     1.159090274, 0.841823173, 0.999086552, 0.750501702, 0.314013512, 0.046461070,
  "Serpine1", 1.247543716, 0.463980368, 1.288475916, 0.956122022, 0.235232651, 0.255086422,
  "Pdgfrb",   1.062398283, 1.133299865, 0.804301852, 1.548528617, 0.482314985, 0.603679252,
  "Itga11",   1.112777392, 0.644891116, 1.242331492, 0,           0,           0,
  "Fap",      2.377298495, 0,           0.622701505, 0,           0.744892799, 0,
  "Acta2",    1.164861613, 1.166285596, 0.668852791, 0.769237780, 0,           0.107011504,
  "Gdf15",    0.618905170, 1.328162687, 1.052932143, 0.922948254, 1.386448447, 2.205953226
)


# =============================================================================
# Part 5: Welch t-tests
# =============================================================================

run_candidate_stats <- function(
    expr_table
) {
  
  long_df <- expr_table %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "mouse_id",
      values_to = "relative_expr"
    ) %>%
    dplyr::mutate(
      condition = dplyr::case_when(
        grepl(
          "^DEN",
          mouse_id
        ) ~ "WD-DEN",
        grepl(
          "^CCl4",
          mouse_id
        ) ~ "FAT-MASH",
        TRUE ~ NA_character_
      ),
      condition = factor(
        condition,
        levels = c(
          "WD-DEN",
          "FAT-MASH"
        )
      )
    )
  
  stats <- long_df %>%
    dplyr::group_by(gene) %>%
    dplyr::group_modify(
      ~{
        
        dat <- .x
        
        den_vals <- dat$relative_expr[
          dat$condition == "WD-DEN"
        ]
        
        fat_vals <- dat$relative_expr[
          dat$condition == "FAT-MASH"
        ]
        
        # Cd163 contained zero expression in all six macrophage samples.
        if (
          length(unique(dat$relative_expr)) == 1
        ) {
          
          return(
            tibble::tibble(
              mean_WD_DEN = mean(den_vals),
              mean_FAT_MASH = mean(fat_vals),
              log2FC_WD_DEN_vs_FAT_MASH = 0,
              welch_p = NA_real_
            )
          )
        }
        
        welch <- t.test(
          relative_expr ~ condition,
          data = dat,
          var.equal = FALSE
        )
        
        tibble::tibble(
          mean_WD_DEN = mean(
            den_vals
          ),
          mean_FAT_MASH = mean(
            fat_vals
          ),
          log2FC_WD_DEN_vs_FAT_MASH =
            log2(
              (
                mean(den_vals) + 0.01
              ) /
                (
                  mean(fat_vals) + 0.01
                )
            ),
          welch_p = welch$p.value
        )
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        log2FC_WD_DEN_vs_FAT_MASH > 0 ~
          "WD-DEN higher",
        log2FC_WD_DEN_vs_FAT_MASH < 0 ~
          "FAT-MASH higher",
        TRUE ~
          "No difference"
      )
    ) %>%
    dplyr::arrange(
      welch_p
    )
  
  list(
    long = long_df,
    stats = stats
  )
}


hep_results <- run_candidate_stats(
  hep_rel_expr
)

endo_results <- run_candidate_stats(
  endo_rel_expr
)

macro_results <- run_candidate_stats(
  macro_rel_expr
)

hsc_results <- run_candidate_stats(
  hsc_rel_expr
)


hep_results$stats
endo_results$stats
macro_results$stats
hsc_results$stats


# =============================================================================
# Part 6: plotting helper
# =============================================================================

plot_candidate_genes <- function(
    expr_table,
    stats_table
) {
  
  gene_order <- expr_table$gene
  
  plot_long <- expr_table %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "mouse_id",
      values_to = "relative_expr"
    ) %>%
    dplyr::mutate(
      condition = dplyr::case_when(
        grepl(
          "^CCl4",
          mouse_id
        ) ~ "CCl4_model",
        grepl(
          "^DEN",
          mouse_id
        ) ~ "DEN_model"
      ),
      condition = factor(
        condition,
        levels = c(
          "CCl4_model",
          "DEN_model"
        )
      ),
      gene = factor(
        gene,
        levels = gene_order
      )
    )
  
  plot_df <- plot_long %>%
    dplyr::group_by(
      gene,
      condition
    ) %>%
    dplyr::summarise(
      mean_expr = mean(
        relative_expr
      ),
      sem_expr =
        sd(relative_expr) /
        sqrt(
          dplyr::n()
        ),
      .groups = "drop"
    )
  
  # Statistical values were calculated in R but significance annotations
  # were added during final figure assembly rather than by geom_text().
  pval_df <- stats_table %>%
    dplyr::mutate(
      gene = factor(
        gene,
        levels = gene_order
      ),
      p_label = dplyr::case_when(
        is.na(welch_p) ~ "",
        welch_p < 0.001 ~ "<0.001",
        welch_p < 0.05 ~ sprintf(
          "%.3f",
          welch_p
        ),
        welch_p < 0.2 ~ sprintf(
          "%.3f",
          welch_p
        ),
        TRUE ~ ""
      )
    )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = gene,
      y = mean_expr,
      fill = condition
    )
  ) +
    ggplot2::geom_col(
      position =
        ggplot2::position_dodge(
          width = 0.72
        ),
      width = 0.65,
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_expr - sem_expr,
        ymax = mean_expr + sem_expr
      ),
      position =
        ggplot2::position_dodge(
          width = 0.72
        ),
      width = 0.22,
      linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.35,
      color = "gray35"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "CCl4_model" = "#9EC5DC",
        "DEN_model" = "#EEC464"
      ),
      labels = c(
        "FAT-MASH",
        "WD-DEN"
      )
    ) +
    ggplot2::scale_y_continuous(
      expand =
        ggplot2::expansion(
          mult = c(
            0,
            0.15
          )
        )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Relative expression (normalized to FAT-MASH)",
      fill = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 12
    ) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          face = "bold",
          color = "black",
          size = 12
        ),
      axis.text.y =
        ggplot2::element_text(
          color = "black"
        ),
      axis.title.y =
        ggplot2::element_text(
          face = "bold"
        ),
      axis.line =
        ggplot2::element_line(
          linewidth = 0.5
        ),
      axis.ticks =
        ggplot2::element_line(
          linewidth = 0.5
        ),
      plot.margin =
        ggplot2::margin(
          10,
          10,
          45,
          10
        ),
      axis.line.y.right =
        ggplot2::element_blank(),
      axis.ticks.y.right =
        ggplot2::element_blank(),
      axis.text.y.right =
        ggplot2::element_blank()
    )
  
  list(
    plot = p,
    plot_data = plot_df,
    p_values = pval_df
  )
}


# =============================================================================
# Part 7: Figure 6C-F plots
# =============================================================================

hep_plot <- plot_candidate_genes(
  hep_rel_expr,
  hep_results$stats
)

endo_plot <- plot_candidate_genes(
  endo_rel_expr,
  endo_results$stats
)

macro_plot <- plot_candidate_genes(
  macro_rel_expr,
  macro_results$stats
)

hsc_plot <- plot_candidate_genes(
  hsc_rel_expr,
  hsc_results$stats
)


p_hep <- hep_plot$plot
p_endo <- endo_plot$plot
p_macro <- macro_plot$plot
p_hsc <- hsc_plot$plot


p_hep
p_endo
p_macro
p_hsc
