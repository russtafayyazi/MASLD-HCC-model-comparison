# Phenotypic comparisons between FAT-MASH and HOT-MASH
# Figure 1

library(dplyr)

# -----------------------------------------------------------------------------
# Histological scores
# -----------------------------------------------------------------------------

histology <- data.frame(
  condition = c(rep("HOT-MASH", 12), rep("FAT-MASH", 11)),
  steatosis = c(
    3,3,3,3,3,3,3,3,3,3,3,2,
    1,2,2,2,3,2,3,2,2,1,3
  ),
  ballooning = c(
    2,2,2,2,2,2,2,2,2,2,2,1,
    0,1,1,1,1,1,1,1,2,0,1
  ),
  inflammation = c(
    2,2,2,1,2,2,2,2,1,2,2,1,
    3,2,3,2,3,2,2,2,2,2,2
  ),
  NAS = c(
    7,7,7,6,7,7,7,7,6,7,7,4,
    4,5,6,5,7,5,6,5,6,3,6
  )
)

histology$condition <- factor(
  histology$condition,
  levels = c("FAT-MASH", "HOT-MASH")
)

score_vars <- c("steatosis", "ballooning", "inflammation", "NAS")

histology_results <- lapply(score_vars, function(v) {

  test <- wilcox.test(
    as.formula(paste(v, "~ condition")),
    data = histology,
    exact = FALSE
  )

  data.frame(
    score = v,
    FAT_MASH_median = median(
      histology[histology$condition == "FAT-MASH", v]
    ),
    HOT_MASH_median = median(
      histology[histology$condition == "HOT-MASH", v]
    ),
    FAT_MASH_IQR = IQR(
      histology[histology$condition == "FAT-MASH", v]
    ),
    HOT_MASH_IQR = IQR(
      histology[histology$condition == "HOT-MASH", v]
    ),
    W = unname(test$statistic),
    p_value = test$p.value
  )
}) %>%
  bind_rows()

# Benjamini-Hochberg correction across the four histological measures
histology_results$BH_adjusted_p <- p.adjust(
  histology_results$p_value,
  method = "BH"
)

histology_results


# -----------------------------------------------------------------------------
# Total lesion burden
# -----------------------------------------------------------------------------

lesions <- data.frame(
  condition = c(rep("HOT-MASH", 12), rep("FAT-MASH", 11)),
  total_lesions = c(
    15,23,7,7,23,24,12,16,9,5,5,10,
    7,2,9,2,3,2,2,1,5,0,0
  )
)

lesion_test <- wilcox.test(
  total_lesions ~ condition,
  data = lesions,
  exact = FALSE
)

lesion_test
lesion_test$p.value


# -----------------------------------------------------------------------------
# Picrosirius red-positive area
# -----------------------------------------------------------------------------

psr <- data.frame(
  condition = c(rep("HOT-MASH", 10), rep("FAT-MASH", 11)),
  PSR = c(
    1.94767, 2.0622, 0.933657, 0.288883, 1.194013,
    2.84545, 0.5142, 0.59364, 0.158375, 0.742517,

    8.451790833, 5.0203246, 15.2457671, 5.3091864,
    6.1523284, 5.522603806, 6.664539436, 6.777238846,
    4.400279421, 2.925964537, 4.713565642
  )
)

psr_test <- wilcox.test(
  PSR ~ condition,
  data = psr,
  exact = FALSE
)

psr_test
psr_test$p.value


# -----------------------------------------------------------------------------
# Serum ALT
# -----------------------------------------------------------------------------

alt <- data.frame(
  condition = c(rep("HOT-MASH", 12), rep("FAT-MASH", 11)),
  ALT = c(
    10.11510537, 11.23782539, 21.80272711, 25.44005667,
    17.83601913, 5.719851248, 8.117584558, 14.3261909,
    7.285284222, 16.46892155, 15.38161856, 2.670444484,

    11.46229508, 13.03606557, 29.35081967, 7.790163934,
    22.37377049, 11.09508197, 18.70163934, 23.79016393,
    11.98688525, 5.114754098, 19.69836066
  )
)

alt_test <- wilcox.test(
  ALT ~ condition,
  data = alt,
  exact = FALSE
)

alt_test
alt_test$p.value
