# Hoshida HCC subtype scoring in tumor-enriched hepatocytes
# Supports Figure 3B-D.
#
# Two complementary analyses are performed:
#   1. Pseudobulk singscore using VST-normalized tumor-cell expression
#      mapped to human ortholog symbols.
#   2. Single-cell UCell scoring using manually mapped mouse orthologs
#      of the Hoshida S1, S2, and S3 signatures.
#
# Internal condition nomenclature retained from the original analysis:
#   CCl4_model = FAT-MASH
#   DEN_model  = HOT-MASH
#   RT_chow    = chow control
#
# Pseudobulk input:
#   tumorCells.VSTnorm.humanOrths.txt

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(singscore)
library(pheatmap)
library(UCell)


# =============================================================================
# Part 1: pseudobulk Hoshida subtype scoring
# =============================================================================

# -----------------------------------------------------------------------------
# Load prepared VST-normalized tumor-cell pseudobulk matrix
# -----------------------------------------------------------------------------

tumorCells.VSTnorm.humanOrths <- read.delim(
  "data/tumorCells.VSTnorm.humanOrths.txt",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

head(tumorCells.VSTnorm.humanOrths)

stopifnot(
  all(c("gene", "ortholog_name") %in%
        colnames(tumorCells.VSTnorm.humanOrths))
)

expr_mat <- tumorCells.VSTnorm.humanOrths %>%
  select(-gene) %>%
  column_to_rownames("ortholog_name") %>%
  as.matrix()

dim(expr_mat)


# -----------------------------------------------------------------------------
# Hoshida human subtype signatures
# -----------------------------------------------------------------------------

S1_genes <- c(
  "ACP5","ACTA2","ADAM15","ADAM8","ADAM9","AEBP1","AIF1","AKT3",
  "ALDOA","ALOX5AP","ANXA1","ANXA4","ANXA5","AP2S1","AQP1","ARF5",
  "ARHGDIB","ARPC1B","ARPC2","ASAH1","ATP1B3","ATP6AP1","ATP6V0B",
  "ATP6V1B2","ATP6V1F","BCL2A1","BLVRA","C1QB","C3AR1","CAPZA1",
  "CBFB","CCL3","CCL5","CCN1","CCN2","CCND2","CCR7","CD151","CD37",
  "CD3D","CD47","CD48","CD53","CD74","CD8A","CDC20","CDC25B","CDH11",
  "CDK2AP1","CELF2","CHN1","COL11A1","COL15A1","COL1A2","COL3A1",
  "COL4A1","COL4A2","COL5A2","COL6A1","COL6A2","CORO1A","CRABP2",
  "CRIP2","CSRP3","CTSC","CTSS","CXCL1","CXCR4","CYBA","CYBB","CYFIP1",
  "CYP1B1","DAB2","DCTN2","DDR1","DDR2","DDX11","DGKA","DGKZ","DNM2",
  "DPYSL2","DUSP5","DUT","EFEMP1","EFNB1","F13A1","FBN1","FBRS","FCGBP",
  "FCGR2A","FGL2","FHL3","FLNA","FUT4","FYB1","GEM","GLIPR1","GNAI2",
  "GNS","GPNMB","GRN","GSTP1","GUCY1A1","GYPC","HCLS1","HEXA","HIF1A",
  "HK1","HLA-DMA","HLA-DOA","HLA-DPB1","HLA-DQA2","HLA-DQB1","ID3",
  "IER3","IFI16","IFI30","IGFBP5","IGKC","IGLL1","IKBKE","IL15RA",
  "IL2RB","IL2RG","IL7R","IQGAP1","IRF1","ITGB2","ITPR3","KLC1","KLF5",
  "LAMB1","LAPTM5","LCP1","LDHB","LGALS1","LGALS3BP","LGALS9","LGMN",
  "LHFPL2","LITAF","LMO4","LOX","LSP1","LTBP2","LTBP3","LTF","LUM",
  "LYN","M6PR","MAP1B","ME1","MFAP1","MGP","MPHOSPH6","MSN","MTHFD2",
  "MYCBP2","NBL1","NPC2","NSMAF","OAZ1","PAK1","PAM","PAPSS1","PCLAF",
  "PEA15","PFN1","PGK1","PIM2","PKD2","PKMYT1","PKN1","PLAUR","PLD3",
  "PLP2","PNP","POLD3","POSTN","PPIC","PPP1CB","PPP4C","PPP4R1","PRKD2",
  "PRMT5","PROCR","PSMD2","PTPRC","PYGB","QSOX1","RAB31","RALB","RALGDS",
  "RCC1","RHOA","RIN2","RIT1","RNASE1","RNASE6","RPA2","RSU1","S100A10",
  "S100A11","S100A13","SLA","SLBP","SLC1A5","SLC2A1","SLC2A5","SLC39A6",
  "SLC7A5","SMAD2","SMARCD1","SPAG8","SRGN","SRI","STK38","STX3","TAGLN",
  "TAX1BP3","TCF4","TFF3","THY1","TIMP2","TMSB4X","TNFRSF1B","TP53BP1",
  "TPM2","TRAF3","TRAF5","TRIP10","TSPAN3","TUBA4A","VCAN","ZNF384"
)

S2_genes <- c(
  "ABCB10","ABCD3","ACP1","ADD3","AFP","AHCY","ARHGAP35","ARID3A",
  "ATF2","ATM","ATP2B1","ATP2B2","ATP5PB","ATXN10","BCAM","BCLAF1",
  "BRD3","BTG3","CASC3","CD46","CDK6","CHKA","CLK2","COL2A1","CPD",
  "CSE1L","CSNK2A1","CSNK2A2","CTNNB1","CUL4A","CXADR","DDX1","DDX18",
  "DEK","EIF4A2","EIF4B","ENPP1","EP300","ERBB3","FBL","FGFR3","FGFR4",
  "FLNB","GBF1","GCN1","GLUD1","GNAI1","GPC3","GTF2I","GTF3C2","H1-0",
  "HELZ","HMGCR","HNRNPA2B1","HNRNPC","HNRNPU","IDI1","IGF2","IGF2R",
  "ITIH2","KLF3","LBR","MAPK6","MEST","NCOA4","NET1","NR2C1","NR5A2",
  "NREP","NT5E","NUP153","PEG3","PHF3","PHKA2","PIGC","PLXNB1","PNN",
  "POFUT1","PPARG","PPP2R1A","PRDX3","PTOV1","RAB4A","RBM39","RPL24",
  "RPL27","RPL31","RPS19","RPS24","RPS25","RPS27","RPS5","RRP1B",
  "SEPHS1","SLC6A2","SLC6A5","SMARCA1","SMARCC1","SNRPE","SNTB1",
  "SREBF2","SSB","SUMO1","SUZ12","TARBP1","TBCE","TFIP11","TIA1",
  "TIAL1","TM9SF4","TP53BP2","TPR","TRIM26","TTC3","UBE2K"
)

S3_genes <- c(
  "ABCB4","ABCC6","ABHD2","ACAA2","ACADM","ACADS","ACADSB","ACADVL",
  "ACO1","ACOX1","ACOX2","ACSL1","ACY1","ADA2","ADH4","ADH6","ADK",
  "AGL","AGXT","AKR1C1","ALAS1","ALDH1A1","ALDH1B1","ALDH2","ALDH3A2",
  "ALDH4A1","ALDH6A1","ALDH7A1","ALDOB","ALPL","AMFR","AMT","ANXA6",
  "AOC1","APCS","APOA1","APOC2","APOC4","APOH","AQP7","ARG1","ARHGEF12",
  "ARSA","ASCL1","ASGR1","ASGR2","ASL","ASS1","ATOX1","ATP5F1D","ATP5PF",
  "AZGP1","BAAT","BDH1","BHMT","BLOC1S1","BLVRB","BPHL","BTD","C1R",
  "C1S","C4A","C4BPA","C8B","CA2","CAT","CBR1","CD14","CD302","CD81",
  "CES1","CFB","CFH","CGREF1","CNGA1","COL18A1","COX5B","CP","CPA3",
  "CPA4","CPB2","CPS1","CRABP1","CRYAA","CRYM","CSTB","CTH","CTSO",
  "CXCL2","CYB5A","CYFIP2","CYP21A2","CYP27A1","CYP2C9","CYP2J2",
  "CYP3A7","DAO","DCAF8","DECR1","DNASE1L3","DPAGT1","DRG2","ECHS1",
  "ECI1","EDNRB","EGFR","EHHADH","EMP2","EPAS1","EPHX1","ETS2","F11",
  "F2","F5","FAH","FANCA","FGB","FGG","FH","FKBP2","FLT4","FMO4","FOXO1",
  "FXR2","GCH1","GCHFR","GCKR","GGH","GHR","GJB1","GLYAT","GOT2","GPT",
  "GPX2","GPX3","GSTA2","GSTO1","GSTZ1","HAAO","HADH","HGD","HMGCS2",
  "HMOX2","HPD","HRG","HSD17B10","HSD17B4","ICAM3","IDH2","IDH3A",
  "IFIT1","IGF1","IL13RA1","IL32","IL6R","IMPA1","INSR","IQGAP2","ISG15",
  "ITIH1","ITIH3","ITIH4","ITPR2","IVD","KCNJ8","KLKB1","KMO","KNG1",
  "LCAT","LONP1","LPIN1","LPIN2","MAOA","MAOB","MAPRE3","MGST2","MME",
  "MMUT","MSMO1","MT2A","MTHFD1","MTHFS","MYLK","MYO1E","NDUFV2","NFIB",
  "NFIC","NFKBIA","NHERF2","NNMT","NRG1","PAH","PAPSS2","PCCA","PCCB",
  "PCK1","PCK2","PDK4","PGM1","PGRMC1","PIK3R1","PKLR","PLA2G2A","PLCG2",
  "PLG","PLGLB2","PNPLA4","POLD4","PON3","PPP2R1B","PROS1","PTGR1","PTS",
  "QDPR","RARRES2","RBP5","RGN","RHOB","RIDA","RNASE4","SBDS","SDC1",
  "SDHB","SDS","SELENBP1","SELENOP","SERPINA3","SERPINA6","SERPINC1",
  "SERPING1","SHB","SHMT1","SLC10A1","SLC16A2","SLC23A1","SLC23A2",
  "SLC2A2","SLC35D1","SLC6A1","SLC6A12","SLC7A2","SLCO2A1","SLPI",
  "SMARCA2","SOAT1","SOD1","SOD2","SORL1","SPAM1","SPARCL1","SRD5A1",
  "SREBF1","SULT2A1","TCEA2","TDO2","TGFBR3","TINAGL1","TJP2","TMBIM6",
  "TMOD1","TOB1","TPMT","TST","UQCRB","VSIG2","ZNF160"
)


# -----------------------------------------------------------------------------
# Pseudobulk singscore
# -----------------------------------------------------------------------------

# Report representation of each subtype signature in the input matrix
length(S1_genes)
length(S2_genes)
length(S3_genes)

length(intersect(S1_genes, rownames(expr_mat)))
length(intersect(S2_genes, rownames(expr_mat)))
length(intersect(S3_genes, rownames(expr_mat)))

rankData <- rankGenes(expr_mat)

S1_score <- simpleScore(
  rankData,
  upSet = S1_genes
)

S2_score <- simpleScore(
  rankData,
  upSet = S2_genes
)

S3_score <- simpleScore(
  rankData,
  upSet = S3_genes
)

pseudobulk_scores <- data.frame(
  sample = colnames(expr_mat),
  S1 = S1_score$TotalScore,
  S2 = S2_score$TotalScore,
  S3 = S3_score$TotalScore
)

# Recover model identity from the original sample names
pseudobulk_scores$condition <- case_when(
  grepl("^CCl4_model_", pseudobulk_scores$sample) ~ "CCl4_model",
  grepl("^DEN_model_", pseudobulk_scores$sample) ~ "DEN_model",
  TRUE ~ NA_character_
)

pseudobulk_scores$condition <- factor(
  pseudobulk_scores$condition,
  levels = c("CCl4_model", "DEN_model")
)

pseudobulk_scores


# Figure 3B: pseudobulk Hoshida subtype heatmap
score_mat <- as.matrix(
  pseudobulk_scores[, c("S1", "S2", "S3")]
)

rownames(score_mat) <- pseudobulk_scores$sample

pheatmap(
  score_mat,
  scale = "none",
  cluster_cols = FALSE,
  treeheight_row = 15
)


# =============================================================================
# Part 2: single-cell Hoshida UCell scoring
# =============================================================================

# -----------------------------------------------------------------------------
# Load primary annotated atlas
# -----------------------------------------------------------------------------

master_integrated <- readRDS(
  "master_integrated_primaryAnnotated.rds"
)

DefaultAssay(master_integrated) <- "RNA"

# Remove assays not required for UCell analysis to reduce memory usage
if ("SCT" %in% Assays(master_integrated)) {
  master_integrated[["SCT"]] <- NULL
}

if ("integrated" %in% Assays(master_integrated)) {
  master_integrated[["integrated"]] <- NULL
}

gc()

Idents(master_integrated) <- "seurat_clusters"


# -----------------------------------------------------------------------------
# Mouse orthologs of Hoshida subtype signatures
# -----------------------------------------------------------------------------
#
# These vectors preserve the mouse ortholog mappings used in the original
# analysis. "N/A" entries represent human signature genes for which no mouse
# feature was included in the mapping; they are automatically removed below
# when intersected with the Seurat feature names.

S1_mouse <- c(
  "Acp5","Acta2","Adam15","Adam8","Adam9","Aebp1","Aif1","Akt3",
  "Aldoa","Alox5ap","Anxa1","Anxa4","Anxa5","Ap2s1","Aqp1","Arf5",
  "Arhgdib","Arpc1b","Arpc2","Asah1","Atp1b3","Atp6ap1","Atp6v0b",
  "Atp6v1b2","Atp6v1f","Bcl2a1c","Blvra","C1qb","C3ar1","Capza1b",
  "Cbfb","Ccl3","Ccl5","Ccn1","Ccn2","Ccnd2","Ccr7","Cd151","Cd37",
  "Cd3d","Cd47","Cd48","Cd53","Cd74","Cd8a","Cdc20","Cdc25b","Cdh11",
  "Cdk2ap1","Celf2","Chn1","Col11a1","Col15a1","Col1a2","Col3a1",
  "Col4a1","Col4a2","Col5a2","Col6a1","Col6a2","Coro1a","Crabp2",
  "Crip2","Csrp3","Ctsc","Ctss","N/A","Cxcr4","Cyba","Cybb","Cyfip1",
  "Cyp1b1","Dab2","Dctn2","Ddr1","Ddr2","Ddx11","Dgka","Dgkz","Dnm2",
  "Dpysl2","Dusp5","Dut","Efemp1","Efnb1","F13a1","Fbn1","Fbrs","N/A",
  "Fcgr2b","Fgl2","Fhl3","Flna","Fut4","Fyb1","Gem","Glipr1","Gnai2",
  "Gns","Gpnmb","Grn","Gstp2","Gucy1a1","Gypc","Hcls1","Hexa","Hif1a",
  "Hk1","H2-DMa","H2-Oa","N/A","N/A","H2-Ab1","Id3","Ier3","Ifi202b",
  "Ifi30","Igfbp5","Igkc","Iglc2","Ikbke","Il15ra","Il2rb","Il2rg",
  "Il7r","Iqgap1","Irf1","Itgb2l","Itpr3","Klc1","Klf5","Lamb1",
  "Laptm5","Lcp1","Ldhb","Lgals1","Lgals3bp","N/A","Lgmn","Lhfpl2",
  "Litaf","Lmo4","Lox","Lsp1","Ltbp2","Ltbp3","Ltf","Lum","Lyn","M6pr",
  "Map1b","Me1","Mfap1b","Mgp","Mphosph6","Msn","Mthfd2","Mycbp2",
  "Nbl1","Npc2","Nsmaf","Oaz1","Pak1","Pam","Papss1","Pclaf","Pea15a",
  "Pfn1","Pgk1","Pim2","Pkd2","Pkmyt1","Pkn1","Plaur","Pld3","Plp2",
  "Pnp2","N/A","Postn","Ppic","Ppp1cb","Ppp4c","Ppp4r1","Prkd2","Prmt5",
  "Procr","Psmd2","Ptprc","Pygb","Qsox1","Rab31","Ralb","Ralgds","Rcc1",
  "Rhoa","Rin2","Rit1","Rnase1","Rnase6","Rpa2","Rsu1","S100a10",
  "S100a11","S100a13","Sla","Slbp","Slc1a5","Slc2a1","Slc2a5","Slc39a6",
  "Slc7a5","Smad2","Smarcd1","Spag8","Srgn","Sri","Stk38","Stx3","Tagln",
  "Tax1bp3","Tcf4","Tff3","Thy1","Timp2","Tmsb4x","Tnfrsf1b","Trp53bp1",
  "Tpm2","Traf3","Traf5","Trip10","Tspan3","Tuba4a","Vcan","Zfp384"
)

S2_mouse <- c(
  "Abcb10","Abcd3","Acp1","Add3","Afp","Ahcy","Arhgap35","Arid3a",
  "Atf2","Atm","Atp2b1","Atp2b2","Atp5pb","Atxn10","Bcam","Bclaf1",
  "Brd3","Btg3","Casc3","Cd46","Cdk6","Chka","Clk2","Col2a1","Cpd",
  "Cse1l","Csnk2a1","Csnk2a2","Ctnnb1","Cul4a","Gm1123","Ddx1","Ddx18",
  "Dek","Eif4a2","Eif4b","Enpp1","Ep300","Erbb3","Fbl","Fgfr3","Fgfr4",
  "Flnb","Gbf1","Gcn1","Glud1","Gnai1","Gpc3","Gtf2i","Gtf3c2","H1f0",
  "Helz","Hmgcr","Hnrnpa2b1","Hnrnpc","Hnrnpu","Idi1","Igf2","Igf2r",
  "Itih2","Klf3","Lbr","Mapk6","Mest","Ncoa4","Net1","Nr2c1","Nr5a2",
  "Nrep","Nt5e","Nup153","Peg3","Phf3","Phka2","Pigc","Plxnb1","Pnn",
  "Pofut1","Pparg","Ppp2r1a","Prdx3","Ptov1","Rab4a","Rbm39","Rpl24",
  "Rpl27","N/A","Rps19","N/A","Rps25","Rps27","Rps5","Rrp1b","Sephs1",
  "Slc6a2","Slc6a5","Smarca1","Smarcc1","Snrpe","Sntb1","Srebf2","Ssb",
  "Sumo1","Suz12","Tarbp1","Tbce","Tfip11","Tia1","Tial1","Tm9sf4",
  "Trp53bp2","Tpr","Trim26","Ttc3","Ube2k"
)

S3_mouse <- c(
  "Abcb4","Abcc6","Abhd2","Acaa2","Acadm","Acads","Acadsb","Acadvl",
  "Aco1","Acox1","Acox2","Acsl1","Acy1","N/A","Adh4","Adh6b","Adk",
  "Agl","Agxt","N/A","Alas1","Aldh1a7","Aldh1b1","Aldh2","Aldh3a2",
  "Aldh4a1","Aldh6a1","Aldh7a1","Aldob","Alpl","Amfr","Amt","Anxa6",
  "Aoc1","Apcs","Apoa1","Apoc2","Apoc4","Apoh","Aqp7","Arg1","Arhgef12",
  "Arsa","Ascl1","Asgr1","Asgr2","Asl","Ass1","Atox1","Atp5f1d","Atp5pf",
  "Azgp1","Baat","Bdh1","Bhmt1b","Bloc1s1","Blvrb","Bphl","Btd","C1ra",
  "C1s1","C4a","C4bp","C8b","Car2","Cat","Cbr1","Cd14","Cd302","Cd81",
  "Ces1d","Cfb","N/A","Cgref1","Cnga1","Col18a1","Cox5b","Cp","Cpa3",
  "Cpa4","Cpb2","Cps1","Crabp1","Cryaa","Crym","Cstb","Cth","Ctso","N/A",
  "Cyb5a","Cyfip2","Cyp21a1","Cyp27a1","N/A","Cyp2j9","Cyp3a25","Dao",
  "Dcaf8","Decr1","Dnase1l3","Dpagt1","Drg2","Echs1","Eci1","Ednrb",
  "Egfr","Ehhadh","Emp2","Epas1","Ephx1","Ets2","F11","F2","F5","Fah",
  "Fanca","Fgb","Fgg","Fh1","N/A","Flt4","Fmo4","Foxo1","Fxr2","Gch1",
  "Gchfr","Gckr","Ggh","Ghr","Gjb1","Glyat","Got2","Gpt","Gpx2","Gpx3",
  "N/A","Gsto1","Gstz1","Haao","Hadh","Hgd","Hmgcs2","Hmox2","Hpd","Hrg",
  "Hsd17b10","Hsd17b4","N/A","Idh2","Idh3a","N/A","Igf1","Il13ra1",
  "N/A","Il6ra","Impa1","Insr","Iqgap2","Isg15","Itih1","Itih3","Itih4",
  "Itpr2","Ivd","Kcnj8","Klkb1","Kmo","Kng1","Lcat","Lonp1","Lpin1",
  "Lpin2","Maoa","Maob","Mapre3","Mgst2","Mme","Mmut","Msmo1","N/A",
  "Mthfd1","Mthfs","Mylk","Myo1e","Ndufv2","Nfib","Nfic","Nfkbia",
  "Nherf2","Nnmt","Nrg1","Pah","Papss2","Pcca","Pccb","Pck1","Pck2",
  "Pdk4","Pgm1","Pgrmc1","Pik3r1","Pklr","Pla2g2a","Plcg2","Plg","Plg",
  "N/A","Pold4","Pon3","Ppp2r1b","Pros1","Ptgr1","Pts","Qdpr","Rarres2",
  "N/A","Rgn","Rhob","Rida","Rnase4","Sbds","Sdc1","Sdhb","Sds","Selenbp1",
  "Selenop","Serpina3n","Serpina6","Serpinc1","Serping1","Shb","Shmt1",
  "Slc10a1","Slc16a2","Slc23a1","Slc23a2","Slc2a2","Slc35d1","Slc6a1",
  "Slc6a12","Slc7a2","Slco2a1","Slpi","Smarca2","Soat1","Sod1","Sod2",
  "Sorl1","Hyal5","Sparcl1","Srd5a1","Srebf1","Sult2a8","Tcea2","Tdo2",
  "Tgfbr3","Tinagl1","Tjp2","Tmbim6","Tmod1","Tob1","Tpmt","Tst","Uqcrb",
  "Vsig2","Zfp160"
)


# -----------------------------------------------------------------------------
# Restrict signatures to genes present in the RNA assay
# -----------------------------------------------------------------------------

obj_genes <- unique(
  rownames(master_integrated[["RNA"]])
)

S1_mouse_present <- intersect(
  S1_mouse,
  obj_genes
)

S2_mouse_present <- intersect(
  S2_mouse,
  obj_genes
)

S3_mouse_present <- intersect(
  S3_mouse,
  obj_genes
)

c(
  S1_total = length(S1_mouse),
  S1_present = length(S1_mouse_present),
  S2_total = length(S2_mouse),
  S2_present = length(S2_mouse_present),
  S3_total = length(S3_mouse),
  S3_present = length(S3_mouse_present)
)


# -----------------------------------------------------------------------------
# Calculate Hoshida UCell scores
# -----------------------------------------------------------------------------

hoshida_signatures <- list(
  Hoshida_S1 = S1_mouse_present,
  Hoshida_S2 = S2_mouse_present,
  Hoshida_S3 = S3_mouse_present
)

master_integrated <- AddModuleScore_UCell(
  obj = master_integrated,
  features = hoshida_signatures,
  assay = "RNA",
  slot = "counts",
  name = NULL
)

grep(
  "^Hoshida_",
  colnames(master_integrated@meta.data),
  value = TRUE
)


# -----------------------------------------------------------------------------
# Restrict downstream Hoshida analysis to tumor-enriched cluster 13
# -----------------------------------------------------------------------------

cl13 <- subset(
  master_integrated,
  idents = "13"
)

# Cluster 13 also contains chow cells in the integrated atlas.
# Hoshida model comparisons use only FAT-MASH and HOT-MASH mice.
cl13_models <- subset(
  cl13,
  subset = condition %in% c(
    "CCl4_model",
    "DEN_model"
  )
)

cl13_models$condition <- factor(
  cl13_models$condition,
  levels = c(
    "CCl4_model",
    "DEN_model"
  )
)


# -----------------------------------------------------------------------------
# Figure 3C: single-cell Hoshida subtype scores
# -----------------------------------------------------------------------------

VlnPlot(
  cl13_models,
  features = c(
    "Hoshida_S1",
    "Hoshida_S2",
    "Hoshida_S3"
  ),
  group.by = "condition",
  pt.size = 0.15,
  combine = TRUE
) &
  scale_fill_manual(
    values = c(
      "CCl4_model" = "#9EC5DC",
      "DEN_model" = "#EEC464"
    ),
    labels = c(
      "CCl4_model" = "FAT-MASH",
      "DEN_model" = "HOT-MASH"
    )
  )


# -----------------------------------------------------------------------------
# Figure 3D: mouse-level mean UCell scores
# -----------------------------------------------------------------------------

cl13_mouse_summary <- cl13_models@meta.data %>%
  group_by(
    sample_id,
    condition
  ) %>%
  summarise(
    mean_S1 = mean(
      Hoshida_S1,
      na.rm = TRUE
    ),
    mean_S2 = mean(
      Hoshida_S2,
      na.rm = TRUE
    ),
    mean_S3 = mean(
      Hoshida_S3,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

cl13_mouse_summary

plot_df <- cl13_mouse_summary %>%
  select(
    sample_id,
    condition,
    mean_S1,
    mean_S2,
    mean_S3
  ) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "subtype",
    values_to = "score"
  ) %>%
  mutate(
    subtype = sub(
      "^mean_",
      "",
      subtype
    ),
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
    y = score
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
  facet_wrap(
    ~ subtype,
    scales = "free_y"
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean UCell score"
  )
