# Data

 `tumorCells.VSTnorm.humanOrths.txt`

VST-normalized pseudobulk expression matrix for tumor-enriched hepatocytes from the six FAT-MASH and HOT-MASH mice used in the tumor analyses.

The file contains the original mouse gene symbol, the mapped human ortholog symbol, and one VST-normalized expression column for each mouse.

This matrix was generated during the original analysis and is used for Hoshida subtype scoring in `scripts/04_tumor_programs_hoshida.R` and human MASLD-HCC signature scoring in `scripts/06_human_MASLD_HCC_signature.R`.

Scripts use repository-relative paths where included data files are referenced.
