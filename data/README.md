`tumorCells.VSTnorm.humanOrths.txt`

VST-normalized pseudobulk expression matrix for tumor-enriched hepatocytes from the six FAT-MASH and HOT-MASH mice used for tumor analyses.

The matrix contains:
- original mouse gene symbols (`gene`)
- mapped human ortholog symbols (`ortholog_name`)
- one VST-normalized expression column per mouse

This matrix was generated during the original analysis and was used as the input for Hoshida subtype scoring in `scripts/04_tumor_programs_hoshida.R`.

Scripts in this repository use paths relative to the repository root unless otherwise noted.
