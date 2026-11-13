# Type II interferon response at the maternal-fetal interface dictates vaccine efficacy against congenital Zika virus infection

This repository contains the R scripts used to perform the human transcriptomic re-analyses and generate Fig. 3i and Fig. 3j of the manuscript.

## Repository Structure

```text
.
├── README.md
├── 01_human_mfi_interferon_receptors.R   # scRNA-seq receptor expression DotPlot (Vento-Tormo et al., 2018)
└── 02_gse234128_htr8_ifn_pathways.R      # Bulk RNA-seq Type I/II IFN pathway scoring (GSE234128)
```

## System Requirements & Dependencies

All analyses were performed in R v4.4.3 on macOS (tested with `Seurat v5.3.0` and `GEOquery v2.74.0`). No non-standard hardware is required.

### CRAN Packages
```r
install.packages(c(
  "Seurat",
  "data.table",
  "Matrix",
  "ggplot2",
  "ggpubr",
  "dplyr",
  "stringr",
  "msigdbr",
  "BiocManager"
))
```

### Bioconductor Packages
```r
BiocManager::install(c("GEOquery", "org.Hs.eg.db"))
```

## Usage & Figure Mapping

### 1. Interferon Receptor Expression at the Human Maternal–Fetal Interface (Fig. 3i)
- **Script:** `01_human_mfi_interferon_receptors.R`
- **Input Data:** Filtered 10x Genomics scRNA-seq count matrix (`targets_filtered.tsv`) and cell metadata (`meta.tsv`) from the human first-trimester maternal–fetal interface dataset (Vento-Tormo et al., *Nature* 2018; ArrayExpress accession `E-MTAB-6701`). Place both files in the repository root directory before running.
- **Run:**
  ```bash
  Rscript 01_human_mfi_interferon_receptors.R
  ```
- **Outputs:** `Figure_Human_MFI_Receptors_DotPlot.pdf` and `Figure_Human_MFI_Receptors_DotPlot.png`

### 2. Type I and Type II IFN Pathway Activity in ZIKV-Infected HTR8 Trophoblasts (Fig. 3j)
- **Script:** `02_gse234128_htr8_ifn_pathways.R`
- **Input Data:** Automatically downloaded from NCBI GEO (`GSE234128`) and MSigDB (`msigdbr`) at runtime (requires an active internet connection).
- **Run:**
  ```bash
  Rscript 02_gse234128_htr8_ifn_pathways.R
  ```
- **Outputs:** `Figure_GSE234128_HTR8_Pathways.pdf` and `Figure_GSE234128_HTR8_Pathways.png`

## License & Citation

This code is released under the **MIT License**. If you use these scripts, please cite our manuscript and the original public datasets:
- Vento-Tormo, R. et al. Single-cell reconstruction of the early maternal–fetal interface in humans. *Nature* **563**, 347–353 (2018).
- NCBI Gene Expression Omnibus accession **GSE234128**.