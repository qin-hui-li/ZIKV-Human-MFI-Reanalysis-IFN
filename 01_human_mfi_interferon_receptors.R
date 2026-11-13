# ==============================================================================
# Title:       Re-analysis of Interferon Receptor Expression at the Human
#              Maternal–Fetal Interface (Vento-Tormo et al., 2018)
# Description: Loads filtered 10x scRNA-seq count matrix and cell metadata,
#              maps annotated cell types to major anatomical compartments
#              (Placental Villi, Extravillous Trophoblasts, Decidua Basalis),
#              performs log-normalization, and generates a publication DotPlot
#              for Type I and Type II interferon pathway components.
# Inputs:      - targets_filtered.tsv (gene-by-cell count matrix)
#              - meta.tsv (cell-level metadata table)
# Outputs:     - Figure_Human_MFI_Receptors_DotPlot.pdf
#              - Figure_Human_MFI_Receptors_DotPlot.png
# ==============================================================================

# --- Load Required Packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

# --- Configuration & Paths ----------------------------------------------------
# Set `base_dir` to the directory containing `targets_filtered.tsv` and `meta.tsv`.
# Defaults to the current working directory for portability across systems.
base_dir <- "."

counts_file <- file.path(base_dir, "targets_filtered.tsv")
meta_file   <- file.path(base_dir, "meta.tsv")

if (!file.exists(counts_file)) stop("Counts matrix not found: ", counts_file)
if (!file.exists(meta_file))   stop("Metadata file not found: ", meta_file)

# ==============================================================================
# 1. Read and Format Filtered Expression Matrix
# ==============================================================================
message("Loading matrix from targets_filtered.tsv...")
expr <- fread(counts_file, header = TRUE, data.table = FALSE)

# Standardize gene symbols (uppercase) and strip residual quotes/whitespace
gene_names    <- toupper(trimws(gsub('"', "", as.character(expr[, 1]))))
cell_barcodes <- trimws(gsub('"', "", colnames(expr)[-1]))

# Convert expression values to a numeric matrix
expr_data <- expr[, -1]
expr_mat  <- matrix(
  as.numeric(as.matrix(expr_data)),
  nrow = length(gene_names),
  ncol = length(cell_barcodes)
)

rownames(expr_mat) <- gene_names
colnames(expr_mat) <- cell_barcodes

# ==============================================================================
# 2. Read Metadata and Align Cell Barcodes
# ==============================================================================
meta <- fread(meta_file, header = TRUE, data.table = FALSE)

cell_id_col    <- if ("cell" %in% colnames(meta)) "cell" else colnames(meta)[1]
rownames(meta) <- trimws(gsub('"', "", as.character(meta[[cell_id_col]])))

# Harmonize barcode delimiters ('.' vs '-') if R altered column headers on import
if (length(intersect(colnames(expr_mat), rownames(meta))) == 0) {
  colnames(expr_mat) <- gsub("\\.", "-", colnames(expr_mat))
}

common_cells <- intersect(colnames(expr_mat), rownames(meta))
if (length(common_cells) == 0) {
  stop("No overlapping cell barcodes found between count matrix and metadata.")
}

expr_mat <- expr_mat[, common_cells, drop = FALSE]
meta     <- meta[common_cells, , drop = FALSE]

# ==============================================================================
# 3. Create Seurat Object and Map Anatomical Compartments
# ==============================================================================
expr_sparse <- as(expr_mat, "dgCMatrix")
seurat_obj  <- CreateSeuratObject(counts = expr_sparse, meta.data = meta)

# Map Vento-Tormo et al. (2018) cell-type annotations to anatomical compartments
vento_to_layer <- c(
  # Placental Villi (fetal villous core, villous trophoblasts, and Hofbauer cells)
  "SCT"  = "Placental Villi",
  "VCT"  = "Placental Villi",
  "HB"   = "Placental Villi",
  "fEPI" = "Placental Villi",
  
  # Invasive Extravillous Trophoblasts
  "EVT"  = "Extravillous Trophoblasts",
  
  # Decidua Basalis (maternal stromal, perivascular, and decidual immune cells)
  "dS1"  = "Decidua Basalis",
  "dS2"  = "Decidua Basalis",
  "dS3"  = "Decidua Basalis",
  "dM1"  = "Decidua Basalis",
  "dM2"  = "Decidua Basalis",
  "dT"   = "Decidua Basalis",
  "dNK1" = "Decidua Basalis",
  "dNK2" = "Decidua Basalis",
  "dNK3" = "Decidua Basalis",
  "PV"   = "Decidua Basalis"
)

# Identify the metadata column containing the reference cell-type labels
anno_col <- NULL
for (col in colnames(seurat_obj@meta.data)) {
  if (any(as.character(seurat_obj@meta.data[[col]]) %in% names(vento_to_layer))) {
    anno_col <- col
    break
  }
}

if (is.null(anno_col)) {
  stop("Could not identify a metadata column matching Vento-Tormo cell annotations.")
}

seurat_obj$UPU_Layer <- unname(vento_to_layer[as.character(seurat_obj@meta.data[[anno_col]])])
seurat_obj$UPU_Layer[is.na(seurat_obj$UPU_Layer)] <- "Other"

# Subset to target anatomical compartments and log-normalize counts
layer_order <- c("Placental Villi", "Extravillous Trophoblasts", "Decidua Basalis")
seurat_sub  <- subset(seurat_obj, subset = UPU_Layer %in% layer_order)
seurat_sub$UPU_Layer <- factor(seurat_sub$UPU_Layer, levels = layer_order)
Idents(seurat_sub)   <- "UPU_Layer"

seurat_sub <- NormalizeData(
  seurat_sub,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ==============================================================================
# 4. Generate and Export Receptor Expression DotPlot
# ==============================================================================
target_genes    <- c("IFNAR1", "IFNAR2", "STAT2", "IFNGR1", "IFNGR2", "STAT1")
available_genes <- intersect(target_genes, rownames(seurat_sub))

message("Available genes for plotting: ", paste(available_genes, collapse = ", "))

p_dotplot <- DotPlot(
  seurat_sub,
  features = available_genes,
  group.by = "UPU_Layer",
  cols     = c("lightgrey", "#990000"),
  dot.scale = 8,
  scale.by = "radius"
) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle  = 45,
      hjust  = 1,
      vjust  = 1,
      face   = "italic",
      color  = "black",
      size   = 12
    ),
    axis.text.y = element_text(
      face  = "bold",
      color = "black",
      size  = 12
    ),
    axis.title = element_blank(),
    plot.title = element_text(
      face  = "bold",
      hjust = 0.5,
      size  = 14
    )
  ) +
  labs(title = "Human Maternal–Fetal Interface (Re-analysis of Vento-Tormo et al., 2018)")

ggsave(
  filename = file.path(base_dir, "Figure_Human_MFI_Receptors_DotPlot.pdf"),
  plot     = p_dotplot,
  width    = 7.5,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = file.path(base_dir, "Figure_Human_MFI_Receptors_DotPlot.png"),
  plot     = p_dotplot,
  width    = 7.5,
  height   = 4.5,
  dpi      = 300
)

# Print package versions for reproducibility verification
sessionInfo()