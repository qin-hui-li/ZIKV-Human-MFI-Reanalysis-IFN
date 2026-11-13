# ==============================================================================
# Title:       Type I and Type II Interferon Pathway Activity in ZIKV-Infected
#              Human HTR8 Trophoblasts (Re-analysis of GSE234128)
# Description: Fetches MSigDB Hallmark interferon gene sets, downloads sample
#              metadata and bulk RNA-seq expression data for GEO accession
#              GSE234128, subsets HTR8 trophoblast samples (Control vs. ZIKV),
#              computes mean log2 pathway expression scores, and exports a
#              two-panel comparison figure with Student's t-test p-values.
# Inputs:      - Downloaded automatically via GEOquery (GSE234128) and msigdbr
# Outputs:     - Figure_GSE234128_HTR8_Pathways.pdf
#              - Figure_GSE234128_HTR8_Pathways.png
# ==============================================================================

# --- Load Required Packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(msigdbr)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(stringr)
})

# Set seed for reproducible horizontal jitter in boxplots
set.seed(42)

# --- Configuration & Paths ----------------------------------------------------
# Defaults to the current working directory for portability across systems.
base_dir <- "."
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Fetch MSigDB Hallmark Interferon Gene Sets
# ==============================================================================
message("1/5: Loading MSigDB Hallmark Gene Sets...")
hs_hallmarks <- msigdbr(species = "Homo sapiens", collection = "H")

gsea_t1 <- hs_hallmarks %>%
  filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene_symbol)

gsea_t2 <- hs_hallmarks %>%
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  pull(gene_symbol)

# ==============================================================================
# 2. Download GSE234128 Metadata (SOFT Format) & Filter HTR8 Samples
# ==============================================================================
message("2/5: Fetching GSE234128 metadata using SOFT format...")
gse_soft <- getGEO("GSE234128", GSEMatrix = FALSE, destdir = base_dir)
gsms     <- GSMList(gse_soft)

pheno_df <- data.frame(
  GSM             = names(gsms),
  Title           = sapply(gsms, function(x) Meta(x)$title),
  Characteristics = sapply(gsms, function(x) paste(Meta(x)$characteristics_ch1, collapse = " | ")),
  stringsAsFactors = FALSE
)

# Subset to HTR8 trophoblast samples and assign experimental condition labels
pheno_htr8 <- pheno_df %>%
  filter(
    str_detect(Title, regex("HTR8|HTR-8|HTR_8|Trophoblast", ignore_case = TRUE)) |
      str_detect(Characteristics, regex("HTR8|HTR-8|HTR_8|Trophoblast", ignore_case = TRUE)) |
      !str_detect(Title, regex("U251", ignore_case = TRUE))
  ) %>%
  mutate(
    Group = case_when(
      str_detect(Title, regex("ZIKV|Zika", ignore_case = TRUE)) |
        str_detect(Characteristics, regex("ZIKV|Zika", ignore_case = TRUE)) ~ "ZIKV",
      TRUE ~ "Control (NF)"
    )
  )

pheno_htr8$Group <- factor(pheno_htr8$Group, levels = c("Control (NF)", "ZIKV"))

cat("\n--- DETECTED HTR8 SAMPLES ---\n")
print(table(pheno_htr8$Group))
cat("-----------------------------\n\n")

# ==============================================================================
# 3. Download & Parse Supplementary Expression Matrix
# ==============================================================================
message("3/5: Downloading supplementary expression matrix...")
supp_files <- getGEOSuppFiles("GSE234128", baseDir = base_dir)
supp_dir   <- file.path(base_dir, "GSE234128")
all_supp   <- list.files(supp_dir, full.names = TRUE, recursive = TRUE)

matrix_file <- all_supp[grep("\\.txt|\\.tsv|\\.csv|\\.gz", all_supp)][1]
if (is.na(matrix_file) || !file.exists(matrix_file)) {
  stop("No supplementary expression matrix found in: ", supp_dir)
}
message("Loading expression matrix: ", basename(matrix_file))

if (grepl("\\.csv$", matrix_file, ignore.case = TRUE)) {
  expr_raw <- read.csv(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  expr_raw <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
}

if (is.character(expr_raw[[1]]) || is.factor(expr_raw[[1]])) {
  gene_ids <- trimws(gsub('"', "", as.character(expr_raw[[1]])))
  num_df   <- expr_raw[, -1, drop = FALSE]
} else {
  gene_ids <- rownames(expr_raw)
  num_df   <- expr_raw
}

expr_mat <- apply(num_df, 2, function(x) suppressWarnings(as.numeric(as.character(x))))
rownames(expr_mat) <- gene_ids
colnames(expr_mat) <- trimws(gsub('"', "", colnames(expr_mat)))
colnames(expr_mat) <- gsub("^X", "", colnames(expr_mat))

# ==============================================================================
# 4. Map Gene Identifiers, Align Samples, and Log2-Transform
# ==============================================================================
if (any(grepl("^ENSG", rownames(expr_mat)))) {
  message("Ensembl IDs detected! Converting to Gene Symbols...")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package 'org.Hs.eg.db' is required for Ensembl ID mapping. Install via BiocManager::install('org.Hs.eg.db').")
  }
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  
  clean_ensembl <- gsub("\\..*", "", rownames(expr_mat))
  symbols <- mapIds(
    org.Hs.eg.db,
    keys      = clean_ensembl,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  valid <- !is.na(symbols)
  expr_mat <- expr_mat[valid, , drop = FALSE]
  rownames(expr_mat) <- toupper(symbols[valid])
} else {
  rownames(expr_mat) <- toupper(rownames(expr_mat))
}

# Match expression matrix columns to filtered HTR8 metadata
common_samples <- intersect(colnames(expr_mat), pheno_htr8$GSM)

if (length(common_samples) > 0) {
  pheno_htr8 <- pheno_htr8 %>% filter(GSM %in% common_samples)
  expr_mat   <- expr_mat[, pheno_htr8$GSM, drop = FALSE]
} else {
  pheno_htr8$MatchCol <- NA
  for (i in seq_len(nrow(pheno_htr8))) {
    hit <- colnames(expr_mat)[grep(pheno_htr8$Title[i], colnames(expr_mat), fixed = TRUE)]
    if (length(hit) > 0) pheno_htr8$MatchCol[i] <- hit[1]
  }
  pheno_htr8 <- pheno_htr8 %>% filter(!is.na(MatchCol))
  expr_mat   <- expr_mat[, pheno_htr8$MatchCol, drop = FALSE]
}

# Log2-transform counts if matrix appears to be on a linear scale
if (max(expr_mat, na.rm = TRUE) > 100) {
  message("Log2 transforming expression values...")
  expr_mat <- log2(expr_mat + 1)
}

# ==============================================================================
# 5. Calculate Pathway Scores & Export Publication Figures
# ==============================================================================
message("4/5: Calculating pathway scores...")

avail_t1 <- intersect(gsea_t1, rownames(expr_mat))
avail_t2 <- intersect(gsea_t2, rownames(expr_mat))

message("Available Type I IFN genes: ", length(avail_t1))
message("Available Type II IFN genes: ", length(avail_t2))

pheno_htr8$TypeI_IFN_Score  <- colMeans(expr_mat[avail_t1, , drop = FALSE], na.rm = TRUE)
pheno_htr8$TypeII_IFN_Score <- colMeans(expr_mat[avail_t2, , drop = FALSE], na.rm = TRUE)

message("5/5: Generating publication plots...")
group_colors <- c("Control (NF)" = "#2B5C8F", "ZIKV" = "#B20000")

tight_theme <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(face = "bold", color = "black"),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.margin     = margin(t = 5, r = 5, b = 5, l = 5)
  )

p_t1 <- ggplot(pheno_htr8, aes(x = Group, y = TypeI_IFN_Score, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.12) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.2) +
  geom_jitter(width = 0.08, size = 2.5, color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  tight_theme +
  labs(
    title = "Type I IFN (MSigDB IFN-Alpha)\nHTR8 Trophoblasts (GSE234128)",
    y     = "Mean Pathway Expression (Log2)"
  ) +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.35)

p_t2 <- ggplot(pheno_htr8, aes(x = Group, y = TypeII_IFN_Score, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.12) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.2) +
  geom_jitter(width = 0.08, size = 2.5, color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  tight_theme +
  labs(
    title = "Type II IFN (MSigDB IFN-Gamma)\nHTR8 Trophoblasts (GSE234128)",
    y     = NULL
  ) +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.35)

combined_fig <- ggarrange(p_t1, p_t2, ncol = 2, labels = c("A", "B"))

ggsave(
  filename = file.path(base_dir, "Figure_GSE234128_HTR8_Pathways.pdf"),
  plot     = combined_fig,
  width    = 6,
  height   = 4,
  dpi      = 300
)

ggsave(
  filename = file.path(base_dir, "Figure_GSE234128_HTR8_Pathways.png"),
  plot     = combined_fig,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# Print package versions for reproducibility verification
sessionInfo()