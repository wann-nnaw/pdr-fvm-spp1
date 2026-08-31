#!/usr/bin/env Rscript

# SPP1 upstream transcription-factor analysis in macrophages
# Run from the DR project root:
#   Rscript codes_GitHub_fixed/10_SPP1_transcription_factor_analysis.R
#
# Strategy:
#   1) export macrophage counts and metadata for pySCENIC;
#   2) if a pySCENIC loom is present, analyse regulon AUC by subtype/SPP1 group;
#   3) compare candidate TF expression between SPP1-high and SPP1-low cells;
#   4) correlate TF and SPP1 expression in sample x subtype pseudobulks;
#   5) integrate regulon activity (when available), DE and correlation evidence.
# Association is not proof of direct transcriptional regulation; candidates should
# be validated by TF perturbation and ChIP/CUT&Tag or promoter-reporter assays.

set.seed(20260705)
options(stringsAsFactors = FALSE)

required_packages <- c(
  "Seurat", "Matrix", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork",
  "pheatmap"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1. Paths and analysis settings
# -----------------------------------------------------------------------------
input_rda <- "rda/UMAP_annotated.rda"
out_dir <- "results/SPP1_transcription_factor"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rda)) {
  stop("Input not found: ", input_rda, ". Run from the DR project root.")
}

macrophage_types <- c("LAM", "IAM", "OSM", "TRM")
min_cells_per_pseudobulk <- 20L
export_pyscenic_input <- TRUE
scenic_loom_file <- file.path(out_dir, "SPP1_macrophage_pySCENIC.loom")

# Candidate panel based on known osteopontin/inflammatory macrophage biology.
# The mechanism labels are priors to interpret, not results inferred from this data.
tf_prior <- tibble::tribble(
  ~TF,      ~prior_pathway,
  "CEBPB",  "myeloid activation / inflammatory enhancer",
  "JUN",    "AP-1",
  "FOS",    "AP-1",
  "JUNB",   "AP-1",
  "FOSL2",  "AP-1",
  "RELA",   "canonical NF-kB",
  "NFKB1",  "canonical NF-kB",
  "STAT3",  "cytokine-JAK/STAT",
  "STAT1",  "interferon-JAK/STAT",
  "HIF1A",  "hypoxia / glycolytic activation",
  "SP1",    "SP/KLF promoter regulation",
  "RUNX2",  "osteopontin promoter regulation",
  "ETS1",   "ETS inflammatory regulation",
  "SPI1",   "myeloid lineage (PU.1)",
  "IRF8",   "myeloid / interferon regulation",
  "PPARG",  "lipid-associated macrophage state",
  "TFEB",   "lysosome / phagosome program",
  "MITF",   "MiT/TFE lysosomal program",
  "KLF4",   "macrophage polarization",
  "MAFB",   "macrophage identity"
)

subtype_cols <- c(LAM = "#3A7CA5", IAM = "#B5502A", OSM = "#4D8B6F", TRM = "#B08C3A")
group_cols <- c("SPP1-low" = "#77A9B8", "SPP1-high" = "#C56A49")

theme_paper <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "grey35", hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), p,
         width = width, height = height, dpi = 600, bg = "white")
}

save_heatmap <- function(mat, annotation_col, stem, width = 8, height = 7) {
  grDevices::pdf(file.path(out_dir, paste0(stem, ".pdf")), width, height)
  pheatmap::pheatmap(
    mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
    annotation_col = annotation_col, border_color = NA,
    color = colorRampPalette(c("#356FA1", "white", "#B9473D"))(101),
    fontsize = 16, fontsize_row = 14, fontsize_col = 16
  )
  grDevices::dev.off()
  grDevices::png(
    file.path(out_dir, paste0(stem, ".png")), width = width,
    height = height, units = "in", res = 600, bg = "white"
  )
  pheatmap::pheatmap(
    mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
    annotation_col = annotation_col, border_color = NA,
    color = colorRampPalette(c("#356FA1", "white", "#B9473D"))(101),
    fontsize = 16, fontsize_row = 14, fontsize_col = 16
  )
  grDevices::dev.off()
}

safe_cor_test <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 4L || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(c(rho = NA_real_, p_value = NA_real_))
  }
  z <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(z$estimate), p_value = z$p.value)
}

# -----------------------------------------------------------------------------
# 2. Load macrophages and define within-subtype SPP1-high/low cells
# -----------------------------------------------------------------------------
load(input_rda)
if (!exists("combined_anno") || !inherits(combined_anno, "Seurat")) {
  stop("The input RDA must contain a Seurat object named combined_anno.")
}
needed_meta <- c("sample", "cell_type")
if (!all(needed_meta %in% colnames(combined_anno@meta.data))) {
  stop("Missing metadata column(s): ",
       paste(setdiff(needed_meta, colnames(combined_anno@meta.data)), collapse = ", "))
}

DefaultAssay(combined_anno) <- "SCT"
mac <- subset(combined_anno, subset = cell_type %in% macrophage_types)
genes_available <- rownames(mac[["SCT"]])
tf_prior <- tf_prior %>% filter(TF %in% genes_available)
analysis_genes <- unique(c("SPP1", tf_prior$TF))
if (!"SPP1" %in% genes_available) stop("SPP1 is absent from the SCT assay.")

expr <- FetchData(mac, vars = analysis_genes)
meta <- mac@meta.data %>%
  mutate(cell = rownames(.), SPP1_expr = expr[, "SPP1"]) %>%
  group_by(cell_type) %>%
  mutate(
    subtype_SPP1_median = median(SPP1_expr, na.rm = TRUE),
    SPP1_group = ifelse(SPP1_expr > subtype_SPP1_median, "SPP1-high", "SPP1-low")
  ) %>%
  ungroup()
mac$SPP1_group <- factor(meta$SPP1_group, levels = c("SPP1-low", "SPP1-high"))
mac$cell_type <- factor(mac$cell_type, levels = macrophage_types)

# pySCENIC input: cells x genes raw counts, matching the reference workflow.
# Low-information genes are removed to keep the input compact and reproducible.
if (export_pyscenic_input) {
  scenic_input_file <- file.path(out_dir, "pySCENIC_macrophage_counts.csv.gz")
  scenic_meta_file <- file.path(out_dir, "pySCENIC_macrophage_metadata.csv")
  if (!file.exists(scenic_input_file)) {
    rna_counts <- GetAssayData(mac, assay = "RNA", layer = "counts")
    keep_scenic_gene <- Matrix::rowSums(rna_counts > 0) >= 3 &
      Matrix::rowSums(rna_counts) >= 10
    scenic_counts <- as.matrix(t(rna_counts[keep_scenic_gene, , drop = FALSE]))
    gz_con <- gzfile(scenic_input_file, open = "wt")
    write.csv(scenic_counts, gz_con, quote = FALSE)
    close(gz_con)
    rm(rna_counts, scenic_counts)
    invisible(gc())
  }
  write.csv(
    mac@meta.data %>%
      mutate(CellID = rownames(.)) %>%
      select(CellID, sample, cell_type, SPP1_group),
    scenic_meta_file, row.names = FALSE
  )
}

# The integrated object contains one SCT model per sample. Recalculate the SCT
# counts to a common minimum median UMI before differential testing.
mac <- PrepSCTFindMarkers(mac, assay = "SCT", verbose = FALSE)

write.csv(
  meta %>% count(cell_type, SPP1_group, name = "cell_count"),
  file.path(out_dir, "SPP1_group_cell_counts.csv"), row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Differential TF expression: SPP1-high versus SPP1-low, adjusted by subtype
# -----------------------------------------------------------------------------
Idents(mac) <- mac$SPP1_group
de_all <- FindMarkers(
  mac,
  ident.1 = "SPP1-high",
  ident.2 = "SPP1-low",
  features = tf_prior$TF,
  assay = "SCT",
  slot = "data",
  test.use = "LR",
  recorrect_umi = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  latent.vars = "cell_type"
)
de_all$TF <- rownames(de_all)
de_all <- as_tibble(de_all) %>%
  rename_with(~ sub("avg_log2FC", "avg_log2FC_high_vs_low", .x)) %>%
  left_join(tf_prior, by = "TF")
write.csv(de_all, file.path(out_dir, "SPP1_high_vs_low_candidate_TF_DE.csv"), row.names = FALSE)

de_by_subtype <- lapply(macrophage_types, function(ct) {
  obj <- subset(mac, subset = cell_type == ct)
  obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
  Idents(obj) <- obj$SPP1_group
  if (length(unique(Idents(obj))) < 2L) return(NULL)
  z <- FindMarkers(
    obj, ident.1 = "SPP1-high", ident.2 = "SPP1-low",
    features = tf_prior$TF, assay = "SCT", slot = "data",
    min.pct = 0, logfc.threshold = 0, test.use = "wilcox",
    recorrect_umi = FALSE
  )
  z$TF <- rownames(z)
  as_tibble(z) %>% mutate(cell_type = ct)
}) %>% bind_rows() %>%
  rename_with(~ sub("avg_log2FC", "avg_log2FC_high_vs_low", .x))
write.csv(de_by_subtype, file.path(out_dir, "SPP1_high_vs_low_TF_DE_by_subtype.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Sample x subtype pseudobulk correlation (reduces single-cell pseudoreplication)
# -----------------------------------------------------------------------------
pb_id <- interaction(mac$sample, mac$cell_type, sep = "__", drop = TRUE)
cell_counts <- table(pb_id)
keep_pb <- names(cell_counts)[cell_counts >= min_cells_per_pseudobulk]
keep_cells <- colnames(mac)[pb_id %in% keep_pb]
mac_pb <- subset(mac, cells = keep_cells)
mac_pb$pseudobulk_id <- droplevels(pb_id[pb_id %in% keep_pb])

pb <- AggregateExpression(
  mac_pb,
  assays = "SCT",
  features = analysis_genes,
  group.by = "pseudobulk_id",
  slot = "data",
  return.seurat = FALSE,
  verbose = FALSE
)$SCT
pb_log <- log1p(as.matrix(pb))

count_lookup <- setNames(as.integer(cell_counts), gsub("_", "-", names(cell_counts)))
pb_meta <- tibble(pseudobulk_id = colnames(pb_log)) %>%
  separate(pseudobulk_id, into = c("sample", "cell_type"), sep = "--", remove = FALSE) %>%
  mutate(cell_count = unname(count_lookup[pseudobulk_id]))

cor_results <- lapply(tf_prior$TF, function(tf) {
  z <- safe_cor_test(pb_log[tf, ], pb_log["SPP1", ])
  tibble(TF = tf, rho = z[["rho"]], p_value = z[["p_value"]])
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  left_join(tf_prior, by = "TF")
write.csv(cor_results, file.path(out_dir, "SPP1_TF_pseudobulk_correlations.csv"), row.names = FALSE)

pb_export <- as.data.frame(t(pb_log)) %>%
  tibble::rownames_to_column("pseudobulk_id") %>%
  left_join(pb_meta, by = "pseudobulk_id")
write.csv(pb_export, file.path(out_dir, "sample_subtype_pseudobulk_expression.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Evidence integration and candidate ranking
# -----------------------------------------------------------------------------
ranked <- de_all %>%
  select(TF, avg_log2FC_high_vs_low, p_val_adj, prior_pathway) %>%
  full_join(cor_results %>% select(TF, rho, cor_p_adj = p_adj), by = "TF") %>%
  mutate(
    de_score = pmax(avg_log2FC_high_vs_low, 0) * pmin(-log10(p_val_adj + 1e-300), 20),
    cor_score = pmax(rho, 0) * pmin(-log10(cor_p_adj + 1e-300), 20),
    # Percentile integration remains defined even if one evidence component is
    # constant (for example, no TF passes the multiple-testing cutoff).
    integrated_score = percent_rank(de_score) + percent_rank(cor_score),
    evidence = case_when(
      avg_log2FC_high_vs_low > 0.1 & p_val_adj < 0.05 & rho > 0.3 & cor_p_adj < 0.1 ~ "DE + correlation",
      avg_log2FC_high_vs_low > 0.1 & p_val_adj < 0.05 ~ "DE only",
      rho > 0.3 & cor_p_adj < 0.1 ~ "correlation only",
      TRUE ~ "weak in this dataset"
    )
  ) %>%
  arrange(desc(integrated_score))
write.csv(ranked, file.path(out_dir, "SPP1_candidate_TF_integrated_ranking.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Visualizations
# -----------------------------------------------------------------------------
p_rank <- ranked %>%
  mutate(TF = factor(TF, levels = rev(TF))) %>%
  ggplot(aes(TF, integrated_score, fill = evidence)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c(
    "DE + correlation" = "#A83E32", "DE only" = "#D78B55",
    "correlation only" = "#477FA3", "weak in this dataset" = "#C7C7C7"
  )) +
  labs(
    title = "Candidate transcription factors associated with SPP1",
    subtitle = "Integrated SPP1-high differential expression and pseudobulk correlation",
    x = NULL, y = "Integrated evidence score", fill = "Evidence"
  ) + theme_paper()
save_plot(p_rank, "01_SPP1_candidate_TF_ranking", 7.2, 5.6)

p_cor <- cor_results %>%
  ggplot(aes(rho, -log10(p_adj + 1e-300), label = TF)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey55") +
  geom_vline(xintercept = 0, color = "grey75") +
  geom_point(aes(color = rho), size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_gradient2(low = "#3A7CA5", mid = "grey85", high = "#B5502A", midpoint = 0) +
  labs(
    title = "TF-SPP1 pseudobulk association",
    subtitle = paste0("Spearman correlation across sample x subtype units (n = ", ncol(pb_log), ")"),
    x = "Spearman rho", y = "-log10(BH-adjusted P)", color = "rho"
  ) + theme_paper()
save_plot(p_cor, "02_TF_SPP1_pseudobulk_correlation", 7.2, 5.2)

top_tfs <- head(ranked$TF[is.finite(ranked$integrated_score)], 8)
heat_df <- pb_export %>%
  select(pseudobulk_id, sample, cell_type, all_of(c("SPP1", top_tfs))) %>%
  pivot_longer(all_of(c("SPP1", top_tfs)), names_to = "gene", values_to = "expression") %>%
  group_by(gene) %>%
  mutate(z = as.numeric(scale(expression))) %>%
  ungroup()

p_heat <- ggplot(heat_df, aes(pseudobulk_id, gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#356FA1", mid = "white", high = "#B9473D", midpoint = 0) +
  facet_grid(~ cell_type, scales = "free_x", space = "free_x") +
  labs(
    title = "SPP1 and top candidate TF expression",
    subtitle = "Gene-wise z-score in sample x macrophage-subtype pseudobulks",
    x = "Pseudobulk", y = NULL, fill = "z-score"
  ) +
  theme_paper(9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_plot(p_heat, "03_SPP1_top_TF_pseudobulk_heatmap", 9.5, 4.8)

violin_df <- FetchData(mac, vars = top_tfs) %>%
  mutate(SPP1_group = mac$SPP1_group, cell_type = mac$cell_type) %>%
  pivot_longer(all_of(top_tfs), names_to = "TF", values_to = "expression")
p_violin <- ggplot(violin_df, aes(SPP1_group, expression, fill = SPP1_group)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  stat_summary(fun = median, geom = "point", size = 0.8, color = "black") +
  facet_grid(TF ~ cell_type, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  labs(
    title = "Top candidate TFs in SPP1-high and SPP1-low macrophages",
    subtitle = "SPP1 groups are defined by the median within each macrophage subtype",
    x = NULL, y = "SCT normalized expression", fill = NULL
  ) + theme_paper(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top",
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_plot(p_violin, "04_top_TF_expression_by_SPP1_group", 10, 10)

# -----------------------------------------------------------------------------
# 7. Optional pySCENIC regulon-AUC analysis
# -----------------------------------------------------------------------------
scenic_status <- paste0(
  "not run: place the pySCENIC loom at ", scenic_loom_file,
  " and install SCopeLoomR"
)
if (file.exists(scenic_loom_file)) {
  if (!requireNamespace("SCopeLoomR", quietly = TRUE)) {
    warning(
      "SCENIC loom found but SCopeLoomR is unavailable. Install it and rerun: ",
      scenic_loom_file
    )
    scenic_status <- "loom found; analysis skipped because SCopeLoomR is unavailable"
  } else {
    loom <- SCopeLoomR::open_loom(scenic_loom_file, mode = "r")
    on.exit(try(SCopeLoomR::close_loom(loom), silent = TRUE), add = TRUE)
    scenic_cell_id <- loom[["col_attrs/CellID"]][]
    regulon_auc <- loom[["col_attrs/RegulonsAUC"]][]
    rownames(regulon_auc) <- scenic_cell_id
    colnames(regulon_auc) <- gsub("[()]", "", colnames(regulon_auc))

    common_cells <- intersect(colnames(mac), rownames(regulon_auc))
    if (length(common_cells) < 100) {
      stop("Fewer than 100 cells overlap between the Seurat object and SCENIC loom.")
    }
    regulon_auc <- regulon_auc[common_cells, , drop = FALSE]
    scenic_meta <- mac@meta.data[common_cells, , drop = FALSE] %>%
      mutate(CellID = rownames(.))

    mean_auc_subtype <- as.data.frame(regulon_auc) %>%
      mutate(cell_type = scenic_meta$cell_type) %>%
      group_by(cell_type) %>%
      summarise(across(where(is.numeric), mean), .groups = "drop")
    auc_mat <- t(as.matrix(mean_auc_subtype[, -1, drop = FALSE]))
    colnames(auc_mat) <- mean_auc_subtype$cell_type
    regulon_variance <- apply(auc_mat, 1, var)
    top_regulons <- names(sort(regulon_variance, decreasing = TRUE))[
      seq_len(min(20, length(regulon_variance)))
    ]
    annotation_col <- data.frame(
      Subtype = factor(colnames(auc_mat), levels = macrophage_types),
      row.names = colnames(auc_mat)
    )
    save_heatmap(
      auc_mat[top_regulons, , drop = FALSE], annotation_col,
      "05_SCENIC_top_regulon_heatmap", width = 7.2, height = 7
    )

    spp1_vec <- FetchData(mac, vars = "SPP1")[common_cells, 1]
    scenic_cor <- lapply(colnames(regulon_auc), function(regulon) {
      z <- safe_cor_test(regulon_auc[, regulon], spp1_vec)
      tibble(regulon = regulon, rho = z[["rho"]], p_value = z[["p_value"]])
    }) %>% bind_rows() %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(desc(rho))
    write.csv(
      scenic_cor, file.path(out_dir, "SCENIC_regulon_SPP1_correlations.csv"),
      row.names = FALSE
    )

    auc_group <- as.data.frame(regulon_auc) %>%
      mutate(SPP1_group = scenic_meta$SPP1_group) %>%
      pivot_longer(-SPP1_group, names_to = "regulon", values_to = "AUC") %>%
      group_by(regulon) %>%
      summarise(
        mean_high = mean(AUC[SPP1_group == "SPP1-high"]),
        mean_low = mean(AUC[SPP1_group == "SPP1-low"]),
        delta_AUC = mean_high - mean_low,
        p_value = wilcox.test(
          AUC[SPP1_group == "SPP1-high"], AUC[SPP1_group == "SPP1-low"]
        )$p.value,
        .groups = "drop"
      ) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      left_join(scenic_cor %>% select(regulon, rho), by = "regulon") %>%
      arrange(desc(delta_AUC))
    write.csv(
      auc_group, file.path(out_dir, "SCENIC_regulon_SPP1_high_vs_low.csv"),
      row.names = FALSE
    )

    scenic_top <- auc_group %>%
      filter(delta_AUC > 0) %>%
      slice_max(order_by = abs(delta_AUC * -log10(p_adj + 1e-300)), n = 20)
    p_scenic <- ggplot(
      scenic_top,
      aes(delta_AUC, -log10(p_adj + 1e-300), label = regulon, colour = rho)
    ) +
      geom_vline(xintercept = 0, colour = "grey70") +
      geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey55") +
      geom_point(size = 2.6) +
      geom_text_repel(size = 2.8, max.overlaps = Inf) +
      scale_colour_gradient2(
        low = "#356FA1", mid = "grey85", high = "#B9473D", midpoint = 0
      ) +
      labs(
        title = "SCENIC regulons associated with the SPP1-high state",
        subtitle = "Activity difference and correlation with SPP1 expression",
        x = "Mean regulon AUC difference (SPP1-high minus SPP1-low)",
        y = "-log10(BH-adjusted P)", colour = "SPP1 rho"
      ) + theme_paper()
    save_plot(p_scenic, "06_SCENIC_SPP1_regulon_activity", 7.2, 5.4)
    scenic_status <- paste0(
      "completed: ", length(common_cells), " cells and ", ncol(regulon_auc),
      " regulons"
    )
  }
}
writeLines(scenic_status, file.path(out_dir, "SCENIC_status.txt"))

summary_lines <- c(
  "SPP1 transcription-factor analysis summary",
  paste0("Macrophages analysed: ", ncol(mac)),
  paste0("Pseudobulk units retained (>=", min_cells_per_pseudobulk, " cells): ", ncol(pb_log)),
  paste0("Candidate TFs tested: ", nrow(tf_prior)),
  paste0("SCENIC: ", scenic_status),
  "",
  "Top integrated candidates:",
  paste0(seq_len(min(10, nrow(ranked))), ". ", head(ranked$TF, 10),
         " | score=", round(head(ranked$integrated_score, 10), 2),
         " | evidence=", head(ranked$evidence, 10)),
  "",
  "Interpretation: candidates supported by both SPP1-high differential expression and",
  "sample-by-subtype pseudobulk correlation are prioritized. These associations do not",
  "establish direct TF binding or causality. Validate with perturbation plus ChIP/CUT&Tag."
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.txt"))

message("Completed. Results written to: ", normalizePath(out_dir))

# -----------------------------------------------------------------------------
# 11. Unified Figure 7 visualizations (panels A-E and supplementary 07-08)
# -----------------------------------------------------------------------------
auc_vis <- read.csv(
  file.path(out_dir, "SCENIC_regulon_AUC.csv.gz"),
  row.names = 1, check.names = FALSE
)
meta_vis <- read.csv(
  file.path(out_dir, "pySCENIC_macrophage_metadata.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
rownames(meta_vis) <- meta_vis$CellID
pt_obj <- readRDS("results/pseudotime/macrophage_with_pseudotime.rds")
vis_cells <- Reduce(intersect, list(rownames(auc_vis), rownames(meta_vis), colnames(pt_obj)))
auc_vis <- as.matrix(auc_vis[vis_cells, , drop = FALSE])
meta_vis <- meta_vis[vis_cells, , drop = FALSE]
pt_obj <- subset(pt_obj, cells = vis_cells)
auc_vis <- auc_vis[colnames(pt_obj), , drop = FALSE]
meta_vis <- meta_vis[colnames(pt_obj), , drop = FALSE]

if (nrow(auc_vis) != 3157 || ncol(auc_vis) != 292) {
  warning(
    "Expected 3157 cells and 292 regulons; observed ",
    nrow(auc_vis), " cells and ", ncol(auc_vis), " regulons."
  )
}

clean_reg <- function(x) sub("\\(\\+\\)$", "", x)
reg_lookup <- setNames(colnames(auc_vis), clean_reg(colnames(auc_vis)))
selected_names <- c("ATF5", "HOXB5", "JUN", "CEBPB", "FOS", "SPI1", "CEBPA", "STAT1")
selected_regs <- unname(reg_lookup[intersect(selected_names, names(reg_lookup))])
selected_regs <- head(selected_regs, 8)
figure7_cols <- c(IAM = "#B9473D", LAM = "#D99532", OSM = "#3976A8", TRM = "#4B9B73")

theme_figure7 <- function(base_size = 19, title_size = 22, subtitle_size = 13) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.45, colour = "black"),
      axis.text = element_text(size = base_size - 2, colour = "black"),
      axis.title = element_text(size = base_size, face = "bold", colour = "black"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = subtitle_size, hjust = 0.5, colour = "grey35"),
      legend.title = element_text(size = base_size - 3, face = "bold"),
      legend.text = element_text(size = base_size - 4),
      strip.background = element_rect(fill = "grey94", colour = "grey80", linewidth = 0.3),
      strip.text = element_text(size = base_size - 1, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(9, 11, 9, 11)
    )
}

save_figure7 <- function(p, stem, width, height) {
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")), p,
    width = width, height = height, units = "in", bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".png")), p,
    width = width, height = height, units = "in", dpi = 600,
    bg = "white", limitsize = FALSE
  )
}

# Figure 7D: top-variable regulon heatmap.
mean_auc_vis <- as.data.frame(auc_vis) %>%
  mutate(cell_type = meta_vis$cell_type) %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")
heat_vis <- t(as.matrix(mean_auc_vis[, -1, drop = FALSE]))
colnames(heat_vis) <- mean_auc_vis$cell_type
top_vis <- names(sort(apply(heat_vis, 1, var), decreasing = TRUE))[1:30]
heat_vis <- heat_vis[top_vis, , drop = FALSE]
rownames(heat_vis) <- clean_reg(rownames(heat_vis))
ann_vis <- data.frame(Subtype = colnames(heat_vis), row.names = colnames(heat_vis))
draw_figure7_heatmap <- function() {
  heatmap_plot <- pheatmap::pheatmap(
    heat_vis, scale = "row", cluster_cols = FALSE,
    annotation_col = ann_vis, annotation_names_col = FALSE,
    annotation_colors = list(Subtype = figure7_cols),
    border_color = NA,
    color = colorRampPalette(c("#356FA1", "white", "#B9473D"))(101),
    fontsize = 14, fontsize_row = 11.5, fontsize_col = 14,
    main = "Macrophage subtype-specific\nregulon activity",
    treeheight_row = 30, silent = TRUE
  )
  main_index <- which(heatmap_plot$gtable$layout$name == "main")
  if (length(main_index) == 1) {
    heatmap_plot$gtable$grobs[[main_index]]$gp <- grid::gpar(
      fontsize = 22, fontface = "bold", fontfamily = "sans"
    )
  }
  heatmap_plot
}
pdf(file.path(out_dir, "05_SCENIC_top_regulon_heatmap.pdf"), 7.2, 6.4)
grid::grid.newpage(); grid::grid.draw(draw_figure7_heatmap()$gtable); dev.off()
png(file.path(out_dir, "05_SCENIC_top_regulon_heatmap.png"),
    7.2, 6.4, units = "in", res = 600, bg = "white")
grid::grid.newpage(); grid::grid.draw(draw_figure7_heatmap()$gtable); dev.off()

# Figure 7B: SPP1-high regulon association.
DefaultAssay(pt_obj) <- "SCT"
spp1_vis <- FetchData(pt_obj, vars = "SPP1", layer = "data")[rownames(auc_vis), 1]
cor_vis <- lapply(colnames(auc_vis), function(reg) {
  z <- suppressWarnings(cor.test(auc_vis[, reg], spp1_vis, method = "spearman", exact = FALSE))
  data.frame(regulon = reg, rho = unname(z$estimate), p_value = z$p.value)
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, "BH")) %>% arrange(desc(rho))
group_vis <- as.data.frame(auc_vis) %>%
  mutate(SPP1_group = meta_vis$SPP1_group) %>%
  pivot_longer(-SPP1_group, names_to = "regulon", values_to = "AUC") %>%
  group_by(regulon) %>%
  summarise(
    mean_high = mean(AUC[SPP1_group == "SPP1-high"]),
    mean_low = mean(AUC[SPP1_group == "SPP1-low"]),
    delta_AUC = mean_high - mean_low,
    p_value = wilcox.test(AUC[SPP1_group == "SPP1-high"], AUC[SPP1_group == "SPP1-low"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, "BH")) %>%
  left_join(cor_vis %>% select(regulon, rho), by = "regulon")
write.csv(cor_vis, file.path(out_dir, "SCENIC_regulon_SPP1_correlations.csv"), row.names = FALSE)
write.csv(group_vis, file.path(out_dir, "SCENIC_regulon_SPP1_high_vs_low.csv"), row.names = FALSE)
label_vis <- group_vis %>%
  mutate(score = abs(delta_AUC) * -log10(p_adj + 1e-300)) %>%
  slice_max(score, n = 25)
p06_vis <- ggplot(group_vis, aes(delta_AUC, -log10(p_adj + 1e-300), colour = rho)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_hline(yintercept = -log10(.05), linetype = 2, colour = "grey50") +
  geom_point(alpha = .68, size = 3.1) +
  geom_text_repel(
    data = label_vis, aes(label = clean_reg(regulon)), size = 5,
    max.overlaps = Inf, min.segment.length = 0, segment.size = .35
  ) +
  scale_colour_gradient2(low = "#356FA1", mid = "grey85", high = "#B9473D", midpoint = 0) +
  labs(
    title = "SCENIC regulons associated with\nthe SPP1-high macrophage state",
    subtitle = paste0(nrow(auc_vis), " cells; ", ncol(auc_vis), " motif-pruned regulons"),
    x = "Mean AUC difference\n(SPP1-high minus SPP1-low)",
    y = "-log10(BH-adjusted P)", colour = "SPP1 rho"
  ) + theme_figure7(18, 22, 14)
save_figure7(p06_vis, "06_SCENIC_SPP1_regulon_activity", 7.5, 5.8)

# Supplementary 07: regulon AUC on t-SNE.
emb_vis <- as.data.frame(Embeddings(pt_obj, "tsne")); emb_vis$CellID <- rownames(emb_vis)
proj_vis <- cbind(emb_vis, as.data.frame(auc_vis[emb_vis$CellID, selected_regs, drop = FALSE])) %>%
  pivot_longer(all_of(selected_regs), names_to = "Regulon", values_to = "AUC") %>%
  mutate(Regulon = clean_reg(Regulon))
p07_vis <- ggplot(proj_vis, aes(tSNE_1, tSNE_2, colour = AUC)) +
  geom_point(size = .32, alpha = .8) + facet_wrap(~Regulon, ncol = 4) +
  scale_colour_viridis_c(option = "magma", trans = "sqrt") + coord_equal() +
  labs(title = "Regulon activity landscape in macrophages", colour = "AUC") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = .5),
        strip.text = element_text(size = 13, face = "bold"))
save_figure7(p07_vis, "07_SCENIC_regulon_activity_tSNE", 10, 5.8)

# Supplementary 08: SPP1-high/low regulon distributions.
long08_vis <- as.data.frame(auc_vis[, selected_regs, drop = FALSE]) %>%
  mutate(SPP1_group = meta_vis[rownames(.), "SPP1_group"]) %>%
  pivot_longer(all_of(selected_regs), names_to = "Regulon", values_to = "AUC") %>%
  mutate(Regulon = clean_reg(Regulon))
p08_vis <- ggplot(long08_vis, aes(SPP1_group, AUC, fill = SPP1_group)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = .8) +
  geom_boxplot(width = .14, outlier.shape = NA, fill = "white", linewidth = .3) +
  facet_wrap(~Regulon, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("SPP1-low" = "#3976A8", "SPP1-high" = "#B9473D")) +
  labs(title = "Regulon activity in SPP1-high and SPP1-low macrophages", x = NULL, y = "AUCell AUC") +
  theme_figure7(13, 18, 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_figure7(p08_vis, "08_SCENIC_regulon_activity_by_SPP1_group", 10, 5.8)

# Figure 7C: subtype regulon dot plot.
q75_vis <- apply(auc_vis[, selected_regs, drop = FALSE], 2, quantile, probs = .75)
dot_vis <- lapply(selected_regs, function(reg) {
  data.frame(cell_type = meta_vis$cell_type, Regulon = clean_reg(reg),
             AUC = auc_vis[, reg], high = auc_vis[, reg] >= q75_vis[reg])
}) %>% bind_rows() %>%
  group_by(cell_type, Regulon) %>%
  summarise(mean_AUC = mean(AUC), high_AUC_fraction = mean(high), .groups = "drop") %>%
  group_by(Regulon) %>% mutate(z_mean_AUC = as.numeric(scale(mean_AUC))) %>% ungroup()
p09_vis <- ggplot(dot_vis, aes(cell_type, Regulon)) +
  geom_point(aes(size = high_AUC_fraction, colour = z_mean_AUC)) +
  scale_colour_gradient2(low = "#356FA1", mid = "white", high = "#B9473D", midpoint = 0) +
  scale_size(range = c(3, 11), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Macrophage subtype-specific\nregulon activity",
       subtitle = "Colour: subtype mean AUC (row z-score)\nSize: cells above global 75th percentile",
       x = NULL, y = NULL, colour = "Mean AUC\n(row z-score)", size = "High-AUC\nfraction") +
  theme_figure7(19, 22, 11) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = .35),
    plot.subtitle = element_text(size = 11, hjust = .5, colour = "grey35",
                                 lineheight = 1.05, margin = margin(b = 7))
  )
save_figure7(p09_vis, "09_SCENIC_regulon_dotplot", 7.2, 6.4)
write.csv(dot_vis, file.path(out_dir, "SCENIC_regulon_subtype_summary.csv"), row.names = FALSE)

# Figure 7A: regulon specificity score.
jsd_vis <- function(p, q) {
  p <- p / sum(p); q <- q / sum(q); m <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  (kl(p, m) + kl(q, m)) / 2
}
rss_vis <- expand.grid(Regulon = colnames(auc_vis), cell_type = unique(meta_vis$cell_type),
                       stringsAsFactors = FALSE) %>%
  rowwise() %>% mutate(RSS = {
    p <- auc_vis[, Regulon] + 1e-12
    q <- as.numeric(meta_vis$cell_type == cell_type) + 1e-12
    1 - sqrt(jsd_vis(p, q))
  }) %>% ungroup() %>%
  group_by(cell_type) %>% arrange(desc(RSS), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>% ungroup()
top_rss_vis <- rss_vis %>% filter(rank <= 5) %>%
  mutate(label = reorder(paste0(clean_reg(Regulon), " - ", cell_type), RSS))
p10_vis <- ggplot(top_rss_vis, aes(RSS, label, colour = cell_type)) +
  geom_segment(aes(x = min(RSS) - .01, xend = RSS, yend = label),
               linewidth = .7, colour = "grey78") +
  geom_point(size = 4.6) + scale_colour_manual(values = figure7_cols) +
  labs(title = "Subtype-specific regulons\nidentified by RSS",
       subtitle = "Top five regulons per subtype; higher RSS indicates greater specificity",
       x = "Regulon specificity score (RSS)", y = NULL, colour = "Subtype") +
  theme_figure7(19, 25, 13)
save_figure7(p10_vis, "10_SCENIC_regulon_specificity_RSS", 8.2, 6.4)
write.csv(rss_vis, file.path(out_dir, "SCENIC_regulon_specificity_scores.csv"), row.names = FALSE)

# Figure 7E: regulon activity along pseudotime.
trend_regs_vis <- selected_regs[seq_len(min(6, length(selected_regs)))]
pt_meta_vis <- pt_obj@meta.data[colnames(pt_obj), c("Pseudotime", "cell_type"), drop = FALSE]
trend_vis <- as.data.frame(auc_vis[, trend_regs_vis, drop = FALSE]) %>%
  mutate(Pseudotime = pt_meta_vis[rownames(.), "Pseudotime"],
         cell_type = pt_meta_vis[rownames(.), "cell_type"]) %>%
  pivot_longer(all_of(trend_regs_vis), names_to = "Regulon", values_to = "AUC") %>%
  mutate(Regulon = clean_reg(Regulon))
p11_vis <- ggplot(trend_vis, aes(Pseudotime, AUC, colour = cell_type)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, linewidth = 1.15, alpha = .12) +
  facet_wrap(~Regulon, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = figure7_cols) +
  labs(title = "Regulon activity dynamics\nalong macrophage pseudotime",
       subtitle = "Smoothed AUCell AUC trends stratified by macrophage subtype",
       x = "Monocle2 pseudotime", y = "AUCell AUC", colour = "Subtype") +
  theme_figure7(13, 17, 12) +
  theme(panel.spacing = grid::unit(10, "pt"), legend.position = "right")
save_figure7(p11_vis, "11_SCENIC_regulon_pseudotime_trends", 10, 6.8)

writeLines(
  c("Unified Figure 7 visualization completed",
    paste0("cells=", nrow(auc_vis)), paste0("regulons=", ncol(auc_vis)),
    "05-11 are generated by 10_SPP1_transcription_factor_analysis.R"),
  file.path(out_dir, "SCENIC_extended_visualization_summary.txt")
)
message("Figure 7 visualizations completed from the unified script.")

# -----------------------------------------------------------------------------
# 审稿补充分析：SCENIC可重复性记录与regulon效应量范围
# -----------------------------------------------------------------------------

review_scenic_dir <- "results/reviewer_supplement/SCENIC_effect_size"
dir.create(review_scenic_dir, recursive = TRUE, showWarnings = FALSE)

review_auc_file <- file.path(out_dir, "SCENIC_regulon_AUC.csv.gz")
review_regulon_file <- file.path(out_dir, "SCENIC_regulons.csv")
review_group_file <- file.path(out_dir, "SCENIC_regulon_SPP1_high_vs_low.csv")
review_cor_file <- file.path(out_dir, "SCENIC_regulon_SPP1_correlations.csv")
review_fallback_file <- file.path(out_dir, "GRNBoost_regulon_AUC.csv.gz")

review_auc <- read.csv(review_auc_file, row.names = 1, check.names = FALSE)
review_group <- read.csv(review_group_file, check.names = FALSE)
review_cor <- read.csv(review_cor_file, check.names = FALSE)

review_effect <- review_group %>%
  left_join(
    review_cor %>% select(regulon, rho_from_correlation_file = rho, correlation_p_adj = p_adj),
    by = "regulon"
  ) %>%
  mutate(
    abs_delta_AUC = abs(delta_AUC),
    abs_rho = abs(rho),
    effect_direction = case_when(
      delta_AUC > 0 ~ "higher_in_SPP1_high",
      delta_AUC < 0 ~ "higher_in_SPP1_low",
      TRUE ~ "no_delta"
    )
  )

review_effect_summary <- data.frame(
  cells = nrow(review_auc),
  motif_pruned_regulons = ncol(review_auc),
  delta_AUC_min = min(review_effect$delta_AUC, na.rm = TRUE),
  delta_AUC_max = max(review_effect$delta_AUC, na.rm = TRUE),
  abs_delta_AUC_median = median(review_effect$abs_delta_AUC, na.rm = TRUE),
  rho_min = min(review_effect$rho, na.rm = TRUE),
  rho_max = max(review_effect$rho, na.rm = TRUE),
  abs_rho_median = median(review_effect$abs_rho, na.rm = TRUE),
  BH_significant_regulons = sum(review_effect$p_adj < 0.05, na.rm = TRUE),
  GRNBoost_AUCell_fallback_used_for_Figure7 = FALSE,
  fallback_file_present = file.exists(review_fallback_file),
  stringsAsFactors = FALSE
)

review_representative <- review_effect %>%
  arrange(desc(abs_rho), desc(abs_delta_AUC)) %>%
  select(regulon, delta_AUC, p_adj, rho, abs_delta_AUC, abs_rho, effect_direction) %>%
  head(30)

write.csv(review_effect, file.path(review_scenic_dir, "SCENIC_regulon_effect_sizes_all.csv"), row.names = FALSE)
write.csv(review_effect_summary, file.path(review_scenic_dir, "SCENIC_regulon_effect_size_summary.csv"), row.names = FALSE)
write.csv(review_representative, file.path(review_scenic_dir, "SCENIC_representative_regulon_effects.csv"), row.names = FALSE)

p_review_effect <- ggplot(
  review_effect,
  aes(x = delta_AUC, y = -log10(p_adj), colour = rho)
) +
  geom_point(alpha = 0.75, size = 1.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  scale_colour_gradient2(low = "#356FA1", mid = "grey90", high = "#B9473D", midpoint = 0) +
  labs(
    title = "SCENIC regulon association with SPP1-high macrophages",
    subtitle = "Colour indicates Spearman correlation with SPP1 expression",
    x = "Mean AUCell AUC difference (SPP1-high minus SPP1-low)",
    y = "-log10(BH-adjusted P)",
    colour = "Spearman rho"
  ) +
  theme_classic(base_size = 11)

p_review_rho_delta <- ggplot(
  review_effect,
  aes(x = rho, y = delta_AUC, colour = p_adj < 0.05)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_point(alpha = 0.75, size = 1.7) +
  scale_colour_manual(values = c("TRUE" = "#B9473D", "FALSE" = "#8EA7B8")) +
  labs(
    title = "Regulon effect-size range",
    x = "Spearman rho with SPP1 expression",
    y = "Delta AUCell AUC",
    colour = "BH < 0.05"
  ) +
  theme_classic(base_size = 11)

ggsave(file.path(review_scenic_dir, "SCENIC_regulon_effect_size_volcano.pdf"),
       p_review_effect, width = 6.6, height = 5.2, bg = "white")
ggsave(file.path(review_scenic_dir, "SCENIC_regulon_effect_size_volcano.png"),
       p_review_effect, width = 6.6, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(review_scenic_dir, "SCENIC_rho_delta_effect_size.pdf"),
       p_review_rho_delta, width = 5.8, height = 4.8, bg = "white")
ggsave(file.path(review_scenic_dir, "SCENIC_rho_delta_effect_size.png"),
       p_review_rho_delta, width = 5.8, height = 4.8, dpi = 300, bg = "white")

writeLines(
  c(
    "Reviewer supplement: SCENIC reproducibility and effect sizes",
    paste0("motif_pruned_AUC_file=", review_auc_file),
    paste0("regulon_file=", review_regulon_file),
    paste0("cells=", nrow(review_auc)),
    paste0("motif_pruned_regulons=", ncol(review_auc)),
    paste0(
      "delta_AUC_range=",
      signif(review_effect_summary$delta_AUC_min, 4), " to ",
      signif(review_effect_summary$delta_AUC_max, 4)
    ),
    paste0(
      "rho_range=",
      signif(review_effect_summary$rho_min, 4), " to ",
      signif(review_effect_summary$rho_max, 4)
    ),
    "Figure 7 uses pySCENIC GRNBoost2 + cisTarget motif-pruned regulons and AUCell AUC.",
    "The GRNBoost-AUCell fallback output was not used for Figure 7 or manuscript conclusions."
  ),
  file.path(review_scenic_dir, "SCENIC_reproducibility_effect_size_summary.txt")
)
