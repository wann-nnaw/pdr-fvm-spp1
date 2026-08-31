library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(forcats)

# ==============================================================================
# GO/KEGG enrichment for macrophage subtype up-regulated DEGs
# ==============================================================================

deg_base_dir <- "results/macrophage_DEG_GSEA"
go_kegg_out_dir <- "results/go_kegg"

comparisons <- c("LAM_vs_others", "IAM_vs_others", "OSM_vs_others", "TRM_vs_others")

deg_files <- file.path(
  deg_base_dir,
  comparisons,
  paste0("degs_", comparisons, "_annotated.csv")
)

names(deg_files) <- comparisons

go_ontology_colors <- c(
  "BP" = "#4E79A7",
  "MF" = "#59A14F",
  "CC" = "#E15759"
)

kegg_gradient_colors <- c(
  "#6BAED6",
  "#74C69D",
  "#F2C14E",
  "#F28E2B",
  "#D95F59"
)

# ==============================================================================
# 1. Nature-style theme  （from go_kegg_fixed.R）
# ==============================================================================

theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid        = element_blank(),
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(color = "black", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size),
      plot.title        = element_text(size = base_size + 2, face = "bold",
                                       hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = base_size - 1, color = "grey45",
                                       hjust = 0.5, margin = margin(b = 6)),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(9, "pt"),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 1, color = "black"),
      legend.position   = "right",
      legend.margin     = margin(0, 0, 0, 4),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.margin       = margin(8, 8, 8, 8)
    )
}

theme_pathway_publication <- function(base_size = 12) {
  theme_nature(base_size = base_size) +
    theme(
      axis.text.x  = element_text(size = 16, face = "bold", color = "black"),
      axis.text.y  = element_text(size = 15, face = "bold", color = "black"),
      axis.title.x = element_text(size = 18, face = "bold", color = "black",
                                   margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, face = "bold", color = "black"),
      plot.title   = element_text(size = 24, face = "bold", hjust = 0.5,
                                   margin = margin(b = 18)),
      plot.subtitle = element_blank(),
      strip.text   = element_text(size = 17, face = "bold", color = "black",
                                   hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold", color = "black"),
      legend.text  = element_text(size = 12, face = "bold", color = "black"),
      legend.key.width  = unit(0.45, "cm"),
      legend.key.height = unit(0.45, "cm"),
      plot.margin  = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# 2. Helper utilities  （from go_kegg_fixed.R）
# ==============================================================================

# Robust "numerator/denominator" → numeric parser
parse_gene_ratio <- function(x) {
  x <- as.character(x)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_real_)
    parts <- strsplit(value, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_real_)
    num <- suppressWarnings(as.numeric(parts[1]))
    den <- suppressWarnings(as.numeric(parts[2]))
    if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
    num / den
  }, numeric(1))
}

# Output dimensions
pathway_plot_width  <- 12
pathway_plot_height <- 8

save_pathway_plot <- function(filename, plot, width, height, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white")
}

# ==============================================================================
# 3. Main enrichment + plot function
# ==============================================================================

run_go_kegg_one <- function(deg_file, comparison_name, out_dir,
                            padj_cutoff  = 0.05,
                            logfc_cutoff = 0.5,
                            go_top_n     = 5,
                            kegg_top_n   = 15) {

  comparison_out_dir <- file.path(out_dir, comparison_name)
  dir.create(comparison_out_dir, recursive = TRUE, showWarnings = FALSE)

  degs <- read.csv(deg_file, stringsAsFactors = FALSE)

  sig_degs <- degs %>%
    filter(p_val_adj < padj_cutoff, avg_log2FC > logfc_cutoff) %>%
    arrange(p_val_adj, desc(avg_log2FC))

  write.csv(
    sig_degs,
    file.path(comparison_out_dir, paste0("sig_genes_", comparison_name, ".csv")),
    row.names = FALSE
  )

  sig_genes <- unique(sig_degs$gene)

  if (length(sig_genes) == 0) {
    empty_df <- data.frame()
    write.csv(empty_df, file.path(comparison_out_dir, paste0("GO_results_",     comparison_name, ".csv")), row.names = FALSE)
    write.csv(empty_df, file.path(comparison_out_dir, paste0("GO_filtered_",    comparison_name, ".csv")), row.names = FALSE)
    write.csv(empty_df, file.path(comparison_out_dir, paste0("KEGG_results_",   comparison_name, ".csv")), row.names = FALSE)
    write.csv(empty_df, file.path(comparison_out_dir, paste0("KEGG_filtered_",  comparison_name, ".csv")), row.names = FALSE)
    return(list(sig_degs = sig_degs, ego = NULL, kk = NULL))
  }

  sig_genes_entrez <- bitr(
    sig_genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )

  cluster_gene <- unique(sig_genes_entrez$ENTREZID)

  # ── GO enrichment ──────────────────────────────────────────────────────────
  ego <- enrichGO(
    gene          = cluster_gene,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE,
    pool          = TRUE
  )

  ego_results <- as.data.frame(ego)

  ego_filtered <- ego_results %>%
    filter(p.adjust < 0.05, Count >= 5, Count <= 500) %>%
    arrange(p.adjust)

  # ── KEGG enrichment ────────────────────────────────────────────────────────
  kk <- enrichKEGG(
    gene          = cluster_gene,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2
  )

  kk_results <- as.data.frame(kk)

  kk_filtered <- kk_results %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust)

  # ── Save results ───────────────────────────────────────────────────────────
  save(
    sig_degs, sig_genes_entrez,
    ego, ego_results, ego_filtered,
    kk, kk_results, kk_filtered,
    file = file.path(comparison_out_dir, paste0("GO_KEGG_results_", comparison_name, ".rda"))
  )

  write.csv(ego_results, file.path(comparison_out_dir, paste0("GO_results_",    comparison_name, ".csv")), row.names = FALSE)
  write.csv(ego_filtered,file.path(comparison_out_dir, paste0("GO_filtered_",   comparison_name, ".csv")), row.names = FALSE)
  write.csv(kk_results,  file.path(comparison_out_dir, paste0("KEGG_results_",  comparison_name, ".csv")), row.names = FALSE)
  write.csv(kk_filtered, file.path(comparison_out_dir, paste0("KEGG_filtered_", comparison_name, ".csv")), row.names = FALSE)

  # ── GO barplot ─────────────────────────────────────────────────────────────
  # x = GeneRatio, color = -log10(p.adjust) gradient
  if (nrow(ego_filtered) > 0) {
    ego_df <- ego_filtered %>%
      group_by(ONTOLOGY) %>%
      slice_min(p.adjust, n = go_top_n, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        GeneRatio_num  = parse_gene_ratio(GeneRatio),
        neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
      ) %>%
      arrange(ONTOLOGY, GeneRatio_num)

    go_barplot <- ggplot(
      ego_df,
      aes(
        x = GeneRatio_num,
        y = fct_reorder(Description, GeneRatio_num)
      )
    ) +
      geom_col(aes(fill = neg_log10_padj), width = 0.72, show.legend = TRUE) +
      scale_fill_gradientn(
        colors = kegg_gradient_colors,
        name   = "-log10(p.adjust)"
      ) +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
      labs(
        title = paste0("Top GO Terms: ", gsub("_", " ", comparison_name)),
        x     = "Gene Ratio",
        y     = NULL
      ) +
      theme_pathway_publication(base_size = 12) +
      theme(
        panel.spacing = unit(0.9, "lines"),
        legend.position   = "right",
        legend.key.height = unit(1, "cm"),
        plot.title = element_text(size = 26, face = "bold",margin = margin(b = 18)),
        
        ## x/y 轴刻度文字大小
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        
        ## x/y 轴标题大小
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        
        ## 图例标题大小，比如 pval、Commun. Prob.
        legend.title = element_text(size = 16, face = "bold"),
        
        ## 图例内容大小，比如 3、min、max
        legend.text = element_text(size = 16),
        legend.margin = margin(b = 8)
      )
    save_pathway_plot(
      file.path(comparison_out_dir, paste0("GO_Top5_Barplot_", comparison_name, ".pdf")),
      go_barplot, width = pathway_plot_width, height = pathway_plot_height
    )
    save_pathway_plot(
      file.path(comparison_out_dir, paste0("GO_Top5_Barplot_", comparison_name, ".png")),
      go_barplot, width = pathway_plot_width, height = pathway_plot_height
    )
  }

  # ── KEGG barplot ───────────────────────────────────────────────────────────
  # x = GeneRatio, color = -log10(p.adjust) gradient
  if (nrow(kk_filtered) > 0) {
    kegg_df <- kk_filtered %>%
      slice_min(p.adjust, n = kegg_top_n, with_ties = FALSE) %>%
      mutate(
        Description_short = gsub(" - Homo sapiens \\(human\\)$", "", Description),
        Description_short = gsub(" pathway$", "", Description_short),
        GeneRatio_num     = parse_gene_ratio(GeneRatio),
        neg_log10_padj    = -log10(pmax(p.adjust, .Machine$double.xmin))
      )

    kegg_barplot <- ggplot(
      kegg_df,
      aes(
        x = GeneRatio_num,
        y = fct_reorder(Description_short, GeneRatio_num)
      )
    ) +
      geom_col(aes(fill = neg_log10_padj), width = 0.72, show.legend = TRUE) +
      scale_fill_gradientn(
        colors = kegg_gradient_colors,
        name   = "-log10(p.adjust)"
      ) +
      labs(
        title = paste0("Top KEGG Pathways: ", gsub("_", " ", comparison_name)),
        x     = "Gene Ratio",
        y     = NULL
      ) +
      theme_pathway_publication(base_size = 12) +
      theme(
        legend.position   = "right",
        legend.key.height = unit(1, "cm"),
        plot.title = element_text(size = 26, face = "bold",margin = margin(b = 18)),
        
        ## x/y 轴刻度文字大小
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        
        ## x/y 轴标题大小
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        
        ## 图例标题大小，比如 pval、Commun. Prob.
        legend.title = element_text(size = 16, face = "bold"),
        
        ## 图例内容大小，比如 3、min、max
        legend.text = element_text(size = 16),
        legend.margin = margin(b = 8)
        
      )

    save_pathway_plot(
      file.path(comparison_out_dir, paste0("KEGG_TopPathways_Barplot_", comparison_name, ".pdf")),
      kegg_barplot, width = pathway_plot_width, height = pathway_plot_height
    )
    save_pathway_plot(
      file.path(comparison_out_dir, paste0("KEGG_TopPathways_Barplot_", comparison_name, ".png")),
      kegg_barplot, width = pathway_plot_width, height = pathway_plot_height
    )
  }

  list(
    sig_degs        = sig_degs,
    sig_genes_entrez = sig_genes_entrez,
    ego             = ego,
    ego_results     = ego_results,
    ego_filtered    = ego_filtered,
    kk              = kk,
    kk_results      = kk_results,
    kk_filtered     = kk_filtered
  )
}

# ==============================================================================
# 运行四组 GO/KEGG 分析
# ==============================================================================

go_kegg_results <- lapply(
  names(deg_files),
  function(x) run_go_kegg_one(
    deg_file        = deg_files[x],
    comparison_name = x,
    out_dir         = go_kegg_out_dir
  )
)

names(go_kegg_results) <- names(deg_files)
