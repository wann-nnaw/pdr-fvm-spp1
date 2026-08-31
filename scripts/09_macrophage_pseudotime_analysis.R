#!/usr/bin/env Rscript

# SPP1 macrophage pseudotime analysis (Monocle 2 / DDRTree)
# Run this script from the DR project root:
#   Rscript codes_GitHub_fixed/09_macrophage_pseudotime_analysis.R

set.seed(20260705)
options(future.globals.maxSize = 2 * 1024^3)

project_library <- file.path(getwd(), ".r-lib")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

# -----------------------------------------------------------------------------
# 1. Project paths and dependency checks
# -----------------------------------------------------------------------------
input_rda <- "rda/UMAP_annotated.rda"
out_dir <- "results/pseudotime"

if (!file.exists("DR.Rproj") || !file.exists(input_rda)) {
  stop(
    "Run this script from the DR project root. Expected files: DR.Rproj and ",
    input_rda,
    call. = FALSE
  )
}

required_packages <- c(
  "Seurat", "monocle", "Biobase", "Matrix", "dplyr", "ggplot2",
  "patchwork", "viridis", "FNN"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install_hint <- if ("monocle" %in% missing_packages) {
    paste0(
      "\nInstall Monocle 2 with:\n",
      "  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
      "install.packages(\"BiocManager\")\n",
      "  BiocManager::install(\"monocle\")"
    )
  } else {
    ""
  }
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    install_hint,
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(Biobase)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

# Monocle 2 calls igraph::graph.dfs(neimode = ...), an argument removed from
# recent igraph releases. Patch only this legacy wrapper while preserving the
# return field (father) expected by Monocle's DDRTree ordering code.
if (packageVersion("igraph") >= "1.3.0") {
  graph_dfs_monocle2 <- function(
      graph, root, mode = c("out", "in", "all", "total"),
      unreachable = TRUE, order = TRUE, order.out = FALSE, father = FALSE,
      dist = FALSE, in.callback = NULL, out.callback = NULL, extra = NULL,
      rho = parent.frame(), neimode) {
    if (!missing(neimode)) {
      mode <- neimode
    }
    result <- igraph::dfs(
      graph = graph,
      root = root,
      mode = mode,
      unreachable = unreachable,
      order = order,
      order.out = order.out,
      parent = father,
      dist = dist,
      in.callback = in.callback,
      out.callback = out.callback,
      extra = extra,
      rho = rho
    )
    result$father <- result$parent
    result
  }
  assignInNamespace("graph.dfs", graph_dfs_monocle2, ns = "igraph")
  monocle_imports <- parent.env(getNamespace("monocle"))
  unlockBinding("graph.dfs", monocle_imports)
  assign("graph.dfs", graph_dfs_monocle2, envir = monocle_imports)
  lockBinding("graph.dfs", monocle_imports)
}

# igraph >= 2.1 removed the vertex-sequence helper nei(). Keep Monocle's
# project2MST logic, replacing only nei() with the stable neighbors() API.
if (packageVersion("igraph") >= "2.1.0") {
  project2mst_monocle2 <- function(cds, Projection_Method) {
    dp_mst <- minSpanningTree(cds)
    Z <- reducedDimS(cds)
    Y <- reducedDimK(cds)
    cds <- findNearestPointOnMST(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex_names <- colnames(Y)[closest_vertex]
    closest_vertex_df <- as.matrix(closest_vertex)
    row.names(closest_vertex_df) <- row.names(closest_vertex)
    tip_leaves <- names(which(igraph::degree(dp_mst) == 1))

    if (!is.function(Projection_Method)) {
      P <- Y[, closest_vertex]
    } else {
      P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
      for (i in seq_along(closest_vertex)) {
        neighbors <- names(igraph::neighbors(
          dp_mst, closest_vertex_names[i], mode = "all"
        ))
        projection <- NULL
        distance <- NULL
        Z_i <- Z[, i]
        for (neighbor in neighbors) {
          if (closest_vertex_names[i] %in% tip_leaves) {
            tmp <- projPointOnLine(
              Z_i, Y[, c(closest_vertex_names[i], neighbor)]
            )
          } else {
            tmp <- Projection_Method(
              Z_i, Y[, c(closest_vertex_names[i], neighbor)]
            )
          }
          projection <- rbind(projection, tmp)
          distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
        }
        if (!methods::is(projection, "matrix")) {
          projection <- as.matrix(projection)
        }
        P[, i] <- projection[which(distance == min(distance))[1], ]
      }
    }

    colnames(P) <- colnames(Z)

    # Build a connected sparse Euclidean neighbor graph. This avoids
    # Monocle's fully connected n^2 graph and is robust to collinear points.
    coords <- t(P)
    cell_names <- rownames(coords)
    n_cells <- nrow(coords)
    k_neighbors <- min(30L, n_cells - 1L)

    repeat {
      knn <- FNN::get.knn(coords, k = k_neighbors, algorithm = "kd_tree")
      neighbor_edges <- data.frame(
        from = rep(cell_names, each = k_neighbors),
        to = cell_names[as.vector(t(knn$nn.index))],
        weight = as.vector(t(knn$nn.dist)),
        stringsAsFactors = FALSE
      )
      edge_key <- paste(
        pmin(neighbor_edges$from, neighbor_edges$to),
        pmax(neighbor_edges$from, neighbor_edges$to),
        sep = "|"
      )
      neighbor_edges <- neighbor_edges[!duplicated(edge_key), , drop = FALSE]
      neighbor_graph <- igraph::graph_from_data_frame(
        neighbor_edges,
        directed = FALSE,
        vertices = data.frame(name = cell_names, stringsAsFactors = FALSE)
      )
      if (igraph::components(neighbor_graph)$no == 1) {
        break
      }
      if (k_neighbors >= n_cells - 1L) {
        stop("Unable to construct a connected projected-cell neighbor graph.")
      }
      k_neighbors <- min(k_neighbors * 2L, n_cells - 1L)
    }

    positive_weights <- neighbor_edges$weight[neighbor_edges$weight > 0]
    if (length(positive_weights) == 0) {
      stop("Projected trajectory coordinates have no positive distances.")
    }
    min_dist <- min(positive_weights)
    igraph::E(neighbor_graph)$weight <-
      igraph::E(neighbor_graph)$weight + min_dist
    dp_mst <- igraph::mst(
      neighbor_graph, weights = igraph::E(neighbor_graph)$weight
    )

    # Monocle's ordering code reads distances only for MST-adjacent cells.
    # Store those distances in the required matrix slot without computing a
    # second all-pairs distance object.
    dp <- matrix(
      0,
      nrow = n_cells,
      ncol = n_cells,
      dimnames = list(cell_names, cell_names)
    )
    mst_edges <- igraph::as_edgelist(dp_mst, names = TRUE)
    mst_weights <- igraph::E(dp_mst)$weight
    dp[cbind(mst_edges[, 1], mst_edges[, 2])] <- mst_weights
    dp[cbind(mst_edges[, 2], mst_edges[, 1])] <- mst_weights
    cellPairwiseDistances(cds) <- dp
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
    cds
  }
  environment(project2mst_monocle2) <- getNamespace("monocle")
  assignInNamespace("project2MST", project2mst_monocle2, ns = "monocle")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subtype_levels <- c("IAM", "LAM", "OSM", "TRM")
# Genes shown as mechanistic readouts along the trajectory. They combine SPP1,
# macrophage-state markers and upstream TF candidates tested in the TF workflow.
dynamic_genes_requested <- c(
  "SPP1", "GPNMB", "C3", "FCGR2B", "CLEC5A", "CD81",
  "JUN", "FOS", "JUNB", "CEBPB", "NFKB1", "RELA", "KLF4", "MAFB"
)
subtype_colors <- c(
  IAM = "#B65C4A",
  LAM = "#2F6F73",
  OSM = "#6E5F8F",
  TRM = "#C49A3A"
)

theme_pseudotime <- function(base_size = 22) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.45, colour = "black"),
      axis.text = element_text(size = base_size - 2, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black", face = "bold"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 17),
      legend.key.height = grid::unit(15, "pt"),
      legend.key.width = grid::unit(20, "pt"),
      legend.spacing.x = grid::unit(8, "pt"),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(
        size = base_size + 2, face = "bold", hjust = 0.5,
        margin = margin(b = 10)
      ),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_plot_pair <- function(plot, filename, width, height, dpi = 600) {
  ggsave(
    file.path(out_dir, paste0(filename, ".pdf")), plot,
    width = width, height = height, units = "in", bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(filename, ".png")), plot,
    width = width, height = height, units = "in", dpi = dpi, bg = "white"
  )
}

# -----------------------------------------------------------------------------
# 2. Load and validate the four macrophage subtypes
# -----------------------------------------------------------------------------
load(input_rda)
if (!exists("combined_anno") || !inherits(combined_anno, "Seurat")) {
  stop("The input RDA must contain a Seurat object named combined_anno.", call. = FALSE)
}
if (!"cell_type" %in% colnames(combined_anno[[]])) {
  stop("combined_anno metadata does not contain cell_type.", call. = FALSE)
}
if (!"sample" %in% colnames(combined_anno[[]])) {
  stop("combined_anno metadata does not contain sample.", call. = FALSE)
}
if (!"RNA" %in% Assays(combined_anno)) {
  stop("combined_anno does not contain an RNA assay.", call. = FALSE)
}

observed_subtypes <- unique(as.character(combined_anno$cell_type))
missing_subtypes <- setdiff(subtype_levels, observed_subtypes)
if (length(missing_subtypes) > 0) {
  stop(
    "Missing macrophage subtype(s): ", paste(missing_subtypes, collapse = ", "),
    call. = FALSE
  )
}

macrophage <- subset(combined_anno, subset = cell_type %in% subtype_levels)
macrophage$cell_type <- factor(as.character(macrophage$cell_type), levels = subtype_levels)
Idents(macrophage) <- macrophage$cell_type

subtype_counts <- table(macrophage$cell_type)
if (any(subtype_counts == 0)) {
  stop("At least one macrophage subtype contains no cells.", call. = FALSE)
}

rna_counts <- GetAssayData(macrophage, assay = "RNA", layer = "counts")
if (nrow(rna_counts) == 0 || ncol(rna_counts) == 0) {
  stop("The RNA counts layer is empty.", call. = FALSE)
}
if (!"SPP1" %in% rownames(rna_counts)) {
  stop("SPP1 is absent from the RNA counts matrix.", call. = FALSE)
}
if (all(rna_counts["SPP1", ] == 0)) {
  stop("SPP1 has zero counts in all macrophage cells.", call. = FALSE)
}

message(
  "Using ", ncol(macrophage), " macrophages: ",
  paste(names(subtype_counts), as.integer(subtype_counts), sep = "=", collapse = ", ")
)

# -----------------------------------------------------------------------------
# 3. Identify subtype markers for Monocle ordering genes
# -----------------------------------------------------------------------------
if (!"SCT" %in% Assays(macrophage)) {
  stop("The SCT assay required for marker selection is absent.", call. = FALSE)
}
DefaultAssay(macrophage) <- "SCT"
macrophage <- PrepSCTFindMarkers(macrophage, assay = "SCT", verbose = FALSE)

markers <- FindAllMarkers(
  macrophage,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  return.thresh = 0.05,
  verbose = FALSE
)

required_marker_columns <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
if (!all(required_marker_columns %in% colnames(markers))) {
  stop(
    "FindAllMarkers output is missing required columns: ",
    paste(setdiff(required_marker_columns, colnames(markers)), collapse = ", "),
    call. = FALSE
  )
}

top_markers <- markers %>%
  filter(
    p_val_adj < 0.05,
    is.finite(avg_log2FC),
    gene %in% rownames(rna_counts)
  ) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup()

marker_counts <- table(factor(top_markers$cluster, levels = subtype_levels))
if (any(marker_counts < 20)) {
  stop(
    "Too few significant ordering markers for: ",
    paste(names(marker_counts)[marker_counts < 20], collapse = ", "),
    ". At least 20 markers per subtype are required.",
    call. = FALSE
  )
}

ordering_genes <- unique(top_markers$gene)
if (length(ordering_genes) < 50) {
  stop("Fewer than 50 unique ordering genes were identified.", call. = FALSE)
}

write.csv(
  as.data.frame(top_markers),
  file.path(out_dir, "ordering_genes.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Build the Monocle 2 trajectory and orient it from the TRM-enriched state
# -----------------------------------------------------------------------------
cell_metadata <- macrophage[[]]
cell_metadata$cell_type <- factor(
  as.character(cell_metadata$cell_type), levels = subtype_levels
)
cell_metadata$cell_id <- rownames(cell_metadata)

dynamic_genes <- intersect(dynamic_genes_requested, rownames(rna_counts))
trajectory_genes <- unique(c(ordering_genes, dynamic_genes))
trajectory_counts <- rna_counts[trajectory_genes, , drop = FALSE]
feature_metadata <- data.frame(
  gene_short_name = rownames(trajectory_counts),
  row.names = rownames(trajectory_counts),
  stringsAsFactors = FALSE
)

pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = feature_metadata)

cds <- newCellDataSet(
  trajectory_counts,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

# The trajectory needs only its ordering genes and SPP1. Release the full
# Seurat object before DDRTree builds dense cell-cell distance matrices.
rm(
  combined_anno, macrophage, rna_counts, trajectory_counts,
  markers, top_markers
)
invisible(gc())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree",
  norm_method = "log",
  pseudo_expr = 1,
  verbose = FALSE
)
cds <- orderCells(cds)

state_subtype_counts <- as.data.frame(
  table(
    State = as.character(pData(cds)$State),
    cell_type = as.character(pData(cds)$cell_type)
  ),
  stringsAsFactors = FALSE
)
state_summary <- state_subtype_counts %>%
  group_by(State) %>%
  mutate(
    state_total = sum(Freq),
    subtype_fraction = ifelse(state_total > 0, Freq / state_total, 0)
  ) %>%
  ungroup()

trm_state_summary <- state_summary %>%
  filter(cell_type == "TRM") %>%
  arrange(desc(subtype_fraction), desc(Freq), as.numeric(State))

if (nrow(trm_state_summary) == 0 || max(trm_state_summary$Freq) == 0) {
  stop("No TRM cells were assigned to a DDRTree state.", call. = FALSE)
}

root_state <- as.character(trm_state_summary$State[[1]])
cds <- orderCells(cds, root_state = root_state)

write.csv(
  state_summary,
  file.path(out_dir, "state_subtype_composition.csv"),
  row.names = FALSE
)
writeLines(
  c(
    paste0("root_state=", root_state),
    paste0(
      "TRM_fraction=",
      signif(trm_state_summary$subtype_fraction[[1]], 5)
    ),
    paste0("TRM_cells=", trm_state_summary$Freq[[1]])
  ),
  file.path(out_dir, "root_state_summary.txt")
)

# -----------------------------------------------------------------------------
# 5. Return pseudotime to Seurat and save source data
# -----------------------------------------------------------------------------
trajectory_metadata <- pData(cds) %>%
  as.data.frame() %>%
  mutate(
    cell_id = rownames(.),
    Pseudotime = as.numeric(Pseudotime),
    State = as.character(State),
    cell_type = factor(as.character(cell_type), levels = subtype_levels)
  )
rownames(trajectory_metadata) <- trajectory_metadata$cell_id

if (all(!is.finite(trajectory_metadata$Pseudotime))) {
  stop("All pseudotime values are non-finite.", call. = FALSE)
}
if (any(!is.finite(trajectory_metadata$Pseudotime))) {
  stop("Some macrophage cells have non-finite pseudotime values.", call. = FALSE)
}

size_factors <- as.numeric(pData(cds)$Size_Factor)
if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Monocle generated invalid size factors.", call. = FALSE)
}
spp1_counts <- as.numeric(exprs(cds)["SPP1", ])
spp1_normalized <- log1p(spp1_counts / size_factors)
pData(cds)$SPP1_expression <- spp1_normalized

trajectory_metadata$SPP1_counts <- spp1_counts
trajectory_metadata$SPP1_expression <- spp1_normalized
trajectory_metadata$is_root_state <- trajectory_metadata$State == root_state

# Reload the annotated Seurat object only after trajectory inference, then
# write the aligned Monocle metadata back to the four macrophage subtypes.
load(input_rda)
macrophage <- subset(combined_anno, subset = cell_type %in% subtype_levels)
macrophage$cell_type <- factor(
  as.character(macrophage$cell_type), levels = subtype_levels
)
rm(combined_anno)
invisible(gc())

macrophage$Pseudotime <- trajectory_metadata[colnames(macrophage), "Pseudotime"]
macrophage$Monocle_State <- trajectory_metadata[colnames(macrophage), "State"]
macrophage$SPP1_monocle_expression <- trajectory_metadata[
  colnames(macrophage), "SPP1_expression"
]

saveRDS(cds, file.path(out_dir, "macrophage_monocle2_cds.rds"))
saveRDS(macrophage, file.path(out_dir, "macrophage_with_pseudotime.rds"))
write.csv(
  trajectory_metadata,
  file.path(out_dir, "macrophage_pseudotime_metadata.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Trajectory overview and SPP1 trajectory
# -----------------------------------------------------------------------------
trajectory_theme <- theme_pseudotime(base_size = 22) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "center",
    legend.margin = margin(0, 0, 5, 0),
    aspect.ratio = 1
  )

p_time <- plot_cell_trajectory(
  cds, color_by = "Pseudotime", show_backbone = TRUE, cell_size = 0.72
) +
  scale_color_viridis_c(option = "C", end = 0.95, name = "Pseudotime") +
  labs(title = "Macrophage pseudotime") +
  trajectory_theme

p_state <- plot_cell_trajectory(
  cds, color_by = "State", show_backbone = TRUE, cell_size = 0.72
) +
  labs(title = paste0("DDRTree state (root: ", root_state, ")")) +
  trajectory_theme

p_subtype <- plot_cell_trajectory(
  cds, color_by = "cell_type", show_backbone = TRUE, cell_size = 0.72
) +
  scale_color_manual(values = subtype_colors, drop = FALSE, name = "Subtype") +
  labs(title = "Macrophage subtype") +
  trajectory_theme

p_sample <- plot_cell_trajectory(
  cds, color_by = "sample", show_backbone = TRUE, cell_size = 0.72
) +
  labs(title = "Sample") +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, title.position = "left")
  ) +
  trajectory_theme +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.width = grid::unit(11, "pt"),
    legend.spacing.x = grid::unit(4, "pt")
  )

# Export the four trajectory views separately. Keeping each legend in its own
# canvas prevents long subtype/sample legends from overlapping or being clipped.
save_plot_pair(
  p_time, "trajectory_pseudotime", width = 7.2, height = 6.6
)
save_plot_pair(
  p_state, "trajectory_state", width = 7.2, height = 6.6
)
save_plot_pair(
  p_subtype, "trajectory_macrophage_subtype", width = 7.2, height = 6.6
)
save_plot_pair(
  p_sample, "trajectory_sample", width = 7.2, height = 6.6
)

p_spp1_trajectory <- plot_cell_trajectory(
  cds,
  color_by = "SPP1_expression",
  show_backbone = TRUE,
  cell_size = 0.78
) +
  scale_color_viridis_c(
    option = "B", end = 0.95, name = "SPP1\nlog-normalized"
  ) +
  labs(title = "SPP1 expression along\nmacrophage trajectory") +
  trajectory_theme

save_plot_pair(
  p_spp1_trajectory, "SPP1_trajectory", width = 7.2, height = 6.6
)

# -----------------------------------------------------------------------------
# 7. SPP1 trend and subtype pseudotime distributions
# -----------------------------------------------------------------------------
plot_df <- trajectory_metadata %>%
  select(
    cell_id, sample, cell_type, Pseudotime, State,
    SPP1_counts, SPP1_expression, is_root_state
  )

p_spp1_trend <- ggplot(
  plot_df,
  aes(Pseudotime, SPP1_expression, colour = cell_type)
) +
  geom_point(size = 0.38, alpha = 0.09, stroke = 0) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    linewidth = 1.15,
    alpha = 0.16
  ) +
  scale_colour_manual(values = subtype_colors, drop = FALSE) +
  labs(
    title = "SPP1 expression across\nmacrophage pseudotime",
    x = "Pseudotime",
    y = "SPP1 expression (log-normalized counts)",
    colour = "Subtype"
  ) +
  theme_pseudotime(base_size = 22) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    aspect.ratio = 1,
    plot.title = element_text(
      size = 24, face = "bold", hjust = 0.5,
      margin = margin(b = 8)
    ),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.width = grid::unit(14, "pt"),
    legend.spacing.x = grid::unit(5, "pt")
  )

save_plot_pair(
  p_spp1_trend, "SPP1_pseudotime_trend", width = 7.2, height = 6.6
)

p_subtype_distribution <- ggplot(
  plot_df,
  aes(cell_type, Pseudotime, fill = cell_type)
) +
  geom_violin(width = 0.85, trim = FALSE, alpha = 0.80, colour = "black", linewidth = 0.3) +
  geom_boxplot(
    width = 0.16, outlier.shape = NA, fill = "white",
    colour = "black", linewidth = 0.35
  ) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  labs(
    title = "Pseudotime distribution by macrophage subtype",
    x = NULL,
    y = "Pseudotime"
  ) +
  theme_pseudotime(base_size = 10) +
  theme(legend.position = "none")

save_plot_pair(
  p_subtype_distribution,
  "subtype_pseudotime_distribution",
  width = 6.2,
  height = 5.2
)

p_subtype_density <- ggplot(
  plot_df,
  aes(Pseudotime, colour = cell_type, fill = cell_type)
) +
  geom_density(linewidth = 0.85, alpha = 0.16, adjust = 1) +
  scale_colour_manual(values = subtype_colors, drop = FALSE) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  labs(
    title = "Macrophage subtype density along pseudotime",
    x = "Pseudotime",
    y = "Density",
    colour = "Subtype",
    fill = "Subtype"
  ) +
  theme_pseudotime(base_size = 10) +
  theme(legend.position = "top")

save_plot_pair(
  p_subtype_density, "subtype_pseudotime_density", width = 7.2, height = 5
)

# -----------------------------------------------------------------------------
# 8. Multi-gene dynamics, State composition and branch-dependent expression
# -----------------------------------------------------------------------------
dynamic_counts <- exprs(cds)[dynamic_genes, , drop = FALSE]
dynamic_long <- as.data.frame(t(dynamic_counts)) %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(Size_Factor = size_factors) %>%
  pivot_longer(
    cols = all_of(dynamic_genes), names_to = "gene", values_to = "counts"
  ) %>%
  mutate(expression = log1p(counts / Size_Factor)) %>%
  left_join(
    trajectory_metadata %>% select(cell_id, Pseudotime, State, cell_type),
    by = "cell_id"
  )
write.csv(
  dynamic_long,
  file.path(out_dir, "dynamic_gene_pseudotime_source_data.csv"),
  row.names = FALSE
)

p_dynamic <- ggplot(dynamic_long, aes(Pseudotime, expression)) +
  geom_point(aes(colour = cell_type), size = 0.22, alpha = 0.08, stroke = 0) +
  geom_smooth(
    aes(colour = cell_type), method = "gam", formula = y ~ s(x, bs = "cs"),
    se = FALSE, linewidth = 0.65
  ) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = subtype_colors, drop = FALSE) +
  labs(
    title = "Macrophage-state genes and candidate TFs along pseudotime",
    subtitle = "Curves are subtype-specific GAM smooths; points represent cells",
    x = "Pseudotime", y = "Log-normalized expression", colour = "Subtype"
  ) +
  theme_pseudotime(base_size = 8) +
  theme(legend.position = "top", strip.background = element_blank())
save_plot_pair(p_dynamic, "dynamic_genes_pseudotime", width = 10.5, height = 8.2)

state_plot_df <- state_summary %>%
  mutate(
    State = factor(State, levels = sort(unique(as.numeric(State)))),
    cell_type = factor(cell_type, levels = subtype_levels)
  )
p_state_composition <- ggplot(
  state_plot_df, aes(State, subtype_fraction * 100, fill = cell_type)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Macrophage subtype composition of trajectory states",
    subtitle = paste0("State ", root_state, " was selected as the TRM-enriched root"),
    x = "DDRTree State", y = "Cells within State (%)", fill = "Subtype"
  ) +
  theme_pseudotime(base_size = 9) +
  theme(legend.position = "top")
save_plot_pair(p_state_composition, "state_subtype_composition", width = 6.8, height = 4.8)

# Monocle branch heatmap is generated only when a valid branch point exists.
# Failure here must not invalidate the principal trajectory analysis.
branch_points <- tryCatch({
  branch_nodes <- names(which(igraph::degree(minSpanningTree(cds)) > 2))
  seq_along(branch_nodes)
}, error = function(e) integer())
branch_genes <- intersect(dynamic_genes, rownames(cds))
branch_status <- "not generated: no valid branch point"
if (length(branch_points) > 0 && length(branch_genes) >= 2) {
  branch_pdf <- file.path(out_dir, "dynamic_genes_branch_heatmap.pdf")
  branch_status <- tryCatch({
    grDevices::pdf(branch_pdf, width = 7.2, height = 6.2)
    plot_genes_branched_heatmap(
      cds[branch_genes, ], branch_point = branch_points[[1]],
      branch_labels = c("Cell fate 1", "Cell fate 2"),
      cluster_rows = TRUE, num_clusters = min(4, length(branch_genes)),
      show_rownames = TRUE, scale_min = -2, scale_max = 2, cores = 1
    )
    grDevices::dev.off()
    "generated"
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    if (file.exists(branch_pdf)) unlink(branch_pdf)
    paste0("not generated: ", conditionMessage(e))
  })
}
writeLines(
  paste0("branch_heatmap=", branch_status),
  file.path(out_dir, "branch_heatmap_status.txt")
)

# -----------------------------------------------------------------------------
# 9. Macrophage-subtype marker dynamics for the Figure 6 lower panel
# -----------------------------------------------------------------------------
marker_sets <- list(
  LAM = c("LGALS3", "ACP5", "PLA2G7", "GPNMB", "CTSD"),
  IAM = c("CCL4L2", "CCL3L1", "CCL3", "CCL4", "CD83"),
  OSM = c("TXNRD1", "SCD", "GCLM", "SQSTM1", "CTSB"),
  TRM = c("STAB1", "MAF", "DAB2", "F13A1", "MRC1")
)
all_markers <- unique(unlist(marker_sets, use.names = FALSE))
missing_marker_genes <- setdiff(all_markers, rownames(macrophage))
if (length(missing_marker_genes) > 0) {
  stop(
    "Marker gene(s) absent from the Seurat object: ",
    paste(missing_marker_genes, collapse = ", "),
    call. = FALSE
  )
}

DefaultAssay(macrophage) <- "SCT"
marker_expression_df <- FetchData(
  macrophage,
  vars = c("Pseudotime", "cell_type", all_markers),
  layer = "data"
) %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(cell_type = as.character(cell_type)) %>%
  pivot_longer(
    cols = all_of(all_markers), names_to = "gene", values_to = "expression"
  )

marker_key <- tibble(
  marker_subtype = rep(names(marker_sets), lengths(marker_sets)),
  gene = unlist(marker_sets, use.names = FALSE),
  marker_order = sequence(lengths(marker_sets))
)

marker_plot_df <- marker_expression_df %>%
  inner_join(marker_key, by = "gene") %>%
  filter(
    cell_type == marker_subtype,
    is.finite(Pseudotime),
    is.finite(expression)
  ) %>%
  mutate(
    marker_subtype = factor(marker_subtype, levels = names(marker_sets)),
    gene = factor(gene, levels = unlist(marker_sets, use.names = FALSE))
  )

if (nrow(marker_plot_df) == 0) {
  stop("No cells remained after matching marker sets to macrophage subtypes.")
}

write.csv(
  marker_plot_df %>% mutate(gene = as.character(gene)),
  file.path(out_dir, "macrophage_marker_pseudotime_source_data.csv"),
  row.names = FALSE
)

theme_marker_pseudotime <- theme_classic(base_size = 12.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.text = element_text(size = 11.5, colour = "black"),
    axis.title = element_text(size = 12.5, face = "bold", colour = "black"),
    strip.background = element_rect(
      fill = "#F2F2F2", colour = "#BDBDBD", linewidth = 0.4
    ),
    strip.text = element_text(size = 12, face = "bold", colour = "black"),
    plot.title = element_text(
      size = 14, face = "bold", hjust = 0,
      margin = margin(b = 4)
    ),
    plot.margin = margin(7, 8, 7, 8),
    panel.spacing.x = grid::unit(9, "pt")
  )

make_marker_subtype_row <- function(subtype) {
  genes <- marker_sets[[subtype]]
  row_df <- marker_plot_df %>%
    filter(marker_subtype == subtype) %>%
    mutate(gene = factor(as.character(gene), levels = genes))

  ggplot(row_df, aes(Pseudotime, expression)) +
    geom_point(
      colour = subtype_colors[[subtype]], alpha = 0.065,
      size = 0.24, stroke = 0
    ) +
    geom_smooth(
      method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE,
      colour = subtype_colors[[subtype]], fill = subtype_colors[[subtype]],
      linewidth = 1.05, alpha = 0.16
    ) +
    facet_wrap(~gene, nrow = 1, scales = "free_y", drop = FALSE) +
    labs(
      title = subtype,
      x = if (subtype == tail(names(marker_sets), 1)) "Pseudotime" else NULL,
      y = "Log-normalized\nexpression"
    ) +
    theme_marker_pseudotime +
    theme(
      axis.text.x = if (subtype == tail(names(marker_sets), 1)) {
        element_text(size = 11.5)
      } else {
        element_blank()
      },
      axis.ticks.x = if (subtype == tail(names(marker_sets), 1)) {
        element_line(linewidth = 0.4)
      } else {
        element_blank()
      }
    )
}

marker_row_plots <- lapply(names(marker_sets), make_marker_subtype_row)
marker_pseudotime_plot <- wrap_plots(marker_row_plots, ncol = 1) +
  plot_layout(heights = rep(1, length(marker_row_plots))) +
  plot_annotation(
    title = "Macrophage subtype marker dynamics along pseudotime",
    theme = theme(
      plot.title = element_text(
        size = 19, face = "bold", hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

save_plot_pair(
  marker_pseudotime_plot,
  "macrophage_marker_pseudotime_4x5",
  width = 11.2,
  height = 10.6
)

# -----------------------------------------------------------------------------
# 10. Final output checks
# -----------------------------------------------------------------------------
expected_files <- c(
  "macrophage_monocle2_cds.rds",
  "macrophage_with_pseudotime.rds",
  "macrophage_pseudotime_metadata.csv",
  "ordering_genes.csv",
  "trajectory_pseudotime.pdf",
  "trajectory_pseudotime.png",
  "trajectory_state.pdf",
  "trajectory_state.png",
  "trajectory_macrophage_subtype.pdf",
  "trajectory_macrophage_subtype.png",
  "trajectory_sample.pdf",
  "trajectory_sample.png",
  "SPP1_trajectory.pdf",
  "SPP1_trajectory.png",
  "SPP1_pseudotime_trend.pdf",
  "SPP1_pseudotime_trend.png",
  "subtype_pseudotime_distribution.pdf",
  "subtype_pseudotime_distribution.png",
  "subtype_pseudotime_density.pdf",
  "subtype_pseudotime_density.png",
  "dynamic_genes_pseudotime.pdf",
  "dynamic_genes_pseudotime.png",
  "state_subtype_composition.pdf",
  "state_subtype_composition.png",
  "dynamic_gene_pseudotime_source_data.csv",
  "branch_heatmap_status.txt",
  "macrophage_marker_pseudotime_4x5.pdf",
  "macrophage_marker_pseudotime_4x5.png",
  "macrophage_marker_pseudotime_source_data.csv"
)
expected_paths <- file.path(out_dir, expected_files)
missing_outputs <- expected_paths[!file.exists(expected_paths)]
empty_outputs <- expected_paths[file.exists(expected_paths) & file.info(expected_paths)$size <= 0]

if (length(missing_outputs) > 0 || length(empty_outputs) > 0) {
  stop(
    "Output validation failed. Missing: ",
    paste(basename(missing_outputs), collapse = ", "),
    "; empty: ", paste(basename(empty_outputs), collapse = ", "),
    call. = FALSE
  )
}

message(
  "Pseudotime analysis completed. Root state: ", root_state,
  ". Results: ", out_dir
)

# -----------------------------------------------------------------------------
# 审稿补充分析：Monocle 2根节点选择敏感性
# -----------------------------------------------------------------------------

review_pt_dir <- "results/reviewer_supplement/pseudotime_root_sensitivity"
dir.create(review_pt_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cds")) {
  cds <- readRDS(file.path(out_dir, "macrophage_monocle2_cds.rds"))
}
if (!exists("macrophage_pt")) {
  macrophage_pt <- readRDS(file.path(out_dir, "macrophage_with_pseudotime.rds"))
}

review_pd <- pData(cds)
review_pd$cell <- rownames(review_pd)
review_pd$State <- as.character(review_pd$State)
review_pd$cell_type <- as.character(review_pd$cell_type)

review_state_comp <- as.data.frame(
  table(review_pd$State, review_pd$cell_type),
  stringsAsFactors = FALSE
)
colnames(review_state_comp) <- c("State", "cell_type", "Freq")
review_state_comp <- review_state_comp %>%
  group_by(State) %>%
  mutate(
    state_total = sum(Freq),
    subtype_fraction = Freq / state_total
  ) %>%
  ungroup()
write.csv(
  review_state_comp,
  file.path(review_pt_dir, "root_sensitivity_state_subtype_composition.csv"),
  row.names = FALSE
)

review_state_root <- review_state_comp %>%
  filter(cell_type == "TRM") %>%
  arrange(desc(subtype_fraction), desc(Freq)) %>%
  slice(1) %>%
  mutate(root_reason = "highest_TRM_fraction")
write.csv(
  review_state_root,
  file.path(review_pt_dir, "root_state_data_driven_selection.csv"),
  row.names = FALSE
)

review_expr <- FetchData(macrophage_pt, vars = "SPP1")[, 1]
names(review_expr) <- colnames(macrophage_pt)
review_original_pt <- review_pd$Pseudotime
names(review_original_pt) <- review_pd$cell
review_states <- sort(unique(review_pd$State))

review_cell_embedding <- as.data.frame(t(reducedDimS(cds)))
colnames(review_cell_embedding) <- c("DDRTree_1", "DDRTree_2")
review_cell_embedding$cell <- rownames(review_cell_embedding)
review_cell_embedding$State <- review_pd[review_cell_embedding$cell, "State"]
review_state_centroids <- review_cell_embedding %>%
  group_by(State) %>%
  summarise(
    DDRTree_1 = mean(DDRTree_1),
    DDRTree_2 = mean(DDRTree_2),
    .groups = "drop"
  )

review_reorder_one <- function(state_id) {
  root_xy <- review_state_centroids %>% filter(State == state_id)
  state_dist <- review_state_centroids %>%
    mutate(
      root_state = state_id,
      state_root_distance = sqrt(
        (DDRTree_1 - root_xy$DDRTree_1[1])^2 +
          (DDRTree_2 - root_xy$DDRTree_2[1])^2
      )
    ) %>%
    select(State, state_root_distance)
  approx_df <- review_pd %>%
    select(cell, State, cell_type) %>%
    left_join(state_dist, by = "State") %>%
    mutate(
      within_state_rank = ave(
        review_original_pt[cell], State,
        FUN = function(x) rank(x, ties.method = "average") / length(x)
      ),
      approximate_pseudotime = state_root_distance + within_state_rank * 1e-3
    )
  shared <- intersect(names(review_original_pt), approx_df$cell)
  approx_pt <- approx_df$approximate_pseudotime[match(shared, approx_df$cell)]
  subtype_medians <- approx_df %>%
    group_by(cell_type) %>%
    summarise(
      median_pseudotime = median(approximate_pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_pseudotime)
  data.frame(
    root_state = state_id,
    status = "state-centroid approximation",
    spearman_vs_state4 = suppressWarnings(cor(
      review_original_pt[shared], approx_pt, method = "spearman",
      use = "complete.obs"
    )),
    spp1_lm_slope = unname(coef(lm(review_expr[shared] ~ approx_pt))[2]),
    subtype_median_order = paste(subtype_medians$cell_type, collapse = " > "),
    stringsAsFactors = FALSE
  )
}

review_root_sensitivity <- lapply(review_states, review_reorder_one) %>% bind_rows()
write.csv(
  review_root_sensitivity,
  file.path(review_pt_dir, "root_state_pseudotime_sensitivity.csv"),
  row.names = FALSE
)

p_review_state_comp <- ggplot(
  review_state_comp,
  aes(x = State, y = subtype_fraction, fill = cell_type)
) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "DDRTree state composition used for root selection",
    x = "DDRTree state",
    y = "Subtype fraction",
    fill = "Subtype"
  ) +
  theme_classic(base_size = 11)

p_review_root_cor <- ggplot(
  review_root_sensitivity,
  aes(x = root_state, y = spearman_vs_state4, fill = root_state == "4")
) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("TRUE" = "#B9473D", "FALSE" = "#8EA7B8"), guide = "none") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Pseudotime rank concordance under alternative roots",
    x = "Alternative root state",
    y = "Spearman rho versus state 4 root"
  ) +
  theme_classic(base_size = 11)

ggsave(file.path(review_pt_dir, "root_state_subtype_composition.pdf"),
       p_review_state_comp, width = 6.5, height = 4.8, bg = "white")
ggsave(file.path(review_pt_dir, "root_state_subtype_composition.png"),
       p_review_state_comp, width = 6.5, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(review_pt_dir, "root_state_spearman_sensitivity.pdf"),
       p_review_root_cor, width = 5.8, height = 4.2, bg = "white")
ggsave(file.path(review_pt_dir, "root_state_spearman_sensitivity.png"),
       p_review_root_cor, width = 5.8, height = 4.2, dpi = 300, bg = "white")

writeLines(
  c(
    "Reviewer supplement: Monocle 2 root-state sensitivity",
    paste0("selected_root_state=", review_state_root$State[1]),
    paste0("selected_root_TRM_fraction=", signif(review_state_root$subtype_fraction[1], 5)),
    paste0("selected_root_TRM_cells=", review_state_root$Freq[1]),
    "Alternative-root analyses are sensitivity checks for transcriptional ordering, not lineage validation."
  ),
  file.path(review_pt_dir, "pseudotime_root_sensitivity_summary.txt")
)
