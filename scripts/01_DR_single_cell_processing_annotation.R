
library(Matrix)
library(data.table)
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(GEOquery)
library(readxl)
library(cowplot)
library(scales)
library(RColorBrewer)

# =============================================================================
# 单细胞数据处理及细胞注释
# =============================================================================

dir.create("results/Marker", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tsne", showWarnings = FALSE, recursive = TRUE)
dir.create("results/QC", showWarnings = FALSE, recursive = TRUE)
dir.create("rda", showWarnings = FALSE)

# =============================================================================
# 1. 下载并读取单细胞数据
# =============================================================================

# (1) GSE165784 (解压后排除GSM5049904_RRD-ERM1_matrix.tsv.gz) # GPL20795
path_1 <- "./data/GSE165784_data/"

# 获取所有tsv文件
files_1 <- list.files(path_1, pattern = "tsv.gz", full.names = TRUE)
seurat_list_1 <- list()

for (f in files_1) {
  # 提取样本名（前10位）
  sample_name <- substr(basename(f), 1, 10)
  
  # 读取矩阵
  mat_1 <- fread(f, data.table = FALSE)
  rownames(mat_1) <- mat_1[[1]]
  mat_1 <- as.matrix(mat_1[, -1])
  
  # 创建Seurat对象
  obj_1 <- CreateSeuratObject(counts = mat_1, min.cells = 3, min.features = 200)
  obj_1$sample <- sample_name
  
  # 存入列表，列表名即为样本名
  seurat_list_1[[sample_name]] <- obj_1
}

# GSM5690478
mat_478 <- Read10X("./data/GSE165784_data/PDR-FM-0609_matrix_10X/")
obj_478 <- CreateSeuratObject(mat_478)
obj_478$sample <- "GSM5690478"

# GSM5690479
mat_479 <- Read10X("./data/GSE165784_data/PDR-ERM-210630_matrix_10X/")
obj_479 <- CreateSeuratObject(mat_479)
obj_479$sample <- "GSM5690479"

seurat_list_1[["GSM5690478"]] <- obj_478
seurat_list_1[["GSM5690479"]] <- obj_479

# =============================================================================
# 2. 质控
# =============================================================================

# 添加线粒体基因表达占比
seurat_list_1 <- lapply(seurat_list_1, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  return(obj)
})

# QC函数
do_QC <- function(obj, name = "") {
  cat("Processing:", name, "- cells before QC:", ncol(obj), "\n")
  obj <- subset(obj,
                subset = nCount_RNA > 1000 &
                  nCount_RNA < 50000 &
                  nFeature_RNA > 200 &
                  nFeature_RNA < 6000 &
                  percent.mt < 15)
  
  cat("Processing:", name, "- cells after QC:", ncol(obj), "\n\n")
  return(obj)
}

seurat_QC_1 <- mapply(do_QC, obj = seurat_list_1, name = names(seurat_list_1), SIMPLIFY = FALSE)

# 标准化SCTransform（仅分析GSE165784，保留GSE165784内部样本间SCT整合）
all_objects <- seurat_QC_1
all_objects <- lapply(all_objects, SCTransform, verbose = FALSE)

save(all_objects, file = "rda/SCTransform.rda")

# 筛选高变基因
set.seed(13)
features <- SelectIntegrationFeatures(all_objects, nfeatures = 3000)

# 整合数据
all_objects <- PrepSCTIntegration(all_objects, anchor.features = features) # 准备整合
anchors <- FindIntegrationAnchors(all_objects, normalization.method = "SCT", anchor.features = features) # 找锚点
combined <- IntegrateData(anchors, normalization.method = "SCT") # 整合

# 准备差异表达
combined_1 <- PrepSCTFindMarkers(combined)
DefaultAssay(combined_1) <- "integrated" # 后续使用 integrated assay

# PCA
combined_PCA <- RunPCA(combined_1, npcs = 50, verbose = FALSE)
ElbowPlot(combined_PCA)   

# 细胞聚类
combined_Neigh <- FindNeighbors(combined_PCA, dims = 1:10)
combined_cluster <- FindClusters(combined_Neigh, resolution = 0.5)

# 可视化
combined_tsne <- RunTSNE(combined_cluster, dims = 1:10)
table(combined_tsne$seurat_clusters)
DimPlot(combined_tsne, reduction = "tsne", label = TRUE)

save(combined_tsne, file = "rda/tsne.rda")

# load("rda/tsne.rda")

# =============================================================================
# 3. 细胞注释
# =============================================================================

# 3.1 FindAllMarkers
DefaultAssay(combined_tsne) <- "SCT"
markers <- FindAllMarkers(combined_tsne, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(markers, "results/tsne/All_cluster_markers.csv")
write.csv(top10, "results/tsne/Top10_markers_each_cluster.csv")

DotPlot(combined_tsne, features = unique(top10$gene)) + RotatedAxis()
top10$gene

# 3.2 手动注释(PMID: 40069725)
# 合并后的细胞类型顺序、marker 与配色作为所有下游图形的唯一配置。
cell_type_order <- c(
  "Macrophage", "Monocyte", "Microglia", "DC",
  "Endothelial", "Fibroblast", "Pericyte", "T cells"
)

marker_panel <- list(
  "Macrophage"  = c("CD68", "CD163", "MRC1", "MSR1", "C1QA", "C1QB", "C1QC", "APOE"),
  "Monocyte"    = c("FCN1", "S100A8", "S100A9", "S100A12", "CD300E"),
  "Microglia"   = c("P2RY12", "CX3CR1", "TMEM119"),
  "DC"          = c("NAPSB", "FCER1A", "CD1C"),
  "Endothelial" = c("CLDN5", "VWF", "SPARCL1"),
  "Fibroblast"  = c("LUM", "THBS2", "DCN"),
  "Pericyte"    = c("PDGFRB", "RGS5", "CSPG4", "MCAM", "NOTCH3"),
  "T cells"     = c("TRBC2", "LTB", "CD2", "CD3D", "CD3E")
)

# 巨噬细胞亚群版 marker
macrophage_subtype_marker_panel <- list(
  "LAM" = c("LGALS3", "ACP5", "PLA2G7", "GPNMB", "CTSD"),
  "IAM" = c("CCL4L2", "CCL3L1", "CCL3", "CCL4", "CD83"),
  "OSM" = c("TXNRD1", "SCD", "GCLM", "SQSTM1", "CTSB"),
  "TRM" = c("STAB1", "MAF", "DAB2", "F13A1", "MRC1", "CD163")
)
macrophage_subtype_marker_panel <- lapply(
  macrophage_subtype_marker_panel,
  function(x) unique(x[x %in% rownames(combined_tsne)])
)

# 保留数据中实际存在的基因，并在细胞类型间按首次出现去重。
marker_seen <- character(0)
marker_panel <- lapply(marker_panel, function(x) {
  x <- unique(x[x %in% rownames(combined_tsne)])
  x <- x[!x %in% marker_seen]
  marker_seen <<- c(marker_seen, x)
  x
})
marker_panel <- marker_panel[lengths(marker_panel) > 0]
marker_genes <- unique(unlist(marker_panel, use.names = FALSE))

cell_type_colors <- c(
  "Macrophage"  = "#3C78A8",
  "Monocyte"    = "#E39C37",
  "Microglia"   = "#C85C5C",
  "DC"          = "#9B72B0",
  "Endothelial" = "#668F3D",
  "Fibroblast"  = "#3F8F72",
  "Pericyte"    = "#A66F58",
  "T cells"     = "#268C92"
)

cell_type_subtype_colors <- c(
  "LAM" = "#2F6F9F",
  "IAM" = "#D95F4C",
  "OSM" = "#E69F00",
  "TRM" = "#4C956C",
  "Monocyte" = "#E39C37",
  "Microglia" = "#9B72B0",
  "Fibroblast" = "#3F8F72",
  "T cells" = "#268C92",
  "Pericyte" = "#A66F58",
  "Endothelial" = "#668F3D",
  "DC" = "#C77CFF"
)

markers_final <- marker_genes

combined_anno <- RenameIdents(
  combined_tsne,
  "0" = "Macrophage",
  "1" = "Macrophage",
  "2" = "Macrophage",
  "3" = "Macrophage",
  "4" = "Monocyte",
  "5" = "Monocyte",
  "6" = "Microglia",
  "7" = "Fibroblast",
  "8" = "T cells",
  "9" = "Pericyte",
  "10" = "Endothelial",
  "11" = "T cells",
  "12" = "DC",
  "13" = "Unknown",
  "14" = "DC"
)

combined_anno$cell_type <- Idents(combined_anno)

# 在同一对象中保留第二套注释：0–3 分别为 LAM/IAM/OSM/TRM，
# 其他 cluster 与大类注释保持一致。
cell_type_subtype_map <- c(
  "0" = "LAM",
  "1" = "IAM",
  "2" = "OSM",
  "3" = "TRM",
  "4" = "Monocyte",
  "5" = "Monocyte",
  "6" = "Microglia",
  "7" = "Fibroblast",
  "8" = "T cells",
  "9" = "Pericyte",
  "10" = "Endothelial",
  "11" = "T cells",
  "12" = "DC",
  "13" = "Unknown",
  "14" = "DC"
)

cell_type_subtype_order <- c(
  "LAM", "IAM", "OSM", "TRM", "Monocyte", "Microglia",
  "Fibroblast", "T cells", "Pericyte", "Endothelial", "DC"
)

cell_type_subtype_full_map <- c(
  "LAM" = "Lipid-associated macrophages",
  "IAM" = "Inflammatory-associated macrophages",
  "OSM" = "Oxidative stress-response macrophages",
  "TRM" = "Tissue-resident macrophages",
  "DC" = "Dendritic cells"
)

cell_type_subtype_full_order <- c(
  "Lipid-associated macrophages",
  "Inflammatory-associated macrophages",
  "Oxidative stress-response macrophages",
  "Tissue-resident macrophages",
  "Monocyte", "Microglia", "Fibroblast", "T cells",
  "Pericyte", "Endothelial", "Dendritic cells"
)

combined_anno$cell_type_subtype <- unname(
  cell_type_subtype_map[as.character(combined_anno$seurat_clusters)]
)
combined_anno$cell_type_subtype_full <- combined_anno$cell_type_subtype
subtype_full_idx <- combined_anno$cell_type_subtype %in% names(cell_type_subtype_full_map)
combined_anno$cell_type_subtype_full[subtype_full_idx] <- unname(
  cell_type_subtype_full_map[combined_anno$cell_type_subtype[subtype_full_idx]]
)

combined_anno <- subset(combined_anno, subset = cell_type != "Unknown")
combined_anno$cell_type <- factor(
  as.character(combined_anno$cell_type),
  levels = cell_type_order
)
combined_anno$cell_type_subtype <- factor(
  combined_anno$cell_type_subtype,
  levels = cell_type_subtype_order
)
combined_anno$cell_type_subtype_full <- factor(
  combined_anno$cell_type_subtype_full,
  levels = cell_type_subtype_full_order
)
Idents(combined_anno) <- "cell_type"
save(combined_anno,file = "rda/UMAP_annotated.rda")
# load("rda/UMAP_annotated.rda")

# =============================================================================
# 可视化
# =============================================================================

# =============================================================================
# 1.  GLOBAL THEME
# =============================================================================

theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_rect(color = "black", fill = NA,
                                       linewidth = 0.5),
      panel.grid        = element_blank(),
      
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(color = "black", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size),
      
      # 标题居中
      plot.title    = element_text(size   = base_size + 2, face  = "bold",
                                   hjust  = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(size   = base_size - 1, color = "grey45",
                                   hjust  = 0.5, margin = margin(b = 6)),
      
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(9, "pt"),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 1, color = "black"),
      legend.position   = "right",
      legend.margin     = margin(0, 0, 0, 4),
      
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold",
                                       hjust = 0.5),
      plot.margin       = margin(8, 8, 8, 8)
    )
}

# 以 p_tsne_before 为基准的统一可视化字体与图例规格。
theme_tsne_reference <- function(x_angle = 0, x_hjust = 0.5,
                                 x_face = "bold") {
  theme_nature() +
    theme(
      legend.justification = "center",
      plot.title = element_text(
        size = 28, face = "bold", hjust = 0.5,
        margin = margin(b = 14)
      ),
      plot.subtitle = element_text(
        size = 16, face = "bold", hjust = 0.5,
        margin = margin(b = 10)
      ),
      axis.text.x = element_text(
        size = 18, face = x_face, color = "black",
        angle = x_angle, hjust = x_hjust, vjust = 1
      ),
      axis.text.y = element_text(size = 18, face = "bold", color = "black"),
      axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 14)),
      axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 14)),
      legend.text = element_text(size = 20, face = "bold"),
      legend.title = element_text(
        size = 24, face = "bold", lineheight = 0.9,
        margin = margin(b = 12), color = "black"
      ),
      plot.margin = margin(16, 16, 16, 16)
    )
}


# =============================================================================
# 2.  tsne  ·  注释前（按 cluster 上色）
# =============================================================================

n_cl <- length(levels(combined_tsne))

# 用与细胞类型同系但更深的色调覆盖 0–12 个 cluster
cluster_pal <- colorRampPalette(
  c("#4A86C8",   # 蓝
    "#5FA898",   # 青
    "#D96B5A",   # 砖红
    "#8A6AAD",   # 紫
    "#E09500",   # 金
    "#D96560",   # 玫红
    "#4E9B7A",   # 绿
    "#7EB3D4",   # 浅蓝
    "#C4856E",   # 暖棕
    "#6B8E23",   # 橄榄绿
    "#C77CFF",   # 淡紫
    "#00A6A6",   # 蓝绿
    "#E6B8A2",   # 肉桂粉
    "#7F7F7F")   # 中灰
)(n_cl)
names(cluster_pal) <- levels(combined_tsne)

p_tsne_before <-
  DimPlot(
    combined_tsne,
    reduction  = "tsne",
    label      = TRUE,
    label.size = 10,
    pt.size    = 2,
    alpha     = 1, 
    cols       = cluster_pal,
    repel      = TRUE
  ) +
  labs(
    title    = "Cell clusters (pre-annotation)",
    subtitle = paste0(format(ncol(combined_tsne), big.mark = ","),
                      " cells  ·  tSNE"),
    x = "tSNE 1", y = "tSNE 2"
  ) +
  theme_tsne_reference() +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4),
                              title = "Cluster"))

p_tsne_before

ggsave("results/tsne/Before_cell_annotation_tsne.pdf",
       p_tsne_before, width = 10, height = 8,
       device = "pdf",dpi = 500)
ggsave("results/tsne/Before_cell_annotation_tsne.png",
       p_tsne_before, width = 10, height = 8, dpi = 300)

# =============================================================================
# 3.  DotPlot  ·  注释前
# =============================================================================

markers_final1 <- marker_panel

dotplot_gradient <- c(
  "#2F6F9F",  
  "#8DB9D3", 
  "#F6F0E8",  
  "#D96B5A"   
)

p_dot_before <-
  DotPlot(combined_tsne, features = markers_final1, dot.scale = 6) +
  scale_color_gradientn(colours = dotplot_gradient,
                        name    = "Scaled average\nexpression",
                        limits  = c(-2, 2), oob = squish) +
  scale_size_continuous(range = c(0.3, 6), name = "Percent\nexpressed") +
  labs(title = "Marker gene expression (pre-annotation)", x = NULL, y = "Cluster") +
  # 注释前气泡图独立字体设置：不与 tSNE 或注释后气泡图共用。
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 34, face = "bold", hjust = 0.5,
      margin = margin(b = 28), color = "black"
    ),
    # D 图不显示 DotPlot 自动生成的分面标题。
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold.italic", color = "black",
                               angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold", color = "black",
                                margin = margin(r = 12)),
    legend.title = element_text(size = 17, face = "bold", color = "black"),
    legend.text = element_text(size = 15, face = "bold", color = "black"),
    legend.key.size = unit(12, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.4),
    panel.spacing.x = unit(2.5, "mm"),
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = unit(12, "pt"),
    plot.margin = margin(16, 16, 16, 16)
  ) +
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(42, "mm"), barheight = unit(4, "mm")),
    size = guide_legend(title.position = "top", title.hjust = 0.5,
                        nrow = 1, byrow = TRUE)
  )

p_dot_before

ggsave("results/Marker/Marker_Gene_Cluster_before_final.pdf",
       p_dot_before, width = 12, height = 10,
       device = "pdf",dpi = 500)
ggsave("results/Marker/Marker_Gene_Cluster_before_final.png",
       p_dot_before, width = 12, height = 10, dpi = 300)

# =============================================================================
# 4.  tsne  ·  注释后
# =============================================================================

p_tsne_after <-
  DimPlot(
    combined_anno,
    reduction  = "tsne",
    group.by   = "cell_type",
    label      = TRUE,
    label.size = 5.5,
    pt.size    = 2,
    alpha     = 1, 
    cols       = cell_type_colors,
    repel      = TRUE
  ) +
  labs(
    title    = "Cell-type annotation",
    subtitle = paste0(format(ncol(combined_anno), big.mark = ","),
                      " cells  ·  tSNE"),
    x = "tSNE 1", y = "tSNE 2"
  ) +
  theme_tsne_reference() +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4),
                              title = "Cell type"))


p_tsne_after

ggsave("results/tsne/After_cell_annotation_tsne.pdf",
       p_tsne_after, width = 10, height = 8,
       device = "pdf",dpi = 500)
ggsave("results/tsne/After_cell_annotation_tsne.png",
       p_tsne_after, width = 10, height = 8, dpi = 300)

# 巨噬细胞亚群与其他细胞的整体 t-SNE
p_tsne_subtype <- DimPlot(
  combined_anno,
  reduction = "tsne",
  group.by = "cell_type_subtype",
  label = TRUE,
  label.size = 5,
  pt.size = 2,
  cols = cell_type_subtype_colors,
  repel = TRUE
) +
  labs(
    title = "Macrophage subtypes and other cell types",
    subtitle = paste0(format(ncol(combined_anno), big.mark = ","),
                      " cells  ·  tSNE"),
    x = "tSNE 1",
    y = "tSNE 2"
  ) +
  theme_tsne_reference() +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 4),
      title = "Cell type"
    )
  )

p_tsne_subtype

ggsave(
  "results/tsne/Macrophage_subtypes_other_cells_tsne.pdf",
  p_tsne_subtype,
  width = 11,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/tsne/Macrophage_subtypes_other_cells_tsne.png",
  p_tsne_subtype,
  width = 11,
  height = 8,
  dpi = 500
)

# 全部细胞重新聚类并进行 UMAP 降维
DefaultAssay(combined_anno) <- "SCT"
dir.create("results/SPP1_FeaturePlot", recursive = TRUE)

combined_umap_recluster <- FindNeighbors(
  combined_anno,
  reduction = "pca",
  dims = 1:10,
  graph.name = c("umap_nn", "umap_snn")
)

set.seed(13)
combined_umap_recluster <- FindClusters(
  combined_umap_recluster,
  graph.name = "umap_snn",
  resolution = 0.5
)

set.seed(13)
combined_umap_recluster <- RunUMAP(
  combined_umap_recluster,
  reduction = "pca",
  dims = 1:10,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)

# UMAP 聚类图与细胞类型注释图
p_umap_cluster <- DimPlot(
  combined_umap_recluster,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 7,
  pt.size = 2,
  cols = cluster_pal,
  repel = TRUE
) +
  labs(
    title = "Cell clusters",
    subtitle = paste0(format(ncol(combined_umap_recluster), big.mark = ","),
                      " cells  ·  UMAP"),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_tsne_reference() +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 4),
      title = "Cluster"
    )
  )

p_umap_annotated <- DimPlot(
  combined_umap_recluster,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  label.size = 5.5,
  pt.size = 2,
  cols = cell_type_colors,
  repel = TRUE
) +
  labs(
    title = "Cell-type annotation",
    subtitle = paste0(format(ncol(combined_umap_recluster), big.mark = ","),
                      " cells  ·  UMAP"),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_tsne_reference() +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 4),
      title = "Cell type"
    )
  )

p_umap_cluster
p_umap_annotated

ggsave(
  "results/tsne/UMAP_recluster_all_cells.pdf",
  p_umap_cluster,
  width = 10,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/tsne/UMAP_recluster_all_cells.png",
  p_umap_cluster,
  width = 10,
  height = 8,
  dpi = 500
)

ggsave(
  "results/tsne/UMAP_cell_type_annotation.pdf",
  p_umap_annotated,
  width = 10,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/tsne/UMAP_cell_type_annotation.png",
  p_umap_annotated,
  width = 10,
  height = 8,
  dpi = 500
)

# UMAP重聚类对应的marker气泡图
p_dot_umap_recluster <- DotPlot(
  combined_umap_recluster,
  features = markers_final1,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  scale_color_gradientn(
    colours = dotplot_gradient,
    name = "Scaled average\nexpression",
    limits = c(-2, 2),
    oob = squish
  ) +
  scale_size_continuous(
    range = c(0.3, 6),
    name = "Percent\nexpressed"
  ) +
  labs(
    title = "Marker gene expression · UMAP clusters",
    x = NULL,
    y = "UMAP cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 30,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 28)
    ),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(
      size = 13,
      face = "bold.italic",
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      margin = margin(r = 12)
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(42, "mm"),
      barheight = unit(4, "mm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  )

p_dot_umap_recluster

ggsave(
  "results/Marker/Marker_Gene_UMAP_recluster.pdf",
  p_dot_umap_recluster,
  width = 12,
  height = 10,
  device = "pdf"
)

ggsave(
  "results/Marker/Marker_Gene_UMAP_recluster.png",
  p_dot_umap_recluster,
  width = 12,
  height = 10,
  dpi = 500
)

# UMAP注释后按细胞类型展示marker表达
p_dot_umap_annotated <- DotPlot(
  combined_umap_recluster,
  features = markers_final1,
  group.by = "cell_type",
  dot.scale = 6
) +
  scale_color_gradientn(
    colours = dotplot_gradient,
    name = "Scaled average\nexpression",
    limits = c(-2, 2),
    oob = squish
  ) +
  scale_size_continuous(
    range = c(0.3, 6),
    name = "Percent\nexpressed"
  ) +
  labs(
    title = "Marker gene expression · annotated UMAP",
    x = NULL,
    y = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 30,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 28)
    ),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(
      size = 13,
      face = "bold.italic",
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      margin = margin(r = 12)
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(42, "mm"),
      barheight = unit(4, "mm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  )

p_dot_umap_annotated

ggsave(
  "results/Marker/Marker_Gene_UMAP_annotated.pdf",
  p_dot_umap_annotated,
  width = 12,
  height = 10,
  device = "pdf"
)

ggsave(
  "results/Marker/Marker_Gene_UMAP_annotated.png",
  p_dot_umap_annotated,
  width = 12,
  height = 10,
  dpi = 500
)

p_spp1_tsne <- FeaturePlot(
  combined_anno,
  features = "SPP1",
  reduction = "tsne",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  pt.size = 1.2,
  cols = c("grey92", "#D73027"),
  combine = FALSE
)[[1]] +
  labs(
    title = "SPP1 expression · t-SNE",
    subtitle = paste0(format(ncol(combined_anno), big.mark = ","),
                      " cells  ·  tSNE"),
    x = "tSNE 1",
    y = "tSNE 2",
    color = "SPP1 expression"
  ) +
  theme_tsne_reference() +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.box.just = "left"
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(3, "mm"),
      barheight = unit(35, "mm")
    )
  )

p_spp1_umap <- FeaturePlot(
  combined_umap_recluster,
  features = "SPP1",
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  pt.size = 1.2,
  cols = c("grey92", "#D73027"),
  combine = FALSE
)[[1]] +
  labs(
    title = "SPP1 expression · UMAP",
    subtitle = paste0(format(ncol(combined_umap_recluster), big.mark = ","),
                      " cells  ·  UMAP"),
    x = "UMAP 1",
    y = "UMAP 2",
    color = "SPP1 expression"
  ) +
  theme_tsne_reference() +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.box.just = "left"
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(3, "mm"),
      barheight = unit(35, "mm")
    )
  )

p_spp1_tsne
p_spp1_umap

# 在全部细胞的 FeaturePlot 中圈出巨噬细胞
tsne_spp1_coords <- as.data.frame(Embeddings(combined_anno, "tsne"))
tsne_spp1_coords$cell_type <- as.character(combined_anno$cell_type)
tsne_spp1_coords$cell_type_subtype <- as.character(
  combined_anno$cell_type_subtype
)
tsne_spp1_coords$SPP1 <- FetchData(
  combined_anno,
  vars = "SPP1",
  layer = "data"
)$SPP1
macrophage_tsne_coords <- tsne_spp1_coords %>%
  filter(cell_type == "Macrophage")

umap_spp1_coords <- as.data.frame(
  Embeddings(combined_umap_recluster, "umap")
)
umap_spp1_coords$cell_type <- as.character(
  combined_umap_recluster$cell_type
)
macrophage_umap_coords <- umap_spp1_coords %>%
  filter(cell_type == "Macrophage")

# 重复用于整体巨噬细胞及4个亚型：仅保留空间上的最大连通主体。
keep_main_spatial_component <- function(data, x_col, y_col, eps) {
  coordinates <- as.matrix(data[, c(x_col, y_col)])
  nearest <- RANN::nn2(
    coordinates,
    k = min(20, nrow(data))
  )

  from <- rep(
    seq_len(nrow(data)),
    each = ncol(nearest$nn.idx) - 1
  )
  to <- as.vector(t(nearest$nn.idx[, -1, drop = FALSE]))
  distance <- as.vector(t(nearest$nn.dists[, -1, drop = FALSE]))
  retained_edges <- distance <= eps

  spatial_graph <- igraph::make_empty_graph(
    n = nrow(data),
    directed = FALSE
  )
  spatial_graph <- igraph::add_edges(
    spatial_graph,
    as.vector(
      t(cbind(from[retained_edges], to[retained_edges]))
    )
  )
  component_id <- igraph::components(spatial_graph)$membership
  main_component <- as.integer(
    names(which.max(table(component_id)))
  )

  data[component_id == main_component, ]
}

macrophage_tsne_hull_coords <- keep_main_spatial_component(
  macrophage_tsne_coords,
  x_col = "tSNE_1",
  y_col = "tSNE_2",
  eps = 1.5
)

macrophage_umap_hull_coords <- keep_main_spatial_component(
  macrophage_umap_coords,
  x_col = "UMAP_1",
  y_col = "UMAP_2",
  eps = 0.2
)

macrophage_subtype_hull_coords <- lapply(
  c("LAM", "IAM", "OSM", "TRM"),
  function(subtype) {
    subtype_coords <- macrophage_tsne_coords %>%
      filter(cell_type_subtype == subtype)
    keep_main_spatial_component(
      subtype_coords,
      x_col = "tSNE_1",
      y_col = "tSNE_2",
      eps = 2.5
    )
  }
) %>%
  bind_rows()

# Figure 3H：圈出IAM和OSM的主要空间区域，沿用Figure 3G的轮廓样式。
iam_osm_tsne_hull_coords <- macrophage_subtype_hull_coords %>%
  filter(cell_type_subtype %in% c("IAM", "OSM"))

p_spp1_tsne_macrophage <- p_spp1_tsne +
  ggforce::geom_mark_hull(
    data = macrophage_tsne_hull_coords,
    aes(
      x = tSNE_1,
      y = tSNE_2
    ),
    inherit.aes = FALSE,
    color = "#245B84",
    fill = NA,
    linewidth = 1.2,
    linetype = "dashed",
    concavity = 4,
    expand = unit(0.8, "mm"),
    radius = unit(2, "mm"),
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -47, y = -41,
    label = "Macrophage",
    color = "#245B84",
    size = 7.5,
    fontface = "bold",
    hjust = 0
  )

p_spp1_umap_macrophage <- p_spp1_umap +
  ggforce::geom_mark_hull(
    data = macrophage_umap_hull_coords,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    inherit.aes = FALSE,
    color = "#245B84",
    fill = NA,
    linewidth = 1.2,
    linetype = "dashed",
    concavity = 4,
    expand = unit(2.5, "mm"),
    radius = unit(3, "mm"),
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 4.6, y = 6.4,
    label = "Macrophage",
    color = "#245B84",
    size = 7.5,
    fontface = "bold",
    hjust = 0.5
  )

# 使用巨噬细胞分析脚本中的原配色分别圈出4个亚型
macrophage_subtype_pal <- c(
  "LAM" = "#2F6F73",
  "IAM" = "#B65C4A",
  "OSM" = "#6E5F8F",
  "TRM" = "#C49A3A"
)

p_spp1_tsne_subtypes <- ggplot(
  tsne_spp1_coords,
  aes(x = tSNE_1, y = tSNE_2)
) +
  geom_point(
    aes(fill = SPP1),
    shape = 21,
    color = "transparent",
    size = 1.2
  ) +
  scale_fill_gradientn(
    colours = c("grey92", "#D73027"),
    name = "SPP1 expression",
    limits = quantile(tsne_spp1_coords$SPP1, c(0.05, 0.95)),
    oob = squish
  ) +
  labs(
    title = "SPP1 expression · macrophage subtypes",
    subtitle = paste0(format(ncol(combined_anno), big.mark = ","),
                      " cells  ·  tSNE"),
    x = "tSNE 1",
    y = "tSNE 2"
  ) +
  theme_tsne_reference() +
  ggforce::geom_mark_hull(
    data = macrophage_subtype_hull_coords,
    aes(
      x = tSNE_1,
      y = tSNE_2,
      color = cell_type_subtype,
      label = cell_type_subtype,
      group = cell_type_subtype
    ),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 1.2,
    linetype = "dashed",
    concavity = 4,
    expand = unit(2.5, "mm"),
    radius = unit(3, "mm"),
    label.fontsize = 11,
    label.fill = "white",
    con.linetype = "dashed",
    key_glyph = "path",
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = macrophage_subtype_pal,
    name = "Macrophage subtype"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = "dashed",
        linewidth = 1.5
      )
    )
  )

p_spp1_tsne_macrophage
p_spp1_umap_macrophage
p_spp1_tsne_subtypes

ggsave(
  "results/SPP1_FeaturePlot/SPP1_FeaturePlot_tsne.pdf",
  p_spp1_tsne_macrophage,
  width = 10,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/SPP1_FeaturePlot/SPP1_FeaturePlot_tsne.png",
  p_spp1_tsne_macrophage,
  width = 10,
  height = 8,
  dpi = 500
)

ggsave(
  "results/SPP1_FeaturePlot/SPP1_FeaturePlot_umap_recluster.pdf",
  p_spp1_umap_macrophage,
  width = 10,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/SPP1_FeaturePlot/SPP1_FeaturePlot_umap_recluster.png",
  p_spp1_umap_macrophage,
  width = 10,
  height = 8,
  dpi = 500
)

ggsave(
  "results/SPP1_FeaturePlot/SPP1_tsne_macrophage_subtypes.pdf",
  p_spp1_tsne_subtypes,
  width = 11,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/SPP1_FeaturePlot/SPP1_tsne_macrophage_subtypes.png",
  p_spp1_tsne_subtypes,
  width = 11,
  height = 8,
  dpi = 500
)

# 提取巨噬细胞，并按照SPP1归一化表达中位数定义SPP1-high
macrophage_spp1 <- subset(
  combined_anno,
  subset = cell_type == "Macrophage"
)
DefaultAssay(macrophage_spp1) <- "SCT"

macrophage_spp1_expression <- FetchData(
  macrophage_spp1,
  vars = "SPP1",
  layer = "data"
)$SPP1
macrophage_spp1_median <- median(macrophage_spp1_expression)
macrophage_spp1$SPP1_group <- ifelse(
  macrophage_spp1_expression > macrophage_spp1_median,
  "SPP1-high",
  "SPP1-low"
)

macrophage_high_coords <- as.data.frame(
  Embeddings(macrophage_spp1, "tsne")
)
macrophage_high_coords$SPP1_group <- macrophage_spp1$SPP1_group
macrophage_high_coords <- macrophage_high_coords %>%
  filter(SPP1_group == "SPP1-high")

p_macrophage_spp1_high <- FeaturePlot(
  macrophage_spp1,
  features = "SPP1",
  reduction = "tsne",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  pt.size = 1.5,
  cols = c("grey92", "#D73027"),
  combine = FALSE
)[[1]] +
  geom_point(
    data = macrophage_high_coords,
    aes(
      x = tSNE_1,
      y = tSNE_2,
      shape = SPP1_group
    ),
    inherit.aes = FALSE,
    size = 1.8,
    stroke = 0.35,
    color = "black",
    fill = NA
  ) +
  ggforce::geom_mark_hull(
    data = iam_osm_tsne_hull_coords,
    aes(
      x = tSNE_1,
      y = tSNE_2
    ),
    inherit.aes = FALSE,
    color = "#245B84",
    fill = NA,
    linewidth = 1.2,
    linetype = "dashed",
    concavity = 4,
    expand = unit(2.5, "mm"),
    radius = unit(3, "mm"),
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -38, y = -20,
    label = "IAM and OSM",
    color = "#245B84",
    size = 7.5,
    fontface = "bold",
    hjust = 0
  ) +
  scale_shape_manual(
    values = c("SPP1-high" = 21),
    name = paste0(
      "Median cutoff = ",
      round(macrophage_spp1_median, 2)
    ),
    guide = guide_legend(
      override.aes = list(
        color = "black",
        fill = NA,
        linetype = 0,
        alpha = 1,
        size = 4
      )
    )
  ) +
  labs(
    title = "SPP1 expression in macrophages",
    subtitle = "SPP1-high: normalized expression > median",
    x = "tSNE 1",
    y = "tSNE 2",
    color = "SPP1 expression"
  ) +
  theme_tsne_reference() +
  theme(
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.box.just = "left"
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(3, "mm"),
      barheight = unit(35, "mm")
    )
  )

p_macrophage_spp1_high

ggsave(
  "results/SPP1_FeaturePlot/Macrophage_SPP1_high_tsne.pdf",
  p_macrophage_spp1_high,
  width = 10,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/SPP1_FeaturePlot/Macrophage_SPP1_high_tsne.png",
  p_macrophage_spp1_high,
  width = 10,
  height = 8,
  dpi = 500
)

p_spp1_feature <- p_spp1_tsne_macrophage +
  p_spp1_umap_macrophage +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

p_spp1_feature

ggsave(
  "results/tsne/SPP1_FeaturePlot_all_cells.pdf",
  p_spp1_feature,
  width = 18,
  height = 8,
  device = "pdf"
)

ggsave(
  "results/tsne/SPP1_FeaturePlot_all_cells.png",
  p_spp1_feature,
  width = 18,
  height = 8,
  dpi = 500
)

# =============================================================================
# 5.  DotPlot  ·  注释后
# =============================================================================

p_dot_after <-
  DotPlot(
    combined_anno,
    features = markers_final1,
    group.by = "cell_type",
    dot.scale = 6
  ) +
  scale_color_gradientn(colours = dotplot_gradient,
                        name    = "Scaled average\nexpression",
                        limits  = c(-2, 2), oob = squish) +
  scale_size_continuous(range = c(0.3, 6), name = "Percent\nexpressed") +
  labs(title = "Marker gene expression by cell type", x = NULL, y = NULL) +
  # 注释后气泡图独立字体设置：分组标题与细胞类型标签适当放大。
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 30, face = "bold", hjust = 0.5,
      margin = margin(b = 28), color = "black"
    ),
    # E 图不显示 DotPlot 自动生成的分面标题。
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 13, face = "bold.italic", color = "black",
                               angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold", color = "black",
                                margin = margin(r = 10)),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    legend.text = element_text(size = 14, face = "bold", color = "black"),
    legend.key.size = unit(12, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.4),
    panel.spacing.x = unit(2.5, "mm"),
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = unit(12, "pt"),
    plot.margin = margin(14, 14, 14, 14)
  ) +
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(42, "mm"), barheight = unit(4, "mm")),
    size = guide_legend(title.position = "top", title.hjust = 0.5,
                        nrow = 1, byrow = TRUE)
  )

p_dot_after

ggsave("results/Marker/Marker_Gene_Cluster_after.pdf",
       p_dot_after, width = 12, height = 10,
       device = "pdf",dpi = 500)
ggsave("results/Marker/Marker_Gene_Cluster_after.png",
       p_dot_after, width = 12, height = 10, dpi = 300)


# =============================================================================
# 6. 堆叠柱状图  ·  细胞类型比例  —  顶刊配色版
# =============================================================================
# 样本中各细胞分布堆叠图
meta <- combined_anno@meta.data

# 统计每个 Sample 中每种 cell type 的数量
cell_counts <- meta %>%
  group_by(sample, cell_type) %>%
  summarise(count = n(), .groups = "drop")

# 计算每个 Sample 中的总细胞数
total_counts <- meta %>%
  group_by(sample) %>%
  summarise(total = n(), .groups = "drop")

# 合并并计算比例
cell_props <- cell_counts %>%
  left_join(total_counts, by = "sample") %>%
  mutate(
    prop = count / total,
    percent = prop * 100
  )

write.csv(
  cell_props,
  "results/tsne/celltype_proportion_by_sample.csv",
  row.names = FALSE
)

# 统计各细胞类型在全部样本已注释细胞中的比例
cell_props_all_samples <- meta %>%
  count(cell_type, name = "count") %>%
  mutate(
    total = sum(count),
    prop = count / total,
    percent = prop * 100
  ) %>%
  arrange(desc(count))

write.csv(
  cell_props_all_samples,
  "results/tsne/celltype_proportion_all_samples.csv",
  row.names = FALSE
)

# 柱状图与 tSNE 使用同一组细胞类型配色和顺序。
cell_type_colors_soft <- cell_type_colors
ct_order <- cell_type_order

# 替换柱状图中的配色
p_bar <- ggplot(cell_props, aes(x = sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.72,
           color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cell_type_colors_soft,
                    breaks = ct_order) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Cell type composition per sample",
       x = NULL, y = "Cell proportion", fill = "Cell type") +
  theme_tsne_reference(x_angle = 45, x_hjust = 1) +
  guides(fill = guide_legend(ncol = 1,
                             keywidth  = unit(10, "pt"),
                             keyheight = unit(10, "pt")))

p_bar

ggsave("results/tsne/celltype_proportion_barplot.pdf",
       p_bar, width = 10, height = 8,
       device = "pdf",dpi = 500)
ggsave("results/tsne/celltype_proportion_barplot.png",
       p_bar, width = 10, height = 8, dpi = 300)


# =============================================================================
# 7.  每个样本中各细胞类型的细胞数目（绝对计数）
#     使用元数据中的 "sample" 列
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# 提取元数据并检查必要列 ----
meta <- combined_anno@meta.data
# 计算每个样本 × 细胞类型的细胞数量 ----
cell_counts <- meta %>%
  group_by(sample, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(sample, cell_type, fill = list(count = 0)) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order))

# 堆叠条形图（绝对细胞数）----
p_count_bar <- ggplot(cell_counts, aes(x = sample, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7,
           color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cell_type_colors_soft,   # 使用已有的高饱和配色
                    breaks = ct_order) +
  labs(title = "Number of cells per cell type across samples",
       x = "Sample", y = "Cell count", fill = "Cell type") +
  theme_tsne_reference(x_angle = 45, x_hjust = 1) +
  guides(fill = guide_legend(ncol = 1, reverse = FALSE))

p_count_bar

ggsave("results/tsne/celltype_absolute_count_per_sample_barplot.pdf",
       p_count_bar, width = 10, height = 8,
       device = "pdf", dpi = 500)
ggsave("results/tsne/celltype_absolute_count_per_sample_barplot.png",
       p_count_bar, width = 10, height = 8, dpi = 300)


# =============================================================================
# 8.  featureplot
# =============================================================================
markers_fp <- list(
  "Macrophage"  = c("CD163", "C1QC"),
  "Monocyte"    = c("FCN1", "S100A9"),
  "Microglia"   = c("CX3CR1", "TMEM119"),
  "DC"          = c("FCER1A", "CD1C"),
  "Endothelial" = c("SPARCL1", "VWF"),
  "Fibroblast"  = c("LUM",   "THBS2"),
  "Pericyte"    = c("MCAM",  "RGS5"),
  "T cells"     = c("CD3D", "CD3E")
)


features_ordered <- unlist(markers_fp, use.names = FALSE)
features_ordered <- features_ordered[features_ordered %in% rownames(combined_anno)]

# ── 高表达颜色（高饱和） ──────────────────────────────────────────────────────
celltype_colors <- c(
  "Macrophage"  = "#245B84",
  "Monocyte"    = "#B86516",
  "Microglia"   = "#9F3F46",
  "DC"          = "#75518B",
  "Endothelial" = "#47702A",
  "Fibroblast"  = "#246A50",
  "Pericyte"    = "#7C4D3E",
  "T cells"     = "#146A70"
)


# ── 10k细胞专用调色：低表达接近白色，高表达全饱和，中间快速跳变 ─────────────
make_balanced_palette <- function(high_col) {
  colorRampPalette(c(
    "#D5D8DC",   # ← 低表达：中灰，能看到细胞轮廓但不抢眼
    "#BFC9CA",   # 低中：轻微过渡
    colorspace::lighten(high_col, 0.55),  # 中低：开始带颜色倾向
    colorspace::lighten(high_col, 0.25),  # 中高：明显有色
    high_col,                              # 高表达：全饱和
    colorspace::darken(high_col, 0.12)    # 最高：略深防过曝
  ))(300)
}


feature_to_celltype <- rep(names(markers_fp), lengths(markers_fp))
names(feature_to_celltype) <- unlist(markers_fp, use.names = FALSE)

plot_one_feature <- function(feature) {
  ct  <- feature_to_celltype[[feature]]
  pal <- make_balanced_palette(celltype_colors[[ct]])
  
  FeaturePlot(
    combined_anno,
    features   = feature,
    reduction  = "tsne",
    pt.size    = 0.6,
    alpha      = 0.9,      # ← 从 0.85 略微提高，灰色点更实
    order      = TRUE,
    cols       = pal,
    min.cutoff = "q05",    # ← 从 q10 放宽回 q05，低表达细胞不被截掉
    max.cutoff = "q97",
    raster     = FALSE
  ) +
    labs(title = feature, subtitle = ct) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 18, face = "bold.italic",
        hjust = 0.5, color = "#1A1A2E", margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        size = 14, face = "bold",
        hjust = 0.5, color = celltype_colors[[ct]], margin = margin(b = 3)
      ),
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(8, 4, 8, 4)
    )
}

p_list <- lapply(features_ordered, plot_one_feature)

# ── 拼图 ──────────────────────────────────────────────────────────────────────
p_feature_final <- wrap_plots(p_list, ncol = 4) +
  plot_annotation(
    title   = "",
    theme   = theme(
      plot.title   = element_text(
        size = 28, face = "bold", hjust = 0.5,
        color = "#111111", margin = margin(b = 14)
      ),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )


ggsave(
  "results/Marker/FeaturePlot_markers.pdf",
  p_feature_final,
  width = 14, height = 12,
  device = "pdf"
)

ggsave(
  "results/Marker/FeaturePlot_markers.png",
  p_feature_final,
  width = 14, height = 12,
  dpi = 500, bg = "white"
)

# =============================================================================
# 输出为MEBOCOST分析数据
# =============================================================================
get_assay_matrix <- function(object, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) GetAssayData(object, assay = assay, slot = layer)
  )
}

export_for_mebocost <- function(scRNA, cells, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  counts_matrix <- get_assay_matrix(scRNA, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
  metadata <- scRNA@meta.data[cells, , drop = FALSE]
  
  if (!all(rownames(metadata) == colnames(counts_matrix))) {
    stop("MEBOCOST export failed: metadata rows do not match count matrix columns.", call. = FALSE)
  }
  
  Matrix::writeMM(counts_matrix, file = file.path(output_dir, "matrix.mtx"))
  writeLines(colnames(counts_matrix), file.path(output_dir, "barcodes.tsv"))
  writeLines(rownames(counts_matrix), file.path(output_dir, "features.tsv"))
  write.csv(metadata, file = file.path(output_dir, "metadata.csv"), quote = TRUE, row.names = TRUE)
}

export_for_mebocost(
  scRNA = combined_anno,
  cells = colnames(combined_anno),
  output_dir = "results/export_for_MEBOCOST"
)


# =============================================================================
# 可视化补充：PCA肘部图、碎石图、质控前后图、高变基因选择图
# =============================================================================

library(ggrepel)

dir.create("results/QC",  showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# A.  PCA 肘部图（Elbow Plot）
# =============================================================================

p_elbow <- ElbowPlot(combined_PCA, ndims = 50) +
  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    color = "#D96B5A",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 11,
    y = Inf,
    label = "dim = 10",
    vjust = 1.5,
    hjust = 0,
    color = "#D96B5A",
    size = 4,
    fontface = "bold"
  ) +
  labs(
    title = "PCA Elbow Plot",
    x = "Principal Component",
    y = "Standard Deviation"
  ) +
  theme_nature(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 28,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 14),
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 22,
      face = "bold",
      color = "black",
      margin = margin(t = 14)
    ),
    axis.title.y = element_text(
      size = 22,
      face = "bold",
      color = "black",
      margin = margin(r = 14)
    ),
    
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 10)
  )

p_elbow

ggsave(
  "results/QC/PCA_elbow_plot.pdf",
  p_elbow,
  width = 6.2,
  height = 6.2,
  device = "pdf",
  dpi = 500,
  bg = "white"
)

ggsave(
  "results/QC/PCA_elbow_plot.png",
  p_elbow,
  width = 6.2,
  height = 6.2,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# B.  碎石图（Scree Plot）—— 每个 PC 的方差解释量 + 累积方差
# =============================================================================

pca_sdev <- combined_PCA[["pca"]]@stdev
var_exp  <- pca_sdev^2 / sum(pca_sdev^2) * 100   # 方差解释百分比
cum_var  <- cumsum(var_exp)
n_plot   <- min(50, length(pca_sdev))
max_ve   <- max(var_exp[1:n_plot])                 # 用于双 Y 轴缩放

scree_df <- data.frame(
  PC   = seq_len(n_plot),
  VE   = var_exp[1:n_plot],
  CumV = cum_var[1:n_plot]
)

p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_col(
    aes(y = VE),
    fill = "#4A86C8",
    alpha = 0.80,
    width = 0.7
  ) +
  geom_line(
    aes(y = CumV / 100 * max_ve),
    color = "#D96B5A",
    linewidth = 0.9,
    group = 1
  ) +
  geom_point(
    aes(y = CumV / 100 * max_ve),
    color = "#D96B5A",
    size = 2.2
  ) +
  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 11,
    y = max_ve * 0.97,
    label = "PC = 10",
    hjust = 0,
    size = 4,
    color = "grey35",
    fontface = "bold"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 20, 30, 40, 50)
  ) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    expand = expansion(mult = c(0, 0.05)),
    sec.axis = sec_axis(
      ~ . / max_ve * 100,
      name = "Cumulative Variance (%)"
    )
  ) +
  labs(
    title = "PCA Scree Plot",
    x = "Principal Component"
  ) +
  theme_nature(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 28,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 14),
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 22,
      face = "bold",
      color = "black",
      margin = margin(t = 14)
    ),
    axis.title.y = element_text(
      size = 22,
      face = "bold",
      color = "black",
      margin = margin(r = 14)
    ),
    
    axis.title.y.right = element_text(
      size = 20,
      face = "bold",
      color = "#D96B5A",
      margin = margin(l = 8)
    ),
    axis.text.y.right = element_text(
      size = 20,
      face = "bold",
      color = "#D96B5A"
    ),
    
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 16, 10, 10)
  )

p_scree

ggsave(
  "results/QC/PCA_scree_plot.pdf",
  p_scree,
  width = 6.2,
  height = 6.2,
  device = "pdf",
  dpi = 500,
  bg = "white"
)

ggsave(
  "results/QC/PCA_scree_plot.png",
  p_scree,
  width = 6.2,
  height = 6.2,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# C.  质控前后图（QC Before / After）
# =============================================================================

# ── C1. 整合质控前后元数据 ──────────────────────────────────────────────────
make_meta_df <- function(obj_list, tag) {
  do.call(rbind, lapply(names(obj_list), function(nm) {
    need_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
    m            <- obj_list[[nm]]@meta.data[,
                      intersect(need_cols, colnames(obj_list[[nm]]@meta.data)),
                      drop = FALSE]
    m$sample     <- nm
    m$QC_status  <- tag
    m
  }))
}

meta_before   <- make_meta_df(seurat_list_1, "Before QC")
meta_after    <- make_meta_df(seurat_QC_1,   "After QC")
meta_qc       <- rbind(meta_before, meta_after)
meta_qc$QC_status <- factor(meta_qc$QC_status,
                             levels = c("Before QC", "After QC"))

n_before <- nrow(meta_before)
n_after  <- nrow(meta_after)
qc_fill  <- c("Before QC" = "#B0BEC5", "After QC" = "#4A86C8")

# ── C2. 小提琴图函数（三项 QC 指标）────────────────────────────────────────
vln_qc <- function(dat, y_var, y_lab, hlines = NULL) {
  p <- ggplot(dat, aes(x = QC_status, y = .data[[y_var]], fill = QC_status)) +
    geom_violin(trim = TRUE, scale = "width",
                alpha = 0.85, linewidth = 0.35, color = "white") +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 color = "black", linewidth = 0.40, alpha = 0.90) +
    scale_fill_manual(values = qc_fill, guide = "none") +
    labs(title = y_lab, x = NULL, y = y_lab) +
    theme_nature(base_size = 10) +
    theme(
      plot.title  = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 9)
    )
  if (!is.null(hlines))
    p <- p + geom_hline(yintercept = hlines, linetype = "dashed",
                        color = "#D96B5A", linewidth = 0.5)
  p
}

p_vln_count   <- vln_qc(meta_qc, "nCount_RNA",   "nCount_RNA",   c(1000, 50000))
p_vln_feature <- vln_qc(meta_qc, "nFeature_RNA", "nFeature_RNA", c(200, 6000))
p_vln_mt      <- vln_qc(meta_qc, "percent.mt",   "MT%",          15)

p_qc_violin <- (p_vln_count | p_vln_feature | p_vln_mt) &
  theme_nature(base_size = 11) &
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10),
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      color = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black",
      margin = margin(t = 14)
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      color = "black",
      margin = margin(r = 14)
    ),
    
    legend.title = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    legend.text = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_qc_violin <- p_qc_violin +
  plot_annotation(
    title = "Quality Control: Before vs. After Filtering",
    theme = theme(
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 14),
        color = "black"
      ),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 16, 10, 10)
    )
  )

p_qc_violin

ggsave(
  "results/QC/QC_violin_before_after.pdf",
  p_qc_violin,
  width = 12,
  height = 5.5,
  device = "pdf",
  dpi = 500,
  bg = "white"
)

ggsave(
  "results/QC/QC_violin_before_after.png",
  p_qc_violin,
  width = 12,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 审稿补充分析：整合后样本混合的替代量化评估与UMAP展示
# =============================================================================

review_batch_dir <- "results/reviewer_supplement/batch_mixing"
dir.create(review_batch_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("combined_anno")) {
  load("rda/UMAP_annotated.rda")
}
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!"umap" %in% names(combined_anno@reductions)) {
  combined_anno <- RunUMAP(
    combined_anno,
    reduction = "pca",
    dims = 1:10,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )
}

review_embedding <- Embeddings(combined_anno, "pca")[, 1:10, drop = FALSE]
review_sample <- as.character(combined_anno$sample)
review_k <- min(30L, nrow(review_embedding) - 1L)

if (requireNamespace("FNN", quietly = TRUE)) {
  review_nn <- FNN::get.knn(review_embedding, k = review_k)$nn.index
} else {
  review_dist <- as.matrix(dist(review_embedding))
  diag(review_dist) <- Inf
  review_nn <- t(apply(review_dist, 1, order)[seq_len(review_k), , drop = FALSE])
}

review_entropy <- function(x) {
  p <- prop.table(table(x))
  -sum(as.numeric(p) * log(as.numeric(p))) / log(length(unique(review_sample)))
}

batch_mixing_cell <- data.frame(
  cell = rownames(review_embedding),
  sample = review_sample,
  cell_type = as.character(combined_anno$cell_type),
  cell_type_subtype = if ("cell_type_subtype" %in% colnames(combined_anno@meta.data)) {
    as.character(combined_anno$cell_type_subtype)
  } else {
    as.character(combined_anno$cell_type_full)
  },
  same_sample_neighbor_fraction = vapply(
    seq_len(nrow(review_nn)),
    function(i) mean(review_sample[review_nn[i, ]] == review_sample[i]),
    numeric(1)
  ),
  local_sample_diversity = vapply(
    seq_len(nrow(review_nn)),
    function(i) review_entropy(review_sample[review_nn[i, ]]),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

batch_mixing_sample_summary <- batch_mixing_cell %>%
  group_by(sample) %>%
  summarise(
    cells = n(),
    mean_same_sample_neighbor_fraction = mean(same_sample_neighbor_fraction),
    median_same_sample_neighbor_fraction = median(same_sample_neighbor_fraction),
    mean_local_sample_diversity = mean(local_sample_diversity),
    median_local_sample_diversity = median(local_sample_diversity),
    .groups = "drop"
  )

batch_mixing_overall_summary <- data.frame(
  cells = nrow(batch_mixing_cell),
  samples = length(unique(batch_mixing_cell$sample)),
  k_neighbors = review_k,
  mean_same_sample_neighbor_fraction = mean(batch_mixing_cell$same_sample_neighbor_fraction),
  median_same_sample_neighbor_fraction = median(batch_mixing_cell$same_sample_neighbor_fraction),
  mean_local_sample_diversity = mean(batch_mixing_cell$local_sample_diversity),
  median_local_sample_diversity = median(batch_mixing_cell$local_sample_diversity)
)

write.csv(
  batch_mixing_cell,
  file.path(review_batch_dir, "integrated_neighbor_mixing_per_cell.csv"),
  row.names = FALSE
)
write.csv(
  batch_mixing_sample_summary,
  file.path(review_batch_dir, "integrated_neighbor_mixing_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  batch_mixing_overall_summary,
  file.path(review_batch_dir, "integrated_neighbor_mixing_overall.csv"),
  row.names = FALSE
)

save_review_plot <- function(p, stem, width, height) {
  ggsave(file.path(review_batch_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, bg = "white")
  ggsave(file.path(review_batch_dir, paste0(stem, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
}

p_review_umap_celltype <- DimPlot(
  combined_anno, reduction = "umap", group.by = "cell_type",
  pt.size = 0.45, label = TRUE, repel = TRUE
) + ggtitle("Integrated UMAP by cell type") + theme_classic()
p_review_umap_sample <- DimPlot(
  combined_anno, reduction = "umap", group.by = "sample", pt.size = 0.45
) + ggtitle("Integrated UMAP by sample") + theme_classic()
p_review_umap_macro <- DimPlot(
  subset(combined_anno, subset = cell_type %in% c("LAM", "IAM", "OSM", "TRM")),
  reduction = "umap", group.by = "cell_type", pt.size = 0.65,
  label = TRUE, repel = TRUE
) + ggtitle("Macrophage-state UMAP") + theme_classic()

save_review_plot(p_review_umap_celltype, "review_integrated_UMAP_cell_type", 7, 5.5)
save_review_plot(p_review_umap_sample, "review_integrated_UMAP_sample", 7, 5.5)
save_review_plot(p_review_umap_macro, "review_macrophage_UMAP_subtype", 6, 5)

writeLines(
  c(
    "Reviewer supplement: integrated sample-mixing metrics",
    paste0("cells=", nrow(batch_mixing_cell)),
    paste0("samples=", length(unique(batch_mixing_cell$sample))),
    paste0("k_neighbors=", review_k),
    paste0(
      "mean_same_sample_neighbor_fraction=",
      signif(batch_mixing_overall_summary$mean_same_sample_neighbor_fraction, 4)
    ),
    paste0(
      "mean_local_sample_diversity=",
      signif(batch_mixing_overall_summary$mean_local_sample_diversity, 4)
    ),
    "These are nearest-neighbour sample-mixing metrics, not kBET/LISI."
  ),
  file.path(review_batch_dir, "batch_mixing_summary.txt")
)
