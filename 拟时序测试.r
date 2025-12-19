"拟时序报错"
library(Seurat)
library(monocle3)
library(Matrix)


pbmc.integrated$celltype <- Idents(pbmc.integrated)
table(Idents(pbmc.integrated))
#画个小提琴图，展示每个细胞类型的 nCount_RNA 分布
VlnPlot(pbmc.integrated, features = "nCount_RNA", group.by = "celltype")

pbmc_cd8 <- subset(
  pbmc.integrated,
  idents = c("Naive CD8+ T", "Cytotoxic CD8+ T")
)
pbmc_cd8 <- JoinLayers(
  object = pbmc_cd8,
  assay = "RNA",
  layers = c("counts", "data")
)
counts_mat <- GetAssayData(
  pbmc_cd8,
  assay = "RNA",
  slot = "counts"
)

gene_anno <- data.frame(
  gene_short_name = rownames(counts_mat),
  row.names = rownames(counts_mat)
)

cell_anno <- pbmc_cd8@meta.data

cds_cd8 <- new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata   = cell_anno,
  gene_metadata   = gene_anno
)

#预处理
cds_cd8 <- preprocess_cds(
  cds_cd8,
  num_dim = 30,
  method = "PCA"
)

#umap降维
cds_cd8 <- reduce_dimension(
  cds_cd8,
  reduction_method = "UMAP"
)

#聚类
cds_cd8 <- cluster_cells(
  cds_cd8,
  reduction_method = "UMAP"
)

#用“Naive CD8+ T”作为根
root_cells <- colnames(cds_cd8)[
  colData(cds_cd8)$celltype == "Naive CD8+ T"
]

#学习轨迹图
cds_cd8 <- learn_graph(cds_cd8)

cds_cd8 <- order_cells(
  cds_cd8,
  reduction_method = "UMAP",
  root_cells = root_cells
)

#拟时序图
plot_cells(
  cds_cd8,
  color_cells_by = "pseudotime",
  label_leaves = TRUE,
  label_branch_points = TRUE
)

#细胞类型着色
plot_cells(
  cds_cd8,
  color_cells_by = "celltype"
)