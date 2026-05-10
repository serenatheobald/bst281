library(data.table)
library(Seurat)
library(Matrix)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)

data_dir <- "/Users/serenatheobald/Downloads/"

count_file <- file.path(data_dir, "GSE144482_gene_by_cell_count_matrix.txt.gz")
cell_anno_file <- file.path(data_dir, "GSE144482_cell_annotation.csv.gz")
gene_anno_file <- file.path(data_dir, "GSE144482_gene_annotation.csv.gz")

counts_raw <- Matrix::readMM(gzfile(count_file))

cell_anno <- read.csv(cell_anno_file, header = TRUE, row.names = 1, check.names = FALSE)
gene_anno <- read.csv(gene_anno_file, header = TRUE, row.names = 1, check.names = FALSE)

gene_names <- make.unique(as.character(gene_anno[["gene_short_name"]]))
cell_names <- make.unique(as.character(cell_anno[["barcode"]]))

rownames(counts_raw) <- gene_names
colnames(counts_raw) <- cell_names

runx1 <- CreateSeuratObject(
  counts = counts_raw,
  meta.data = cell_anno,
  project = "GSE144482_Runx1_LSK",
  min.cells = 3,
  min.features = 200
)

runx1 <- NormalizeData(runx1, verbose = FALSE)
runx1 <- FindVariableFeatures(runx1, nfeatures = 2000, verbose = FALSE)
runx1 <- ScaleData(runx1, features = VariableFeatures(runx1), verbose = FALSE)
runx1 <- RunPCA(runx1, features = VariableFeatures(runx1), npcs = 30, verbose = FALSE)
runx1 <- FindNeighbors(runx1, reduction = "pca", dims = 1:20, verbose = FALSE)
runx1 <- FindClusters(runx1, resolution = 0.5, verbose = FALSE)
runx1 <- RunUMAP(runx1, reduction = "pca", dims = 1:20, verbose = FALSE)

print(runx1)

cds <- as.cell_data_set(runx1)

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)

cd34_cells <- rownames(colData(cds))[
  colData(cds)$matched_ctype_refined == "CD34"
]

print(length(cd34_cells))
print(head(cd34_cells))

cds <- order_cells(cds, root_cells = cd34_cells)

p1 <- plot_cells(cds, color_cells_by = "Dataset", label_cell_groups = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)
p3 <- plot_cells(cds, color_cells_by = "matched_ctype_refined", label_cell_groups = FALSE)

ggsave("~/monocle_dataset_plot.png", p1, width = 7, height = 5, dpi = 300)
ggsave("~/monocle_pseudotime_plot.png", p2, width = 7, height = 5, dpi = 300)
ggsave("~/monocle_celltype_plot.png", p3, width = 7, height = 5, dpi = 300)

stat_df <- data.frame(
  cell = rownames(colData(cds)),
  pseudotime = pseudotime(cds),
  Dataset = colData(cds)$Dataset,
  celltype = colData(cds)$matched_ctype_refined,
  cluster = colData(cds)$seurat_clusters
)

stat_df <- stat_df[is.finite(stat_df$pseudotime), ]

pseudotime_test <- wilcox.test(
  pseudotime ~ Dataset,
  data = stat_df
)

celltype_table <- table(
  stat_df$celltype,
  stat_df$Dataset
)

celltype_chisq <- chisq.test(celltype_table)

cluster_table <- table(
  stat_df$cluster,
  stat_df$Dataset
)

cluster_chisq <- chisq.test(cluster_table)

print(pseudotime_test)
print(celltype_chisq)
print(cluster_chisq)

test_summary <- data.frame(
  Test = c(
    "Wilcoxon rank-sum test",
    "Chi-square test: cell type occupancy",
    "Chi-square test: cluster occupancy"
  ),
  Purpose = c(
    "Compare pseudotime distributions between WT and KO cells",
    "Test whether WT and KO differ in annotated cell-type occupancy",
    "Test whether WT and KO differ in Seurat cluster occupancy"
  ),
  Justification = c(
    "Wilcoxon test was used because pseudotime is continuous, non-normal, and often skewed in trajectory analyses",
    "Chi-square test was used because cell-type occupancy is categorical count data across genotype groups",
    "Chi-square test was used because cluster occupancy is categorical count data across genotype groups"
  ),
  Statistic = c(
    paste("W =", as.numeric(pseudotime_test$statistic)),
    paste("X-squared =", round(as.numeric(celltype_chisq$statistic), 2)),
    paste("X-squared =", round(as.numeric(cluster_chisq$statistic), 2))
  ),
  Degrees_of_freedom = c(
    NA,
    as.numeric(celltype_chisq$parameter),
    as.numeric(cluster_chisq$parameter)
  ),
  P_value = c(
    pseudotime_test$p.value,
    celltype_chisq$p.value,
    cluster_chisq$p.value
  ),
  Interpretation = c(
    "Tests whether WT and KO cells differ in their overall relative pseudotime positions",
    "Tests whether WT and KO cells are distributed differently across author-annotated differentiation states",
    "Tests whether WT and KO cells are distributed differently across unsupervised Seurat clusters"
  )
)

print(test_summary)

write.csv(stat_df, "~/monocle_pseudotime_metadata.csv", row.names = FALSE)
write.csv(as.data.frame(celltype_table), "~/monocle_celltype_dataset_table.csv", row.names = FALSE)
write.csv(as.data.frame(cluster_table), "~/monocle_cluster_dataset_table.csv", row.names = FALSE)
write.csv(test_summary, "~/monocle_test_summary.csv", row.names = FALSE)

sink("~/monocle_statistical_tests.txt")

cat("Wilcoxon Test: Pseudotime by Dataset\n\n")
print(pseudotime_test)

cat("\n\nChi-square Test: Cell Type by Dataset\n\n")
print(celltype_chisq)

cat("\n\nChi-square Test: Cluster by Dataset\n\n")
print(cluster_chisq)

cat("\n\nSummary Table\n\n")
print(test_summary)

sink()

saveRDS(runx1, "~/runx1_seurat.rds")
saveRDS(cds, "~/runx1_monocle_cds.rds")

print("DONE")