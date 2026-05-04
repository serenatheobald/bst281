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



cd34_cells <- rownames(colData(cds))[colData(cds)$matched_ctype_refined == "CD34"]

print(length(cd34_cells))
print(head(cd34_cells))

cds <- order_cells(cds, root_cells = cd34_cells)

print(colnames(colData(cds)))

p1 <- plot_cells(cds, color_cells_by = "Dataset", label_cell_groups = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)
p3 <- plot_cells(cds, color_cells_by = "matched_ctype_refined", label_cell_groups = FALSE)

ggsave("~/monocle_dataset_plot.png", p1, width = 7, height = 5, dpi = 300)
ggsave("~/monocle_pseudotime_plot.png", p2, width = 7, height = 5, dpi = 300)
ggsave("~/monocle_celltype_plot.png", p3, width = 7, height = 5, dpi = 300)

saveRDS(runx1, "~/runx1_seurat.rds")
saveRDS(cds, "~/runx1_monocle_cds.rds")

print("DONE")