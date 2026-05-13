README.txt

BST281 — Module 3 Supplementary Data Packet
Serena Theobald
Email: stheobald@hsph.harvard.edu

Project Title:
Runx1 Deficiency Alters Hematopoietic Differentiation-State Occupancy and Trajectory-Associated Transcriptional Programs

Overview:
This supplementary packet contains the data files, R scripts, processed outputs, statistical summaries, and figures used for Module 3 trajectory analysis of the GSE144482 single-cell RNA-seq dataset. The analysis focused on reconstructing hematopoietic differentiation trajectories using Monocle3 to determine how Runx1 deficiency alters developmental-state occupancy, pseudotime progression, and inflammatory lineage-associated transcriptional programs.

Dataset Source:
The single-cell RNA-seq dataset was obtained from the NCBI Gene Expression Omnibus (GEO):

GSE144482
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144482

The dataset contains wild-type (WT) and Runx1 knockout (KO) hematopoietic stem and progenitor cells with accompanying refined cell-type annotations.

Input Data Files:
- GSE144482_gene_by_cell_count_matrix.txt.gz
  Sparse gene-by-cell count matrix used for downstream preprocessing and trajectory reconstruction.

- GSE144482_cell_annotation.csv.gz
  Cell-level metadata including genotype labels and refined hematopoietic cell-type annotations.

- GSE144482_gene_annotation.csv.gz
  Gene annotation metadata associated with the count matrix.

Code Files:
- runx1_monocle3_trajectory_and_gene_analysis copy.R
  Main analysis workflow used for preprocessing, dimensionality reduction, Monocle3 trajectory inference, pseudotime ordering, trajectory-associated gene analysis, and figure generation.

- runx1_monocle3_stats_analysis copy.R
  Statistical analysis workflow used for pseudotime, cluster occupancy, and cell-type occupancy

  

Output Files:
- monocle_trajectory_genes_all.csv
  Complete list of trajectory-associated genes identified using Monocle3 graph_test() with Moran’s I statistics and Benjamini–Hochberg false discovery rate (FDR) correction.

- monocle_top100_trajectory_genes.csv
  Top 100 significantly trajectory-associated genes ranked by adjusted q-value.

- monocle_test_summary.csv
  Summary table of Wilcoxon rank-sum and chi-square statistical tests.

- monocle_pseudotime_metadata.csv
  Cell-level pseudotime assignments and metadata exported from Monocle3.

- monocle_celltype_dataset_table.csv
  Cell-type occupancy table comparing WT and KO distributions.

- monocle_cluster_dataset_table.csv
  Cluster occupancy table comparing WT and KO distributions.

- monocle_statistical_tests.txt
  Full statistical test outputs generated during analysis.

Figure Files:
- monocle_dataset_plot.png
  Shared Monocle3 trajectory colored by genotype.

- monocle_pseudotime_plot.png
  Shared trajectory colored by inferred pseudotime progression.

- monocle_celltype_plot.png
  Trajectory colored by refined hematopoietic cell-type annotations.

- monocle_top_trajectory_genes_pseudotime.png
  Representative trajectory-associated genes plotted across pseudotime.

Software and Packages:
Analyses were performed in R using:
- Seurat v5
- Monocle3 v1.3
- SeuratWrappers
- Matrix
- ggplot2

Analysis Summary:
WT and KO cells were jointly embedded into a shared Monocle3 trajectory landscape. CD34-positive progenitor cells were used as the trajectory root population for pseudotime inference. Trajectory-associated genes were identified using Monocle3 graph_test() with Benjamini–Hochberg FDR correction. Additional statistical analyses included Wilcoxon rank-sum testing for pseudotime comparisons and chi-square testing for cell-type and cluster occupancy differences between WT and KO cells.
