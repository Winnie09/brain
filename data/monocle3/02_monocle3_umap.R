expression_matrix <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/proc/expr_sub.rds')
cell_metadata <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/meta/meta_sub.rds')
rownames(cell_metadata) <- cell_metadata[,1]
gene_annotation <- data.frame(gene = rownames(expression_matrix),
                              order = seq(1, nrow(expression_matrix)),
                              stringsAsFactors = FALSE)
rownames(gene_annotation) <- rownames(expression_matrix)
library(monocle3)
library(ggplot2)
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/monocle3/'
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)

cds <- align_cds(cds, alignment_group = "donor")

cds <- reduce_dimension(cds)

ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")


cds <- learn_graph(cds)

p <- plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
ggsave(paste0(pdir, 'monocle3_umap.png'), p, width = 8, height = 6, dpi = 200)
