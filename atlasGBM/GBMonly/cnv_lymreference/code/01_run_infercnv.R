library(infercnv)
library(here)
setwd(here())
cut <- 0.1
meta <- readRDS('atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
allp <- unique(meta$Sample.ID)
# writeLines(allp, 'atlasGBM/GBMonly/cnv_myelymreference/code/unfinished.txt')

p <- as.character(commandArgs(trailingOnly = T)[[1]])

print(p)
anno <- read.table(paste0('atlasGBM/GBMonly/cnv_lymreference/data/data/36nonNormal_myelym_combined_', p, '_cellanno.txt'), sep = '\t', as.is = TRUE)
dir.create(paste0('atlasGBM/GBMonly/cnv_lymreference/res/',p,'/cutoff', cut, '/output'), recursive = TRUE, showWarnings = FALSE)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0('atlasGBM/GBMonly/cnv_lymreference/data/data/36nonNormal_myelym_combined_',p, '_countmatrix.txt'),
                                    annotations_file=paste0('atlasGBM/GBMonly/cnv_lymreference/data/data/36nonNormal_myelym_combined_', p, '_cellanno.txt'),
                                    delim="\t",
                                    gene_order_file='atlasGBM/GBMonly/cnv_lymreference/data/data/36nonNormal_myelym_combined_gr.txt',
                                    ref_group_names=unique(anno[!grepl('cluster', anno[,2]),2]))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=cut, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0('atlasGBM/GBMonly/cnv_lymreference/res/',p,'/cutoff', cut, '/output'), 
                             cluster_by_groups=FALSE, 
                             denoise=TRUE,
                             HMM=FALSE)

