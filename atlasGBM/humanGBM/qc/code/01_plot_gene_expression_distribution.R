library(here)
here()
setwd(here())
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanGBM/qc/plot/'
mat =  readRDS('./atlasGBM/humanGBM/data/combine_mat.rds')
meta =  readRDS('./atlasGBM/humanGBM/data/meta_allcell.rds')
rownames(meta) <- meta$cell
# selectcell <- meta$cell[!grepl('GBM', meta$study) & !grepl('2014_Pollen_NatBiotech', meta$study)] ####
# mat <- mat[,selectcell]
mat <- mat[, colSums(mat > 0) > 500]
mat <- mat[!grepl('^MT-',rownames(mat)),]
mat <- mat[rowMeans(mat > 0) >= 0.01,]
meta <- meta[colnames(mat), ]
str(mat)

library(ggplot2)
library(RColorBrewer)

samp <- sapply(sort(unique(meta$study)), function(i){
    tmp <- mat[, rownames(meta)[meta$study == i]]
    v = as.vector(tmp)
    set.seed(12345)
    sample(v, 1e5)  
  })
colnames(samp) <- sort(unique(meta$study))
pd <- reshape2::melt(samp)

colnames(pd) <- c('samplecell','study','expression')

pdf(paste0(pdir, 'distribution_sampled_gene_expression.pdf'), width = 7, height = 12)
ggplot(data = pd, aes(x = study, y = expression, fill = study)) +
  # geom_jitter(alpha = 0.2, size = 0.2) +
  geom_violin(scale = 'area', trim = T, alpha = 0.5) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab('sampled gene expression') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, 'Dark2'))(length(unique(pd$study)))) +
  coord_flip()
dev.off()

## only use non-zero
samp <- sapply(sort(unique(meta$study)), function(i){
    tmp <- mat[, rownames(meta)[meta$study == i]]
    v = as.vector(tmp)
    v <- v[v != 0]
    set.seed(12345)
    sample(v, 1e5)  
  })
colnames(samp) <- sort(unique(meta$study))
pd <- reshape2::melt(samp)

colnames(pd) <- c('samplecell','study','expression')

pdf(paste0(pdir, 'distribution_sampled_gene_expression_excluede_zeros.pdf'), width = 7, height = 12)
ggplot(data = pd, aes(x = study, y = expression, fill = study)) +
  geom_violin(scale = 'area', trim = T, alpha = 0.5) +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, 'Dark2'))(length(unique(pd$study)))) +
  ylab('sampled gene expression (exclude zeros)') +
  coord_flip()
dev.off()
