af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/', pattern = '.csv')
tblist <- sapply(af,function(f) {
  d <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/', f), as.is = T)
  rownames(d) <- d[,2]
  ud = unique(d[,1])
  marker = 'NA'
  if (grepl('Amplification', ud) & grepl('>=3', ud)) marker <- '+3'
  if (grepl('Amplification', ud) & grepl('>=4', ud)) marker <- '+4'
  if (grepl('Amplification', ud) & grepl('>=5', ud)) marker <- '+5'
  if (grepl('Amplification', ud) & grepl('>=6', ud)) marker <- '+6'
  if (grepl('Amplification', ud) & grepl('>=7', ud)) marker <- '+7'
  if (grepl('Amplification', ud) & grepl('>=8', ud)) marker <- '+8'
  if (grepl('Deletion: Homozygous', ud)) marker <- '-homo'
  if (grepl('Deletion: Loss of Heterozygosity', ud)) marker <- '-heter'
  colnames(d) <- paste0(colnames(d), marker)
  as.matrix(d[,-c(1:2)])
},simplify = F)

tb <- t(do.call(cbind,tblist))

tb <- tb[rowMeans(tb > 0.1) > 0.1, ]

tb <- tb[apply(tb,1,sd) > 0,]

pb <- prcomp(t(tb), scale=T)$x

library(ggplot2)
library(ggrepel)
library(cowplot)
pdtext <- data.frame(PC1 = pb[,1], PC2 = pb[,2], sample = rownames(pb), stringsAsFactors = F)


source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/plot/exome_sample_pca.pdf', width = 11, height = 12)
ggplot(data = pdtext, aes(x = PC1, y = PC2,label=sample)) + geom_point() + geom_text_repel()
dev.off()

library(pheatmap)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/plot/exome_sample_pattern.pdf', width = 20, height = 6)
pheatmap(t(tb), cluster_cols = F, cluster_rows = F)
dev.off()


