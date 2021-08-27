list = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/nonmaglinant_cells.rds')
cell = unlist(list)
cell = cell[!is.na(cell)]
clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
patient = gsub('_.*', '', names(clu))  

samp <- unique(sub('_.*','',cell))
clusub <- clu[patient %in% samp]
t1 <- table(clusub[names(clusub) %in% cell])
t2 <- table(clusub)
res = sort(t1/t2[names(t1)])
write.csv(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/nonmalignant_cell_proportion_in_clusters.csv')

library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/plot/nonmalignant_proportion_in_clusters.pdf', width = 4, height = 3.5)
pd = data.frame(cluster = paste0('cluster ', names(res)), prop = as.numeric(res), stringsAsFactors = F)
pd[,1] = factor(pd[,1], levels = pd[,1])

ggplot(data = pd, aes(x = cluster, y = prop)) +
  geom_col() +
  ylab('non-malignant cells proportion') +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()
