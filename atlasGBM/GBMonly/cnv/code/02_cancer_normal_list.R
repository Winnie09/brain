pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/dend_cluster/'
library(ape)
list <- list()
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM006/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM006']] <- names(clu)[which(clu %in% c(1:3,9))]  ### select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM009/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM009']] <- names(clu)[which(clu %in% c(1:5,8:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM010/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM010']] <- names(clu)[which(clu %in% c(4:6))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM013/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM013']] <- names(clu)[which(clu %in% c(5,7:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM029/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM029']] <- names(clu)[which(clu %in% c(1:2,7:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM030/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM030']] <- names(clu)[which(clu %in% c(2:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM035/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM035']] <- names(clu)[which(clu %in% setdiff(1:10,c(7,9)))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM045/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM045']] <- names(clu)[which(clu %in% setdiff(1:10,9))]


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM036/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM036']] <- names(clu)[which(clu %in% c(4))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM043/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM043']] <- names(clu)[which(clu %in% c(7))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM048/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM048']] <- names(clu)[which(clu %in% setdiff(1:10, 1))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM050/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM050']] <- names(clu)[which(clu %in% setdiff(1:10, 6:7))]  ## select cancer cells
  
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM052/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM052']] <- names(clu)[which(clu %in% setdiff(1:10, 1:2))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM055/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM055']] <- names(clu)[which(clu %in% setdiff(1:10, 2:4))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM062/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM062']] <- names(clu)[which(clu %in% c(1:5))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM065/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM065']] <- names(clu)[which(clu %in% c(3:7))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM069/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM069']] <- names(clu)[which(clu %in% c(5:10))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM071/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM071']] <- names(clu)[which(clu %in% c(7))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM074/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM074']] <- names(clu)[which(clu %in% setdiff(1:10, 1))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM081/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM081']] <- NULL


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM087/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM087']] <- names(clu)[which(clu %in% c(1:5))]  ## select cancer cells



d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM037/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM037']] <- names(clu)[which(clu %in% c(1))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM049/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM049']] <- names(clu)[which(clu %in% c(1:8))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM051/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM051']] <- names(clu)[which(clu %in% c(1:8))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM054/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM054']] <- names(clu)[which(clu %in% setdiff(1:10,9))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM056/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM056']] <- names(clu)[which(clu %in% c(3:10))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM059/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM059']] <- NULL ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM064/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM064']] <- NULL  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM068/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM068']] <- names(clu)[which(clu %in% c(5:10))]  ## select cancer cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM070/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM070']] <- names(clu)[which(clu %in% c(1:7))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM073/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM073']] <- names(clu)[which(clu %in% setdiff(1:10, 8:9))]  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM075/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM075']] <- NULL  ## select cancer cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM082/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM082']] <- names(clu)[which(clu %in% c(1:3))]  ## select cancer cells

rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/res/'
dir.create(rdir)
saveRDS(list, paste0(rdir, 'cancer_cells.rds'))


