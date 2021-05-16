pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/dend_cluster/'
library(ape)
list <- list()
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM006/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM006.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM006']] <- NA  ### select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM009/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM009.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM009']] <- NA  ### 

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM010/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM010.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM010']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM013/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM013.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM013']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM029/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM029.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM029']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM030/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM030.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM030']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM035/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM035.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM035']] <- NA



d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM036/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM036.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM036']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM037/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM037.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM037']] <- NA  ## select non-malignant cells




d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM043/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM043.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM043']] <- NA ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM045/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM045.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM045']] <- NA


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM048/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM048.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM048']] <- NA ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM049/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM049.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM049']] <- NA  ## select non-malignant cells, mild mutations

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM050/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM050.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM050']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM051/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM051.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM051']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM052/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM052.odf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM052']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM054/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM054.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM054']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM055/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM055.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM055']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM056/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM056.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM056']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM059/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM059.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM059']] <- NA ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM060/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM060.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM060']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM062/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM062.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM062']] <- names(clu)[which(clu %in% c(6))]  ## select non-malignant cells !!!!!!

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM064/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM064.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM064']] <- names(clu)[which(clu %in% c(3,4))]   ## weak !!!!!



d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM065/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM065.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM065']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM066/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM066.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM066']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM068/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM068.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM068']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM069/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM069.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM069']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM070/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM070.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM070']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM071/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM071.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM071']] <- NA  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM073/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM073.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM073']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM074/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM074.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM074']] <- NA ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM075/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM075.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM075']] <- NA  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM081/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM081.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM081']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM082/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM082.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
#list[['GBM082']] <- names(clu)[which(clu %in% c(1:3))]  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM087/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM087.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM087']] <- NA  ## select non-malignant cells


rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/res/'
dir.create(rdir)
saveRDS(list, paste0(rdir, 'nonmaglinant_cells_high_confident.rds'))





