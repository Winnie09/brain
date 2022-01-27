pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/dend_cluster/'
dir.create(pdir)
cludir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/cnvclu/'
library(ape)
ap = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/summary/', pattern = 'GBM')
ap = sub('.png', '', ap)

for (p in setdiff(ap, c('GBM069', 'GBM074'))){
  print(p)
  d <- read.tree(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/', p,'/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
  clu <- cutree(as.hclust(d),10)
  saveRDS(clu, paste0(cludir, p, '.rds'))
  pdf(paste0(pdir, p, '.pdf'), width = 5, height = 7)
  print(plot(d, tip.color=clu, cex = 0.5))
  dev.off()
}
  

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM069/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),15)
saveRDS(clu, paste0(cludir, 'GBM069.rds'))
pdf(paste0(pdir, 'GBM069.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM074/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),20)
saveRDS(clu, paste0(cludir, 'GBM074.rds'))
pdf(paste0(pdir, 'GBM074.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()

###########
list <- list()
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/GBM006/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM006']] <- names(clu)[clu == 1] ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM009/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM009']] <- names(clu)[which(clu %in% c(1:2))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM010/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM010']] <- names(clu)[which(clu %in% c(1))]

## s/p immunotherapy
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM013/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM013']] <- names(clu)[which(clu %in% c(7,9,10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM029/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM029']] <- names(clu)[which(clu %in% c(1,3))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM030/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM030']] <- names(clu)[which(clu %in% c(1:2))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM035/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM035']] <- names(clu)[which(clu %in% c(1:3))]



d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM036/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM036']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM037/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM037']] <- names(clu)[which(clu %in% c(2))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM043/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM043.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM043']] <- names(clu)[which(clu %in% c(2))] 

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM045/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM045']] <- names(clu)[which(clu %in% c(4))]


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM048/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM048']] <- names(clu)[which(clu %in% c(9))]  ## select non-malignant cells


d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM049/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM049']] <- names(clu)[which(clu %in% c(3))]  ## select non-malignant cells, mild mutations
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM050/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
pdf(paste0(pdir, 'GBM050.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM050']] <- names(clu)[which(clu %in% c(6))] 
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM051/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM051']] <- names(clu)[which(clu %in% c(7))]   ## select non-malignant cells
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM052/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM052']] <- names(clu)[which(clu %in% c(1,2))]   ## select non-malignant cells
# 
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM054/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM054']] <- names(clu)[which(clu %in% c(7))]   ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM055/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM055']] <- names(clu)[which(clu %in% c(10))]  ## select non-malignant cells

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM056/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM056']] <- names(clu)[which(clu %in% c(3,4))]

### new
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM057/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM057']] <- names(clu)[which(clu %in% c(6))]  ## select non-malignant cells



d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM059/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM059']] <- names(clu)[which(clu %in% c(8,9))]  ## select non-malignant cells
# 
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM060/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM060']] <- names(clu)[which(clu %in% c(6,7))]  ## select non-malignant cells

## s/p immunotherapy
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM062/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM062']] <- NA

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM064/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM064']] <- NA



d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM065/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM065']] <- NA
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM066/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM066']] <- names(clu)[which(clu %in% c(1,2))]  ## select non-malignant cells
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM068/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM068']] <- names(clu)[which(clu %in% c(5,6))]  ## select non-malignant cells
# 
# ## difficult ! half of cluster 10 seems to be normal
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM069/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),15)
pdf(paste0(pdir, 'GBM069.pdf'), width = 5, height = 7)
plot(d, tip.color=clu, cex = 0.5)
dev.off()
list[['GBM069']] <- names(clu)[which(clu %in% c(14))]  ## select non-malignant cells
# 
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM070/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM070']] <- names(clu)[which(clu %in% c(2))]  ## select non-malignant cells
# 
## s/p immunotherapy
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM071/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM071']] <- NA  ## select non-malignant cells
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM073/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM073']] <- NA  ## select non-malignant cells
# 
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM074/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),20)
list[['GBM074']] <- names(clu)[which(clu %in% c(17:20))]
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM075/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM075']] <- names(clu)[which(clu %in% c(8,9,10))]   
# 
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM081/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM081']] <- names(clu)[which(clu %in% c(1,2))]
# 
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM082/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM082']] <- names(clu)[which(clu %in% c(8))]


## s/p immunotherapy
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/GBM087/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
list[['GBM087']] <- NA  ## This sample is super wiered


rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/'
dir.create(rdir)
saveRDS(list, paste0(rdir, 'nonmaglinant_cells.rds'))



