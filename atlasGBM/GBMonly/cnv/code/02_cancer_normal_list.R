library(ape)
list <- list()
d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM006/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM006']] <- names(clu)[which(clu %in% c(1:3,9))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM009/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM009']] <- names(clu)[which(clu %in% c(1:5,8:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM010/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM010']] <- names(clu)[which(clu %in% c(4:6))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM013/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM013']] <- names(clu)[which(clu %in% c(5,7:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM029/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM029']] <- names(clu)[which(clu %in% c(1:2,7:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM030/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM030']] <- names(clu)[which(clu %in% c(2:10))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM035/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM035']] <- names(clu)[which(clu %in% setdiff(1:10,c(7,9)))]

d <- read.tree('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/GBM045/cutoff0.1/output/infercnv.observations_dendrogram.txt')
clu <- cutree(as.hclust(d),10)
#plot(d,tip.color=clu)
list[['GBM045']] <- names(clu)[which(clu %in% setdiff(1:10,9))]

rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv/res/'
dir.create(rdir)
saveRDS(list, paste0(rdir, 'cancer_normal_list.rds'))
