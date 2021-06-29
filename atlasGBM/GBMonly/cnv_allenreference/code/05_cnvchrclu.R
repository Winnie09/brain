library(ape)
ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/'
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/res/cnvchrclu/'
dir.create(rdir)
p = 'GBM006'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1))]
group2 = names(clu)[which(clu %in% c(2:7, 10))]
group3 = names(clu)[which(clu %in% c(8))]
group4 = names(clu)[which(clu %in% c(9))]
myclu = rep(1:4, c(length(group1), length(group2), length(group3), length(group4)))
names(myclu) =  c(group1, group2, group3, group4)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM009'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:2))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM010'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1))]
group2 = names(clu)[which(clu %in% c(2))]
group3 = names(clu)[which(clu %in% c(3:10))]
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM013'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(2:4))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))



p = 'GBM029'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:2))]
group2 = names(clu)[which(clu %in% c(4))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM030'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:2))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM035'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(2:3))]
group2 = names(clu)[which(clu %in% c(1,4,5))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM036'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(7))]
group2 = names(clu)[which(clu %in% c(8))]
group3 = names(clu)[which(clu %in% c(1:2))]
group4 = names(clu)[which(clu %in% c(3:4))]
group5 = setdiff(names(clu), c(group1, group2, group3, group4))
myclu = rep(1:5, c(length(group1), length(group2), length(group3), length(group4), length(group5)))
names(myclu) =  c(group1, group2, group3, group4, group5)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM037'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(2))]
group2 = names(clu)[which(clu %in% c(1))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM043'  
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(7))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM045'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(6,7))]
group2 = names(clu)[which(clu %in% c(1))]
group3 = names(clu)[which(clu %in% c(8:10))]
group4 = setdiff(names(clu), c(group1, group2, group3))
myclu = rep(1:4, c(length(group1), length(group2), length(group3), length(group4)))
names(myclu) =  c(group1, group2, group3, group4)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM048'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:2))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM049'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1))]
group2 = names(clu)[which(clu %in% c(2:3))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM050'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(7))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM051'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(10))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM052'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:3))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM054'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(9))]
group2 = names(clu)[which(clu %in% c(1:4, 7,8,10))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM055'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
myclu = rep(1, length(clu))
names(myclu) =  names(clu)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM056' 
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(3:6))]
group2 = names(clu)[which(clu %in% c(7:10))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM059'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(8))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM060'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(7:10))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))



p = 'GBM062'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(6))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM064'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:4))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM065'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(2))]
group2 = names(clu)[which(clu %in% c(1,3,4))]
group3 = setdiff(names(clu), c(group1, group2))
myclu = rep(1:3, c(length(group1), length(group2), length(group3)))
names(myclu) =  c(group1, group2, group3)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM066'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(7,8))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM068'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:3))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM069'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:4))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM070'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
myclu = rep(1, length(clu))
names(myclu) =  names(clu)
saveRDS(myclu, paste0(rdir, p, '.rds'))


p = 'GBM071'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
myclu = rep(1, length(clu))
names(myclu) =  names(clu)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM073'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
myclu = rep(1, length(clu))
names(myclu) =  names(clu)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM074'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1:3))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM075'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(1))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM081'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
myclu = rep(1, length(clu))
names(myclu) =  names(clu)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM082'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(6:7))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))

p = 'GBM087'
d <- read.tree(paste0(ddir,p, '/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
clu <- cutree(as.hclust(d),10)
group1 = names(clu)[which(clu %in% c(5))]
group2 = setdiff(names(clu), group1)
myclu = rep(1:2, c(length(group1), length(group2)))
names(myclu) =  c(group1, group2)
saveRDS(myclu, paste0(rdir, p, '.rds'))



