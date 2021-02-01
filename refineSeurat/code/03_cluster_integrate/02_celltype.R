clu1 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu1.rds')
clu2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/clu2.rds')

fullmeta = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/old/20200212/result/full/meta.rds')
meta1 = fullmeta[match(names(clu1), fullmeta$cell), 'celltype_super']
meta2 = fullmeta[match(names(clu2), fullmeta$cell), 'celltype_super']

ct1 <- sapply(1:max(clu1), function(c) {
    tab = table(meta1[match(c,clu1)])
  names(tab)[which.max(tab)]
})

ct2 <- sapply(1:max(clu2), function(c) {
    tab = table(meta2[match(c,clu2)])
  names(tab)[which.max(tab)]
})
names(ct1) = paste0('clu',1:max(clu1))
names(ct2) = paste0('clu',1:max(clu2))
saveRDS(ct1,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/ct1.rds')
saveRDS(ct2,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/ct2.rds')

