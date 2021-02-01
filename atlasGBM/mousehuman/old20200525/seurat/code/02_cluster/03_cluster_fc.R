setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result')
key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'atlasOnly_homolog_sub'
# key = 'atlas_GBM_homolog_sub'
library(Seurat)
d = readRDS(paste0('./',key,'/integrated.rds'))
count = as.matrix(d@assays$RNA@counts)
dlist = list()
dlist[['human']] = readRDS(paste0('./',key,'/cluster_mean/species/human.rds'))
dlist[['mouse']] = readRDS(paste0('./',key,'/cluster_mean/species/mouse.rds'))
a <- sapply(names(dlist), function(n){
  print(n)
  d = dlist[[n]]
  library(preprocessCore)
dn <- dimnames(d)
d <- normalize.quantiles(d)
dimnames(d) <- dn
gl <- sapply(colnames(d), function(i){
  v = setdiff(colnames(d), i)
  tmpmat <- sapply(v, function(j){
    tmp = d[, colnames(d) == i] - d[, colnames(d) == j]    
  })
  u = sort(rowMeans(tmpmat), decreasing = T)[1:100]
  paste0(names(u), '(', round(u,3), ')')
})
#colnames(gl) = paste0('clu',colnames(gl))
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/',key,'/cluster_fc/species/'), showWarnings = F, recursive = T)  
saveRDS(gl,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/',key,'/cluster_fc/species/',n,'.rds'))  
gl = gl[,order(as.numeric(gsub('clu','',colnames(gl))))]
write.csv(gl, quote=F, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/',key,'/cluster_fc/species/',n,'.csv'))

})

