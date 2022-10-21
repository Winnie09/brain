d <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/paired.rds')
s <- sub(':.*','',colnames(d))
m <- sapply(unique(s),function(i) {
  id <- which(s==i)
  rowMeans(d[,id,drop=F],na.rm=T)
})
saveRDS(m,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/pb.rds')

