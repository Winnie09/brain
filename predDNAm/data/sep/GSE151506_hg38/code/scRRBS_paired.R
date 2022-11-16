r <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/unpaired.rds')
p <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/pairing/pair.rds')
r <- r[,p[,2]]
colnames(r) <- p$newcellname
r <- r[rowMeans(is.na(r)) < 1,]
sp <- do.call(rbind,strsplit(rownames(r),'-'))
rownames(r) <- paste0(sp[,1],'_',as.numeric(sp[,2]))
saveRDS(r,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/paired.rds')
