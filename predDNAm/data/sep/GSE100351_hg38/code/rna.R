source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
library(data.table)
rf <- list.files('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/raw/',pattern = 'rna.RPKM.gz')
cf <- list.files('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/raw/',pattern = 'cpgmeth.bed.gz')

names(rf) <- sub('GSM[0-9]*_','',sub("_rna.RPKM.gz",'',rf))
names(cf) <- sub('GSM[0-9]*_','',sub("_cpgmeth.bed.gz",'',cf))
int <- intersect(names(rf),names(cf))
# m <- fread('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/meta/meta.txt',data.table=F)
# dim(m[m$Description%in%int,])
rf <- rf[int]
cf <- cf[int]

rna <- sapply(rf,function(f) {
  d <- fread(paste0('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/raw/',f),data.table=F)
  tapply(d$RPKM,list(d[,1]),sum)
  #tapply(d$RPKM/d$transcript_length*1000,list(d[,1]),sum)
},simplify = F)
rna <- do.call(cbind,rna)
tpm <- rpkm2tpm(rna)
tpm <- log2(tpm + 1)
saveRDS(tpm,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/proc/tpm.rds')

