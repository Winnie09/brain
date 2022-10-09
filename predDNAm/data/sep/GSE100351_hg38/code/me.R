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

me <- sapply(cf,function(f) {
  d <- fread(paste0('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/raw/',f),data.table=F)
  v <- d[,5]/1000
  names(v) <- paste0(d[,1],'_',d[,2]+1)
  v
},simplify = F)
ag <- unique(unlist(sapply(me,names)))
mm <- matrix(NA,nrow=length(ag),ncol=length(me),dimnames = list(ag,names(me)))
for (i in names(me))
  mm[names(me[[i]]),i] <- me[[i]]
saveRDS(mm,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/proc/me.rds')

