ddir2 <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/raw/nygc/CSV/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'
allf = list.files(ddir2, pattern = 'MGH')

library(data.table)
## read in all samples
alld <- lapply(allf, function(f) {
  print(f)
  me <- fread(paste0(ddir2, f),data.table=F)
  rownames(me) <- me[, 1]
  me <- me[, -1]
  me <- as.matrix(me)
  me <- matrix(as.numeric(me),nrow=nrow(me),dimnames = dimnames(me))
  me
})

a <- unique(unlist(sapply(alld,rownames)))
b <- unlist(sapply(alld,colnames))
m <- matrix(NA,nrow=length(a),ncol=length(b),dimnames = list(a,b))
for (i in 1:length(alld)) {
  m[rownames(alld[[i]]),colnames(alld[[i]])] <- alld[[i]]
}

s <- sub('-.*','',rownames(m))
m <- m[nchar(s) <= 5,]
colnames(m) <- sub('.R1.fastq_bismark_bt2_pe.bismark.cov','',colnames(m))
saveRDS(m,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/unpaired.rds')

