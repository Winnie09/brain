library(data.table)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_zhong_nature/raw/expr')
d <- lapply(af,function(f) {
  d <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_zhong_nature/raw/expr/',f),data.table=F)
  row.names(d) <- d[,1]
  d <- as.matrix(d[,-1])
})
d <- do.call(cbind,d)
saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_zhong_nature/proc/expr.rds')

m <- read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_zhong_nature/raw/Sample.txt',as.is=T,header=T,sep='\t')
m <- m[,c('Title','Characteristic')]
m[,1] <- sub('_[0-9]*$','',m[,1])
ms <- t(sapply(m[,2],function(i) strsplit(i,'; ')[[1]],USE.NAMES = F))
ms <- sub('.*: ','',ms)
m <- cbind(m[,1],ms[,-4])
colnames(m) <- c('exp','Location','Time','Gender')
m <- unique(m)
mv <- rep(NA,ncol(d))
for (i in m[,'exp']) {
  mv[grep(i,colnames(d))] <- i
}
sm <- m[match(mv,m[,1]),]
sm[,1] <- colnames(d)
colnames(sm)[1] <- 'Cell'
saveRDS(sm,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_zhong_nature/proc/meta.rds')
