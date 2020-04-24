library(data.table)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/raw')
data <- sapply(af,function(f) {
  d <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/raw/',f),data.table = F)
  row.names(d) <- d[,1]
  as.matrix(d[,-1])
})
data <- do.call(cbind,data)
region <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/meta/region.csv',header=F)
regv <- sub(' *$','',paste(region[,3],region[,1]))
names(regv) <- region[,2]

cn <- paste0('2018_Fan_CellResearch:cell_',1:ncol(data))
cnm <- t(sapply(colnames(data),function(i) strsplit(i,'_')[[1]][1:2],USE.NAMES = F))
meta <- data.frame(cell=cn,location=regv[cnm[,1]],time=paste0('mid-gestation embryo week ',sub('W.*','',cnm[,2])),gender=sub('.*W','',cnm[,2]),species='human')
colnames(data) <- cn

saveRDS(data,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')
saveRDS(meta,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/meta.rds')

