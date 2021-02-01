af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Darmanis_PNAS/raw')
res <- sapply(af,function(f) {
  d <- read.table(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Darmanis_PNAS/raw/',f),header=F)
  v <- d[,2]
  names(v) <- d[,1]
  v
})
library(GEOsearch)
meta <- SampleDetail('GSE67835')
samp <- sub('_.*','',colnames(res))
meta <- meta[match(samp,meta[,2]),]
mm <- t(sapply(meta$Characteristic,function(i) strsplit(i,'; ')[[1]],USE.NAMES = F))
# cn <- paste0('2015_Darmanis_PNAS:cell_',1:ncol(res))
m <- data.frame(cell=colnames(res),celltype=sub('cell type: ','',mm[,2]),time=sub('age: ','',mm[,3]),species='human',location=sub('tissue: ','',mm[,1]),stringsAsFactors = F)  
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Darmanis_PNAS/proc/count.rds')
saveRDS(m,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Darmanis_PNAS/proc/meta.rds')

