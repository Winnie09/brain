d <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/mutation.csv',as.is=T)
d <- d[grep('^GBM',d$SampleID),]
d$SampleID <- sub('WES','',sub('T.*','',sub('PT.*','',d$SampleID)))
d <- d[,c('SampleID','Gene.Sym')]
d <- unique(d)
d <- d[grepl('[A-Z]',d[,2]),]
library(reshape2)
p <- dcast(data.frame(d,value=1),'SampleID~Gene.Sym')
rownames(p) <- p[,1]
p <- p[,-1]
p <- as.matrix(p)
p[is.na(p)] <- 0
saveRDS(p,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/mutation.rds')
