#key = 'referenceVec_allenDEG'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/')
library(ggplot2)
res <-  sapply(unique(meta$study), function(i){
  tmp = u[which(meta$study==i),]
  sum(tmp[,1]>5)/nrow(tmp)   
})

pd <- data.frame(method=names(res),prop=res,type=ifelse(grepl('^GBM',names(res)),'GBM','Atlas'),stringsAsFactors = F)
pd$method <- factor(pd$method,levels=names(sort(res)))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/referenceVec_atlasDEG/umap/umap0.5.pdf',height=6,width=6)
ggplot(pd,aes(x=method,y=prop,fill=type)) + geom_bar(stat='identity') + theme_classic() + coord_flip()+
  ylab('prop.cells.umap1.>0.5')
dev.off()
