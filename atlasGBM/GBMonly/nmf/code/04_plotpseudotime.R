factorscore <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')

library(umap)
set.seed(12345)
u <- umap(factorscore)$layout

pd <- do.call(rbind,lapply(1:ncol(factorscore),function(i) {
  data.frame(UMAP1=u[,1],UMAP2=u[,2],score=factorscore[,i],id=i)
}))
pd$score[pd$score < quantile(pd$score,0.05)] <- quantile(pd$score,0.05)
pd$score[pd$score > quantile(pd$score,0.95)] <- quantile(pd$score,0.95)
library(ggplot2)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/pseudotime/umap.pdf'),width=20,height=5)
ggplot(pd,aes(x=UMAP1,y=UMAP2,col=score)) + geom_point(size=0.1) + theme_classic() + facet_wrap(~id,nrow=1) + scale_color_gradientn(colours = rainbow(15)[c(11:1,15)])
dev.off()

cc <- factorscore[,2] > 0

library(ggplot2)
plotfunc <- function(p1,p2) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/pseudotime/',p1,'_',p2,'.pdf'))
  print(ggplot(data.frame(p1=factorscore[,p1],p2=factorscore[,p2],cc=cc),aes(x=p1,y=p2,col=cc)) + geom_point() + theme_classic())
  dev.off()  
}
plotfunc(1,2)
plotfunc(1,3)
plotfunc(1,4)
plotfunc(2,3)
plotfunc(2,4)
plotfunc(3,4)

