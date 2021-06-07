library(ggplot2)
library(RColorBrewer)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/GBM/atlas/data/pb.rds')
m <- readRDS('/home-4/zji4@jhu.edu/scratch/GBM/atlas/data/meta.rds')
m <- m[m$Pathology=='GBM' & m$Treatment=='Untreated',]
d <- sapply(d,function(i) i[,colnames(i) %in% m$Sample.ID,drop=F])
d <- d[sapply(d,ncol)==nrow(m)]
for (i in names(d)) rownames(d[[i]]) <- paste0(i,':',rownames(d[[i]]))
d <- do.call(rbind,d)
d <- d[rowMeans(d > 0.1) > 0.1,]
cv <- apply(d,1,sd)/rowMeans(d)
d <- d[cv > 0.5,]
pr <- prcomp(t(d),scale. = T)$x
saveRDS(pr, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pbpt/res/pca_expr.rds')
cv <- brewer.pal(4,'Set1')
names(cv) <- sort(unique(m$Location))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pbpt/plot/pca_expr.pdf',width=5,height=4)
ggplot(data.frame(PC1=pr[,1],PC2=pr[,2],location=m[match(rownames(pr),m$Sample.ID),'Location']),aes(x=PC1,y=PC2,col=location)) + geom_point() + theme_classic() + scale_color_manual(values=cv)
dev.off()

