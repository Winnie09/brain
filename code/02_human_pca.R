data = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/genebycell.rds')
batch <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/batch.rds')
batch = rep(names(batch), batch)
suppressMessages(library(scran))
## highly variable genes
# cn <- colnames(data)
fit <- trendVar(data)
decomp <- decomposeVar(data,fit)
gs <- row.names(decomp)[decomp[,'total'] > decomp[,'tech']]
data <- data[gs,]  ###  [1:13662, 1:299]
scmd <- sapply(1:length(unique(batch)),function(sp) {
  paste0("data[gs,batch==unique(batch)[",sp,"]]")
})

pr = prcomp(t(data),scale. = T)
saveRDS(pr[c('sdev','x')],'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/pca.rds')

dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/plot/merge_human/',recursive = T,showWarnings = F)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/plot/merge_human/pca_statistic.pdf',width=12,height=3)
par(mfrow=c(1,4))
results = decomp
plot(results$mean, results$total,xlab='Mean normalized log-expression per gene',ylab='Variance of the normalized log-expression per gene')
o <- order(results$mean)
lines(results$mean[o], results$tech[o], col="red", lwd=2)
plot(results$mean, results$bio,xlab='Mean normalized log-expression per gene',ylab='Biological component of the variance')
plot(results$total, results$tech,xlab='Variance of the normalized log-expression per gene',ylab='Technical component of the variance')
abline(a=0,b=1,col='red')
plot(pr$sdev[1:100],xlab='number of PC', ylab='sd',main='',pch=20)
abline(v=30,col='red')
dev.off()

