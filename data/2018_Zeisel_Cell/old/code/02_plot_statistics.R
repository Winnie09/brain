m = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/processed/count/mat.rds')
meta = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/processed/meta/cell.rds')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/plot/statistics_before_filter.pdf',width=6,height = 6)
par(mfrow=c(2,2))
v = colSums(m>0)
hist(v,col='grey',breaks=50,main='',xlab='num.genes.detected',ylab='num.cells')

v = 100*rowMeans(m>0)
hist(v,col='grey',breaks=50,main='',xlab='%cells.expressed',ylab='num.genes')

v = colSums(m)
hist(v,col='grey',breaks=50,main='',xlab='sum.UMI.count',ylab='num.cells')

v = 100*colSums(m[grepl('^mt-',rownames(m)),])/ colSums(m)
hist(v,col='grey',breaks=50,main='',xlab='%mitochondrial.count',ylab='num.cells')

dev.off()
