setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/')
plotfunc<- function(m,title){
  v = colSums(m>0)
  hist(v,col='grey',breaks=50,main=title,xlab='num.genes.detected',ylab='num.cells')
  
  v = 100*rowMeans(m>0)
  hist(v,col='grey',breaks=50,main=title,xlab='%cells.expressed',ylab='num.genes')
  
  v = colSums(m)
  hist(v,col='grey',breaks=50,main=title,xlab='sum.UMI.count',ylab='num.cells')
  
  v = 100*colSums(m[grepl('^MT-',rownames(m)),])/ colSums(m)
  hist(v,col='grey',breaks=50,main=title,xlab='%mitochondrial.count',ylab='num.cells')
}

pdf('./plot/statistics_before_filter.pdf',width=12,height = 4)
par(mfrow=c(3,4))

m = readRDS('./processed/count/CerebellarHem.rds')
plotfunc(m,'CerebellarHem')

m = readRDS('./processed/count/FrontalCorte.rds')
plotfunc(m,'FrontalCorte')

m = readRDS('./processed/count/VisualCortex.rds')
plotfunc(m,'VisualCortex')

dev.off()
