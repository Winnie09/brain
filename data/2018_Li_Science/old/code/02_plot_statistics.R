setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/')
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

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/plot/statistics_before_filter.pdf',width=10,height = 5)
par(mfrow=c(2,4))

m = readRDS('./processed/count/Sestan.fetalHuman.Psychencode.rds')
meta = readRDS('./processed/meta/Sestan.fetalHuman.Psychencode.rds')
plotfunc(m,'fetal')

m = readRDS('./processed/count/Sestan.adultHumanNuclei.Psychencode.rds')
meta = readRDS('./processed/meta/Sestan.adultHumanNuclei.Psychencode.rds')
plotfunc(m,'adult')

dev.off()


