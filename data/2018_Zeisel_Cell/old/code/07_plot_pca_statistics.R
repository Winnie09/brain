setwd("/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Zeisel_Cell/")
pr = readRDS('./result/pca.rds')
v = sort(pr$sdev,decreasing = T)[1:100]
pdf('./plot/pca_statistics.pdf',width=3,height=3)
plot(v~seq(1,100),pch=20,main='',xlab='num.PCs',ylab='sdev')
dev.off()

