setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/')
m = readRDS('./processed/count/mat.rds')
m = m[,colSums(m>0)>=1000]

v = colSums(m[grepl('^mt-',rownames(m)),])/colSums(m)
m = m[,v<0.1]
m = m[!grepl('^mt-',rownames(m)),]
dir.create("./processed/qc/",showWarnings = F,recursive = T)
saveRDS(m,'./processed/qc/mat.rds')

