setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/')
m = readRDS('./processed/count/Sestan.fetalHuman.Psychencode.rds')
m = m[,colSums(m>0)>=1000]
v = colSums(m[grepl('^MT-',rownames(m)),])/colSums(m)
m = m[,v<0.1]
m = m[!grepl('^MT-',rownames(m)),]
dir.create("./processed/qc/",showWarnings = F,recursive = T)
saveRDS(m,'./processed/qc/Sestan.fetalHuman.Psychencode.rds')


m = readRDS('./processed/count/Sestan.adultHumanNuclei.Psychencode.rds')
m = m[,colSums(m>0)>=1000]
v = colSums(m[grepl('^MT-',rownames(m)),])/colSums(m)
m = m[,v<0.1]
m = m[!grepl('^MT-',rownames(m)),]
dir.create("./processed/qc/",showWarnings = F,recursive = T)
saveRDS(m,'./processed/qc/Sestan.adultHumanNuclei.Psychencode.rds')

