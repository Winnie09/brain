
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/')
af = list.files('./processed/count/')
for (f in af){
  m = readRDS(paste0('./processed/count/',f))
  m = m[,colSums(m>0)>=500]
  v = colSums(m[grepl('^MT-',rownames(m)),])/colSums(m)
  m = m[,v<0.1]
  m = m[!grepl('^MT-',rownames(m)),]
  dir.create("./processed/qc/",showWarnings = F,recursive = T)
  saveRDS(m,paste0('./processed/qc/',f))
  
}


