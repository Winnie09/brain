setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')
d = list()
af1 = list.files('./2018_Lake_NatBiotech/processed/scran/')
for (f in af1){
  d[[f]] = readRDS(paste0('./2018_Lake_NatBiotech/processed/scran/',f))  
}
af2 = list.files('./2018_Li_Science/processed/scran/')
af2 = "Sestan.adultHumanNuclei.Psychencode.rds"
for (f in af2){
  d[[f]] = readRDS(paste0('./2018_Li_Science/processed/scran/',f)) 
}

ag = sapply(d, function(i)rownames(i))

g = intersect(ag[[1]],ag[[2]])
for (i in 3:length(ag)){
  g = intersect(g,ag[[i]])
}

df = sapply(d,function(i){
  i[g,]
})
df = do.call(cbind,df)
saveRDS(df,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/genebycell.rds')

ct = sapply(d,ncol)
af2 = sub('Sestan.','',af2)
af2 = sub('.Psychencode','',af2)
names(ct) = sub('.rds','',c(paste0('2018_Lake_NatBiotech;',af1),paste0('2018_Li_Science;',af2)))
saveRDS(ct,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/batch.rds')
