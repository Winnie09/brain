load('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/raw/Sestan.adultHumanNuclei.Psychencode.Rdata')
rownames(umi2) = sub('.*\\|','',rownames(umi2))
saveRDS(umi2,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/processed/count/Sestan.adultHumanNuclei.Psychencode.rds')

saveRDS(meta2,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/processed/meta/Sestan.adultHumanNuclei.Psychencode.rds')

rm(list=ls())
load('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/raw/Sestan.fetalHuman.Psychencode.Rdata')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
count2 = as.matrix(count2)
m <- rmDupGenesNameOnly(count2, sub('.*\\|','',rownames(count2)))
saveRDS(m,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/processed/count/Sestan.fetalHuman.Psychencode.rds')
saveRDS(meta2,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/processed/meta/Sestan.fetalHuman.Psychencode.rds')

