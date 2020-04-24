setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/')
af = list.files('./processed/scran/')
for (f in af){
  cnt = readRDS(paste0('./processed/scran/',f))
  cnt = log2(cnt+1)
  dir.create('./processed/log2/',showWarnings = F, recursive = T)
  saveRDS(cnt,file=paste0('./processed/log2/',f))  
}

