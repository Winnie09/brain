rm(list=ls())
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/res/twogroup/'
num = sapply(list.files(ddir, pattern = '.rds'), function(f){
  print(f)
  Res <- readRDS(paste0(ddir, f))
  print(names(Res))
  stat = Res$statistics
  sum(stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05)
})
names(num) = sub('.rds', '', names(num))
num = num[order(num)]
saveRDS(num, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/plot/twogroup/perf/num_diffgene.rds')

library(ggplot2)
pd = data.frame(numdiffgene = num, comparison = names(num), stringsAsFactors = F)
pd[,2] = factor(pd[,2], levels = names(num))

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/plot/twogroup/perf/num_diffgene.pdf', width = 3.2, height = 3.5)
ggplot(data=pd, aes(x = comparison, y = numdiffgene)) +
  geom_bar(stat = 'identity')+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('comparison') + ylab('number of XDE genes')
dev.off()

