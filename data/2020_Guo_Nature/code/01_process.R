library(data.table)
d <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/meta/HCL_Fig1_cell_Info.csv',data.table = F)
d[,1] <- paste0('2020_Guo_Nature:',d[,1])
e <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/proc/expr.rds')
d <- d[match(colnames(e),d[,1]),]

d <- d[,c(1,4,6,7)]
colnames(d) <- c('cell','time','donor','celltype')
d$study <- '2020_Guo_Nature'
saveRDS(d, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/meta/meta.rds')

ct <- sub('.*:', '', colnames(e))
ct <- sub('_.*','', ct)

id <- which((ct %in% c('HESC', 'AdultCerebellum', 'FetalBrain')))
e.sub <- e[, id]
d.sub <- d[id, ]
saveRDS(d.sub, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/meta/meta_sub.rds')
saveRDS(e.sub, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/proc/expr_sub.rds')



