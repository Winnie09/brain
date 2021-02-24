library(here)
setwd(here())
f <- as.character(commandArgs(trailingOnly = T)[[1]][1])
print(f)

ddir <- paste0('atlasGBM/mouseGBM/pseudotime/data/testvar/', f)
rdir <- paste0('atlasGBM/mouseGBM/pseudotime/res/testvar/', f)
dir.create(rdir, showWarnings = F, recursive = T)

expr <- readRDS(paste0(ddir, '/expr.rds'))
pseudotime <-  readRDS(paste0(ddir, '/pseudotime.rds'))
cellanno <- readRDS(paste0(ddir, '/cellanno.rds'))
design <- readRDS(paste0(ddir, '/design.rds'))

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design = design, type = 'Variable')
saveRDS(res, paste0(rdir, '/res.rds'))



