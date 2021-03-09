library(here)
setwd(here())
# fd <- 'pt_astrocyte'
# f <- 'Occipital_Frontal'
fd <- as.character(commandArgs(trailingOnly = T)[[1]])
f <- as.character(commandArgs(trailingOnly = T)[[2]])
print(f)
ddir <- paste0('atlasGBM/mouseGBM/pseudotime/data/testvar/', fd, '/', f)
rdir <- paste0('atlasGBM/mouseGBM/pseudotime/res/testvar_saver_chisq/', fd, '/', f)
dir.create(rdir, showWarnings = F, recursive = T)

expr <- readRDS(paste0(ddir, '/expr.rds'))
expr <- expr[rowMeans(expr > 0) > 0.01, ]
pseudotime <-  readRDS(paste0(ddir, '/pseudotime.rds'))
cellanno <- readRDS(paste0(ddir, '/cellanno.rds'))
design <- readRDS(paste0(ddir, '/design.rds'))

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design = design, test.type = 'Variable', ncores = 12, test.method = 'chisq')
saveRDS(res, paste0(rdir, '/res.rds'))


