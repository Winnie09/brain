library(NMF)
library(parallel)
num.nmf <- 4
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmfnew/res/', num.nmf, '/')
dir.create(rdir, showWarnings = F, recursive = T)
dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
meta <- meta[meta[, 'Pathology'] == 'GBM' & meta[, 'Tumor.Grade'] == 'IV' & meta[, 'Treatment'] == 'Untreated', ]
# i = unique(meta[, 'Location'])[1]
# allp <- unique(meta[meta[,'Location']==i, 'Sample.ID'])

# ==========================
# only use maglinant cells !!
# =========================
nonmalignant <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/res/nonmalignant_cells_high_confident.rds')
nonmalignant <- unlist(nonmalignant)
nonmalignant <- nonmalignant[!is.na(nonmalignant)]

# ==============
# gene filtering 
# ==============
i <- as.numeric(commandArgs(trailingOnly = T)[[1]])
print(i)
p <- unique(meta[, 'Sample.ID'])[i]
print(p)
keep <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/filtergenelist/keep.rds')

# ==================================
# center expr for each gene, and nmf
# ==================================
mat <- dlist[[p]][keep,]
mat <- mat[, setdiff(colnames(mat), nonmalignant)]
mat <- mat[rowSums(mat) > 0, ]
mat <- mat - rowMeans(mat) ## centered expression
mat[mat < 0]  <-  0  ## convert negatives to zero
res <- nmf(mat, num.nmf, seed = 12345) ## nmf 
saveRDS(res, paste0(rdir, 'nmf_', p, '.rds'))


