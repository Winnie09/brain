library(NMF)
library(parallel)
num.nmf <- 4
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmfoptimr/res/nmfres/')
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
allp <- unique(meta[, 'Sample.ID'])
gene <- lapply(allp, function(p){
  mattmp <- dlist[[p]]          
  rownames(mattmp[rowMeans(mattmp>0) > 0.01, ])
})
gene2 <- table(unlist(gene))
gene2 <- names(gene2)[gene2==length(allp)]
gene2 <- gene2[!grepl('^MT-',gene2)]
gene2 <- gene2[!grepl('^RPL|^RPS',gene2)]

# ==================================
# center expr for each gene, and nmf
# ==================================
i <- as.numeric(commandArgs(trailingOnly = T)[[1]])
print(i)
p <- unique(meta[, 'Sample.ID'])[i]
print(p)
af <- list.files(rdir)
af <- sub('.rds', '', sub('nmf_', '', af))
if (p %in% af) break 
mat <- dlist[[p]][gene2,]
mat <- mat[, setdiff(colnames(mat), nonmalignant)]
mat <- mat[rowSums(mat) > 0, ]
mat <- mat - rowMeans(mat) ## centered expression
mat[mat < 0]  <-  0  ## convert negatives to zero
# res <- nmf(mat, num.nmf, seed = 12345) ## manually set num.nmf = 4
if(requireNamespace("Biobase", quietly=TRUE)){
  # perform 10 runs for each value of r in range 2:6
  if(requireNamespace("Biobase", quietly=TRUE)){
    pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmfoptimr/plot/estim.r_', p, '.pdf'))
    plot(estim.r)
    dev.off()
  }
  estim.r <- nmf(mat, 2:6, nrun=30, seed=123456)
}
saveRDS(estim.r, paste0(rdir, 'estim.r_', p, '.rds'))

