library(Matrix)
suppressMessages(library(Seurat))
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellorigin/compareToHumanCortexAtlas/plot/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellorigin/compareToHumanCortexAtlas/res/'

### read in atlas
atlas <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/pseudobulk/data/humancortexatlas_pb.rds')
# colnames(atlas)[colnames(atlas) == 'Intermediate progenitor cell (IPC)'] = 'IPC'
# colnames(atlas)[grepl('Unclear', colnames(atlas))] <- 'Other non-neuronal'
atlas <- atlas[rowSums(atlas >= 1) >= 1,]

### read in GBM samples
gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
table(as.numeric(gbm@active.ident))
normcell = names(gbm@active.ident)[as.numeric(gbm@active.ident) %in% c(26, 21, 18, 16, 13)]
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
meta <- meta[,-1]
# meta <- meta[meta[,'Pathology'] == 'GBM' & meta[,'Treatment'] == 'Untreated' & meta[,'Tumor.Grade'] == 'IV', -1]
meta <- unique(meta)
meta <- meta[,apply(meta,2,function(i) length(unique(i))) > 1]
rownames(meta) <- meta[,1]
meta <- meta[,-1]
# meta <- gbm@meta.data
# str(meta)
# 
# meta <- meta[!duplicated(meta), ]
gbm <- gbm$RNA@data
# normcell <- unlist(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/res/nonmalignant_cells_high_confident.rds'))

# normcell <- normcell[!is.na(normcell)]
gbm <- gbm[,!colnames(gbm) %in% normcell]
gbm <- gbm[rowMeans(gbm > 0) > 0.01,]

### intersect genes between atlas and GBM
int <- intersect(rownames(atlas),rownames(gbm))
gbm <- gbm[int,]
atlas <- atlas[int, ]

### calculate correlation 
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

set.seed(12345)
clunum <- 100
gc <- kmeans(scalematrix(atlas),clunum)$cluster
atlasgc <- t(sapply(1:clunum,function(i) colMeans(atlas[gc==i,])))
gbmgc <- t(sapply(1:clunum,function(i) colMeans(gbm[gc==i,])))
saveRDS(gbmgc, paste0(rdir, 'gbmgc.rds'))

gbmgc <- scalematrix(gbmgc)
atlasgc <- scalematrix(atlasgc)
cm <- scalematrix(t(gbmgc)) %*% t(scalematrix(t(atlasgc)))/(nrow(gbmgc)-1)            
saveRDS(cm, paste0(rdir, 'cm.rds'))

samp <- sub('_.*','',rownames(cm))
library(pheatmap)
meta$selected <- as.numeric(meta[,'Pathology'] == 'GBM' & meta[,'Treatment'] == 'Untreated' & meta[,'Tumor.Grade'] == 'IV')
for (cutoff in c(0.2,0.3,0.5,0.6)) {
  ccm <- sapply(unique(samp),function(i) {
    colMeans((cm > cutoff)[samp==i,])
  })
  pdf(paste0(pdir, 'full_',cutoff,'.pdf'), width = 12, height = 17)
  pheatmap(ccm,annotation_col = meta)
  dev.off()
}

meta <- meta[meta[,'Pathology'] == 'GBM' & meta[,'Treatment'] == 'Untreated' & meta[,'Tumor.Grade'] == 'IV',]
for (cutoff in c(0.2,0.3,0.5,0.6)) {
  print(cutoff)
  ccm <- sapply(unique(samp),function(i) {
    colMeans((cm > cutoff)[samp==i,])
  })[,rownames(meta)]
  pdf(paste0(pdir, 'selected_',cutoff,'.pdf'), width = 12, height = 17)
  pheatmap(ccm,annotation_col = meta)
  dev.off()
}

## exclude GBM cluster 17 (microglia), 13, 18 (endothelial)
