library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/mouseGBM/integrate/harmony/seuratGene/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
## select nonNormal GBM samples and do gene filtering
mat =  readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/proc/tumor/all/matrix/count.rds')
str(mat)
allp <- sub('_.*', '', colnames(mat))
table(allp)
id <- which(!allp %in% c('GBM076', 'GBM077', 'GBM086', 'GBM089', 'GBM090', 'GBM094', 'GBM047'))
mat <- mat[,id]
str(mat)
allp <- allp[id]
saveRDS(mat, 'atlasGBM/GBMonly/data/36nonNormal_combined_count_mat.rds')
mat <- log2(t(t(mat)/(colSums(mat) / 1e6) + 1))
str(mat)
summary(colSums(mat))
saveRDS(mat, 'atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_mat.rds')
mat <- mat[!grepl('^MT-',rownames(mat)),]
str(mat)
mat <- mat[rowMeans(mat > 0) >= 0.01,]
str(mat)
nfeatures <- 2000
brain.ls <- list()
for (i in unique(allp)){
  brain.ls[[i]] <- mat[,allp==i]
}

mouseatlas <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/data/combine_mat.rds')
mousemeta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/data/meta_allcell.rds')
identical(mousemeta$cell, colnames(mouseatlas))
mousemeta$mydonor <- paste0(mousemeta$study, ';',mousemeta$time)
for (i in unique(mousemeta$mydonor)){
  brain.ls[[i]] <- mouseatlas[, mousemeta$mydonor == i]
}
allp <- c(as.character(allp), as.character(mousemeta$mydonor))

intgn <- sapply(brain.ls,row.names)
intgn <- table(unlist(intgn))
intgn <- names(intgn)[intgn==length(brain.ls)]
for (i in 1:length(brain.ls)) {
  brain.ls[[i]] <- brain.ls[[i]][intgn,]
  colnames(brain.ls[[i]]) = paste0(names(brain.ls)[i], ';', colnames(brain.ls[[i]]))
}
brain <- do.call(cbind, brain.ls)

brain <- CreateSeuratObject(counts = brain, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- data.frame(brain@meta.data, mydonor = allp, normalTumor = ifelse(grepl('GBM', allp), 'GBM', 'mouseatlas'))
brain.ls <- SplitObject(brain,split.by = 'mydonor')

for(i in 1:length(brain.ls)){
    brain.ls[[i]] <- CreateSeuratObject(brain.ls[[i]])
    brain.ls[[i]] <- FindVariableFeatures(brain.ls[[i]],selection.method = "vst", nfeatures = nfeatures)
    print(i)
}

integrated.features <- SelectIntegrationFeatures(brain.ls, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
saveRDS(integrated.features,paste0(rdir,nfeatures,'features.rds'))

brain <- ScaleData(brain,features = integrated.features,verbose = FALSE)
brain <- RunPCA(brain,npcs = 30,features = integrated.features, verbose = FALSE)
brain <- RunHarmony(brain, c("normalTumor", "mydonor"))
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony_tmp.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))

print(sessionInfo())




