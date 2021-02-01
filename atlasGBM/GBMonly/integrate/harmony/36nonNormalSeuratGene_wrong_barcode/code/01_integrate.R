library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
mat =  readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/proc/tumor/all/matrix/count.rds')
str(mat)
allp <- sub('_.*', '', colnames(mat))
table(allp)
id <- which(!allp %in% c('GBM076', 'GBM077', 'GBM086', 'GBM089', 'GBM090', 'GBM094', 'GBM047'))
mat <- mat[,id]
str(mat)
mat <- mat[, colSums(mat > 0) > 500]
mat <- log2(t(t(mat)/(colSums(mat) / 1e6) + 1))
str(mat)
summary(colSums(mat))
mat <- mat[!grepl('^MT-',rownames(mat)),]
str(mat)
mat <- mat[rowMeans(mat > 0) >= 0.01,]
str(mat)
nfeatures <- 2000

brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- data.frame(brain@meta.data, patient = allp[id])
brain.ls <- SplitObject(brain,split.by = 'patient')
print(length(brain.ls))

for(i in 1:length(brain.ls)){
    brain.ls[[i]] <- NormalizeData(brain.ls[[i]], verbose = FALSE)
    brain.ls[[i]] <- FindVariableFeatures(brain.ls[[i]],selection.method = "vst", nfeatures = nfeatures)
    print(i)
}

integrated.features <- SelectIntegrationFeatures(brain.ls, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
saveRDS(integrated.features,paste0(rdir,nfeatures,'features.rds'))

brain <- ScaleData(brain,features = integrated.features,verbose = FALSE)
brain <- RunPCA(brain,npcs = 30,features = integrated.features, verbose = FALSE)
brain <- RunHarmony(brain, c("patient"))
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony_tmp.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))

print(sessionInfo())


