library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/humanatlas/integrate/harmony/seuratGene/res/'  ##### 
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
mat =  readRDS('./atlasGBM/humanGBM/data/combine_mat.rds')
meta =  readRDS('./atlasGBM/humanGBM/data/meta_allcell.rds')
rownames(meta) <- meta$cell
selectcell <- meta$cell[!grepl('GBM', meta$study) & !grepl('2014_Pollen_NatBiotech', meta$study)]
mat <- mat[,selectcell]
mat <- mat[, colSums(mat > 0) > 500]
mat <- mat[rowMeans(mat > 0) >= 0.01,]
mat <- mat[!grepl('^MT-',rownames(mat)),]
meta <- meta[colnames(mat), ]
nfeatures <- 2000
brain <- CreateSeuratObject(counts = mat, project = "atlasGBM") ## input is log-normalized => no more Normalized().
brain@meta.data <- cbind(brain@meta.data, meta)
brain@meta.data$mydonor <- paste0(meta$study, ';',meta$gender,';', meta$age,';', meta$location, ';',meta$time, ';',meta$donor)
brain.ls <- SplitObject(brain,split.by = 'mydonor')
print(length(brain.ls))
for(i in 1:length(brain.ls)){
    brain.ls[[i]] <- FindVariableFeatures(brain.ls[[i]],selection.method = "vst", nfeatures = nfeatures)
    print(i)
}

integrated.features <- SelectIntegrationFeatures(brain.ls, nfeatures = nfeatures)
print(sort(integrated.features))
print(length(integrated.features))
#saveRDS(integrated.features,paste0(path,'newhvg-hvg_',npcs,'pcs-',nfeatures,'features.rds'))

brain <- ScaleData(brain,features = integrated.features,verbose = FALSE)
brain <- RunPCA(brain,npcs = 30,features = integrated.features, verbose = FALSE)
brain <- RunHarmony(brain, c("study", "mydonor"))
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony_tmp.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony.rds'))
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap_embeddings.rds'))



