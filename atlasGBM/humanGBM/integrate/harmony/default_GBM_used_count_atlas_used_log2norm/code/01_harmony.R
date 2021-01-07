library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/human/integrate/harmony/default/res/'
mat =  readRDS('./atlasGBM/human/data/combine_mat.rds')
study = sub(';.*', '', colnames(mat))
brain <- CreateSeuratObject(counts = mat, project = "brain", min.cells = 5) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = brain@var.genes, npcs = 20, verbose = FALSE)
brain@meta.data$stim <- study
brain@meta.data$study <- study
# saveRDS(brain,paste0(rdir, 'humanAtlas_GBM_harmony.rds'))

brain <- RunHarmony(brain, "study")

pca_embeddings <- Embeddings(brain, 'pca')
saveRDS(pca_embeddings, paste0(rdir, 'pca_embeddings.rds'))

harmony_embeddings <- Embeddings(brain, 'harmony')
saveRDS(harmony_embeddings, paste0(rdir, 'harmony_embeddings.rds'))

brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
umap <- brain@reductions$umap@cell.embeddings
saveRDS(umap, paste0(rdir, 'umap.rds'))

saveRDS(brain,paste0(rdir, 'humanAtlas_GBM_harmony.rds'))
