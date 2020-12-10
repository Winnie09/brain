library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')     
rdir = './atlasGBM/humanatlas/integrate/harmony/default/res/'
dir.create(rdir, recursive = T)
mat =  readRDS('./atlasGBM/humanGBM/data/combine_mat.rds')
meta =  readRDS('./atlasGBM/humanGBM/data/meta_allcell.rds')
rownames(meta) <- meta$cell
selectcell <- meta$cell[!grepl('GBM', meta$study) & !grepl('2014_Pollen_NatBiotech', meta$study)]
mat <- mat[,selectcell]
mat <- mat[, colSums(mat > 0) > 500]
mat <- mat[rowMeans(mat > 0) >= 0.01,]
mat <- mat[!grepl('^MT-',rownames(mat)),]
meta <- meta[colnames(mat), ]

study = sub(';.*', '', colnames(mat))
brain <- CreateSeuratObject(counts = mat, project = "brain", min.cells = 5) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = brain@var.genes, npcs = 20, verbose = FALSE)

brain@meta.data$donor <- paste0(meta$study, ';',meta$gender,';', meta$age,';', meta$location, ';',meta$time, ';',meta$donor)
brain@meta.data$study <- study
brain <- RunHarmony(brain, c("study", "donor"))
saveRDS(brain,paste0(rdir, 'humanAtlas_harmony_object.rds'))

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

