library(Seurat)
library(here)
setwd(here())

## prepare reference: select allen brain 10x cells 
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/countfunc.R')
atlas <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/count.rds'), platform='10x', log2normalize = FALSE, species='human')

colnames(atlas) <- gsub('allen:', 'allen_human_10x:', colnames(atlas))

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/meta_allcell.rds')
meta[,1] <- sub('allen_human_10x;', '', meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[colnames(atlas), ]

selcell <- lapply(unique(meta[,7]), function(ct){
  tmp <- rownames(meta[which(meta[,7] == ct), ])
  if (length(tmp) > 500){
    tmp[1:500]
  } else {
    tmp
  }
})
selcell <- unlist(selcell)
atlas <- atlas[, selcell] ##
colnames(atlas) <- paste0('atlas_',meta[selcell,7],':',1:ncol(atlas))

## prepare query: select GBM cells
brain <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
clu <- paste0('cluster ', as.character(Idents(brain)))
names(clu) <- colnames(brain@assays$RNA@counts)
patient = brain@meta.data[,'Sample.ID']
rm(brain)

expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_count_mat.rds')
int <- intersect(rownames(expr), rownames(atlas))
str(int)
expr <- expr[int, ]
atlas <- atlas[int, ]

expr <- cbind(expr,atlas)
rm(atlas)
expr <- CreateSeuratObject(counts = expr, project = "expr")
expr <- NormalizeData(expr)
expr <- FindVariableFeatures(expr)
expr <- ScaleData(expr)
expr <- RunPCA(expr)
expr <- FindNeighbors(expr, dims = 1:10)
expr <- FindClusters(expr, resolution = 0.5)
expr <- RunUMAP(expr, dims = 1:10)

expr@meta.data$orig.ident <- sub('GBM.*','GBM',expr@meta.data$orig.ident)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cluster_allenreference/res/nointer.pdf')
DimPlot(expr,group.by='orig.ident')
dev.off()


saveRDS(expr,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cluster_allenreference/res/nointer.rds')
