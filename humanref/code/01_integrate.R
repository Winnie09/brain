library(Seurat)
reference <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/data/alltogether_cortex_v3.rds')
#Idents(reference) <- 'orig.ident'
meta = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/data/meta.tsv', header = T, sep = '\t')
rownames(meta) <- meta[,1]
reference@meta.data <- cbind(reference@meta.data, meta[rownames(reference@meta.data), 4:6])
reference <- RunPCA(reference)
reference <- RunUMAP(reference,dims=1:50,return.model=TRUE)
saveRDS(reference, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/reference.rds')
# $ : chr [1:19368] "FO538757.2" "AP006222.2" "RP4-669L17.10" "RP5-857K21.4" ...
# $ : chr [1:58145]

target <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
#Idents(target) <- 'orig.ident'

anchors <- FindTransferAnchors(
  reference = reference,
  query = target,
  k.anchor=50,
  reference.reduction='pca',
  dims = 1:50
)

target <- MapQuery(
  anchorset = anchors,
  query = target,
  reference = reference,
  reference.reduction='pca',
  reduction.model = "umap"
)
saveRDS(target,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/integrate.rds')

