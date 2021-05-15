library(Seurat)
set.seed(12345)
reference <- readRDS('/home-4/zji4@jhu.edu/scratch/GBM/atlas/humanref/data/alltogether_cortex_v3.rds')
meta <- read.table('/home-4/zji4@jhu.edu/scratch/GBM/atlas/humanref/data/allcortex_metadata.tsv',sep='\t',header=T,row.names=1)
for (i in colnames(meta))
  reference@meta.data[,i] <- meta[rownames(reference@meta.data),i]
#Idents(reference) <- 'orig.ident'
reference <- RunPCA(reference)
reference <- RunUMAP(reference,dims=1:50,return.model=TRUE)
saveRDS(reference,file='/home-4/zji4@jhu.edu/scratch/GBM/atlas/humanref/res/reference.rds')

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
  refdata = list(celltype = "Cell.Type"),
  reference.reduction='pca',
  reduction.model = "umap"
)
saveRDS(target,file='/home-4/zji4@jhu.edu/scratch/GBM/atlas/humanref/res/integrate.rds')

