library(Seurat)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
r <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/humanAtlas_harmony.rds')

anchors <- FindTransferAnchors(reference = r, query = d,dims = 1:30, reference.reduction = "pca",k.filter=NA)

d <- TransferData(anchorset = anchors, reference = r, query = d,refdata = list(celltype = "seurat_clusters"))
d <- IntegrateEmbeddings(anchorset = anchors, reference = r,query = d, new.reduction.name = "ref.pca")
d <- ProjectUMAP(query = d, query.reduction = "ref.pca", reference = r, reference.reduction = "pca", reduction.model = "umap")

saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project.rds')
