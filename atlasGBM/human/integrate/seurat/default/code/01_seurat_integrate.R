rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/human/integrate/seurat/default/res/'
library(Seurat)
dlist = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/human/data/dlist.rds')
for (i in names(dlist)) {
  print(i)
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]],project=i))
}
d.anchors <- FindIntegrationAnchors(object.list = dlist)
d.integrated <- IntegrateData(anchorset = d.anchors)
saveRDS(d.integrated, paste0(rdir, 'human_atlas_GBM_Seurat.rds'))
