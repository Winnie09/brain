list = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/nonmaglinant_cells.rds')
cell = unlist(list)
clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
patient = gsub('_.*', '', names(clu))  

samp <- unique(sub('_.*','',cell))
clusub <- clu[patient %in% samp]
t1 <- table(clusub[names(clusub) %in% cell])
t2 <- table(clusub)
t1/t2[names(t1)]


