clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')

normalclu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/normalclu.rds')

normalcell = clu[clu %in% normalclu]
cancercell = clu[!clu %in% normalclu]
str(normalcell)
str(cancercell)

saveRDS(normalcell, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/normalcell_final.rds')
saveRDS(cancercell, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/cancercell_final.rds')
