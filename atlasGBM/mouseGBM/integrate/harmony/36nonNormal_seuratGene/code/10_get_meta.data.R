res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/atlasGBM_harmony_with_ct.rds')
meta <- res@meta.data
saveRDS(meta, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/meta.data.rds')

