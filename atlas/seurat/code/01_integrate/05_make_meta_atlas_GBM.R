setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
fullmeta = readRDS('./brain/atlas/result/atlasOnly_homolog_sub/meta.rds')
u = readRDS(paste0('./brain/atlas/result/atlas_GBM_homolog_sub/umap/umap.rds'))
cell = rownames(u)[grepl('GBM',rownames(u))]
study = gsub('_.*','', cell)
addmeta = data.frame(cell = cell,
                     study = study,
                     celltype_super = 'unknown',
                     species = 'human',
                     time_super = 'adult',stringsAsFactors = F)
meta = rbind(fullmeta, addmeta)
meta = meta[match(rownames(u), as.character(meta$cell)), ]
saveRDS(meta, './brain/atlas/result/atlas_GBM_homolog_sub/meta.rds')
