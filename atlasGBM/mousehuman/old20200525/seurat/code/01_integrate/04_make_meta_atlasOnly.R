setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
fullmeta = readRDS('./brain/atlas/old/20200212/result/full/meta.rds')
u = readRDS('./brain/atlas/result/atlasOnly_homolog_sub/umap/umap.rds')
s = gsub(':.*','', rownames(u))
meta = data.frame(cell = rownames(u), study = sub(':.*','',rownames(u)))
meta = cbind(meta, 
             celltype_super = fullmeta[match(meta$cell, fullmeta$cell),'celltype_super'],
             species = fullmeta[match(meta$cell, fullmeta$cell),'species'],
             time_super = fullmeta[match(meta$cell, fullmeta$cell),'time_super'])

meta[grepl('2017_Nowakowski_Science',meta$study),'species'] = 'human'
meta[grepl('2016_Tasic_NatNeuro',meta$study),'species'] = 'mouse'
meta[grepl('2015_Zeisel_Science',meta$study),'species'] = 'mouse'
meta[grepl('2017_Nowakowski_Science',meta$study),'time_super'] = 'developing'
meta[grepl('2016_Tasic_NatNeuro',meta$study),'time_super'] = 'adult'
meta[grepl('2015_Zeisel_Science',meta$study),'time_super'] = 'adult'
meta[is.na(meta$species),'species'] <- 'mouse'
meta[is.na(meta$celltype_super),'celltype_super'] <- 'unknown'
saveRDS(meta, './brain/atlas/result/atlasOnly_homolog_sub/meta.rds')
