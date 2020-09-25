setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result')
fullmeta = readRDS('./atlasOnly_homolog_sub/meta.rds')
u = readRDS(paste0('./atlas_GBM_homolog_sub/umap/umap.rds'))
cell = rownames(u)[grepl('GBM',rownames(u))]
study = gsub('_.*','', cell)
addmeta = data.frame(cell = cell,
                     study = study,
                     celltype_super = 'unknown',
                     species = 'human',
                     time_super = 'adult',stringsAsFactors = F)
meta = rbind(fullmeta, addmeta)
meta = meta[match(rownames(u), as.character(meta$cell)), ]
meta$study_super = sapply(as.character(meta$study), function(i) ifelse(grepl('GBM',i),'GBM',i))

for (i in c("study",'study_super','celltype_super','species','time_super')){
  meta[is.na(meta[,i]),i] = 'unknown'
  tab = table(meta[,i])
  meta[,i] = paste0(meta[,i],'(',tab[match(meta[,i],names(tab))],')')    
  if (i == 'study'|i== 'study_super')
  for (j in 1:nrow(meta)){
    if (grepl('2018_Rosenberg',meta[j,i])){
      meta[j,i] = sub(')', paste0('/156049)'),meta[j,i] )
    } else if (grepl('2018_Zeisel',meta[j,i])) {
      meta[j,i] = sub(')', paste0('/492949)'),meta[j,i] )
    } else if (grepl('allen',meta[j,i])){
      meta[j,i] = sub(')', paste0('/47509)'),meta[j,i] )
    }
  }
}
saveRDS(meta, './atlas_GBM_homolog_sub/meta.rds')
