library(Seurat)
ref.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')
reference.vec = c('2018_Li_Science:fetal','2018_Li_Science:adult','2018_Lake_NatBiotech:VisualCortex','2018_Lake_NatBiotech:FrontalCortex','2018_Lake_NatBiotech:CerebellarHem','2018_Fan_CellResearch','2016_La:iPSMoleculeCounts','2016_La:ESMoleculeCounts','2016_La:EmbryoMoleculeCounts')
### make meta
meta <- sapply(c(reference.vec,'allen'),function(i){
  print(i)
  a = sub(':.*','',i)
  b = sub('.*:','',i)
  if (a != b){
    meta = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/',a,'/proc/meta/',b,'.rds'))
  } else if(a=='allen'){
    meta = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/',a,'/data/proc/meta.rds'))    
  } else {
    meta = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/',a,'/proc/meta.rds'))
  }
  tmpmet <- sapply(c('cell','celltype','location','species','time','gender'),function(type){
    if (type %in% colnames(meta)){
      tmp = as.character(meta[,type] )
    } else {
      tmp = rep('unknown', nrow(meta))
    }
    
  })
  tmpmet
})
meta = data.frame(do.call(rbind,meta),stringsAsFactors = F)
meta$celltype[meta$celltype == 'oligodendrocytes precursor cells'] <- 'oligodendrocyte precursor cells'

meta$celltype_super <- sapply(meta$celltype,function(i){
  if (grepl('embryonic stem cell', i)){
    'ESC'
  } else if (grepl('embryonic', i)){
    'embryonic'
  } else if (grepl('iPSC',i)){
    'iPSC'
  } else if (grepl('neu',i) | grepl('neural',i)){
    'neuronal'
  } else if (grepl('progenitor',i)){
    'progenitor'
  } else if (grepl('unknown',i)){
    'unknown'
  } else {
    i
  }
})


meta$location_super <- sapply(meta$location,function(i){
  if (grepl('Frontal',i)){
    'frontal'
  } else if (grepl('Inferior',i)){
    'inferior'
  } else if (grepl('Parietal', i)){
    'parietal'
  } else if (grepl('Temporal',i)){
    'temporal'
  } else if (grepl('unkown',i)){
    'unknown'
  } else {
    i
  }
})

meta$time_super <- unlist(sapply(meta$time,function(i){
  if (i == 'hESC differentiated day_0'){
    return ('hESC')
  } else if (grepl('post',i)){
    return('post')
  } else if (grepl('embryo',i) | grepl('hESC',i)){
    return('embryo')
  } else if (grepl('adult',i)){
    return('adult')
  } else if (grepl('iPSC',i)){
    return('iPSC')
  } else if (!is.null(i)){
    return(i)
  } else {
    return('unknown')
  }
}))
meta$study = sub(':.*','',meta$cell)
rownames(meta) = meta[,'cell']
ref.integrated@meta.data = cbind(ref.integrated@meta.data[,1:3],meta[rownames(ref.integrated@meta.data),])
## pca, umap
ref.integrated <- ScaleData(ref.integrated, verbose = FALSE)
ref.integrated <- RunPCA(ref.integrated, npcs = 30, verbose = FALSE)
ref.integrated <- RunUMAP(ref.integrated, reduction = "pca", dims = 1:30)
saveRDS(ref.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')
