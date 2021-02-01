setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/old/20200212')
library(Seurat)
tmp = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/meta.rds') ## 240295
tmp = tmp[,!colnames(tmp)%in%c('nFeature_RNA','nCount_RNA','orig.ident')]
meta_rosen = readRDS('./brain/data/2018_Rosenberg_Science/proc/meta.rds')  ## 156049
meta_rosen$time <- 'developing'
meta_zeisel <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/meta.rds') ## 160796
colnames(meta_zeisel)[2] <- 'location'
meta_zeisel$time <- 'adult'
meta_zeisel = meta_zeisel[,colnames(meta_rosen)]

addmeta = rbind(meta_zeisel, meta_rosen)
addmeta$study = sub(':.*','',addmeta$cell)
addmeta$species = 'mouse'
addmeta$gender = 'unknown'  ## 316845

addmeta$celltype[addmeta$celltype == 'oligodendrocytes precursor cells'] <- 'oligodendrocyte precursor cells'
addmeta$celltype_super <- sapply(addmeta$celltype,function(i){
  if (grepl('astrocytes', i) | grepl('Astrocytes',i)){
    'astrocytes'
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


addmeta$location_super <- sapply(addmeta$location,function(i){
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

addmeta$time_super <- sapply(addmeta$time,function(i){
  v  = c('post','adult','hESC','iPSC','unknown','developing')
  for (sv in v){
    if (grepl(sv,i))  return(sv)
  }
})
tmp$species <- 'human'
fullmeta = rbind(tmp, addmeta[,colnames(tmp)])
ct = fullmeta$celltype
ct[ct=="vascular and leptomeningeal cells"] <- "vascular and leptomeningeal cell"
ct[ct=="oligodendrocytes"] <- "oligodendrocyte"
ct[ct=="oligodendrocyte precursor cells"] <- "oligodendrocyte precursor cell"
ct[ct=="oligodendroctye precursor cells"] <- "oligodendrocyte precursor cell"
ct[ct=="Neurons"] <- "neuron"
ct[ct=="ependyma"] <- "ependymal"
ct[ct=="Ependymal"]  <- "ependymal" 
ct[ct=="endothelia cells"] <- "endothelial cell"
ct[ct=="Astrocytes"] <- "astrocyte"
ct[ct=="astrocytes"] <- "astrocyte"
ct[ct=="Vascular"] <- "vascular"
id = setdiff(grep('cell',ct),grep('stem cell', ct))
ct[id] <- sub(' cell','', ct[id])
ct[id] <- sub(' cells','', ct[id])
fullmeta$celltype = ct

fullmeta$celltype_super <- sapply(fullmeta$celltype,function(i){
  if (grepl('embryonic stem cell', i)){
    'ESC'
  } else if (grepl('embryonic ', i)){
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
# saveRDS(fullmeta,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/old/20200212/result/full/meta.rds')
saveRDS(fullmeta,'./brain/atlas/res/full/meta.rds')
