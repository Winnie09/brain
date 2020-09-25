setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')
library(Seurat)
tmp = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/meta.rds') ## 240295
tmp = tmp[,!colnames(tmp)%in%c('nFeature_RNA','nCount_RNA','orig.ident')]

m <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/human/data/combine_mat.rds')
ap <- sub(';.*', '', colnames(m))
length(table(ap))
acell <- sub('.*;', '', colnames(m))

m1 <- data.frame(cell = colnames(m)[!grepl('GBM', ap)], 
                 study = ap[!grepl('GBM', ap)],
                 stringsAsFactors = FALSE)
m2 <- data.frame(cell = colnames(m)[grepl('GBM', ap)], 
                 study = ap[grepl('GBM', ap)],
                 stringsAsFactors = FALSE)

metalist <- list()
for (s in unique(m1$study)){
  print(s)
  f = 'meta.rds'
  m1.tmp <- m1[m1$study == s, ]
  if (f %in% list.files(paste0( s, '/proc/'))){
    tmp <- readRDS(paste0( s, '/proc/',f))
  } else {
    if (f %in% list.files(paste0( s, '/proc/meta'))){
      tmp <- readRDS(paste0( s, '/proc/meta/',f))
    } else {
      af <- list.files(paste0( s, '/proc/meta'))  
      if (length(af) >= 1){
        af <- af[!af%in% c('meta.rds')]
        tmpp <- lapply(af, function(f){
          a <- readRDS(paste0( s, '/proc/meta/',f))
          a <- data.frame(apply(a, 2, as.character), stringsAsFactors = FALSE)
        })
        tmp <- Reduce(function(...) merge(..., all=T), tmpp)
      }
    } 
  }
  
  colnames(tmp) <- tolower(colnames(tmp))
  tmp[,'cell'] <- as.character(tmp[,'cell'])
  tmp <- data.frame(apply(tmp, 2, as.character), stringsAsFactors = FALSE)
  colnames(tmp)[which(colnames(tmp) == 'sex')] <- 'gender'
  metatmp <- cbind(m1.tmp, tmp[match(sub('.*;', '',m1.tmp$cell), tmp[,'cell']), ])
  colnames(metatmp)[3] <- 'cell.shortname'
  metalist[[s]] <- metatmp
}


meta1 <- reshape::merge_recurse(metalist)

meta.gbm <- read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/doc/metas_20200329.csv', header = T, sep = ',', as.is = TRUE)

meta2 <- cbind(m2,meta.gbm[match(m2$study, meta.gbm$Sample.ID), seq(1, ncol(meta.gbm))])
colnames(meta2)[1:9] <- tolower(colnames(meta2)[1:9])
colnames(meta2)[which(colnames(meta2) == 'sex')] <- 'gender'
meta <- reshape::merge_recurse(list(meta1, meta2))
  
# ------------
meta$celltype[meta$celltype == 'oligodendrocytes precursor cells'] <- 'oligodendrocyte precursor cells'
meta$celltype_super <- sapply(meta$celltype,function(i){
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

meta$time_super <- sapply(meta$time,function(i){
  v  = c('post','adult','hESC','iPSC','unknown','developing')
  flag = 0
  for (sv in v){
    if (grepl(sv,i))  {
      flag = 1
      return(sv)
    }
  }
  if (!flag) return(i)
})

ct = meta$celltype
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
meta$celltype = ct

meta$celltype_super <- sapply(meta$celltype,function(i){
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

saveRDS(meta, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/human/data/meta_allcell.rds')

