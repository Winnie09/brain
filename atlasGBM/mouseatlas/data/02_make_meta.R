setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')
mat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/data/combine_mat.rds')
ap <- sub(';.*','',colnames(mat))
m <- data.frame(cell = colnames(mat)[!grepl('GBM', ap)], 
                 study = ap[!grepl('GBM', ap)],
                 stringsAsFactors = FALSE)
metalist <- list()
m1 <- readRDS('./2016_La/proc/meta/MouseAdultDAMoleculeCounts.rds')
m2 <- readRDS('./2016_La/proc/meta/MouseEmbryoMoleculeCounts.rds')
metalist[["2016_La"]] <- rbind(m1,m2)
metalist[["2015_Zeisel_Science"]] <- data.frame(cell = sub('.*;','',m[m$study == '2015_Zeisel_Science','cell']), celltype = rep(NA, nrow(m[m$study == '2015_Zeisel_Science',])), time = NA, location = NA, sex = NA,  stringsAsFactors = FALSE)
metalist[["2016_Tasic_NatNeuro"]] <- data.frame(cell = sub('.*;','', m[m$study == '2016_Tasic_NatNeuro','cell']), celltype = rep(NA, nrow(m[m$study == '2016_Tasic_NatNeuro',])), time = NA, location = NA, sex = NA,  stringsAsFactors = FALSE)
metalist[["2018_Rosenberg"]] <- readRDS('./2018_Rosenberg_Science/proc/meta.rds')
metalist[["2018_Zeisel"]] <- readRDS('./2018_Zeisel_Cell/proc/meta.rds')
m1 <- readRDS('./2019_loo_natcom/proc/meta_E14.rds')
m2 <- readRDS('./2019_loo_natcom/proc/meta_P0.rds')
metalist[["2019_loo_natcom"]]  <- rbind(m1,m2)

for (i in names(metalist)){
  metatmp <- metalist[[i]]
  colnames(metatmp) <- tolower(colnames(metatmp))
  metatmp$cell <- paste0(i, ';', metatmp$cell)
  metatmp$study <- i
  colnames(metatmp)[which(colnames(metatmp) == 'gender')] <- 'sex'
  metatmp <- data.frame(apply(metatmp, 2, as.character), stringsAsFactors = FALSE)
  metalist[[i]] <- metatmp
}
meta <- reshape::merge_recurse(metalist)
meta <- meta[match(colnames(mat), meta$cell), ]

# ------------
ct = tolower(meta$celltype)
ct[ct == 'oligodendrocytes precursor cells'] <- 'oligodendrocyte precursor cells'
ct[ct == 'oligodendrocytes'] <- 'oligodendrocyte'
ct[ct=="vascular and leptomeningeal cells"] <- "vascular and leptomeningeal cell"
ct[ct=="oligodendrocytes"] <- "oligodendrocyte"
ct[ct=="oligodendrocyte precursor cells"] <- "oligodendrocyte precursor cell"
ct[ct=="oligodendroctye precursor cells"] <- "oligodendrocyte precursor cell"
ct[ct=="neurons"] <- "neuron"
ct[ct=="ependyma"] <- "ependymal"
ct[ct=="ependymal"]  <- "ependymal" 
ct[ct=="endothelia cells"] <- "endothelial cell"
ct[ct=="astrocytes"] <- "astrocyte"
ct[ct=="astrocytes"] <- "astrocyte"
ct[ct=="vascular"] <- "vascular"
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

meta$time <- sub('embryo E', 'E', meta$time)

meta$time_super <- sapply(meta$time,function(i){
  if (grepl('E', i)){
    return('E')
  } else if (grepl('p', i)){
    return('P')
  } else {
    return(i)
  }
})

meta$species <- 'mouse'
meta$sex <- tolower(meta$sex)
meta$sex[which(meta$sex == 'f')] <- 'female'
meta$sex[which(meta$sex == 'm')] <- 'male'
saveRDS(meta, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/data/meta_allcell.rds')
