### this data generation process is relied on previously generated data in cnv_myelymreference_with_doublet. 
### In cnv_myelymreference_with_doublet, the data is generated basedon the right process, but the original GBM cells were not filtered by doubletFinder. 
### In this step, we further remove the cells found by doubletFinder. 

library(data.table)
library(here)
setwd(here())

af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference_with_doublet/data/data/', pattern =  '_countmatrix.txt')
ap = gsub('36nonNormal_myelym_combined_', '', af)
ap = gsub('_countmatrix.txt', '', ap)

for (p in ap){
  print(p)
  ### cell annotation file
  anno = read.table(paste0('atlasGBM/GBMonly/cnv_myelymreference_with_doublet/data/data/36nonNormal_myelym_combined_', p, '_cellanno.txt'), sep= '\t', as.is = T)
  colnames(anno) = c('cell', 'cluster')
  anno_ref = anno[anno[,2] %in% c('lymph', 'myeloid'), ]
  
  clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
  patient = gsub('_.*', '', names(clu))  
  id.p <- which(patient == p)
  tab <- table(clu[id.p])
  id <- which(clu %in% names(tab)[tab > 2])
  id <- intersect(id, id.p)
  
  anno_gbm <- data.frame(cell = names(clu)[id], cluster =  paste0('cluster ', clu[id]), stringsAsFactors = FALSE)
  anno <- rbind(anno_gbm, anno_ref)
  
  
  ### write gene by cell expression matrix
  expr = fread(file = paste0('atlasGBM/GBMonly/cnv_myelymreference_with_doublet/data/data/36nonNormal_myelym_combined_', p, '_countmatrix.txt'), sep = '\t', stringsAsFactors = F, header = T, data.table = F)
  rn = expr[,1]
  expr = expr[,-1]
  rownames(expr) = rn
  
  ## intersect expression cells and annotation cells
  intcell = intersect(colnames(expr), anno[,1])
  str(intcell)
  # anno_int = anno[anno[,1] %in% intcell, ] ## do not use this because tumor and mye/lym have several (~10) duplicated names. 
  anno_int = anno[match(intcell, anno[,1]), ]  ## in that case, we preseve the tumor samples and ignore just several reference cells.
str(anno_int)
  write.table(anno_int, paste0('atlasGBM/GBMonly/cnv_myelymreference/data/data/36nonNormal_myelym_combined_', p, '_cellanno.txt'), quote=F, sep="\t", row.names = FALSE, col.names = FALSE)  ##
  expr = expr[, anno_int[,1]]
  fwrite(expr, file = paste0('atlasGBM/GBMonly/cnv_myelymreference/data/data/36nonNormal_myelym_combined_', p, '_countmatrix.txt'), quote = F, sep = '\t', row.names = T)
  
  print(identical(ncol(expr), nrow(anno_int)))
  
}



