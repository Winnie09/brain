### this data generation process is relied on previously generated data in cnv_myelymreference_with_doublet. 
### In cnv_myelymreference_with_doublet, the data is generated basedon the right process, but the original GBM cells were not filtered by doubletFinder. 
### In this step, we further remove the cells found by doubletFinder. 

library(data.table)
library(here)
setwd(here())
af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference_with_doublet/data/data/', pattern =  '_countmatrix.txt')
ap = gsub('36nonNormal_myelym_combined_', '', af)
ap = gsub('_countmatrix.txt', '', ap)

clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
patient = gsub('_.*', '', names(clu))  

p_remain = setdiff(unique(patient), ap)   ## patients that need to generate data

### get reference cell anno  
anno = read.table(paste0('atlasGBM/GBMonly/cnv_myelymreference_with_doublet/data/data/36nonNormal_myelym_combined_GBM006_cellanno.txt'), sep= '\t', as.is = T)
colnames(anno) = c('cell', 'cluster')
anno_ref = anno[anno[,2] %in% c('lymph', 'myeloid'), ]
rm(anno)
  
### get reference gene expression
expr = fread(file = paste0('atlasGBM/GBMonly/cnv_myelymreference_with_doublet/data/data/36nonNormal_myelym_combined_GBM006_countmatrix.txt'), sep = '\t', stringsAsFactors = F, header = T, data.table = F)
rn = expr[,1]
expr = expr[,-1]
rownames(expr) = rn
refcell = intersect(colnames(expr), anno[anno[,2] %in% c('lymph', 'myeloid'),1])
refexpr = expr[,refcell]
rm(expr)

### get GBM all cell gene expression
exprall <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_count_mat.rds')
exprall.p = gsub('_.*', '', colnames(exprall))

for (p in p_remain){
  print(p)
  ### cell annotation file (using clu names from seurat that after rm doublets)
  id.p <- which(patient == p)
  tab <- table(clu[id.p])
  id <- which(clu %in% names(tab)[tab > 2])
  id <- intersect(id, id.p)
  anno_gbm <- data.frame(cell = names(clu)[id], cluster =  paste0('cluster ', clu[id]), stringsAsFactors = FALSE)
  anno <- rbind(anno_gbm, anno_ref)
  
  ### gene by cell expression matrix
  expr.p = exprall[rownames(refexpr), exprall.p == p]
  expr = cbind(expr.p, refexpr)
  
  ## intersect expression cells and annotation cells
  intcell = intersect(colnames(expr), anno[,1])
  str(intcell)
  anno_int = anno[match(intcell, anno[,1]), ]
  str(anno_int)
  write.table(anno_int, paste0('atlasGBM/GBMonly/cnv_myelymreference/data/data/36nonNormal_myelym_combined_', p, '_cellanno.txt'), quote=F, sep="\t", row.names = FALSE, col.names = FALSE)  ##
  
  expr_int = expr[, anno_int[,1]]
  dim(expr_int)
  fwrite(expr_int, file = paste0('atlasGBM/GBMonly/cnv_myelymreference/data/data/36nonNormal_myelym_combined_', p, '_countmatrix.txt'), quote = F, sep = '\t', row.names = T)
  print(identical(ncol(expr_int), nrow(anno_int)))
}


