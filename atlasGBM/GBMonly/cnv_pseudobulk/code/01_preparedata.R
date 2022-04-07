library(data.table)
af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', pattern = '_countmatrix.txt')

mat.gbm = lapply(af, function(f){
  print(f)
  cnt = fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', f), data.table = F)
  cellanno = fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', sub('_countmatrix', '_cellanno', f)), data.table = F, header = F)
  rownames(cnt) = cnt[,1]
  cnt = as.matrix(cnt[,-1])
  v = rowSums(cnt[, cellanno[grep('cluster',cellanno[,2]),1]])  
})
mat.gbm = do.call(cbind, mat.gbm)

mat.lym = lapply(af, function(f){
  print(f)
  cnt = fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', f), data.table = F)
  cellanno = fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', sub('_countmatrix', '_cellanno', f)), data.table = F, header = F)
  rownames(cnt) = cnt[,1]
  cnt = as.matrix(cnt[,-1])
  v = rowSums(cnt[, cellanno[grep('lymph',cellanno[,2]),1]])  
})
mat.lym = do.call(cbind, mat.lym)



mat.mye = lapply(af, function(f){
  print(f)
  cnt = fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', f), data.table = F)
  cellanno = fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/data/data/', sub('_countmatrix', '_cellanno', f)), data.table = F, header = F)
  rownames(cnt) = cnt[,1]
  cnt = as.matrix(cnt[,-1])
  v = rowSums(cnt[, cellanno[grep('myeloid',cellanno[,2]),1]])  
})
mat.mye = do.call(cbind, mat.mye)
colnames(mat.gbm) <- colnames(mat.lym) <- colnames(mat.mye) <- gsub('36nonNormal_myelym_combined_', '', gsub('_countmatrix.txt', '', af))

colnames(mat.lym) = paste0(colnames(mat.lym), '_lymph')
colnames(mat.mye) = paste0(colnames(mat.mye), '_myeloid')
mat = cbind(mat.gbm, mat.mye, mat.lym)
mat = as.data.frame(mat)

meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/metas_20200329.csv', as.is = TRUE)
rownames(meta) = meta[,1]
meta = meta[colnames(mat.gbm), ]
meta[grep('Astrocytoma', meta[,'Pathology']),'Pathology'] = 'AST'
meta[grep('Oligodendroglioma', meta[,'Pathology']),'Pathology'] = 'OLI'
meta[grep('immunotherapy', meta[,'Treatment']),'Treatment'] = 'Immunotherapy'

cellanno = data.frame(cell = colnames(mat), 
                      cluster = c(paste0(meta[,'Pathology'], ';', meta[,'Treatment']),
                                  rep('ref_myeloid', ncol(mat.mye)),
                                  rep('ref_lymph', ncol(mat.lym))),
                      stringsAsFactors = FALSE)
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_pseudobulk/data/'
write.table(cellanno, paste0(rdir, 'cellanno.txt'), quote=F, sep="\t", row.names = FALSE, col.names = FALSE)  ##
fwrite(mat, file = paste0(rdir, 'countmatrix.txt'), quote = F, sep = '\t', row.names = TRUE)

