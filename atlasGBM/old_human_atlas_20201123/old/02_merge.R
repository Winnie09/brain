setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/')
mat = readRDS('./proc/scran/norm_genebycell.rds')
meta = read.csv('./raw/sample_annotations.csv',stringsAsFactors = F)
identical(as.character(meta$sample_name), colnames(mat))
ct = meta[,'cell_type_alias_label']
majorCellClass = sub(' .*','', ct)
ag0 = rownames(mat)

af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/processed/scran/')
d <- sapply(af,function(f){
  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/processed/scran/',f))
})
batch1 = paste0('2018_Lake_NatBiotech;',rep(sub('.rds','',af), sapply(d,ncol)))
ag1 = rownames(d[[1]])
for ( i in 2:length(d)){
  ag1 <- intersect(rownames(d[[i]]),ag)  
}


d2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/processed/scran/Sestan.adultHumanNuclei.Psychencode.rds')
ag2 = rownames(d2)
meta2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/processed/meta/Sestan.adultHumanNuclei.Psychencode.rds')
meta2 = meta2[match(colnames(d2),rownames(meta2)),]

ag = intersect(ag0, intersect(ag1,ag2))
d1 = sapply(d,function(subd){
  subd[ag,]
})
d1 = do.call(cbind,d1)


res <- cbind(mat[ag,], d1[ag,],d2[ag,])
batch <- c(rep('allenAtlas',ncol(mat)),  batch1, rep('2018_Li_Science;adult',ncol(d2)))
cellMeta = data.frame(cell = colnames(res), 
                      batch = batch, 
                      subtype = c(ct, colnames(d1), meta2$subtype),
                      ctype = c(majorCellClass, sub('_.*','',colnames(d1)), meta2$ctype)
                      )
saveRDS(res,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/genebycell.rds')
saveRDS(batch,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/cellMeta.rds')
