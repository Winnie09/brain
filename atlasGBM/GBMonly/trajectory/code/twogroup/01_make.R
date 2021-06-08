ord <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pseudotime_3_4.rds')
pt <- 1:length(ord)
names(pt) <- ord
expr <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/imp/saver/combine/log2norm_matrix_qc.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/meta.data.rds')
meta <- meta[,c('orig.ident','Location')]
meta <- unique(meta)
meta <- meta[meta[,1] %in% unique(sub('_.*','',ord)),]
uniloc <- unique(meta[,2])
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

comp <- expand.grid(1:4,1:4)
comp <- comp[comp[,1] > comp[,2],]
compid <- as.numeric(commandArgs(trailingOnly = T)) #####
print(paste0('compid = ', compid))
loc <- uniloc[unlist(comp[compid,])]
selectsamp <- as.character(meta[meta[,2] %in% loc,1])
selectcell <- ord[sub('_.*','',ord) %in% selectsamp]

expr <- expr[,selectcell]
expr <- expr[rowMeans(expr > 0.1) > 0.01,]
cellanno <- data.frame(cell=selectcell,sample=sub('_.*','',selectcell),stringsAsFactors = F,row.names = selectcell)
design=data.frame(intercept=1,contrast=as.numeric(selectsamp %in% meta[meta[,2]==loc[1],1]),row.names = selectsamp)
identical(sort(rownames(design)), sort(unique(cellanno[,2])))
pt <- pt[selectcell]
res <- testpt(expr,test.type='VARIABLE', cellanno=cellanno, pseudotime=pt, design = design,ncores=20)

saveRDS(res,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/res/twogroup/',loc[1],'_',loc[2],'.rds'))

