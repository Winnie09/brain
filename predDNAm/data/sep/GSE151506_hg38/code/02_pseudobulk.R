rm(list=ls())
# ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRNA/'
# pb <- readRDS(paste0(ddir, 'pseudobulk_ge.rds'))

## meta data
mdir <- '/scratch16/hji7/whou10/brain/predDNAm/data/sep/GSE151506_hg38/meta/'
tb1 <- read.csv(paste0(mdir, 'GSE151506_suppltb5_GBM_malignant_classification.csv'))
tb2 <- read.csv(paste0(mdir, 'GSE151506_suppltb5_IDHMutant_malignant_classification.csv'))
selcell <- c(tb1$CellAssignment, tb2$CellAssignment)
names(selcell) <- c(tb1[,12], tb2[,11])
tmpname <- c(tb1[,1], tb2[,2])
tmpname <- tmpname[nchar(selcell) > 0]
selcell <- selcell[nchar(selcell) > 0]
names(selcell) <- sub('P_1', 'P1', names(selcell))
names(selcell) <- sub('P_2', 'P2', names(selcell))
names(selcell) <- sub('\\.', '-', names(selcell))
names(selcell) <- gsub('GBM_', '', names(selcell))


## scRRBS
ddir2 <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/raw/nygc/CSV/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'
allf = list.files(ddir2, pattern = 'MGH')
alls <- substr(allf, 1, 6)
alls <- sub('\\.', '', alls)
names(alls) <- allf

selcelldf <- data.frame(type = selcell, 
                        scRNA_name = tmpname, 
                        scRRBS_name = names(selcell),
                        stringsAsFactors = FALSE)
saveRDS(selcelldf, paste0(rdir, 'cells_with_annotation.rds'))

## read in all samples
alld <- lapply(allf, function(f) {
  # alld <- lapply(allf[1:2], function(f) {  #########################
  print(f)
  me <- read.csv(paste0(ddir2, f))
  rownames(me) <- me[, 1]
  me <- me[, -1]
  me
})
str(alld)

## check colnames
for (i in 1:length(alld)){
  ## no match
  # i = 17
  # i = 16
  # i = 15 # no
  ## 
  print(i)
  me = alld[[i]]
  cn <- sub('\\.R.*','',colnames(me))
  cn <- gsub('\\.','-',cn)
  cn <- gsub('-1-Pool','_P1-Pool',cn)
  cn <- gsub('-2-Pool','_P2-Pool',cn)
  cn <- gsub('-3-Pool','_P3-Pool',cn)
  cn <- gsub('-4-Pool','_P4-Pool',cn)
  cn <- gsub('GBM_','',cn)
  colnames(me) <- cn
  tmp <- names(selcell)[grep(sub('\\..*','',sub('_.*','',allf[i])), names(selcell))]
  int <- intersect(cn, tmp)
  str(int)
}
saveRDS(alld, paste0(rdir, 'cells_scRRBS_each_sample_annnotated.rds'))

## pseudobulk for each sample
pb.allf <- lapply(allf, function(f) {
  # pb.allf <- lapply(allf[1:2], function(f) {  ############################
  print(f)
  me <- read.csv(paste0(ddir2, f))
  rownames(me) <- me[, 1]
  me <- me[, -1]
  ## change column names for scRRBS matrix
  cn <- sub('\\.R.*','',colnames(me))
  cn <- gsub('\\.','-',cn)
  cn <- gsub('-1-Pool','_P1-Pool',cn)
  cn <- gsub('-2-Pool','_P2-Pool',cn)
  cn <- gsub('-3-Pool','_P3-Pool',cn)
  cn <- gsub('-4-Pool','_P4-Pool',cn)
  cn <- gsub('GBM_','',cn)
  colnames(me) <- cn
  ## check whether colnames exist in cells with annotations
  len <- length(intersect(cn, names(selcell)))
  ### if yes
  if (len > 0){
    me.tmp <- me[, intersect(cn, names(selcell))]
    dn <- dimnames(me.tmp)
    me.tmp <- as.numeric(as.matrix(me.tmp))
    me.tmp <- matrix(me.tmp, nrow = nrow(me))
    dimnames(me.tmp) <- dn
    
    ## calculate pseudobulks within this tumor, average each celltype
    ### select the cells with annotations
    selcell.sub <- selcell[colnames(me.tmp)] 
    
    ## average within each celltype
    pb <- sapply(sort(unique(selcell.sub)), function(i){
      v <- rowMeans(me.tmp[,selcell.sub == i, drop = FALSE], na.rm = TRUE)
    })
    str(pb)
    ## colnames(pb) <- paste0(sub('.csv','',f), ';', sort(unique(selcell.sub)))
    pb <- pb[rowMeans(is.na(pb)) < 1,]
  } else {
    return(NULL)
  }
})
names(pb.allf) <- allf
saveRDS(pb.allf, paste0(rdir, 'pseudobulk_scRRBS_each_sample.rds'))

v <- unlist(sapply(pb.allf, ncol))
names(pb.allf)
pb.allf <- pb.allf[names(v)]

alls <- alls[names(pb.allf)]
str(alls)

## pseudobulk for each tumor (multiple samples)
pb.allt <- lapply(sort(unique(alls)), function(s){
  print(s)
  af <- names(alls)[alls == s]
  tmp <- pb.allf[names(pb.allf) %in% af]
  ## because the row names (CpGs) in each sample are different
  ## new a matrix to save the union of all row names
  ## augment the columns to incluce all cell types
  rn = sort(unique(unlist(sapply(tmp, rownames))))
  cn <- sort(unique(selcell))
  for (i in 1:length(tmp)){
    print(i)
    mat <- matrix(NA, nrow = length(rn), ncol = length(cn))
    dimnames(mat) <- list(rn, cn)
    mat[rownames(tmp[[i]]), colnames(tmp[[i]])] <- tmp[[i]]
    tmp[[i]] <- mat
  }
  
  ## calculate pseudobulks within this tumor, average each celltype, na.rm
  mat.all <- sapply(tmp, function(i){
    as.vector(i)
  })
  
  colMeans(mat.all, na.rm = T)
  mat.all <- rowMeans(mat.all, na.rm = TRUE)
  
  mat.all <- matrix(mat.all, nrow = length(rn))
  dimnames(mat.all) <- list(rn, paste0(s, ';', cn))
  mat.all <- mat.all[,!is.na(colMeans(mat.all, na.rm = T))]
})
names(pb.allt) <- sort(unique(alls))
saveRDS(pb.allt, paste0(rdir, 'pseudobulk_scRRBS_each_tumor_list.rds'))


## combine all tumors' pseudobulk in a large matrix
rn = sort(unique(unlist(sapply(pb.allt, rownames))))
cn <- unlist(sapply(pb.allt, colnames))
mat <- matrix(NA, nrow = length(rn), ncol = length(cn))
dimnames(mat) <- list(rn, cn)
for (i in 1:length(pb.allt)){
  mat[rownames(pb.allt[[i]]), colnames(pb.allt[[i]])] <- pb.allt[[i]]
}
saveRDS(mat, paste0(rdir, 'pseudobulk_scRRBS_each_tumor_matrix.rds'))
str(mat) ## [1:16750766, 1:72] 


