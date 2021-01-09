

########3
library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')     
countfunc <- function(mat, platform='fluidigm', log2normalize = TRUE, species='human', log2.only = FALSE, genename.sep = NULL, use.homolog = FALSE, subsample = NULL) {
  if (!is.null(genename.sep)){
    rownames(mat) <- gsub(paste0(genename.sep, '.*'),'',rownames(mat))
  }
  if (log2normalize){
    if (platform == '10x'){
      mat = log2CPM_from_10x_count(mat)
    } else {
      mat = log2TPM_from_fludigm_count(mat)
    } 
  } else if (log2.only){
    mat = log2(mat + 1)
  } else if (log2normalize & log2.only){
    print("log2normalize and log2.only CANNOT be both TRUE! ")
  }
  mat <- mat[!duplicated(rownames(mat)),]
  if (use.homolog){
    library(data.table)
    row.names(mat) <- toupper(row.names(mat))
    homolog = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/homologgene/homologene.data.txt', data.table=F)
    homolog = homolog[homolog[,2] %in% c('9606','10090'),]
    num <- sapply(unique(homolog[,1]),function(i) { 
      sum(homolog[homolog[,1]==i,2]=='9606')==1 & sum(homolog[homolog[,1]==i,2]=='10090')==1
    })
    tar <- unique(homolog[,1])[which(num)]
    
    t1 <- homolog[homolog[,1] %in% tar & homolog[,2]=='9606',]
    t2 <- homolog[homolog[,1] %in% tar & homolog[,2]=='10090',]
    identical(t1[,1],t2[,1])
    gn <- cbind(t1[,4],t2[,4])
    gn[,1] <- toupper(gn[,1])
    gn[,2] <- toupper(gn[,2])

    if (species=='human') {
      mat <- mat[intersect(row.names(mat),gn[,1]),]
    } else if (species == 'mouse'){
      mat <- mat[intersect(row.names(mat),gn[,2]),]
      row.names(mat) <- gn[match(row.names(mat),gn[,2]),1]
    } else {
      print('Mouse or human only! ')
    }
    fullmat <- matrix(0,nrow=nrow(gn),ncol=ncol(mat))
    dimnames(fullmat) <- list(gn[,1], colnames(mat))
    fullmat[rownames(mat), colnames(mat)] <- mat
  } else {
    fullmat = mat
  }
  rm(mat)  
  set.seed(12345)
  if (!is.null(subsample)){
    if (ncol(fullmat) > subsample) fullmat <- fullmat[,sample(1:ncol(fullmat), subsample)]
  }
  fullmat  
}

## load atlas
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
dlist = list()
## 1
a <- countfunc(readRDS('./2018_Li_Science/proc/expr/adult.rds'), platform = '10x', species='human')
b <- countfunc(readRDS('./2018_Li_Science/proc/expr/fetal.rds'), species='human')
int = intersect(rownames(a), rownames(b))
dlist[['2018_Li_Science']] <- cbind(a[int,], b[int,])
## 2
a <- countfunc(readRDS('./2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'), platform = '10x', species='human')
b <- countfunc(readRDS('./2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'), platform = '10x', species='human')
c <- countfunc(readRDS('./2018_Lake_NatBiotech/proc/count/VisualCortex.rds'), platform = '10x', species='human')
int = intersect(intersect(rownames(a), rownames(b)), rownames(c))
dlist[['2018_Lake_NatBiotech']] <- cbind(a[int,], b[int,], c[int,]) 
## 3
dlist[['2018_Fan_CellResearch']] <- countfunc(readRDS('./2018_Fan_CellResearch/proc/log2expr.rds'), log2normalize=FALSE, log2.only = FALSE, species='human') 
## 4
dlist[['allen_human_fluidigm']] <- countfunc(readRDS('./allen/data/proc/count.rds'), species='human')
## 5
dlist[['2017_Nowakowski_Science']] <- countfunc(readRDS('./2017_Nowakowski_Science/proc/log2tpm.rds'),log2normalize=F, species='human')  ## add
## 6
dlist[['2018_zhong_nature']] <- countfunc(readRDS('./2018_zhong_nature/proc/tpm.rds'), platform = 'fluidigm', log2normalize = FALSE, log2.only = TRUE, species = 'human')
## 7
dlist[['2015_Darmanis_PNAS']] <- countfunc(readRDS('./2015_Darmanis_PNAS/proc/count.rds'), log2normalize = TRUE, platform = 'fluidigm')
## 8
dlist[['2014_Pollen_NatBiotech']] <- readRDS('./2014_Pollen_NatBiotech/proc/log2tpm.rds')
## 9
dlist[['allen_human_10x']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/count.rds'), platform='10x', log2normalize = TRUE, species='human')
## 10
expr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/proc/expr.rds')
ct <- sub('_.*', '', colnames(expr))
mat <- expr[, ct %in% c('HESC', 'AdultCerebellum', 'FetalBrain')]
colnames(mat) <- paste0('2020_Guo_Nature:', colnames(mat))
dlist[['2020_Guo_Nature']] <- countfunc(mat, platform = '10x', log2normalize = T, species = 'human')

## gbm
print(Sys.time())
gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix_20200329/count.rds')
print(Sys.time())

## select intersect genes
intgn <- sapply(dlist,row.names)
intgn <- table(unlist(intgn))
intgn <- names(intgn)[intgn==length(dlist)]
intgn <- intersect(rownames(gbm), intgn)

for (i in 1:length(dlist)) {
  dlist[[i]] <- dlist[[i]][intgn,]
  colnames(dlist[[i]]) = paste0(names(dlist)[i], ';', colnames(dlist[[i]]))
}
saveRDS(dlist, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/dlist.rds')

mat = do.call(cbind, dlist)
saveRDS(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/combine_mat.rds')


