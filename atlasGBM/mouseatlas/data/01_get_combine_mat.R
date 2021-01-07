## reference code:/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/integrate_GBM_to_atlas_referenceLi/code
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')     
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
countfunc <- function(mat, platform='fluidigm', log2normalize = TRUE, species='human', log2.only = FALSE, genename.sep = NULL, use.homolog = FALSE, subsample = NULL) {
  if (!is.null(genename.sep)){
    rownames(mat) <- gsub(paste0(genename.sep, '.*'),'',rownames(mat))
  }
  if (log2normalize){
    if (platform == '10x'){
      mat = log2CPM_from_10x_count(mat)
    } else {
      mat = log2TPM_from_fludigm_count(mat, species = species)
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
dlist = list()
## 1
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseEmbryoMoleculeCounts.rds')
str(a)
id.a <- which(colSums(a>0)>500)
a2 <- countfunc(a[,id.a], species='mouse', use.homolog = TRUE)
b = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseAdultDAMoleculeCounts.rds')
str(b)
id.b <- which(colSums(b>0)>500)
b2 <- countfunc(b[,id.b], species='mouse', use.homolog = TRUE)
int <- intersect(rownames(a2), rownames(b2))
dlist[['2016_La']] <- cbind(a2[int, ], b2[int, ])

## 2
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/log2tpm.rds')
str(a)
id <- which(colSums(a>0)>500)
dlist[['2015_Zeisel_Science']] <- countfunc(a[,id], log2normalize=FALSE, species='mouse', platform = 'fluidigm', use.homolog = TRUE)

## 3
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/log2tpm.rds')
id <- which(colSums(a>0)>500)
dlist[['2016_Tasic_NatNeuro']] <- countfunc(a[,id], log2normalize=FALSE, species='mouse', platform = 'fluidigm', use.homolog = TRUE)

## 4
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/count.rds')
str(a)
id <- which(colSums(a>0)>500)
dlist[['2018_Rosenberg']] <- countfunc(a[,id], platform = 'fluidigm', species='mouse', log2normalize=TRUE, use.homolog = TRUE)

## 5
# a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/sub10000mat.rds')
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/mat.rds')
str(a)
id <- colnames(a)[which(colSums(a>0)>500)]
dlist[['2018_Zeisel']] <- countfunc(a[, id], platform = '10x',species='mouse', log2normalize=TRUE, use.homolog = TRUE)

## 6
a <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/proc/count_E14.rds')
str(a)
id.a <- which(colSums(a>0)>500)
b <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/proc/count_P0.rds')
str(b)
id.b <- which(colSums(b>0)>500)
int <- intersect(rownames(a), rownames(b))
m <- cbind(a[int, id.a], b[int, id.b])
dlist[['2019_loo_natcom']] <- countfunc(m, platform = '10x', species = 'mouse', log2normalize=TRUE, use.homolog = TRUE)

## GBM
print(Sys.time())
gbm <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/proc/tumor/all/matrix/count.rds') ## num [1:20330, 1:289919] 
rownames(gbm) <- toupper(rownames(gbm))
print(Sys.time())

## select intersect genes
intgn <- sapply(dlist, row.names)
intgn <- table(unlist(intgn))
intgn <- names(intgn)[intgn==length(dlist)]
print(length(intgn))
intgn <- intersect(rownames(gbm), intgn)
print(length(intgn))
for (i in 1:length(dlist)) {
  dlist[[i]] <- dlist[[i]][intgn,]
  colnames(dlist[[i]]) = paste0(names(dlist)[i], ';', colnames(dlist[[i]]))
}
saveRDS(dlist,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/data/dlist.rds')
mat = do.call(cbind, dlist)
saveRDS(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseatlas/data/combine_mat.rds')


