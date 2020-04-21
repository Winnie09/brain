pca_lmFilter <- function(genebycellmat){
  d = genebycellmat
  d = d[rowSums(d)>0,]
  m = rowMeans(d)
  v = apply(d,1,sd)
  resid = resid(lm(v~bs(m)))
  d = d[names(sort(resid,decreasing = T))[1:1e3], ]
  prcomp(t(d),scale. = T)
}


pca_scranFilter <- function(genebycellmat){
  data = genebycellmat
  data = data[rowSums(data)>0,]
  suppressMessages(library(scran))
  fit <- trendVar(data)
  decomp <- decomposeVar(data,fit)
  gs <- row.names(decomp)[decomp[,'total'] > decomp[,'tech']]
  data <- data[gs,]  ###  [1:13662, 1:299]
  prcomp(t(data),scale. = T)  
}

pca_scranFilter_1k <- function(genebycellmat){
  data = genebycellmat
  data = data[rowSums(data)>0,]
  suppressMessages(library(scran))
  fit <- trendVar(data)
  decomp <- decomposeVar(data,fit)
  diff <- as.numeric(decomp[,'total']) - as.numeric(decomp[,'tech'])
  names(diff) = rownames(diff)
  diff = sort(diff[diff>0],decreasing = T)
  diff <- diff[1:(min(1000, length(diff)))]
  gs = names(diff)
  data <- data[gs,]  ###  [1:13662, 1:299]
  prcomp(t(data),scale. = T)  
}

getMatrixCluster <- function(data.path, result.path, numVariableGene = 1000, num.cluster= ncol(mat)/5, platform = 'none', use.pc = F){
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
  if (platform == '10x'){
    mat <-  log2CPM_from_10x_count(readRDS(data.path))  
  } else if (platform == 'fluidigm'){
    mat <-  log2TPM_from_fludigm_count(readRDS(data.path))  
  } else {
    mat <- readRDS(data.path)  
  }
  row.names(mat) <- toupper(row.names(mat)) 
  mat <- mat[!duplicated(row.names(mat)),]
  mat = mat[rowSums(mat>0)>0.01,]
  if (use.pc == F){
    print('calculating variable genes...')
    gs = findVariableGene(mat,numVariableGene)
    set.seed(12345)
    print('performing k-means clustering...')
    clu = kmeans(t(mat[gs,]), floor(num.cluster))$cluster
  } else {
    print('calculating principal components...')
    pr = PCA(genebycell_mat = mat, save.pca = FALSE, plot.statistics=FALSE, PC_for_columns = TRUE, findVariableGenes = TRUE, numPC = NULL)
    set.seed(12345)
    print('performing k-means clustering...')
    clu = kmeans(pr, floor(num.cluster))$cluster  
  }
  dir.create(result.path, recursive = T, showWarnings = F)
  saveRDS(clu, paste0(result.path,'clu.rds'))
  print('calculating cluster mean...')
  data = sapply(1:max(clu), function(c) rowMeans(mat[,names(clu[clu==c]),drop=F]) )
  colnames(data) = paste0('clu',1:max(clu))
  saveRDS(data, paste0(result.path,'cluster_mat.rds'))
  return(0)
}
