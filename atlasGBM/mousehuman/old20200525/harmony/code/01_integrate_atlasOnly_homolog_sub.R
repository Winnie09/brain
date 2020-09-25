library(Seurat)
library(harmony)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
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

countfunc <- function(mat, platform='fluidigm', normalize = TRUE,species='human', log.transform.only = FALSE) {
  row.names(mat) <- toupper(row.names(mat))
  rownames(mat) <- gsub('_.*','',rownames(mat))
  if (normalize){
    if (platform == '10x'){
      mat = log2CPM_from_10x_count(mat)
    } else {
      mat = log2TPM_from_fludigm_count(mat)
    } 
  }
  if (log.transform.only){
    mat = log2(mat + 1)
  }
  mat <- mat[!duplicated(rownames(mat)),]
  if (species=='human') {
    mat <- mat[intersect(row.names(mat),gn[,1]),]
  } else {
    mat <- mat[intersect(row.names(mat),gn[,2]),]
    row.names(mat) <- gn[match(row.names(mat),gn[,2]),1]
  }
  fullmat <- matrix(0,nrow=nrow(gn),ncol=ncol(mat))
  dimnames(fullmat) <- list(gn[,1], colnames(mat))
  fullmat[rownames(mat), colnames(mat)] <- mat
  rm(mat)
  set.seed(12345)
  if (ncol(fullmat) > 10000) fullmat <- fullmat[,sample(1:ncol(fullmat),10000)]
  fullmat  
}

## load atlas
dlist = list()
print('1')
dlist[['2018_Zeisel']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/sub10000mat.rds'), platform = '10x',species='mouse')
print('2')
dlist[['2018_Rosenberg']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/expr.rds'),species='mouse')
print('3')
a <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds'), platform = '10x', species='human')
b <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds'), species='human')
dlist[['2018_Li_Science']] <- cbind(a,b)
print('4')
a <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'), platform = '10x', species='human')
b <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'), platform = '10x', species='human')
c <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/VisualCortex.rds'), platform = '10x', species='human')
dlist[['2018_Lake_NatBiotech']] <- cbind(a, b, c)
###########
print('5')
dlist[['2018_Fan_CellResearch']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds'), normalize=FALSE, log.transform.only = T, species='human') 
print('6')
dlist[['allen']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds'), species='human')
print('7')
dlist[['2017_Nowakowski_Science']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2017_Nowakowski_Science/proc/log2tpm.rds'),normalize=F, species='human')  ## add
print('8')
dlist[['2016_Tasic_NatNeuro']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/log2tpm.rds'),normalize=F, species='mouse')
print('9')
a <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/EmbryoMoleculeCounts.rds'), species='human')
print('9')
b <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/ESMoleculeCounts.rds'), species='human')
print('10')
c <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseEmbryoMoleculeCounts.rds'), species='mouse')
print('11')
d <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/MouseAdultDAMoleculeCounts.rds'), species='mouse')
dlist[['2016_La']] <- cbind(a,b,c,d)
print('10')
dlist[['2015_Zeisel_Science']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/log2tpm.rds'), normalize=F, species='mouse')
print('11')
###########################################################################33333
study = rep(names(dlist), sapply(dlist,ncol))
dlist <- do.call(cbind,dlist)
brain <- CreateSeuratObject(counts = dlist, project = "brain", min.cells = 5) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = brain@var.genes, npcs = 20, verbose = FALSE)
rm(dlist)
brain@meta.data$stim <- study
brain@meta.data$dataset <- study
saveRDS(brain,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/integrated.rds')

brain <- RunHarmony(brain, "dataset")

pca_embeddings <- Embeddings(brain, 'pca')
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/',showWarnings = F, recursive = T)
saveRDS(pca_embeddings,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/pca_embeddings.rds')

harmony_embeddings <- Embeddings(brain, 'harmony')
saveRDS(harmony_embeddings,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/harmony_embeddings.rds')

print('18')
brain <- RunUMAP(brain, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

uamp <- brain@reductions$umap@cell.embeddings
saveRDS(umap,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/umap.rds')

