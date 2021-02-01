library(Seurat)
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/resource/01_function.R')
## load atlas
dlist = list()
print('1')
dlist[['2018_Fan_CellResearch']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/cluster_mat.rds')
print('2')
dlist[['2018_Lake_NatBiotech:VisualCortex']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/cluster_mat/VisualCortex/cluster_mat.rds')
dlist[['2018_Lake_NatBiotech:CerebellarHem']] <- log2CPM_from_10x_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/CerebellarHem.rds'))
dlist[['2018_Lake_NatBiotech:FrontalCortex']] <- log2CPM_from_10x_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/FrontalCortex.rds'))
print('3')
dlist[['2018_Li_Science:adult']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/cluster_mat/adult/cluster_mat.rds')
dlist[['2018_Li_Science:fetal']] <- log2CPM_from_10x_count(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds'))
print('4')
dlist[['allen']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/cluster_mat/cluster_mat.rds'))
print('5')
###### running >>>>>>>>>>>>>>
dlist[['2018_Rosenberg']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/cluster_mat.rds'))
row.names(dlist[['2018_Rosenberg']]) <- toupper(row.names(dlist[['2018_Rosenberg']]))  ## mouse
dlist[['2018_Rosenberg']] <- dlist[['2018_Rosenberg']][!duplicated(row.names(dlist[['2018_Rosenberg']])),]
print('6')
dlist[['2018_Zeisel']] <- countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/cluster_mat.rds'))
row.names(dlist[['2018_Zeisel']]) <- toupper(row.names(dlist[['2018_Zeisel']]))  ## mouse
dlist[['2018_Zeisel']] <- dlist[['2018_Zeisel']][!duplicated(row.names(dlist[['2018_Zeisel']])),]
##############3 <<<<<<<<<<<<<
print('7')
dlist[['2017_Nowakowski_Science']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2017_Nowakowski_Science/proc/cluster_mat.rds')  ## add
print('8')
dlist[['2016_La:EmbryoMoleculeCounts']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/cluster_mat/EmbryoMoleculeCounts/cluster_mat.rds')
print('9')
dlist[['2016_La:ESMoleculeCounts']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/cluster_mat/ESMoleculeCounts/cluster_mat.rds')
print('10')
dlist[['2016_La:MouseEmbryoMoleculeCounts.rds']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/cluster_mat/MouseEmbryoMoleculeCounts/cluster_mat.rds')
print('11')
dlist[['2016_La:MouseAdultDAMoleculeCounts.rds']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/cluster_mat/MouseAdultDAMoleculeCounts/cluster_mat.rds')
print('12')
dlist[['2016_Tasic_NatNeuro']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/cluster_mat.rds')
row.names(dlist[['2016_Tasic_NatNeuro']]) <- toupper(row.names(dlist[['2016_Tasic_NatNeuro']])) ## mouse
dlist[['2016_Tasic_NatNeuro']] <- dlist[['2016_Tasic_NatNeuro']][!duplicated(row.names(dlist[['2016_Tasic_NatNeuro']])),]
print('13')
dlist[['2015_Zeisel_Science']] <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/cluster_mat.rds') ## add
row.names(dlist[['2015_Zeisel_Science']]) <- toupper(row.names(dlist[['2015_Zeisel_Science']]))  ## mouse
dlist[['2015_Zeisel_Science']] <- dlist[['2015_Zeisel_Science']][!duplicated(row.names(dlist[['2015_Zeisel_Science']])),]
print('14')
# gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/count.rds')
gn <- table(unlist(sapply(dlist,row.names)))
gn <- names(gn)[gn==length(dlist)]
print('16')
homolog = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/homologgene/homologene.data.txt', data.table=F)
humanhomolog = homolog[which(homolog[,2]=='9606'),4]
mousehomolog = toupper(homolog[which(homolog[,2]=='10090'),4])
hlgene = length(intersect(humanhomolog,mousehomolog))
gn = intersect(gn, hlgene)
str(gn)
print('17')
# gn <- intersect(gn,row.names(gbm))
# gbm <- gbm[gn,]
# p = sub('_.*','',colnames(gbm))
for (i in names(dlist)) {
  print('findvariablefeatures...')
  print(i)
  dlist[[i]] <- FindVariableFeatures(CreateSeuratObject(dlist[[i]][gn,],project=i))
}
print('FindIntegrationAnchors...')
# for (sp in unique(p)) {
#   dlist[[sp]] <- FindVariableFeatures(CreateSeuratObject(countfunc(gbm[,p==sp]),project=sp))
# }
# reference.vec = c('2018_Rosenberg','2018_Zeisel','2018_Li_Science:fetal','2018_Li_Science:adult','2018_Lake_NatBiotech:VisualCortex','2018_Lake_NatBiotech:FrontalCortex','2018_Lake_NatBiotech:CerebellarHem','2018_Fan_CellResearch','2016_La:ESMoleculeCounts','2016_La:EmbryoMoleculeCounts')
# d.anchors <- FindIntegrationAnchors(object.list = dlist, reference = grep('GBM',names(dlist),invert = T))
# d.integrated <- IntegrateData(anchorset = d.anchors)
d.anchors <- FindIntegrationAnchors(object.list = dlist)
print('IntegrateData...')
d.integrated <- IntegrateData(anchorset = d.anchors)
print('20')
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlasOnly_homolog_sub_cluster_mat/',showWarnings = F, recursive = T)
saveRDS(d.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlasOnly_homolog_sub_cluster_mat/integrated.rds')
