library(here)
here()
meta <- readRDS(here('atlasGBM','humanatlas','data','meta_allcell.rds'))
dlist <- readRDS(here('atlasGBM','humanatlas','data','dlist.rds'))
dlist <- dlist[!grepl('GBM',names(dlist))]
dlist <- dlist[!grepl('Pollen',names(dlist))]
dlist <- dlist[!grepl('allen_human_fluidigm',names(dlist))]
dlist <- do.call(cbind,dlist)
meta <- meta[!grepl('Pollen', meta$study), ]
meta <- meta[!grepl('GBM', meta$study), ]
meta <- meta[!grepl('allen_human_fluidigm', meta$study), ]

ct <- meta$celltype
names(ct) <- sub('.*;','',meta$cell)
ct <- ct[colnames(dlist)]

study <- meta$study
names(study) <- sub('.*;','',meta$cell)
study <- study[colnames(dlist)]

identical(sub('.*;','',meta$cell), colnames(dlist))
mydonor <- paste0(meta$study, ':',meta$gender,':', meta$age,':', meta$location, ':',meta$time, ':',meta$donor)
# > length(unique(mydonor))
# [1] 76

sct <- table(ct)
sct <- names(sct)[sct > 100]
sct <- setdiff(sct,'unknown')
samp <- paste0(ct,';',mydonor)

d <- sapply(unique(samp),function(ss) {
  if (!grepl('NA',sub(';.*','',ss)))
  rowSums(dlist[,samp==ss, drop=F])
})
if (is.list(d)) {
  d <- d[sapply(d,length)>0]
  d <- do.call(cbind, d)
}

d <- sweep(d,2,colSums(d)/1e6,'/')
d <- log2(d + 1)
library(preprocessCore)
dn <- dimnames(d)
d <- normalize.quantiles(d)
dimnames(d) <- dn

dcn <- sub(';.*','',colnames(d))
ddonor <- sub('.*;','', colnames(d))
eid <- expand.grid(1:length(sct),1:length(sct))
eid <- eid[eid[,1] > eid[,2],]
library(limma)
library(Matrix)


gl <- sapply(1:nrow(eid),function(i) {
  print(i)
  id <- c(which(dcn==sct[eid[i,1]]),which(dcn==sct[eid[i,2]]))
  
  tab <- table(ddonor[id])
  sel <- names(tab)[tab==2]
  id <- id[ddonor[id] %in% sel]
  uidct <- unique(dcn[id])
  dif <- rowMeans(d[,dcn[id]==uidct[1]])-rowMeans(d[,dcn[id]==uidct[2]])
  c(names(head(sort(dif),100)),names(head(sort(dif,decreasing = T),100)))
})
gl <- unique(as.vector(gl))
saveRDS(gl,file=here('atlasGBM','humanatlas','diffgene','res', 'celltype_diffgene_pairwise_difference_includeOnlyAllen10.rds'))


gl <- sapply(1:nrow(eid),function(i) {
  print(i)
  id <- c(which(dcn==sct[eid[i,1]]),which(dcn==sct[eid[i,2]]))
  if (length(id) > 2 & sum(table(ddonor[id])>=2) == length(unique(ddonor[id]))) {
    des <- model.matrix(~dcn[id]+ as.character(ddonor[id]))  
    res <- topTable(eBayes(lmFit(d[,id,drop=FALSE],design=des)),n=nrow(d),coef=2)
    row.names(res)[res[,'adj.P.Val'] < 0.05]
  }
})
len <- sapply(gl, length)
names(len) <- paste0(sct[eid[,1]], ';', sct[eid[,2]])
gl <- unique(unlist(gl))
saveRDS(gl,file=here('atlasGBM','humanatlas','diffgene','res', 'celltype_diffgene_pairwise_limma_includeOnlyAllen10x.rds'))



