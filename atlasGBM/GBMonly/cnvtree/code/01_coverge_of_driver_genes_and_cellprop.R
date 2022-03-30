library(Seurat)
library(ggplot2)
library(scattermore)
library(RColorBrewer)
library(ggrepel)
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/brain')
cnv = readRDS('atlasGBM/GBMonly/cnvsel/cnv/summary/summary.rds')
alls = sub('_.*', '', names(cnv))


## subset 21 GBM untreated
m <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')
m <- unique(m[,-c(1:4)])
rownames(m) <- m$Sample.ID
m <- m[, c(4:7, 10:18)]
sampsel = rownames(m[m$Pathology == 'GBM' & m$Tumor.Grade == 'IV' & m$Treatment == 'Untreated', ])
str(sampsel)
str(alls)
cnv = cnv[alls %in% sampsel]
str(cnv)
alls = sub('_.*', '', names(cnv))



## 
tb = table(cnv, alls)
tb = tb/colSums(tb)
dim(tb)
tb[1:3, 1:3]
tmp = tb > 0.05
summary(rowSums(tmp))

tb2 = tb[which(rowSums(tb > 0.05)>1), ]
str(tb2)

# s= alls[1]


# cnvsel = cnv[alls == s]
# str(cnvsel)
# table(cnvsel)

cnv.uni = rownames(tb2)
len = sapply(cnv.uni, function(i){
  length(strsplit(i,';')[[1]])
})
length(cnv.uni)
str(len)
table(len)

trlist = list()
trlist[[1]] = c("chr10-", "chr10-;chr13-", "chr10-;chr13-;chr7+")
trlist[[2]] = c("chr10-", "chr10-;chr22-", "chr10-;chr22-;chr7+")
trlist[[3]] = c("chr10-", "chr10-;chr7+", "chr10-;chr21+;chr7+")
trlist[[4]] = c("chr10-", "chr10-;chr7+", "chr10-;chr22-;chr7+")
trlist[[5]] = c('chr18+')
trlist[[6]] = c('chr7+', "chr10-;chr7+", "chr10-;chr22-;chr7+")
trlist[[7]] = c('chr13-')
trlist[[8]] = c('chr21+')


## celltype
clu <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000//res/cluster.rds')
libsize = table(sub('_.*', '', names(clu)))


ct <- read.table('doc/GBM_cluster_annotations.csv', as.is = TRUE, sep = ',', header = T)
ct <- ct[, c(1,3)]
ct <- ct[match(clu, as.character(ct[,1])),2]
ct.gbm <- ct
names(ct.gbm) <- names(clu)
str(ct.gbm)
str(cnv)
ct.gbm = ct.gbm[names(cnv)]
df = data.frame(cell = names(ct.gbm), sample = sub('_.*','', names(ct.gbm)), ct = ct.gbm, cnv = cnv,  stringsAsFactors = FALSE)
rownames(df) = df[,1]

## for each cnv pattern, plot, sample by ct. 
tr.tmp = trlist[[1]]
cnv.mark = tr.tmp[1]
m = m[sampsel, !colnames(m) %in% c('Pathology', 'Tumor.Grade','Treatment', 'IDH.status', "CDKN2A..2G.loss")]
m$Sex = ifelse(m$Sex == 'F', 1, 0)
m$MGMT.status = ifelse(m$MGMT.status == 'methylated', 1, 0)
m$EGFR.amplification = ifelse(m$EGFR.amplification == 'amplified', 1, 0)
m$PDGFRA.amplification = ifelse(m$PDGFRA.amplification == 'amplified', 1, 0)
m$TP53.mutation = ifelse(m$TP53.mutation == 'mutated', 1, 0)
m$EGFR.mutation = ifelse(m$EGFR.mutation == 'mutated', 1, 0)
m$PTEN.mutation = ifelse(m$PTEN.mutation == 'mutated', 1, 0)




for (tr.tmp in trlist){
  for (cnv.mark in tr.tmp){
    # cnv.mark = 'chr10-;chr22-'
    print(tr.tmp)
    print(cnv.mark)
    cnv.tmp = cnv[cnv == cnv.mark]
    
    df.tmp = df[names(cnv.tmp),]
    cellprop = table(df.tmp[,2], df.tmp[,3])
    cellprop = cellprop/as.vector(libsize[rownames(cellprop)])
    str(cellprop)
    
    pdir <- paste0('atlasGBM/GBMonly/cnvtree/', '/', cnv.mark, '/')
    dir.create(pdir, showWarnings = FALSE)
    pdf(paste0(pdir, 'hm_cellprop_on_sample_by_ct.pdf'), width = 5, height = 6)
    pheatmap(cellprop)
    dev.off()
    
    ## for each cnv pattern, test sample covariate and ct association
    m.tmp = m[rownames(cellprop), ]
    pval = matrix(NA, nrow = ncol(m.tmp), ncol = ncol(cellprop))
    dimnames(pval) = list(colnames(m.tmp), colnames(cellprop))
    
    for (tmpct in colnames(pval)){
      for (tmpmeta in rownames(pval)){
        print(paste0(tmpct, ';', tmpmeta))
        fstat = summary(lm(cellprop[, tmpct] ~ model.matrix(~m.tmp[, tmpmeta])[,-1]))$fstatistic
        pval[tmpmeta, tmpct] = pf(fstat[1],fstat[2],fstat[3],lower.tail = F) ## tail
      }
    }
    
    fdrmat = matrix(p.adjust(as.vector(pval)), nrow = nrow(pval), ncol = ncol(pval))
    dimnames(fdrmat) = dimnames(pval)
    

    mat_breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 0.3,0.5,1)
    library(viridis)

    pdf(paste0(pdir, 'celltype_meta_association.pdf'), width = 5, height = 4.2)
    pheatmap(fdrmat, color = inferno(length(mat_breaks) - 1), breaks = mat_breaks, show_rownames = T, show_colnames = T)
    dev.off()
    ## for each cnv pattern, plot the pie chart for sample groups of significant covariate
    # sig.id = which(fdrmat < 0.1, arr.ind = T)
    # print(length(sig.id))
    # 
    # fdrmat[sig.id[1], sig.id[2]]
    # 
    # c(rownames(fdrmat)[sig.id[1]], colnames(fdrmat)[sig.id[2]])
    # 
    # str(cellprop)
    # v1 = cellprop[rownames(m)[m$TP53.mutation=='1'], 'OLGs']
    # 
    # v2 = cellprop[intersect(rownames(cellprop),rownames(m)[m$TP53.mutation=='0']), 'OLGs']
    # 
    # pd = rbind(data.frame(cellprop = v1, type = 'mutated'),
    #       data.frame(cellprop = v2, type = 'unmutated'))
  
    # library(ggplot2)  
    # ggplot(data = pd, aes(x = type, y = cellprop, color = type)) + 
    #   geom_boxplot()
  }
}




