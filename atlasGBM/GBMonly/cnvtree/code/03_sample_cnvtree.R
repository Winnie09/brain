d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/summary/summary.rds')
ap = sub('_.*', '', names(d))
library(igraph)

meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/metas_20200329.csv', row.names = 1)
meta = meta[unique(ap),]

for (cutoff in c(0.01, 0.05)){
  print(cutoff)
  for (p in unique(ap)){
    pdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvtree/sampleplot_cutoff', cutoff, '/')
    dir.create(pdir, recursive = T, showWarnings = F)
    print(p)
    dp = d[ap == p]
    tab = table(dp)
    udp = names(tab[tab> length(dp)*cutoff])
    u = unique(unlist(sapply(udp, function(i) strsplit(i, ';'))))
    mat = matrix(0, nrow = length(udp), ncol = length(u))
    rownames(mat) = udp
    colnames(mat) = u
    if (nrow(mat) == 0) next
    for (i in udp){
      v = strsplit(i, ';')
      for (j in v)
        mat[i, j] = 1
    }  
    
    eg <- expand.grid(1:nrow(mat),1:nrow(mat))
    eg <- eg[eg[,1] > eg[,2],]
    if (nrow(eg) == 0) next
    
    id <- which(rowSums(abs(mat[eg[,1],, drop=F]-mat[eg[,2],, drop=F]))==1)
    if (length(id) == 0) next
    
    eg <- eg[id,,drop=F]
    rs <- rowSums(mat)
    egrs <- cbind(rs[eg[,1]],rs[eg[,2]])
    egrev <- eg[,c(2,1)]
    
    egfinal <- t(sapply(1:nrow(egrs),function(i) {
      if (egrs[i,1] < egrs[i,2]) {
        unlist(eg[i,])
      } else {
        unlist(egrev[i,])
      }
    }))
    
    eglist <- data.frame(rownames(mat)[egfinal[,1]],rownames(mat)[egfinal[,2]])
    saveRDS(list(eglist = eglist, mat = mat), paste0(pdir, 'cnvtree_adjacentMatrix_and_cnvPatternMatrix/', p, '_cnvtree_adjacentMatrix_and_cnvPatternMatrix.rds'))
    oranges <- colorRampPalette(c("dark red", "gold"))
    len = rowSums(mat)
    col = oranges(max(len)+1)
    col = col[len+1]
    names(col) = rownames(mat)
    g <- graph_from_data_frame(eglist,vertices=data.frame(vertices=names(col),color=col),directed = T)
    pdf(paste0(pdir, 'cnvtree/', p, '_cnvtree.pdf'), width = 8, height = 8)
    print(plot(g, main = paste0(p, ';', meta[p, 'Pathology'], ';Grade', meta[p, 'Tumor.Grade'], ';', meta[p, 'Treatment'])))
    dev.off()
  }
}



