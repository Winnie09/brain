infer_adjmat <- function(udp){
  ## udp: a vector of cnv profiles, no duplicated. 
  u = unique(unlist(sapply(udp, function(i) strsplit(i, ';'))))
  mat = matrix(0, nrow = length(udp), ncol = length(u))
  rownames(mat) = udp
  colnames(mat) = u
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
  return(list(eglist = eglist, mat = mat))
}



library(igraph)
d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/summary/summary.rds')
ap = sub('_.*', '', names(d))
cutoff = 0.01
pdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvtree/allsampleplot_cutoff', cutoff, '/')
dir.create(pdir, recursive = T, showWarnings = F)

udp <- unique(unlist(sapply(unique(ap), function(p){
  print(p)
  dp = d[ap == p]
  tab = table(dp)
  names(tab[tab> length(dp)*cutoff])    
})))
res = infer_adjmat(udp)
eglist = res[['eglist']]
mat = res[['mat']]
saveRDS(res, paste0(pdir, 'cnvtree_adjacentMatrix_and_cnvPatternMatrix.rds'))
oranges <- colorRampPalette(c("dark red", "gold"))
len = rowSums(mat)
col = oranges(max(len)+1)
col = col[len+1]
names(col) = rownames(mat)
g <- graph_from_data_frame(eglist,vertices=data.frame(vertices=names(col),color=col),directed = T)
pdf(paste0(pdir, 'cnvtree.pdf'), width = 10, height = 10)
print(plot(g, main = paste0('all 36 samples, #CNV=', nrow(mat))))
dev.off()


