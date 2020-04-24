input.seed = as.numeric(commandArgs(trailingOnly = T)[1])
gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/count.rds')
gbm = log2(t(t(gbm)/colSums(gbm)/1e6)+1)
order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/order.rds')
order = order[grepl('GBM', order$sample_name),]
ap = gsub('_.*','', order$sample_name)
t_true = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/t_statistics.rds')
pval = 0

dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3_pm/',input.seed))
gbm.bak = gbm

for (myseed in (input.seed*100) : (input.seed*100 + 99)){
  print(myseed)
  ## <----- randomization 
  set.seed(myseed)
  id = sample(1:nrow(order),nrow(order), replace=T) #######
  gbm = gbm.bak[, order$sample_name[id]]
  ## --------------------->
  res <- sapply(unique(ap), function(s){
    trainData = gbm[, grepl(s, colnames(gbm))]
    trainX = order[grepl(s, order$sample_name[id]),'Pseudotime']
    num.base = 10
    knots = seq(min(order$Pseudotime),max(order$Pseudotime),length.out=num.base+2)[2:(num.base+1)]
    library(splines)
    base = cbind(1,bs(trainX,knots = knots))
    # plot(trainX,base[,1])
    colidx = NULL
    for (ii in seq(2,ncol(base))){
      if (length(unique(base[,ii])) == 1) colidx = c(colidx, ii)
    }
    if (length(colidx)) base = base[,-colidx]
    para = chol2inv(chol(crossprod(base))) %*% t(base) %*% t(trainData) ##
    saveRDS(t(para), paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3_pm/',input.seed,'/spline_coef_', s, '.rds'))
    return(0)
  })
        
  ########## a new file
  meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/doc/Meta_Wenpin_10292019.csv')
  af = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3_pm/',input.seed))
  af = gsub('spline_coef_','',af[grepl('GBM',af)])
  existf = gsub('.rds','',af)
  
  get_coef_list <- function(samples){
    m = list()
    for (j in 1:length(samples)){
      print(j)
      fn = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3_pm/',input.seed,'/spline_coef_',samples[j],'.rds')
  
      coef = readRDS(fn)
      for (i in 1:ncol(coef)){
        if (j == 1) m[[i]] = coef[,i]
        else m[[i]] = cbind(m[[i]], coef[,i])  
      }  

    }
    m  
  }
  
  existf = gsub('.rds','',gsub('spline_coef_','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3_pm/',input.seed))))
  g1 = as.character(unique(meta$Location)[2])
  g2 = as.character(unique(meta$Location)[4])
  ap1 = as.character(meta[as.character(meta$Location)==g1,'Sample.ID'])
  ap1 = intersect(gsub(' ','',ap1), existf)
  ap2 = as.character(meta[as.character(meta$Location)==g2,'Sample.ID'])
  ap2 = intersect(gsub(' ','',ap2), existf)
  ## <----- randomization 
  tmp = c(ap1, ap2)
  set.seed(myseed)
  ap1 = sample(tmp, length(ap1))
  ap2 = setdiff(tmp, ap1)
  ## -------------------->
  m1 = get_coef_list(ap1)
  m2 = get_coef_list(ap2)
  
  t <- sapply(1:length(m1), function(i){
    sapply(rownames(m1[[1]]), function(gene){
      v1 = m1[[i]][gene,]
      v2 = m2[[i]][gene,]
      (mean(v1) - mean(v2))/(sqrt(var(v1)/length(v1) + var(v2)/length(v2)))
    })
  })
  pval <- ((abs(t)>t_true)-0) + pval
  str(pval)
}
saveRDS(pval, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3_pm/',input.seed,'/pvalue.rds'))
  

