meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/doc/Meta_Wenpin_10292019.csv')
af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/')
af = gsub('spline_coef_','',af[grepl('GBM',af)])
existf = gsub('.rds','',af)

get_coef_list <- function(samples){
  m = list()
  for (j in 1:length(samples)){
    fn = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/spline_coef_',samples[j],'.rds')
    if (file.exists(fn)){
        coef = readRDS(fn)
        for (i in 1:ncol(coef)){
          if (j == 1) m[[i]] = coef[,i]
          else m[[i]] = cbind(m[[i]], coef[,i])  
        }  
    }
  }
  m  
}

g1 = as.character(unique(meta$Location)[2])
g2 = as.character(unique(meta$Location)[4])
ap1 = as.character(meta[as.character(meta$Location)==g1,'Sample.ID'])
ap1 = gsub(' ','',ap1)
m1 = get_coef_list(ap1)

ap2 = as.character(meta[as.character(meta$Location)==g2,'Sample.ID'])
ap2 = gsub(' ','',ap2)
m2 = get_coef_list(ap2)

t <- sapply(1:length(m1), function(i){
  sapply(rownames(m1[[1]]), function(gene){
    v1 = m1[[i]][gene,]
    v2 = m2[[i]][gene,]
    (mean(v1) - mean(v2))/(sqrt(var(v1)/length(v1) + var(v2)/length(v2)))
  })
})


saveRDS(t, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/t_statistics.rds')

