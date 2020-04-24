setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result')
df = readRDS('./6_10_3/significant_genes.rds')
order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/order.rds')
order = order[grepl('GBM', order$sample_name),]
ap = gsub('_.*','', order$sample_name)
gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/saver.rds')
gbm = gbm[, order$sample_name]
meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/doc/Meta_Wenpin_10292019.csv')

g1 = as.character(unique(meta$Location)[2])
g2 = as.character(unique(meta$Location)[4])
ap1 = as.character(meta[as.character(meta$Location)==g1,'Sample.ID'])
ap1 = gsub(' ','',ap1)
ap2 = as.character(meta[as.character(meta$Location)==g2,'Sample.ID'])
ap2 = gsub(' ','',ap2)

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
get_plotdata <- function(samples, g){
  pd = NULL
  pd <- sapply(samples, function(i){
    if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/spline_coef_',i,'.rds'))){
      print(i)
      tmp = gbm[g, ap == i,drop=F]
      time = order[match(colnames(tmp),order$sample_name),'Pseudotime']
      coef = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/spline_coef_',i,'.rds'))
      pred = get_spline_fit(trainData=tmp,trainX=time,fit.min = min(order$Pseudotime),fit.max=max(order$Pseudotime))      
      data.frame(time = time, expr = tmp[1,], pred = pred[1,], sample = i, stringsAsFactors = F)  
    }
  })
  pd = pd[!sapply(pd, function(i) is.null(nrow(i)))]
  pd = do.call(rbind,pd)
}

for (i in 1:30){
  g = rownames(df)[df[,2]>0][i]
  pd1 = get_plotdata(ap1,g)
  pd2 = get_plotdata(ap2,g)
  library(ggplot2)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/plot/6_10_3/',g,'_1.pdf'),width=12,height=3)
  print(ggplot(data=pd1) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample))
  dev.off()
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/plot/6_10_3/',g,'_2.pdf'),width=12,height=3)
  print(ggplot(data=pd2) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample))
  dev.off()
}


