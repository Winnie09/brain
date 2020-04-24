library(data.table)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/raw/')
for (f in af) {
  if (f %in% af[c(1,2,3)]) {
    d <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/raw/',f),data.table = F)
    d <- as.matrix(d)
    cv <- d[2,-c(1:2)]
    ctv <- d[3,-c(1:2)]
    timev <- d[4,-c(1:2)]
    gn <- d[-c(1:5),1]
    d <- d[-c(1:5),-c(1:2)]
    d <- matrix(as.numeric(d),nrow=length(gn))
    row.names(d) <- gn
    namev <- paste0('2016_La:',sub('\\..*','',f),':cell_',1:ncol(d))
    colnames(d) <- namev
    e <- d
    if (f==af[1]) {
      ac <- sub('^h','',ctv)
      ac[grep('Rgl',ac)] <- 'radial glia-like cells'
      conv <- c('OMTN'= 'oculomotor and trochlear nucleus', 'Sert'='serotonergic', 'NbM'= 'medial neuroblast', 'NbDA'= 'neuroblast dopaminergic', 'DA0'= 'dopaminergic neurons','DA1'= 'dopaminergic neurons','DA2'= 'dopaminergic neurons', 'RN'= 'red nucleus neurons', 'NbGaba'='GABAergic neuroblasts','Gaba'='GABAergic neurons', 'mNbL1-2'= 'lateral neuroblasts', 'NbML1'='mediolateral neuroblasts','NbML2'='mediolateral neuroblasts','NbML3'='mediolateral neuroblasts','NbML4'='mediolateral neuroblasts','NbML5'='mediolateral neuroblasts', 'NProg'= 'neuronal progenitor', 'ProgFPM'= 'progenitor medial floorplate', 'ProgFPL'='progenitor lateral floorplate', 'ProgM'='progenitor midline', 'ProgBP'='progenitor basal plate', 'Rgl1-3'= 'radial glia-like cells', 'Mgl'= 'microglia', 'Endo'= 'endothelial cells', 'Peric'= 'pericytes', 'Epend'= 'ependymal', 'OPC'= 'oligodendrocyte precursor cells','Unk'='unknown')
      ac[ac %in% names(conv)] <- conv[ac][ac %in% names(conv)]  
      timev <- paste0('embryo ', timev)
    } else if (f==af[2]) {
      ac <- sub('^e','',ctv)
      ac[grep('Rgl',ac)] <- 'embryonic radial glia-like cells'
      ac[grep('Nb',ac)] <- 'embryonic neuroblasts'
      ac[grep('Prog',ac)] <- 'embryonic progenitor'
      ac[grep('SC',ac)] <- 'embryonic stem cell'
      timev <- paste0('hESC differentiated ', timev)
    } else if (f==af[3]) {
      ac <- sub('^i','',ctv)
      ac[grep('DA',ac)] <- 'iPSC dopaminergic neurons'
      ac[grep('MN',ac)] <- 'iPSC motor neurons'
      ac[grep('Nb',ac)] <- 'iPSC neuroblasts'
      ac[grep('Prog',ac)] <- 'iPSC progenitor'
      ac[grep('Rgl',ac)] <- 'iPSC radial glia-like cells'
      ac[grep('RN',ac)] <- 'iPSC red nucleus neurons'
      timev <- paste0('iPSC differentiated ', timev)
    }
    m <- data.frame(cell=namev,rawcelltype=ctv,celltype=ac,time=timev,species='human',location='ventral midbrain',stringsAsFactors = F)  
  } else if (f %in% af[4]) {
    d <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/raw/',f),data.table = F)
    d <- as.matrix(d)
    cv <- d[3,-c(1:2)]
    ctv <- d[4,-c(1:2)]
    gn <- d[-c(1:5),1]
    d <- d[-c(1:5),-c(1:2)]
    d <- matrix(as.numeric(d),nrow=length(gn))
    row.names(d) <- gn
    namev <- paste0('2016_La:',sub('\\..*','',f),':cell_',1:ncol(d))
    colnames(d) <- namev
    e <- d
    ctv[grep('VTA',ctv)] <- 'ventral midbrain ventral tegmental area'
    ctv[grep('SNC',ctv)] <- 'ventral midbrain substantia nigra pars compacta'
    m <- data.frame(cell=namev,rawcelltype='dopaminergic neurons',celltype='dopaminergic neurons',time='Adult',species='mouse',location=ctv,stringsAsFactors = F)  
  } else if (f %in% af[5]) {
    d <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/raw/',f),data.table = F)
    d <- as.matrix(d)
    cv <- d[3,-c(1:2)]
    ctv <- d[4,-c(1:2)]
    gn <- d[-c(1:7),1]
    timev <- d[5,-c(1:2)]
    d <- d[-c(1:7),-c(1:2)]
    d <- matrix(as.numeric(d),nrow=length(gn))
    row.names(d) <- gn
    namev <- paste0('2016_La:',sub('\\..*','',f),':cell_',1:ncol(d))
    colnames(d) <- namev
    e <- d
    ac <- sub('^m','',ctv)
    ac[grep('Rgl',ac)] <- 'radial glia-like cells'
    ac[grep('Gaba',ac)] <- 'GABAergic neurons'
    conv <- c('Epen'='ependymal cells','OMTN'= 'oculomotor and trochlear nucleus', 'Sert'='serotonergic', 'NbM'= 'medial neuroblast', 'NbDA'= 'neuroblast dopaminergic', 'DA0'= 'dopaminergic neurons','DA1'= 'dopaminergic neurons','DA2'= 'dopaminergic neurons', 'RN'= 'red nucleus neurons', 'NbGaba'='GABAergic neuroblasts', 'NbL1'= 'lateral neuroblasts', 'NbL2'= 'lateral neuroblasts', 'NbML1'='mediolateral neuroblasts','NbML2'='mediolateral neuroblasts','NbML3'='mediolateral neuroblasts','NbML4'='mediolateral neuroblasts','NbML5'='mediolateral neuroblasts', 'NProg'= 'neuronal progenitor', 'ProgFPM'= 'progenitor medial floorplate', 'ProgFPL'='progenitor lateral floorplate', 'ProgM'='progenitor midline', 'ProgBP'='progenitor basal plate', 'Rgl1-3'= 'radial glia-like cells', 'Mgl'= 'microglia', 'Endo'= 'endothelial cells', 'Peric'= 'pericytes', 'Epend'= 'ependymal', 'OPC'= 'oligodendrocyte precursor cells','Unk'='unknown')
    ac[ac %in% names(conv)] <- conv[ac][ac %in% names(conv)]  
    m <- data.frame(cell=namev,rawcelltype=ctv,celltype=ac,time=paste0('embryo ',timev),species='mouse',location='ventral midbrain',stringsAsFactors = F)  
  }
  #row.names(e) <- sub('_.*','',row.names(e))
  #e <- e[!duplicated(row.names(e)),]
  saveRDS(e,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/expr/',sub('.*_','',sub('\\..*','',f)),'.rds'))
  saveRDS(m,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_La/proc/meta/',sub('.*_','',sub('\\..*','',f)),'.rds'))
}


