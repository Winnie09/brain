library(here)
setwd(here())
i <- as.numeric(commandArgs(trailingOnly = T)[[1]])
j <- as.numeric(commandArgs(trailingOnly = T)[[2]])
print(i)
print(j)
print(seq(1,10)[i])
print(seq(1,10)[j])

