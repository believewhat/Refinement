library("bnlearn")
library(xlsx)
data=read.dsc("G:/机器学习/代码/sampling/asia.dsc")
sim=rbn(data,100000)
write.csv(sim,file="G:/机器学习/代码/sampling/asia_100000.csv")