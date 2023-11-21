getwd()
setwd("/Users/sundaming/R")
library(ggplot2)
dat <- read.table("Candida albicans.csv",sep=",",header = TRUE)
library(ggpubr)
ggplot(data=dat, aes(x=CA, y=PGD2))+geom_point(color="black")+stat_smooth(method="lm",color="black")+stat_cor(data=dat, method = "spearman")
ggplot(data=dat, aes(x=CA, y=PGD2))+geom_point(color="black")+stat_smooth(method="lm",color="black")

