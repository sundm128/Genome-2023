data<-read.table("pca.txt",header=T,row.names=1)
dim(data)
pca =prcomp(data, scale= TRUE)
summary(pca)
group=factor(c(rep("M",8),rep("MH",8),rep("MC",8),rep("MHC",8)))
library(ggplot2)
pca_reuslt<-as.data.frame(pca$x)
write.csv(pca_reuslt, file = "pca_reuslt")
pca_reuslt<-cbind(pca_reuslt,group)
ggplot(pca_reuslt,aes(x=PC1,y=PC2,color=group))+ geom_point(shape= 16,size=6)+
  stat_ellipse(level = 0.95, show.legend = F) +
  theme_bw() +
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))

  