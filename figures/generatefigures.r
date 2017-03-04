#this should be run after running the regression code

expn = read.csv(paste("../data/", "TMtissueAge RMA JunLi.txt",sep=""),skip=1,header=T,stringsAsFactors = F,sep="\t")
info = read.table(paste("../data/","TM datasheetinfo_tissue_revised.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)

info2 = info[order(info$Cel.Filename),]
library(limma)
library("preprocessCore")

#supplementary figure 1
plot(info2$Age, info2$RIN,col=as.factor(info2$type),pch=19,main="Plotting sample age versus RIN (red: remnant, black: cornea)",cex.main=0.9)

#heatmap
rowstokeep = which(limma_5pcs_pval < 0.0001)
genestoplot = unique(genes[rowstokeep])
redblue100<-rgb(read.table("redblue100.txt",sep='\t',row.names=1,header=T))
toplot = RESIDUALS[rowstokeep,]

heatmap(t(toplot), col=redblue100,labRow=info2$Age, labCol=genestoplot)
#standardize
library(MASS)
toplot2 = scale(t(toplot))
toplot2 = toplot2[order(info2$Age),]
heatmap(toplot2, col=redblue100,labRow=sort(info2$Age), labCol=genestoplot,Rowv = NA)

#panel on right


#figure 2
#creating plot of expression residuals (i.e. expression corrected for covariates) by age


par(mfrow=c(3,3),mar=c(3.8,3.7,4,2))
for(gene in genestoplot){
  rowstoplot = which(genes[rowstokeep]==gene)
  if(length(rowstoplot)> 1){
    rowstoplot = rowstoplot[1]
  }
  plot(as.numeric(RESIDUALS[rowstokeep[rowstoplot],]),info2$Age,xlab="Residuals",ylab="Age (years)",main=gene,cex.main=0.9,pch=19)
  abline(lm(info2$Age ~ as.numeric(RESIDUALS[rowstokeep[rowstoplot],])))
}

library(ggplot2)
plots = list()
i = 1
for(gene in genestoplot){
  rowstoplot = which(genes[rowstokeep]==gene)
  if(length(rowstoplot)> 1){
    rowstoplot = rowstoplot[1]
  }
   datamat = as.data.frame(cbind(RESIDUALS[rowstokeep[rowstoplot],], info2$Age))
   p[[i]] <- ggplot(datamat, aes(x=V1,y=V2)) + geom_point() +geom_smooth(method=lm,se=FALSE) + scale_x_continuous(name="Residuals") + scale_y_continuous(name="Age(yrs)") + ggtitle(gene) +theme_bw()+ theme(plot.title=element_text(size=9),axis.title=element_text(size=8),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
   i = i +1
}
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],ncol=3)
