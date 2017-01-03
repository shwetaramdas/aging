#Limma using 4 covariates

library(limma)
library("preprocessCore")

FOLDER = "C:/Users/sramdas/Dropbox/"

#read in input data
expn = read.csv("C:/Users/sramdas/Dropbox/TM_Li//original_data/TMtissueAge RMA JunLi.txt",skip=1,header=T,stringsAsFactors = F,sep="\t")
info = read.table(paste(FOLDER,"TM_Li/redo_analysis/aug2016/NEWANALYSIS_Nov22_2016/TM datasheetinfo_tissue_revised.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)

#normalize expression data
info2 = info[order(info$Cel.Filename),]
quantile_norm_expn = normalize.quantiles(as.matrix(expn[,7:ncol(expn)]))
boxplot(quantile_norm_expn,cex=0.25)
expn2 = quantile_norm_expn

samples = names(expn[7:length(names(expn))])
info2 = info[order(info$Cel.Filename),]


##getting PCs
pca = prcomp(t(quantile_norm_expn))

#limma with top 10 PCs and sex,RIN,plate, Race and tissue type as covariates
RESIDUALS = matrix(0,nrow=nrow(expn2),ncol=ncol(expn2))
design2 = model.matrix(~1+ info2$Age + info2$RIN + as.factor(info2$array.plate) + as.factor(info2$Race) + as.factor(info2$Sex) + as.factor(info2$type)+ pca$x[,1]+pca$x[,2]+pca$x[,3]+pca$x[,4]+pca$x[,5]+pca$x[,6]+pca$x[,7]+pca$x[,8]+pca$x[,9]+pca$x[,10])
limma_10pcs_pval = c()
limma_10pcs_coef = c()
for(i in 1:nrow(expn2)){
  residuals = summary(lm(as.numeric(expn2[i,]) ~  pca$x[,1] + pca$x[,2] + pca$x[,3] + pca$x[,4] + pca$x[,5] + pca$x[,6] + pca$x[,7] + pca$x[,8] + pca$x[,9] + pca$x[,10]+ as.factor(info2$Sex) + info2$RIN + as.factor(info2$array.plate) + as.factor(info2$type)) + as.factor(info2$Race))$residuals
  print(length(residuals))
  RESIDUALS[i,] = residuals
  
  fit = lmFit(t(expn2[i,]), design2)
  model = eBayes(fit)
  limma_10pcs_pval = c(limma_10pcs_pval, model$p.value[1,2])
  limma_10pcs_coef = c(limma_10pcs_coef, model$coefficients[1,2])
}
limma_10pcs_qval = p.adjust(limma_10pcs_pval,method="fdr")

#if a gene has multiple probes, get the most significant p value
limma_towrite = cbind(expn[,2], limma_10pcs_pval, limma_10pcs_coef)
limma_towrite[,1] = gsub(" ///.*","",limma_towrite[,1])
limma_towrite = limma_towrite[order(limma_towrite[,1], as.numeric(limma_towrite[,2])),]
for(i in nrow(limma_towrite):2){
  if(limma_towrite[i,1] == limma_towrite[i-1,1]){
    limma_towrite = limma_towrite[-i,]
  }
}
limma_towrite = limma_towrite[order(as.numeric(limma_towrite[,2])),]
write.table(limma_towrite, file="limma_10PCs.txt",row.names=F,col.names=F,quote=F,sep="\t")

plot(limma_10pcs_coef,-log10(limma_10pcs_pval),pch=20)
points(limma_10pcs_coef[which(limma_10pcs_pval < 0.05)],-log10(limma_10pcs_pval[which(limma_10pcs_pval < 0.05)]),pch=20,col="red")
