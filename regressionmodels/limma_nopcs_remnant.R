library(limma)
library("preprocessCore")


#Limma using 4 covariates

#read in expression file and infor file containing covariates
expn = read.csv("../data/TMtissueAge RMA JunLi.txt",skip=1,header=T,stringsAsFactors = F,sep="\t")
info = read.table(paste("../data/","TM datasheetinfo_tissue_revised.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)

#sort the info file by the same order in the expn data frame
info2 = info[order(info$Cel.Filename),]

#the second column of the expression matrix contains gene symbols: get unique gene names
expn[,2] = gsub("///.*","",expn[,2])
genes = unique(expn[,2])
genes = gsub("///.*","",genes)
genes = unique(genes)

#remove first 6 annotation columns from the matrix
expn2 = expn[,-c(1:6)]

#for each gene, expression value = average expression across all probes for that gene
expngene = matrix(0, nrow = length(genes), ncol = 93)
j = 1
for(gene in genes){
  rows = which(expn[,2] == gene)
  if(length(rows) > 1){
    expngene[j,] = as.numeric(apply(expn2[rows,],2,mean))
  }else{
    expngene[j,] = as.numeric(expn2[rows,])
  }
  j = j + 1
}

#normalize data
expngene = normalize.quantiles(expngene)

#keeping only 72 remnant samples
expngene = expngene[,which(info2$type == "remnant")]
info2 = info[order(info$Cel.Filename),]
samples = info2$Cel.Filename[info2$type=="remnant"]
info2 = info2[info2$type=="remnant",]


##with no pcs as covariates
pca = prcomp(t(expngene))
apply(pca$x, 2, cor.test, info2$Age)[1:10]
screeplot(pca,72)

design2 = model.matrix(~1+ info2$Age + info2$RIN + as.factor(info2$array.plate) + as.factor(info2$Race) + as.factor(info2$Sex))
limma_nopcs_pval = c()
limma_nopcs_coef = c()
RESIDUALS = matrix(0,nrow=length(genes),ncol=72)
for(i in 1:nrow(expngene)){
  fit = lmFit(t(expngene[i,]), design2)
  model = eBayes(fit)
  limma_nopcs_pval = c(limma_nopcs_pval, model$p.value[1,2])
  limma_nopcs_coef = c(limma_nopcs_coef, model$coefficients[1,2])
  
  RESIDUALS[i,] = summary(lm(expngene[i,] ~ info2$RIN + as.factor(info2$array.plate) + as.factor(info2$Race) + as.factor(info2$Sex) ))$residuals
}
limma_nopcs_qval = p.adjust(limma_nopcs_pval,method="fdr")

limma_towrite = cbind(genes, limma_nopcs_pval, limma_nopcs_coef)
limma_towrite[,1] = gsub(" ///.*","",limma_towrite[,1])
limma_towrite = limma_towrite[order(limma_towrite[,1], as.numeric(limma_towrite[,2])),]
for(i in nrow(limma_towrite):2){
  if(limma_towrite[i,1] == limma_towrite[i-1,1]){
    limma_towrite = limma_towrite[-i,]
  }
}
limma_towrite = limma_towrite[order(as.numeric(as.character(limma_towrite[,2]))),]
write.table(limma_towrite, file="limma_noPCs_remnant.txt",row.names=F,col.names=F,quote=F,sep="\t")

#qq plot
library("qqman")
png("Limma_nopcs_remnant.png")
qq(as.numeric(limma_5pcs_pval),main="Limma, no pcs, remnant")
dev.off()


#volcano plot
png("volcanoplot_limma_nopcs_remnant.png")
plot(limma_5pcs_coef,-log10(limma_5pcs_pval),pch=20,main="Limma, no PCs, remnant")
points(limma_5pcs_coef[which(limma_5pcs_pval < 0.05)],-log10(limma_5pcs_pval[which(limma_5pcs_pval < 0.05)]),pch=20,col="red")
dev.off()


