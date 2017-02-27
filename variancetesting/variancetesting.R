#Limma using 4 covariates
library(limma)
library("preprocessCore")

expn = read.csv("../data//TMtissueAge RMA JunLi.txt",skip=1,header=T,stringsAsFactors = F,sep="\t")
info = read.table(paste("../data/","TM datasheetinfo_tissue_revised.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)

info2 = info[order(info$Cel.Filename),]

expn[,2] = gsub("///.*","",expn[,2])
genes = unique(expn[,2])
genes = gsub("///.*","",genes)
genes = unique(genes)

expn2 = expn[,-c(1:6)]
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

expngene_unnorm = expngene
pca = prcomp(t(expngene))

#keeping only remnant samples
expngene = expngene[,which(info2$type == "remnant")]
info2 = info[order(info$Cel.Filename),]
samples = info2$Cel.Filename[info2$type=="remnant"]
info2 = info2[info2$type=="remnant",]

expngene = normalize.quantiles(expngene)
expngene_unnorm = expngene_unnorm[,which(info2$type == "remnant")]

#categorizing samples into young and old
young = which(info2$Age <=35)
old = which(info2$Age > 50)

vars_young = c()
vars_old = c()
for(i in 1:nrow(expngene)){
  vars_young = c(vars_young,var(expngene_unnorm[i,young]))
  vars_old = c(vars_old,var(expngene_unnorm[i,old]))
}

#old versus young samples
yg = expngene[,young]
og = expngene[,old]

#variances
#do genes have higher variance in old samples?
vars_young = c()
vars_old = c()
p_val = c()
f_val = c()
DIFFS = c()
for(i in 1:nrow(expngene)){
  vars_young = c(vars_young,var(yg[i,]))
  vars_old = c(vars_old,var(og[i,]))
  p_val = c(p_val, var.test(yg[i,],og[i,])[[3]])
  f_val = c(f_val, var.test(yg[i,],og[i,])[[1]])
}
DIFFS = vars_old - vars_young
genes_withvardiff = cbind(genes, p_val, DIFFS)
genes_withvardiff = genes_withvardiff[order(as.numeric(genes_withvardiff[,2])),]
write.table(genes_withvardiff, file="genes_variancebetweenoldandyoung.txt",row.names=F,col.names=F,quote=F,sep="\t")

#levene's test
vars_young = c()
vars_old = c()
p2_val = c()
f2_val = c()
DIFFS = c()
for(i in 1:nrow(expngene)){
  vars_young = c(vars_young,var(yg[i,]))
  vars_old = c(vars_old,var(og[i,]))
  e = c(yg[i,],og[i,])
  g = c(rep(0,length(yg[i,])), rep(1,length(og[i,]))) 
  f2_val = c(f2_val, leveneTest(e,g)[[2]][1])
  p2_val = c(p2_val, leveneTest(e,g)[[3]][1])
}
DIFFS = vars_old - vars_young
genes_withvardiff = cbind(genes, p_val, DIFFS)
genes_withvardiff = genes_withvardiff[order(as.numeric(genes_withvardiff[,2])),]
write.table(genes_withvardiff, file="genes_variancebetweenoldandyoung_levenetest.txt",row.names=F,col.names=F,quote=F,sep="\t")

#difference in variance, randomly
diffs = c()
for(perm in 1:100){
  print(perm)
  y1 = sample(1:72,22,replace=F)
  o1 = sample(c(1:72)[-y1], 33,replace=F)
  vars_y = c()
  vars_o = c()
  p_val = c()
  for(i in 1:nrow(expngene)){
    e1 = expngene[i,y1]
	e2 = expngene[i,o1]
    vars_y = c(vars_y,var(e1))
    vars_o = c(vars_o,var(e2))

  }
  diffs = c(diffs, mean(vars_y) - mean(vars_o))
}

DIFFS = vars_old - vars_young
genes_with_vardiff = cbind(genes, DIFFS)

#random permutations, levene's test 
diffs = c()
f2stats = c()
p2vals = c()
for(perm in 1:100){
  print(perm)
  y1 = sample(1:72,22,replace=F)
  o1 = sample(c(1:72)[-y1], 33,replace=F)
  vars_y = c()
  vars_o = c()
  p_val = c()
  stat = c()
  for(i in 1:nrow(expngene)){
    e = c(expngene[i,y1],expngene[i,o1])
  	g = c(rep(0,22), rep(1,33)) 
  	stat = c(stat, leveneTest(e,g)[[2]][1])
  	p_val = c(p_val, leveneTest(e,g)[[3]][1])	
	}
  diffs = c(diffs, mean(vars_y) - mean(vars_o))
  f2stats = c(f2stats, mean(stat))
  p2vals = c(p2vals, mean(p2vals))
}

#pathways 
# Use devtools to install the package
install.packages("devtools")
library(devtools)
devtools::install_github("stephenturner/annotables")

library(dplyr)
library(annotables)
grch37 = as.data.frame(grch37)

pathways = read.table("../lrpathresults/limma_nopcs_averagingprobespergene_remnant",header=T,stringsAsFactors=F,sep="\t",comment.char="",quote="")


All_pathways_young_cor = c()
All_pathways_old_cor = c()
variances_old = c()
variances_young = c()
tokeep = which(pathways$P.Value < 1)
lrpath_tokeep = pathways[tokeep,]

for(PATHWAY in 1:nrow(lrpath_tokeep)){
  print(PATHWAY)
  siggenes = strsplit(pathways$SigGenes[PATHWAY], split=", ")
  #  install.packages("devtools")
  
  #get gene symbols from ENTREZ gene ids
  g = c()
  for (gene in siggenes[[1]]){
    
    gene = sub('"',"",gene)
    gene = sub(' ', '',gene)
    
    g = c(g, unlist(grch37[which(grch37[,2]==as.numeric(gene)),3]))
    #   print(grch37[which(grch37[,2]==as.numeric(gene)),3])
  }
  g = unique(g)

  rowstokeep = which(genes %in% g) #which genes in our expression data are in this geneset
  if(length(rowstokeep) == 0){
    next;
  }

  #rowstokeep in young versus old
  young_cor = c()
  old_cor = c()
  variance_gene_young = c()
  variance_gene_old = c()
  if(length(rowstokeep)==1){
    next
  }

  for(i in 1:(length(rowstokeep)-1)){
      variance_gene_young = c(variance_gene_young, var(yg[rowstokeep[i],]))
      variance_gene_old = c(variance_gene_old,var(og[rowstokeep[i],]))
      for(j in (i+1):length(rowstokeep)){
        young_cor = c(young_cor, cor(as.numeric(yg[rowstokeep[i],]), as.numeric(yg[rowstokeep[j],])))
        old_cor = c(old_cor, cor(as.numeric(og[rowstokeep[i],]), as.numeric(og[rowstokeep[j],])))
      }
  }
  variance_gene_young = c(variance_gene_young, var(yg[rowstokeep[length(rowstokeep)],]))
  variance_gene_old = c(variance_gene_old,var(og[rowstokeep[length(rowstokeep)],]))
  
  young_cor = median(young_cor)
  old_cor = median(old_cor)
  All_pathways_old_cor = c(All_pathways_old_cor, old_cor)
  All_pathways_young_cor= c(All_pathways_young_cor, young_cor)  
  
  variances_old = c(variances_old, median(variance_gene_old)  )
  variances_young = c(variances_young,median(variance_gene_young))
} 

