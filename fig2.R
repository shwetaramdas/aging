#creating plot of expression residuals (i.e. expression corrected for covariates) by age

rowstokeep = which(limma_10pcs_pval < 0.00005)
genestoplot = unique(expn[rowstokeep,2])

par(mfrow=c(3,4))
for(gene in genestoplot){
  rowstoplot = which(expn[rowstokeep,2]==gene)
  if(length(rowstoplot)> 1){
    rowstoplot = rowstoplot[1]
  }
  plot(as.numeric(RESIDUALS[rowstokeep[rowstoplot],]),info2$Age,xlab="Residuals",ylab="Age",main=gene,cex.main=0.9)
  abline(lm(info2$Age ~ as.numeric(RESIDUALS[rowstokeep[rowstoplot],])))

  
}
