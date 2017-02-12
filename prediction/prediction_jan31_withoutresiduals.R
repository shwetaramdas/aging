#reading in residuals of expressions corrected for sex, RIN, type, race, plate and top 10 PCs
RES = read.table("../TMtissueAge RMA JunLi.txt",header=T,stringsAsFactors=F,sep="\t")

info = read.table("../TM datasheetinfo_tissue.txt",stringsAsFactors=F,header=T,sep="\t")
info2 = info[order(info$Cel.Filename),]

#creating objects in which to store predictions/training sets and pvalues 
PREDICTIONS = matrix(-999, nrow=100, ncol=93)
training = matrix(0, nrow=100, ncol=93)
PVALUES = matrix(0, nrow=100, ncol=nrow(RES))

#run 100 permutations
for(perm in 1:100){
  print(perm)

  #select 70 samples for the training set, remainder for the test set
  train = sample(93,70)
  test = c(1:93)[-which(c(1:93) %in% train)]
  
  #list to store the coefficients and pvalues
  coeffs = c()
  pvals = c()
  training[perm,train] = 1
  
  #for each probe, build a model: RESIDUALS ~ Age, get the pvalues and coefficients for the probe
  for(i in 1:nrow(RES)){
    model = lm(as.numeric(RES[i,train]) ~ info2$Age[train] + info2$RIN[train] + info2$Sex[train] + info2$type[train] + info2$array.plate[train])
    pval = summary(model)[[4]][2,4]
    coeff = summary(model)[[4]][2,1]
    pvals = c(pvals, pval)
    coeffs = c(coeffs, coeff)
  }
  
  CUTOFF = 0.0001
  PVALUES[perm,] = pvals
  
  #pick those rows with pvals < CUTOFF
  rowstokeep = which(pvals < CUTOFF)
 
  #skip prediction in this permutation if there are no genes that pass the cutoff
  if(length(rowstokeep)==0){next;}

  #create an object containing age, residuals of significant genes, rin, 10 PCs, sex, plate, type
  data = as.data.frame(cbind(as.numeric(info2$Age), t(RES[rowstokeep,]), info2$RIN,info2$Sex,info2$type, info2$array.plate))
  data[,1] = as.numeric(data[,1])
  for(colnum in 1:(ncol(data)-4)){
	data[,(colnum+1)] = as.numeric(as.character(data[,(colnum+1)]))
  }
  colnames = c(1:(length(rowstokeep)+5))
  colnames = paste("V",colnames,sep="")
  colnames(data) = colnames
  
  #build a model age ~ residuals of significant genes + covariates
  fit2 = lm(V1 ~ . ,data=data[train,])

  #predict on test data
  PREDICTIONS[perm,test] =  predict(fit2, data[test,])  

}

predictions_average = c()
for(i in 1:93){
  predictions_average = c(predictions_average, mean(PREDICTIONS[which(PREDICTIONS[,i]>=0),i]))
}

plot(predictions_average, info2$Age,xlab="Predicted Ages", ylab="Actual Ages",pch=16)
abline(lm(info2$Age ~ predictions_average))
cor.test(predictions_average, info2$Age)
