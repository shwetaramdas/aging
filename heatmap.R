#Run after running runregression_github.R

#heatmap
rowstokeep = which(limma_10pcs_pval < 0.00005)
hmcol<-brewer.pal(11,"RdBu")
h = heatmap(t(RESIDUALS[rowstokeep,]), col=hmcol,labRow=info2$Age,labCol=expn[rowstokeep,2])

#barplot on right
order = h$rowInd
barplot(info2$Age[order], horiz=T)

#to get figure 1, I paste the two figures obtained above