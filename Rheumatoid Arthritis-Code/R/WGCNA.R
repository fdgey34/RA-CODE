library(WGCNA)
library(stringr)
library(dplyr)
library(ggplot2)

exprData<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\combine.csv', row.names = 1, sep = ',', header = TRUE)
datExpr0<-data.frame(t(exprData))

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
datExpr<-datExpr0

sampleTree<-hclust(dist(datExpr), method = 'average')
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1, cex.main = 2)

powers<-c(c(1:10), seq(from = 12, to = 20, by = 2))
sft<-pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed', verbose = 5)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", 
     type="n", main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     labels = powers, cex = 0.9, col = 'red')
abline(h=0.85,col="red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = 'Soft Threshold (power)', 
     ylab = 'Mean Connectivity', type = 'n', main = paste('Mean connectivity'))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = 'red')

power<-8
cor <- WGCNA::cor
net<-blockwiseModules(datExpr, power = power, maxBlockSize = ncol(datExpr), 
                      TOMType = 'signed', minModuleSize = 30, 
                      reassignThreshold = 0, mergeCutHeight = 0.25,
                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                      saveTOMs=TRUE, corType = 'pearson', 
                      maxPOutliers = 1, loadTOMs=TRUE,
                      saveTOMFileBase = 'GSEcombineTOM',
                      verbose = 3)
cor<-stats::cor
mergedColors<-labels2colors(net$colors)
moduleColors<-mergedColors
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, cex.colorLabels = 1.1, cex.dendroLabels = 1)

MEs<-net$MEs
MEs_col<-MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = 1)
nSelect = 1000
select<-sample(ncol(datExpr), size = nSelect)
selectTOM<-dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
plotDiss = selectTOM^7
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

TOM = TOMsimilarityFromExpr(datExpr, power = 1)
modules = c("turquoise", "blue")
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste('turquoise', ".edges.txt", sep=""),
                               nodeFile = paste('turquoise', ".nodes.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.4,
                               nodeNames = modProbes,                               
                               nodeAttr = moduleColors[inModule])
edges<-cyt$edgeData
node<-cyt$nodeData
rownames(node)<-node$nodeName
node<-node[,-1]
df_na<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\result\\differential analysis\\differential_analysis.csv', row.names = 1)
up<-rownames(df_na)[df_na$threshold == 'up']
down<-rownames(df_na)[df_na$threshold == 'down']
node$threshold<-'non'
node[intersect(up, rownames(node)),3]<-'up'
node[intersect(down, rownames(node)),3]<-'down'
colnames(node)<-c("altName", "module", "threshold")
write.csv(edges[, 1:3], file = 'D:\\work\\Rheumatoid arthritis\\result\\WGCNA\\edges.csv', row.names = FALSE)
write.csv(node[, 2:3], file = 'D:\\work\\Rheumatoid arthritis\\result\\WGCNA\\node.csv', row.names = rownames(node))

traitData <- read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\combine_group.csv', 
                      sep=',', header=T, row.names=1)
sampleName <- rownames(datExpr)
traitData <- traitData[match(sampleName, rownames(traitData)), ]
modTraitCor <- cor(MEs_col, traitData, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nrow(datExpr))
modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=TRUE)
modTraitCor = modTraitCorP$bicor
modTraitP   = modTraitCorP$p

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab.x = 1.4, 
               cex.lab.y = 0.715, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), nrow(datExpr)))
geneModuleMembershipA = bicorAndPvalue(datExpr, MEs_col, robustY=TRUE)
geneModuleMembership = geneModuleMembershipA$bicor
MMPvalue = geneModuleMembershipA$p

geneTraitCor = as.data.frame(cor(datExpr, traitData, use = "p"))
geneTraitP = as.data.frame(corPvalueStudent(
  as.matrix(geneTraitCor), nrow(datExpr)))
geneTraitCorA = bicorAndPvalue(datExpr, traitData, robustY=TRUE)
geneTraitCor = as.data.frame(geneTraitCorA$bicor)
geneTraitP   = as.data.frame(geneTraitCorA$p)

module = "turquoise"
pheno = "status"
modNames = substring(colnames(MEs_col), 3)
# ??ȡ??ע????
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# ??ȡģ???ڵĻ???
moduleGenes = moduleColors == module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
gene_module<-cbind(colnames(datExpr), moduleColors)
write.csv(gene_module, file = 'D:\\work\\Rheumatoid arthritis\\data\\gene_module.csv')
