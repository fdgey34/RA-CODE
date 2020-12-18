library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
targets<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\combine.csv', sep = ',', row.names = 1, header = TRUE, stringsAsFactors = FALSE)

group<-c()
for (j in 1:ncol(targets)) {
  group[j]<-substring(colnames(targets)[j], 1, 2)
}

dge<-DGEList(counts = targets, group = group)
design <- model.matrix(~0+group)
colnames(design)=levels(factor(group))
rownames(design)=colnames(targets)
v<-voom(dge, design, normalize = 'quantile', plot = TRUE)
data<-as.data.frame(v)
#write.csv(data, file = 'D:\\work\\Rheumatoid arthritis\\data\\combine.csv')

contrast.matrix<-makeContrasts(paste0(unique(group),collapse = "-"), levels = design)
fit1<-lmFit(v, design)
fit2<-contrasts.fit(fit1, contrast.matrix)
fit2<-eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

df_ana<-as.data.frame(nrDEG)
df_ana$threshold[df_ana$P.Value < 0.05 & df_ana$logFC > 0.5]<-'up'
df_ana$threshold[df_ana$P.Value < 0.05 & df_ana$logFC < -0.5]<-'down'
df_ana$threshold[df_ana$P.Value >= 0.05 | (df_ana$logFC <= 0.5 & df_ana$logFC >= -0.5)]<-'non'
ggplot(df_ana, aes(x=logFC, y=-log10(P.Value), color=threshold)) + xlab('Log 2 (fold-change)') + ylab('-Log 10 (P-value)') + geom_point(size = 3.2) + scale_color_manual(values = c('#66FF00', '#666666', '#FF3300')) + expand_limits(x = c(-5, 5)) + geom_hline(yintercept = 1.30103, linetype=1, color = 'black', size = 0.5) + geom_vline(xintercept = c(-0.5,0.5), color = 'black', linetype=1, size = 0.5) + ggtitle('Downregulated                         Upregulated') + theme(plot.title = element_text(hjust = 0.4, face='bold', size = 25, family = 'serif'), axis.text.x = element_text(hjust = 0.4, face = 'bold', size = 20, colour = 'black', family = 'serif'), axis.text.y = element_text(hjust = 0.5, face = 'bold', size = 20, colour = 'black', family = 'serif'), axis.title.x = element_text(face = 'bold', size = 20, family = 'serif'), axis.title.y = element_text(face = 'bold', size=20, family = 'serif'), legend.position = 'none')
#write.csv(df_ana, file = 'D:\\work\\Rheumatoid arthritis\\result\\differential analysis\\differential_analysis.csv')

up_regulated_name<-as.character(rownames(df_ana)[order(df_ana$logFC, decreasing = TRUE)][1:20])
down_regulated_name<-as.character(rownames(df_ana)[order(df_ana$logFC)][1:20])
matrix<-as.matrix(rbind(data[up_regulated_name,], data[down_regulated_name,]))
annotation<-HeatmapAnnotation(df = data.frame(Type = group), 
                              col = list(Type = c('RA' = 'orangered', 'HC' = 'darkblue')), 
                              annotation_name_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily = 'serif'), 
                              annotation_legend_param = list(title_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily = 'serif'), labels_gp = gpar(fontsize = 14, fontface = 'bold', fontfamily = 'serif')))
Heatmap(matrix, 
        heatmap_legend_param = list(legend_height = unit(5, "cm"), title_gp = gpar(fontsize = 14, fontface = 'bold', fontfamily = 'serif'), labels_gp = gpar(fontsize = 15, fontface = 'bold', fontfamily = 'serif')), 
        col = colorRampPalette(c('blue', 'black', 'red'))(100), 
        column_names_gp = gpar(fontsize = 1, fontface = 'bold', fontfamily = 'serif'), 
        row_names_gp = gpar(fontsize = 13, fontface = 'bold', fontfamily = 'serif'), 
        name = ' ', 
        show_column_names = FALSE, 
        top_annotation = annotation)

