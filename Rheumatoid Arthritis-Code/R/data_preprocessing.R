library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(Rtsne)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)

GSE93272<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE93272.csv', row.names = 1, header = TRUE, sep = ',')
GSE45291<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE45291.csv', row.names = 1, header = TRUE, sep = ',')
GSE74143<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE74143.csv', row.names = 1, header = TRUE, sep = ',')
GSE65010<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE65010.csv', row.names = 1, header = TRUE, sep = ',')
GSE15573<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE15573.csv', row.names = 1, header = TRUE, sep = ',')
GSE61635<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE61635.csv', row.names = 1, header = TRUE, sep = ',')
GSE65391<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE65391.csv', row.names = 1, header = TRUE, sep = ',')
GSE138458<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE138458.csv', row.names = 1, header = TRUE, sep = ',')
GSE143272<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE143272.csv', row.names = 1, header = TRUE, sep = ',')
GSE113469<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE113469.csv', row.names = 1, header = TRUE, sep = ',')
GSE50772<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE50772.csv', row.names = 1, header = TRUE, sep = ',')
GSE55457<-read.csv(file = 'D:\\work\\Rheumatoid arthritis\\data\\GSE55457.csv', row.names = 1, header = TRUE, sep = ',')

combine_index<-intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(rownames(GSE93272), rownames(GSE45291)), rownames(GSE74143)), rownames(GSE65010)), rownames(GSE15573)), rownames(GSE61635)), rownames(GSE65391)), rownames(GSE138458)), rownames(GSE143272)), rownames(GSE113469)), rownames(GSE50772)), rownames(GSE55457))

GSE93272<-GSE93272[combine_index, ]
GSE45291<-GSE45291[combine_index, ]
GSE74143<-GSE74143[combine_index, ]
GSE65010<-GSE65010[combine_index, ]
GSE15573<-GSE15573[combine_index, ]
GSE61635<-GSE61635[combine_index, ]
GSE65391<-GSE65391[combine_index, ]
GSE138458<-GSE138458[combine_index, ]
GSE143272<-GSE143272[combine_index, ]
GSE113469<-GSE113469[combine_index, ]
GSE50772<-GSE50772[combine_index, ]
GSE55457<-GSE55457[combine_index, ]

data<-cbind(GSE93272, GSE45291, GSE74143, GSE65010, GSE15573, GSE61635, GSE65391, GSE138458, GSE143272, GSE113469, GSE50772)
labels<-c()
for (i in 1:ncol(data)) {
  if (i <= 697){
    labels[i]<-substring(colnames(data)[i], 4, 5)
  }
  if (i > 697){
    labels[i]<-substring(colnames(data)[i], 4, 6)
  }
}

colors<-rainbow(length(unique(labels)))
names(colors)<-unique(labels)
tsne<- Rtsne(t(data), dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
plot(tsne$Y, main="Raw", col=colors[labels], type = 'p', pch = 19, cex = 1.4, cex.axis = 2, cex.main = 2, font.axis = 2, xlab = '', ylab = '')
legend('bottomleft', title = 'Batch', pch = 16, legend = unique(labels), col = colors, ncol = 3, cex = 1.35, text.font = 2, pt.cex = 2.5)

batch<-c()
status<-c()
for (i in 1:ncol(data)) {
  status[i]<-substring(colnames(data)[i], 1, 2)
  if (i <= 697){
    index<-substring(colnames(data)[i], 4, 5)
    if (index == 'B1'){
      batch[i]<-1
    }
    if (index == 'B2'){
      batch[i]<-2
    }
    if (index == 'B3'){
      batch[i]<-3
    }
    if (index == 'B4'){
      batch[i]<-4
    }
    if (index == 'B5'){
      batch[i]<-5
    }
    if (index == 'B6'){
      batch[i]<-6
    }
    if (index == 'B7'){
      batch[i]<-7
    }
    if (index == 'B8'){
      batch[i]<-8
    }
    if (index == 'B9'){
      batch[i]<-9
    }
  }
  if (i > 697){
    index<-substring(colnames(data)[i], 4, 6)
    if (index == 'B10'){
      batch[i]<-10
    }
    if (index == 'B11'){
      batch[i]<-11
    }
  }
}
status<-as.factor(status)
design<-model.matrix(~status)
data_norm<-as.data.frame(removeBatchEffect(data, batch = batch, design = design))

tsne<-Rtsne(t(data_norm), dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
plot(tsne$Y, main="Norm", col=colors[labels], type = 'p', pch = 19, cex = 1.4, cex.axis = 2, cex.main = 2, font.axis = 2, xlab = '', ylab = '')
legend('bottomleft', title = 'Batch', pch = 16, legend = unique(labels), col = colors, ncol = 3, cex = 1.35, text.font = 2, pt.cex = 2.5)

write.csv(data_norm, file = 'D:\\work\\Rheumatoid arthritis\\data\\combine.csv')
write.csv(GSE55457, file = 'D:\\work\\Rheumatoid arthritis\\data\\test.csv')
