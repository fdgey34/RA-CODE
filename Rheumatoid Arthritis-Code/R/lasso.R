library(glmnet)
library(ggplot2)
lasso<-read.csv(file = 'D:/work/Rheumatoid arthritis/data/lasso.csv', row.names = 1)

x<-t(lasso)
y<-c()
for (i in 1:length(colnames(lasso))) {
  if (substr(colnames(lasso)[i], 1, 2) == 'RA'){
    y[i] = 1
  }
  else if(substr(colnames(lasso)[i], 1, 2) == 'HC'){
    y[i] = 0
  }
}

fit<-glmnet(x, y, family = 'binomial')
plot(fit, lwd = 2, xvar = 'lambda')

cv<-cv.glmnet(x, y, family = 'binomial')
plot(cv)

coefficients<-as.matrix(coef(cv, s=cv$lambda.min))
coeffi<-as.data.frame(sort(coefficients[coefficients[, 1]!=0, ][-1], decreasing = TRUE))
colnames(coeffi)<-'weight'
coef<-rownames(coeffi)
weight<-coeffi$weight
ggplot(coeffi, aes(coef, weight)) + coord_flip() + ggtitle('Coefficient Weights') + geom_bar(aes(fill=factor((weight<0)+1)), stat="identity") + theme(plot.title = element_text(hjust = 0.4, face='bold', size = 20, family = 'serif'), axis.text.x = element_text(hjust = 0.4, face = 'bold', size = 20, colour = 'black', family = 'serif'), axis.text.y = element_text(hjust = 0.5, face = 'bold', size = 16, colour = 'black', family = 'serif'), axis.title.x = element_text(face = 'bold', size = 20, family = 'serif'), axis.title.y = element_text(face = 'bold', size=20, family = 'serif'), legend.position = 'none')

train<-lasso[coef, ]
write.csv(train, 'D:/work/Rheumatoid arthritis/data/train.csv')
