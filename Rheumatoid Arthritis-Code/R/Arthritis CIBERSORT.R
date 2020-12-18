
rm(list = ls()) 
library(readxl)
setwd("D:/")
setwd("CIBERSORT/")

#CIBERSORT.GSE55457.csv was retrieved from https://cibersort.stanford.edu/

cibersort <- as.data.frame(read.csv("CIBERSORT.GSE55457.csv"))
cibersort <- cibersort[,-c((ncol(cibersort)-2):ncol(cibersort))]

cibersort$Input.Sample
#replace dots with dashes
for (i in 1:ncol(cibersort)){
  colnames(cibersort)[i] <- gsub("\\."," ",colnames(cibersort)[i])
} 

#Patients and case grouping,GSM1337304--GSM1337313 control;others are case
control_cibersort <-  cibersort[c(1:10), ]
case_cibersort <-  cibersort[c(11:23), ]

#ranking according to specific cell
case <- case_cibersort[order(case_cibersort$`Macrophages M1`, decreasing = FALSE), ]
control <- control_cibersort[order(control_cibersort$`Macrophages M1`, decreasing = FALSE), ]

library(reshape2)
row.names(case) <- case$`Input Sample`
row.names(control) <- control$`Input Sample`
case$`Input Sample` = NULL
control$`Input Sample` = NULL

CellType <- colnames(cibersort)
CellType <- CellType[-1]

case_data_frame  <- data.frame(t(case),CellType)
control_data_frame  <- data.frame(t(control),CellType)

case_data_frame <- melt(case_data_frame, id= "CellType")
control_data_frame <- melt(control_data_frame, id= "CellType")

#group labeling
group <- c(rep("case", nrow(case_data_frame)))
group <- as.data.frame(group)
gcase <- data.frame(case_data_frame,group)
group <- c(rep("control", nrow(control_data_frame)))
group <- as.data.frame(group)
gcontrol <- data.frame(control_data_frame,group)

final <- rbind (gcase,gcontrol)
names(final)[2]='sample_id'

# stacking barplot
library(ggplot2)

png(filename = "Arthritis Cibersort Stack.png", res = 300, width = 5000, height = 3300)
ggplot(data = final, mapping = aes(x = sample_id,fill=CellType, y = value*100, width = 1)) + 
  geom_col (position='stack') +
  labs(x='sample_id',y='Relative Abundance (%)')+
  facet_grid(.~group,scales ="free_x",space = "free_x")+
  scale_y_continuous(expand=c(0, 0))+
  theme (axis.title=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
  )
dev.off()

#dotplot
caselist <- list()
controllist <- list()

for (i in CellType){
  caselist[i]<- subset(case_cibersort, select=i)
  controllist[i]<- subset(control_cibersort, select=i)
}

final_zp = c()
for (i in CellType) {
  vectorcase <- unlist(caselist[[i]])
  vectorcontrol <- unlist(controllist[[i]])
  wilcox1 <- wilcox.test(vectorcase, vectorcontrol, alternative = "two.sided",exact=FALSE)
  final_zp = rbind(final_zp, c(wilcox1$statistic, wilcox1$p.value))
}
final_zp = as.data.frame(final_zp, row.names = CellType)
colnames(final_zp) = c("Z_score", "p_val")

final_zp$Z_score = as.vector(scale(final_zp$Z_score))
final_zp$p_val <- (p.adjust(as.vector(final_zp$p_val),method = "fdr", n=length(final_zp$p_val)))

#-log10 transformation
final_zp$p_val <- -log10(final_zp$p_val)

library(ggrepel)

png(filename = "Arthritis Cibersort Dot.png", res = 300, width = 4500, height = 3300)
ggplot(final_zp, aes(x=Z_score,y=p_val, fill=CellType)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.05, position = "dodge") + 
  scale_fill_hue()+
  geom_hline(aes(yintercept=-log10(0.05)),colour="#990000", linetype="dashed",size=1)+ 
  geom_vline(aes(xintercept=0), colour="#990000", linetype="dashed",size=1)+  
  geom_text_repel(label=rownames(final_zp),size=3)+   
  labs(x="Z score", y = "-log10(q Value)") + 
  theme_classic()
dev.off()

