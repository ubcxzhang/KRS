####plot####
rm(list = ls())
dat <- c("CCLE","GDSC")
dd <- 2
setwd("/Users/caoxiaowen/Desktop/KRS-code-data/")
load(paste0("result/rda files/bias.abs.",dat[dd],"train.rda"))

####main plot####
library(reshape2)
library(ggsignif)
library(ggplot2)

melt(mse5)
mse <- melt(mse5)
table(mse$Var2)
mse$Var2 <- as.factor(c(rep("RS",7),rep("KRS",7),rep("L21",7),rep("KMTrace",7),rep("KBMTL",7)))
mse$Var2 <- factor(mse$Var2,levels=c("KBMTL","KMTrace","L21","RS","KRS"))

mseplot <- ggplot(mse,aes(x=Var2,y=value))+geom_boxplot( )+geom_jitter(width = 0.05)+  labs(x="",y="MSE"
)+geom_hline(yintercept = median(mse[,3][which(mse$Var2=="KRS")]),color = "steelblue")+
  theme_bw()+
  theme(panel.grid.major=element_line(size = 0.25, linetype = 'solid',
                                      colour = "light gray"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())+
  theme(legend.text = element_text( size = 16,face = 'bold'),legend.title = element_text( size = 13,face = 'bold'),
        axis.text.x = element_text(size = 12,face = 'bold',color="black"),axis.text.y = element_text(size = 12,face = 'bold',color="black"),
        axis.title.y = element_text(size = 12,face = 'bold'),
  )

png(file=paste0("figure/train_",dat[dd],"_mse_figure.png"),width=600,height=300,family="Times")
mseplot
dev.off()


####supplymentary plot#####

table(bias.abs$Var2)
if(dd==1){
  drugs <- c("X17.AAG" ,  "GSK.1904529A", "Bosutinib"
             ,"RDEA119",      "Nutlin.3a" ,   "PD.0325901"
             , "Epothilone.B") #CCLE
}else{
  
  drugs <- c("X17.AAG" , "AEW541",    "AZD0530",
             "AZD6244",   "Nutlin.3", "PD.0325901", "Paclitaxel" )} #GDSC as train
method <- c("KBMTL" ,"KMTrace",     "L21",      "RS" )
bias.abs$col.p <- rep(0,length(bias.abs$Var1))

for (ii in 1:7) {
  for (jj in 1:4) {
    
    indices <- which(bias.abs$Var2 == drugs[ii] & bias.abs$Var1 == method[jj])
    ttest <- t.test(bias.abs[indices, 3], mu = 0, alternative = "less")
    if(ttest$p.value < 0.05) { pvalue <- "Significant Improvement"
    }else{ pvalue <- "No Significant Improvement"}
    print(paste0("drug: ",drugs[ii]," method: ", method[jj]," p.value: ", pvalue))
    bias.abs$col.p[indices] <- rep(pvalue, length(indices))
  }
}
table(bias.abs$col.p)
head(bias.abs)
jj <- 4 #krs & rs

bias.abs$col.p <- as.factor(bias.abs$col.p)
color_mapping <- c("No Significant Improvement" = "lightblue", "Significant Improvement" = "red")
color_mapping
p <- ggplot(bias.abs[which(bias.abs$Var1==method[jj]),], aes(x = factor(Var2), y = value,fill = col.p)) +
  geom_boxplot( outlier.shape = NA,alpha = 0.5) +ylim(c(-0.4,0.2))+  scale_fill_manual(name="",values = color_mapping)+
  labs(x = "", y =  expression(abs(Bias[italic("KRS")]) - abs(Bias[italic("RS")])), fill = "p.value") +geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle(paste0("The Difference of Absolute Bias between KRS and RS\n (",dat[dd]," as Training Dataset)"))+
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12))
file_name <- paste("supplementary figure/boxplot_train_",dat[dd],"_bias_krs and rs.png", sep = "")
ggsave(file_name, plot = p, width = 10, height = 6)

#####enrichment analysis plot#####
# The enrichment analysis utilized the genes filtered 
# from the previous section and was run using tools outside of R software, 
# therefore only the code for plotting is provided
####go & kegg####
rm(list = ls())
setwd("/Users/caoxiaowen/Desktop/KRS-code-data/")
load(paste0("result/rda files/enrithmentplot.rda"))
ls() #convert the results of kegg & go into rda file
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(grid)
go1$dataset <- rep("GDSC",length(go1$term.description))
go2$dataset <- rep("CCLE",length(go2$term.description))

go <- rbind(go2[1:20,],go1[1:20,])
class(go$dataset)
go$dataset<-as.factor(go$dataset)
class(go$term.description)
class(go$observed.gene.count)
class(go$false.discovery.rate)
summary(go$false.discovery.rate)
termrank <- go$term.description[order(go$false.discovery.rate)]
go$term.description <- factor(go$term.description,levels=unique(termrank))
goplot <- ggplot(go, aes(x =dataset,y = term.description,fill=false.discovery.rate,color=false.discovery.rate))+ geom_point(aes(size=observed.gene.count),alpha=0.8)+
  labs(x="",title="GO Analysis")+
  #+facet_wrap(~dataset)
  theme_bw()+theme(panel.grid.major=element_line(size = 0.25, linetype = 'solid',
                                                 colour = "light gray"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())+#guides(fill=F)+
  theme(legend.text = element_text( size = 13,face = 'bold'),legend.title = element_text( size = 13,face = 'bold'),axis.text.x = element_text(size = 13,face = 'bold',color="black"),axis.text.y = element_text(size = 18,face = 'bold',color="black"),
        axis.title.x = element_blank()#element_text(size = 10,face = 'bold'),
        ,axis.title.y = element_blank(),#+ggtitle("GO Analysis")
        plot.title = element_text(size = 15,color="black",hjust = 0.5,face = 'bold'))

kegg2$dataset <- rep("GDSC",length(kegg2$term.description))
kegg1$dataset <- rep("CCLE",length(kegg1$term.description))

kegg <- rbind(kegg1[1:20,],kegg2[1:20,])
class(kegg$dataset)
kegg$dataset <- as.factor(kegg$dataset)
class(kegg$term.description)
class(kegg$observed.gene.count)
class(kegg$false.discovery.rate)
summary(kegg$false.discovery.rate)
termrank <- kegg$term.description[order(kegg$false.discovery.rate)]
kegg$term.description <- factor(kegg$term.description,levels=unique(termrank))
keggplot < -ggplot(kegg, aes(x =dataset,y = term.description,fill=false.discovery.rate,color=false.discovery.rate))+ geom_point(aes(size=observed.gene.count),alpha=0.8)+
  labs(x="",title="KEGG Analysis")+
  #+facet_wrap(~dataset)
  theme_bw()+theme(panel.grid.major=element_line(size = 0.25, linetype = 'solid',
                                                 colour = "light gray"),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())+#guides(fill=F)+
  theme(legend.text = element_text( size = 13,face = 'bold'),legend.title = element_text( size = 13,face = 'bold'),
        axis.text.x = element_text(size = 13,face = 'bold',color="black"),axis.text.y = element_text(size = 18,face = 'bold',color="black"),
        axis.title.x = element_blank()#element_text(size = 10,face = 'bold'),
        ,axis.title.y = element_blank(),#+ggtitle("GO Analysis")
        plot.title = element_text(size = 15,color="black",hjust = 0.5,face = 'bold'))+guides(size=guide_legend(order=1))



png(file="figure/enrichment analysis.png",width=1000,height=1000)
plot_grid(
  goplot, keggplot,
  labels = c("a","b"),  label_size = 20,ncol = 1,align = "v")
dev.off()

####cmap####
rm(list = ls())
setwd("/Users/caoxiaowen/Desktop/KRS-code-data/")

library(ggplot2)
library(readr)
library(gridExtra)
library(grid)
data <- read.csv("result/csv files/ccle.csv")
data$MOA <- as.factor(data$MOA)
data$Name <- as.factor(data$Name)

A <- ggplot(data,aes(x=Name,y=MOA,colour=Score)) +
  geom_point(size=6) +
  labs(x="Compound",y="MOA",title= "connectivity map for CCLE")+
  theme(plot.title=element_text(size=15,colour="black",face = 'bold', hjust=0.5))+
  theme(axis.title.y=element_text(size=15,colour="black"))+
  theme(axis.title.x=element_text(size=15,colour="black"))+
  theme(legend.text=element_text(size=14,colour="black"),legend.title=element_text(size=12))+ 
  theme(axis.text.x=element_text(size = 18,face = 'bold',color="black",angle=90, hjust=1, vjust=.5))+
  theme(axis.text.y=element_text(size = 18,face = 'bold',color="black"))+
  theme(panel.grid=element_line(color = "grey80",size = 0.5), panel.background = element_blank())

data <- read.csv("result/csv files/gdsc.csv")
data$MOA <- as.factor(data$MOA)
data$Name <- as.factor(data$Name)
B <- ggplot(data,aes(x=Name,y=MOA,colour=Score)) +
  geom_point(size=6) +
  labs(x="Compound",y="MOA",title= "connectivity map for GDSC")+
  theme(plot.title=element_text(size=15,colour="black", hjust=0.5,face = 'bold'))+
  theme(axis.title.y=element_text(size=15,colour="black"))+
  theme(axis.title.x=element_text(size=15,colour="black"))+
  theme(legend.text=element_text(size=14,colour="black"),legend.title=element_text(size=12))+ 
  theme(axis.text.x=element_text(size = 18,face = 'bold',color="black",angle=90, hjust=1, vjust=.5))+
  theme(axis.text.y=element_text(size = 18,face = 'bold',color="black"))+
  theme(panel.grid=element_line(color = "grey80",size = 0.5), panel.background = element_blank())

png(filename = "figure/cmap.png",width = 1000,height = 900)
plot_grid(
  A, B,
  labels = c("a","b"),  label_size = 20,ncol = 1,align = "v"
)
dev.off()

####note: The PPI network of the hub gene was created using the Cytoscape ####
