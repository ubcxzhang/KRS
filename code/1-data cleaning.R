rm(list=ls())

##set your path
setwd("/Users/caoxiaowen/Desktop/KRS-code-data/")

#####data cleaning#####
###read outcomes(drug response) from CCLE and GDSC
yy1 <- read.table("data/CCLE_DrugResponse.csv")
yy2 <- read.csv("data/GDSC_DrugResponse.csv")
rownames(yy2) <- yy2[,1]
yy2 <- yy2[,-1]
names(yy2) <- sub("_IC_50", "", names(yy2))

###read the predictors(gene expression) data
xx1 <- read.table("data/CCLE_ExpData.csv")
xx2 <- read.csv("data/GDSC_ExpData.csv")
names(xx1) <- sub("_expr", "", names(xx1))
rownames(xx2) <- rownames(yy2)

##Prioritize the removal of rows or columns 
##with the most NA values based on the data from outcomes
while(length(which(is.na(yy1))) > 0){
  na.c <- c()
  na.r <- c()
  for (ii in 1:dim(yy1)[2]) {
    na.c <- c(na.c, length(which(is.na(yy1[,ii]))))
  }
  for (ii in 1:dim(yy1)[1]) {
    na.r <- c(na.r, length(which(is.na(yy1[ii,]))))
  }
  if(max(na.c) > max(na.r)){
    yy1 <- yy1[, -which.max(na.c)]
  }else{
    yy1 <- yy1[-which.max(na.r),]
    xx1 <- xx1[-which.max(na.r),]
  }
  
}
while(length(which(is.na(yy2))) > 0){
  na.c <- c()
  na.r <- c()
  for (ii in 1:dim(yy2)[2]) {
    na.c <- c(na.c, length(which(is.na(yy2[,ii]))))
  }
  for (ii in 1:dim(yy2)[1]) {
    na.r <- c(na.r, length(which(is.na(yy2[ii,]))))
  }
  if(max(na.c) > max(na.r)){
    yy2 <- yy2[, -which.max(na.c)]
  }else{
    yy2 <- yy2[-which.max(na.r),]
    xx2 <- xx2[-which.max(na.r),]
  }
}

dim(xx1)
dim(xx2)
dim(yy1)
dim(yy2)


##select similar drugs based on the reference
yy1 <- yy1[,c("X17.AAG" , "AEW541"  ,  "AZD0530",  "AZD6244"  ,"Nutlin.3" ,"PD.0325901","Paclitaxel")]
yy2 <- yy2[,c("X17.AAG" , "GSK.1904529A" ,"Bosutinib" , "RDEA119"   , "Nutlin.3a"    ,"PD.0325901","Epothilone.B")]

##select overlapping genes from CCLE and GDSC
index.gene1 <-  names(xx1) %in% names(xx2)
table(index.gene1)
xxx1 <- xx1[,index.gene1]
index.gene2 <-  names(xx2) %in% names(xx1)
table(index.gene2)
xxx2 <- xx2[,index.gene2]

xx1<-xxx1
xx2.1<-xxx2 # keep the gene orders consistent
xx2 <- xx2.1[, match(names(xx1), names(xx2.1))] # keep the gene orders consistent

rw<-which(yy1==0) %% dim(yy1)[1] # Remove the samples where the drug reaction is meaningless
which(rw==0)
yy1<-yy1[-rw,]
xx1<-xx1[-rw,] # yy2 is ok
save(xx1,xx2,yy1,yy2,file = "result/rda files/drug7_gene1067.rda")

###some checking

# # Assuming xx1 and xx2 are data frames or matrices
# rownames_xx1 <- rownames(xx1)
# rownames_xx2 <- rownames(xx2)
# 
# # Find overlapping row names
# overlap <- intersect(rownames_xx1, rownames_xx2)
# 
# # Check if there are any overlapping row names
# if (length(overlap) > 0) {
#   print("There are overlapping row names.")
#   print(overlap)
# } else {
#   print("There are no overlapping row names.")
# }
# 
# colnames_xx1 <- colnames(xx1)
# colnames_xx2 <- colnames(xx2)
# 
# # Find overlapping col names
# overlap <- intersect(colnames_xx1, colnames_xx2)
# 
# # Check if there are any overlapping col names
# if (length(overlap) > 0) {
#   print("There are overlapping row names.")
# 
#   print(length(overlap))
# } else {
#   print("There are no overlapping row names.")
# }
