####gene selection####
rm(list = ls())
setwd("/Users/caoxiaowen/Desktop/KRS-code-data/")
load("result/rda files/drug7_gene1067.rda")
xx1 <- as.matrix(xx1)
xx2 <- as.matrix(xx2)

feature.glm1 <- matrix(NA,dim(xx1)[2],dim(yy1)[2])
feature.glm2 <- matrix(NA,dim(xx1)[2],dim(yy1)[2])
for (ii in 1:dim(yy1)[2]) {
  set.seed(ii)
  glmnet.mol <- try(cv.glmnet(xx1,yy1[,ii],family="gaussian"))
  if("try-error" %in% class(glmnet.mol))
  {
    next
  }else{
    coef.fit0 <- coef(glmnet.mol,s=glmnet.mol$lambda.1se)
    coef.fit0 <- coef.fit0[-1,]
    coef.idex <- which(coef.fit0!=0)
    
    feature.glm1[coef.idex,ii] <- 1
    feature.glm2[coef.idex,ii] <- coef.fit0[coef.idex]
    cat(paste0("-",ii))}
}


gene.id <- matrix(100,dim(xx1)[2],2)
for (ii in 1:dim(feature.glm1)[2]) {
  
  gene.id[which(feature.glm1[,ii]!=100),1] <- 1
  gene.id[which(feature.glm2[,ii]!=100),2] <- feature.glm2[which(feature.glm2[,ii]!=100),ii]
  cat(paste0("-",ii))
}

gene.name <- colnames(xx1)
selectedgenexx1 <- gene.name[which(gene.id[,1]==1)]
upgene1 <- gene.name[which(gene.id[,2]>0 & gene.id[,2]!=100)]
downgene1 <- gene.name[which(gene.id[,2]<0)]


feature.glm1 <- matrix(NA,dim(xx1)[2],dim(yy1)[2])
feature.glm2 <- matrix(NA,dim(xx1)[2],dim(yy1)[2])
for (ii in 1:dim(yy2)[2]) {
  set.seed(ii)
  glmnet.mol <- try(cv.glmnet(xx2,yy2[,ii],family="gaussian"))
  if("try-error" %in% class(glmnet.mol))
  {
    next
  }else{
    coef.fit0 <- coef(glmnet.mol,s=glmnet.mol$lambda.1se)
    coef.fit0 <- coef.fit0[-1,]
    coef.idex <- which(coef.fit0!=0)
    
    feature.glm1[coef.idex,ii] <- 1
    feature.glm2[coef.idex,ii] <- coef.fit0[coef.idex]
    cat(paste0("-",ii))}
}


gene.id <- matrix(100,dim(xx1)[2],2)
for (ii in 1:dim(feature.glm1)[2]) {
  
  gene.id[which(feature.glm1[,ii]!=100),1] <- 1
  gene.id[which(feature.glm2[,ii]!=100),2] <- feature.glm2[which(feature.glm2[,ii]!=100),ii]
  cat(paste0("-",ii))
}

gene.name <- colnames(xx1)
selectedgenexx2 <- gene.name[which(gene.id[,1]==1)]
upgene2 <- gene.name[which(gene.id[,2]>0 & gene.id[,2]!=100)]
downgene2 <- gene.name[which(gene.id[,2]<0)]


save(selectedgenexx1,selectedgenexx2, 
     upgene1,upgene2, downgene1, downgene2, file = "result/rda files/selectedgene.rda")
write.csv(selectedgenexx1,file = "result/csv files/xx1gene(up-down).csv")
write.csv(selectedgenexx2,file = "result/csv files/xx2gene(up-down).csv")

##NOTE: use the xx1gene(up-down).csv and xx2gene(up-down).csv on the string website to get the go, kegg and PPI result.
