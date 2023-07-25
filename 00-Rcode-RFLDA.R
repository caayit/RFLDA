###########################################20191226 IRFMDA ################################################### 
###########################################20191226 load R packages ########################################## 
library(openxlsx)
library(randomForest)
library(ROCR)
library(plyr)
##############################################################################################################



###########################################20191226 data preparation  ######################################### 
### L:lncRNA(240*1)
L <- read.xlsx("./01-lncRNAs-240.xlsx", sheet = 1, colNames = FALSE)
### D:diseases(412*1)
D <- read.xlsx("./02-diseases-412.xlsx", sheet = 1, colNames = FALSE)
### M:miRNA(495*1)
M <- read.xlsx("./03-miRNAs-495.xlsx", sheet = 1, colNames = FALSE)
### LL:lncRNA-lncRNA functional similarities (240*240)
LL <- read.xlsx("./04-lncRNA-lncRNA.xlsx", sheet = 1, colNames = FALSE)
### LD:lncRNA-disease associations (240*412)
LD <- read.xlsx("./05-lncRNA-disease.xlsx", sheet = 1, colNames = FALSE)
### MD:miRNA-disease associations (495*412)
MD <- read.xlsx("./06-miRNA-disease.xlsx", sheet = 1, colNames = FALSE)
### DD:disease-disease semantic similarities (412*412)
DD <- read.xlsx("./07-disease-disease.xlsx", sheet = 1, colNames = FALSE)
### LM:lncRNA-miRNA interactions (240*495)
LM <- read.xlsx("./08-lncRNA-miRNA.xlsx", sheet = 1, colNames = FALSE)

##### consturcting sample dataset 
### Represent lncRNA by 1147 features, L1147: (240*1147)
L1147 <- cbind(LL, LM, LD)
### adding lncRNA name column for L1147 
L1148 <- cbind(L[,1], L1147)
### Represent disease by 1147 features, D1147: (412*1147)
### DL: the transfer matrix of LD
DL <- t(LD)
### DM: the transfer matrix of MD
DM <- t(MD)
### D1147: (412*1147)
D1147 <- cbind(DL, DM, DD)
### adding disease name column for D1147
D1148 <- cbind(D[,1], D1147)
### merge L1148 and D1148 to LDALL (98880*(2+1147+1147=2296)) 
LDALL <- merge(x = L1148, y = D1148, by = NULL)

### Adjust column position of LDALL (98880*(2+1147+1147=2296)) 
d1 <- subset(LDALL,select=1)
d2 <- subset(LDALL,select=1149)
d3 <- subset(LDALL,select=c(-1,-1149))
LDALL <-data.frame(d1,d2,d3)
write.xlsx(LDALL, "E:/RFLDA/data/lncRNA-disease-ALL.xlsx", colNames = FALSE)
LDALL <- read.xlsx("E:/RFLDA/data/lncRNA-disease-ALL.xlsx", sheet = 1, colNames = FALSE)

### exclude columns with full zero values,(98880*(3+1952=1955)) 
d1 <- subset(LDALL,select=c(1,2))
d2 <- subset(LDALL,select=c(-1,-2))
d2 <- d2[,which(colSums(d2) > 0)] 
LDExcl0 <- cbind(d1, d2)
write.xlsx(LDExcl0, "E:/RFLDA/data/lncRNA-disease-Excl0.xlsx", colNames = FALSE)
LDExcl0 <- read.xlsx("E:/RFLDA/data/lncRNA-disease-Excl0.xlsx", sheet = 1, colNames = FALSE)

### construct known lncRNA-disease association pairs
xy <- which(LD[,] == 1, arr.ind = TRUE)
LDA <- cbind(L[xy[,1],1],D[xy[,2],1])
write.xlsx(LDA, "E:/RFLDA/data/lncRNA-disease-associations-2697.xlsx", colNames = FALSE)
LDA <- read.xlsx("E:/RFLDA/data/lncRNA-disease-associations-2697.xlsx", sheet = 1, colNames = FALSE)

### set label for each sample, 1 for known LDA, 0 for unknown LDA (98880*1955)
label.fr <- data.frame("label" = c(0))
LDExcl0 <- cbind(label.fr, LDExcl0)
for (i in 1:nrow(LDExcl0))
{
    for (j in 1:nrow(LDA))
    {
        if((LDExcl0[i,2] == LDA[j,1]) && (LDExcl0[i,3] == LDA[j,2]))
        {
           LDExcl0[i,1] <- 1
        }
    }    
}

### regularize samples using maximum-minimum regularization method (98880*1955)
B1=subset(LDExcl0[,], select=-c(label,X1,X2))
center <- sweep(B1, 2, apply(B1, 2, min),'-') 
R <- apply(B1, 2, max) - apply(B1,2,min)
B2<- sweep(center, 2, R, "/")        
B3 <- subset(LDExcl0[,], select=c(label,X1,X2))
LDExcl0 <- cbind(B3, B2)
write.xlsx(LDExcl0, "E:/RFLDA/data/lncRNA-disease-Excl0-regulation.xlsx", colNames = FALSE)
LDExcl0 <- read.xlsx("E:/RFLDA/data/lncRNA-disease-Excl0-regulation.xlsx", sheet = 1, colNames = FALSE)

### consturct positive samples (2697*1955)
PositiveSample <- LDExcl0[LDExcl0[,1]==1, ]
write.xlsx(PositiveSample, "E:/RFLDA/data/PositiveSample.xlsx", colNames = FALSE)
PositiveSample <- read.xlsx("E:/RFLDA/data/PositiveSample.xlsx", sheet = 1, colNames = FALSE)

### consturct unlabeled samples ((98880-2697=96183)*1955)
UnlabeledSample <- LDExcl0[LDExcl0[,1]==0, ]
write.xlsx(UnlabeledSample, "E:/RFLDA/data/UnlabeledSample.xlsx", colNames = FALSE)
UnlabeledSample <- read.xlsx("E:/RFLDA/data/UnlabeledSample.xlsx", sheet = 1, colNames = FALSE)

### consturct negative samples (2697*1955)
set.seed(1234)
sp <- sample(nrow(UnlabeledSample), nrow(PositiveSample), replace = FALSE, prob = NULL)
sps<-sort(sp)
NegativeSample <- UnlabeledSample[sps,]
write.xlsx(NegativeSample, "E:/RFLDA/data/NegativeSample.xlsx", colNames = FALSE)
NegativeSample <- read.xlsx("E:/RFLDA/data/NegativeSample.xlsx", sheet = 1, colNames = FALSE)

### construct training sample set by combining positive and negative samples (5394*1955)
TrainingSample <- rbind(PositiveSample, NegativeSample)
write.xlsx(TrainingSample, "E:/RFLDA/data/TrainingSample.xlsx", colNames = FALSE)
TrainingSample <- read.xlsx("E:/RFLDA/data/TrainingSample.xlsx", sheet = 1, colNames = FALSE)
##############################################################################################################



########################### compute variable importance score using RandonForest #############################
### read training sample set 
TrainingSample <- read.xlsx("E:/RFLDA/data/TrainingSample.xlsx", sheet = 1, colNames = FALSE)
B1 <- subset(TrainingSample[,], select=-(X2:X3))

### compute variable importance score
ComputeFeatureImportance<-function(){
    ### train RandomForest model with parameter mtry=1952/3
    rf=randomForest(X1~.,data = B1,mtry=651, importance = TRUE, ntree=500, na.action=na.omit)
    ### get variable importance score 
    im <- t(round(importance(rf),3))[2,]
    
    ### repeat 9 times
    for(i in 2:10)
    {
        ### train RandomForest model with parameter mtry=1952/3
        rf <- randomForest(X1~.,data = B1,mtry=651, importance = TRUE, ntree=500, na.action=na.omit)
        ### accumulate variable importance score 
        im <- im + t(round(importance(rf),3))[2,]
    }
    
    ### compute average variable importance score of 10 runnings
    im <- im/10
    
    ### sort features by their variable importance score
    fsort <- sort(im,decreasing=TRUE)
    
    ### store feature names and correspongding variable importance score 
    fs <- data.frame(attr(fsort,"names"),fsort)
    write.xlsx(fs, "E:/RFLDA/data/FeatureScore.xlsx", colNames = FALSE)
}
system.time(ComputeFeatureImportance())
##############################################################################################################



############# training RandomForest model with top 50, 100, ..., 1950 features ###############################
### read training sample set consisting of 5394 lncRNA-disease pairs (5394*(3+1952))
B <- read.xlsx("./TrainingSample.xlsx", sheet = 1, colNames = FALSE)
### read variable importance score of each feature
fs <- read.xlsx("./FeatureScore.xlsx", sheet = 1, colNames = FALSE)

ComputeClassificationAccuracy<-function(){
    ### gaccuracy is used to record classification accuracy on differnt training sample subset
    gaccuracy<-matrix(0,nrow=39,ncol=2)

    tt <- 50
    ### nrow(df)=1952
    while (tt < nrow(fs))   
    {
        ### classification accuracy in each fold in 10-fold crossing validataion 
        laccuracy<-c(0,0,0,0,0,0,0,0,0,0)    
        ### average classification accuracy in 10-fold crossing validataion
        lmeanaccuracy<-0

        ### construct training sample subset consisted X1(label) + top tt feature
        ttt <- fs[1:tt,1]
        B1 <- subset(B[,], select=ttt)
        B3 <- subset(B[,], select=X1)
        B4 <- cbind(B3, B1)   

        ### random resampling in 10-fold crossing validation
        ind<- sample(10, nrow(B4), replace = TRUE, prob=c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  0.1, 0.1))

        for(i in 1:10)
        {
            ### train RandomForest model
            rf=randomForest(X1~.,data = B4[ind != i,],mtry=floor(tt/3), importance = TRUE, ntree=500, na.action=na.omit)              
            ### predict using RadomForest model  
            pred <- predict(rf, B4[ind == i,])
         
            ### judging the category of test samples with threshold = 0.5
            for(j in 1:length(pred)) 
            {   if(pred[j] >= 0.5)
                {
                   pred[j] = 1
                }else
               {
                   pred[j] = 0
                }
            }   
    
            ### construct confusion matrix
            t=table(observed = B4[ind==i, "X1"], predicted = pred)
        
            ### compute classification accuracy
            laccuracy[i]=(t[1,1]+t[2,2])/(t[1,1]+t[1,2]+t[2,1]+t[2,2])
            lmeanaccuracy=lmeanaccuracy+laccuracy[i]/10				
        }
    
        gaccuracy[tt/50,1]=tt
        gaccuracy[tt/50,2]=lmeanaccuracy
    
        tt <- tt+50
        print(tt)
    }
    gaccuracy <- data.frame(gaccuracy)
    write.xlsx(gaccuracy, "E:/RFLDA/data/TrainingSample-gaccuracy.xlsx",colNames = FALSE)
}
system.time(ComputeClassificationAccuracy())
##############################################################################################################



########################## 5-fold crossing validation on training sample set with 300 features ###############
### read variable importance score of each feature
fs <- read.xlsx("E:/RFLDA/data/FeatureScore.xlsx", sheet = 1, colNames = FALSE)

### read training sample set consisting of 5394 lncRNA-disease pairs (5394*(3+1952))
B <- read.xlsx("E:/RFLDA/data/TrainingSample.xlsx", sheet = 1, colNames = FALSE)
### extract subset consisting of top 300 featues
tt <- 300
ttt <- fs[1:tt,1]
B1 <- subset(B[,], select=ttt)
B2 <- subset(B[,], select=X1)
### TB is training sample set without column X2（lncRNA name）and X3（disease name）
TB <- cbind(B2, B1)   

### read unlabeld sample set consisting of 96183 lncRNA-disease pairs ((98880-2697=96183)*(3+1952=1955))
BB <- read.xlsx("E:/RFLDA/data/UnlabeledSample.xlsx", sheet = 1, colNames = FALSE)
### extract subset consisting of top 300 featues
tt <- 300
ttt <- fs[1:tt,1]
B1 <- subset(BB[,], select=ttt)
B2 <- subset(BB[,], select=X1)
### NB is unlabeled sample set without column X2（lncRNA name）and X3（disease name）
NB <- cbind(B2, B1)   

### 5-fold crossing validation
FiveFoldCrossingValidation <- function(){
    sumauc <- 0
    sumap <- 0
    PTB <- TB[1:2697,]
    NTB <- TB[2698:5394,] 	

    for(i in 1:4)
    {
	### training sample set
        PTB1 <- PTB[-(((540*(i-1))+1):(540*i)),]
        NTB1 <- NTB[-(((540*(i-1))+1):(540*i)),]
        TrainB <- rbind(PTB1, NTB1)
                
        ### test sample set
        PTB2 <- PTB[(((540*(i-1))+1):(540*i)),]
        TestB <- rbind(PTB2,NB)
        
        ### train RandomForest Model with parameter，try=the number of features（300）/3
        rf=randomForest(X1~.,data = TrainB, mtry=100, importance = TRUE, ntree=500, na.action=na.omit)
  
        ### predict using RandomForest Model
        pred <- predict(rf, TestB)
        
        pred1 <- prediction(pred, TestB$X1)

        ### computing a simple ROC curve (x-axis: fpr, y-axis: tpr)
        ### roc <- performance(pred1, "tpr", "fpr")
        ### plot(roc, main = "ROC chart")
        
	### compute AUC value    
        auc <- performance(pred1, "auc")@y.values
        print(auc)
        sumauc <- sumauc + as.numeric(auc[[1]])

        ### draw ROC precision/recall curve (x-axis: recall, y-axis: precision)
	### perf1 <- performance(pred1, "prec", "rec")
	### plot(perf1)
        
        ### compute AUPR valute 
        prec <- performance(pred1, "prec")@y.values
        rec <- performance(pred1, "rec")@y.values
        ap <-0
        cur_rec <- rec[[1]][2]
        cur_prec <- prec[[1]][2]     
        for (j in 3:length(rec[[1]])) {
          if(prec[[1]][j] >= cur_prec)
          {
            cur_prec = prec[[1]][j]
          }
          if (abs(cur_rec - rec[[1]][j]) > 0) {
           ap = ap + cur_prec * abs(cur_rec - rec[[1]][j])
          }
          cur_rec = rec[[1]][j]
       }
       print(ap)
       sumap <- sumap + ap
    }
      
      i <- 5
      ### training sample set
      PTB1 <- PTB[-(((540*(i-1))+1):2697),]
      NTB1 <- NTB[-(((540*(i-1))+1):2697),]
      TrainB <- rbind(PTB1, NTB1)
        
      ### test sample set
      PTB2 <- PTB[(((540*(i-1))+1):2697),]
      TestB <- rbind(PTB2,NB)
     
      ### train RandomForest Model with parameter，try=the number of features（300）/3
      rf=randomForest(X1~.,data = TrainB, mtry=100, importance = TRUE, ntree=500, na.action=na.omit)
  
      ### predict using RandomForest Model
      pred <- predict(rf, TestB)
        
      pred1 <- prediction(pred, TestB$X1)

      ### draw ROC curve
      ### roc <- performance(pred1, "tpr", "fpr")
      ### plot(roc, main = "ROC chart")

      ### compute AUC valute    
      auc <- performance(pred1, "auc")@y.values
      print(auc)
      sumauc <- sumauc + as.numeric(auc[[1]])
      sumauc <- sumauc/5
      print(sumauc)

      ### draw ROC precision/recall curve (x-axis: recall, y-axis: precision)
      ### perf1 <- performance(pred1, "prec", "rec")
      ### plot(perf1)
      
      ### compute AUPR valute 
      prec <- performance(pred1, "prec")@y.values
      rec <- performance(pred1, "rec")@y.values
      ap <-0
      cur_rec <- rec[[1]][2]
      cur_prec <- prec[[1]][2]     
      for (j in 3:length(rec[[1]])) {
          if(prec[[1]][j] >= cur_prec)
          {
            cur_prec = prec[[1]][j]
          }
          if (abs(cur_rec - rec[[1]][j]) > 0) {
           ap = ap + cur_prec * abs(cur_rec - rec[[1]][j])
          }
          cur_rec = rec[[1]][j]
     }
     print(ap)
     sumap <- sumap + ap
     sumap <- sumap/5
     print(sumap)
}
system.time(FiveFoldCrossingValidation())
##############################################################################################################



################################### predict all lncRNA-disease samples with 300 fearues#######################
### read variable importance score of each feature
fs <- read.xlsx("E:/RFLDA/data/FeatureScore.xlsx", sheet = 1, colNames = FALSE)

### read training sample set consisting of 5394 lncRNA-disease pairs (5394*(3+1952))
B <- read.xlsx("E:/RFLDA/data/TrainingSample.xlsx", sheet = 1, colNames = FALSE)
### extract subset consisting of top 300 featues
tt <- 300
ttt <- fs[1:tt,1]
B1 <- subset(B[,], select=ttt)
B2 <- subset(B[,], select=X1)
### TB is training sample set without column X2（lncRNA name）and X3（disease name）
TB <- cbind(B2, B1)   

### read unlabeld sample set consisting of 96183 lncRNA-disease pairs ((98880-2697=96183)*(3+1952=1955))
BB <- read.xlsx("E:/RFLDA/data/UnlabeledSample.xlsx", sheet = 1, colNames = FALSE)
### extract subset consisting of top 300 featues
tt <- 300
ttt <- fs[1:tt,1]
B1 <- subset(BB[,], select=ttt)
B2 <- subset(BB[,], select=X1)
### NB is unlabeled sample set without column X2（lncRNA name）and X3（disease name）
NB <- cbind(B2, B1)   

### train RandomForest Model with parameter，try=the number of features（300）/3
rf=randomForest(X1~.,data = TB, mtry=100, importance = TRUE, ntree=500, na.action=na.omit)
  
### predict all unknwon samples using RandomForest Model
pred <- predict(rf, NB)
NB1 <- BB[,1:3]
UnlabeledSampleScore <- cbind(NB1, data.frame(pred))
write.xlsx(UnlabeledSampleScore, "E:/RFLDA/data/UnlabeledSampleScore-300-features.xlsx",colNames = FALSE)
################################################################################################################
