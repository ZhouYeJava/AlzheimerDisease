#!/usr/bin/Rscript

###########################
###Prediction Evaluation###
###########################

###Zhou Ye###
###08/04/2014###

rm(list=ls())
setwd("/Users/zhouye/Google Drive/ADPS_PPMI/Program/")
time <- "m24" #m12 or m24
num <- 1000 #bootstrap times

###Process Data###
adni <- read.csv("ADNI_Data.csv", stringsAsFactors=F)
adni <- subset(adni, adni$PHASE==1)
biomarkers <- c("ADAS", "RHIPPO", "MMSE", "RAVLT30", "TMTB", "ATEC", "SVV")
sign <- c(1, -1, -1, -1, 1, -1, 1)
adps <- read.csv("ADPS_bl_7.csv", stringsAsFactors=F)
lr <- read.csv(paste0("LR_", time, ".csv"), stringsAsFactors=F)
bl_data <- subset(adni, adni$VISCODE=="bl")

###Construct Labels###
rid <- sort(unique(adni$RID))
current <- subset(adni, adni$VISCODE=="bl")$CURRENTDIAG
future <- subset(adni, adni$VISCODE==time)$CURRENTDIAG
label <- data.frame(RID=rid, BL=current, TIME=future)
label <- subset(label, !is.na(label$BL) & !is.na(label$TIME))
label$LABEL <- 0
label[label$TIME==3,]$LABEL <- 1
label <- subset(label, label$BL==2)
label$BL <- NULL
label$TIME <- NULL

###Result Data Set###
result <- data.frame(matrix(NA, length(biomarkers)+2, 13))
colnames(result) <- c("Biomarker", "Mean", "Std", "LB", "UB", "ACC", "PPV", "NPV", "SEN", "SPE", "P", "Convert", "NonConvert")

###Classification###
classify <- function(data) {
    data$PREDICT <- NA
    data <- data[order(data$FEATURE),]
    iterator <- cumsum(1-2*data$LABEL)
    bestInd <- which.max(iterator)
    data[1:bestInd,]$PREDICT <- 0
    if (bestInd<nrow(data)) {
        data[(bestInd+1):nrow(data),]$PREDICT <- 1
    }
    cls <- list(accuracy=0, ppv=0, npv=0, sen=0, spe=0)
    cls$accuracy <- length(which(data$PREDICT==data$LABEL))/nrow(data)
    cls$ppv <- length(which(data$PREDICT==1 & data$LABEL==1))/length(which(data$PREDICT==1))
    cls$npv <- length(which(data$PREDICT==0 & data$LABEL==0))/length(which(data$PREDICT==0))
    cls$sen <- length(which(data$PREDICT==1 & data$LABEL==1))/length(which(data$LABEL==1))
    cls$spe <- length(which(data$PREDICT==0 & data$LABEL==0))/length(which(data$LABEL==0))
    return(cls)
}
permutation <- function(data) {
    labelNew <- sample(data$LABEL, length(data$LABEL), replace=F)
    data$LABEL <- labelNew
    return(data)
}

###Single Biomarkers###
for (i in 1:length(biomarkers)) {
    predictor <- data.frame(RID=bl_data$RID, FEATURE=bl_data[,biomarkers[i]]*sign[i])
    predictor <- subset(predictor, !is.na(predictor$FEATURE))
    data_set <- merge(predictor, label, by="RID")
    result[i,"Biomarker"] <- biomarkers[i]
    #bootstrap accuracy
    accuracy <- vector(length=num)
    for (j in 1:num) {
        index <- sample(1:nrow(data_set), nrow(data_set), replace=T)
        temp <- classify(data_set[index,])
        accuracy[j] <- temp$accuracy
    }
    accuracy <- sort(accuracy)
    result[i,"Mean"] <- mean(accuracy)
    result[i,"Std"] <- sd(accuracy)
    result[i,"LB"] <- accuracy[as.integer(num*0.025)]
    result[i,"UB"] <- accuracy[as.integer(num*0.975)]
    #real accuracy
    real <- classify(data_set)
    result[i,"ACC"] <- real$accuracy
    result[i,"PPV"] <- real$ppv
    result[i,"NPV"] <- real$npv
    result[i,"SEN"] <- real$sen
    result[i,"SPE"] <- real$spe
    #permutation test
    p <- 0
    for (j in 1:num)  {
        permute <- classify(permutation(data_set))
        if (permute$accuracy>real$accuracy) {
            p <- p+1
        }
    }
    result[i,"P"] <- p/num
    #statistics
    result[i,"Convert"] <- sum(data_set$LABEL==1)/nrow(data_set)
    result[i,"NonConvert"] <- sum(data_set$LABEL==0)/nrow(data_set)
}

###Logistic Regression###
predictor <- data.frame(RID=lr$RID, FEATURE=lr$bl)
predictor <- subset(predictor, !is.na(predictor$FEATURE))
data_set <- merge(predictor, label, by="RID")
i <- length(biomarkers)+1
result[i,"Biomarker"] <- "LR"
#bootstrap accuracy
accuracy <- vector(length=num)
for (j in 1:num) {
    index <- sample(1:nrow(data_set), nrow(data_set), replace=T)
    temp <- classify(data_set[index,])
    accuracy[j] <- temp$accuracy
}
accuracy <- sort(accuracy)
result[i,"Mean"] <- mean(accuracy)
result[i,"Std"] <- sd(accuracy)
result[i,"LB"] <- accuracy[as.integer(num*0.025)]
result[i,"UB"] <- accuracy[as.integer(num*0.975)]
#real accuracy
real <- classify(data_set)
result[i,"ACC"] <- real$accuracy
result[i,"PPV"] <- real$ppv
result[i,"NPV"] <- real$npv
result[i,"SEN"] <- real$sen
result[i,"SPE"] <- real$spe
#permutation test
p <- 0
for (j in 1:num)  {
    permute <- classify(permutation(data_set))
    if (permute$accuracy>real$accuracy) {
        p <- p+1
    }
}
result[i,"P"] <- p/num
#statistics
result[i,"Convert"] <- sum(data_set$LABEL==1)/nrow(data_set)
result[i,"NonConvert"] <- sum(data_set$LABEL==0)/nrow(data_set)

###ADPS###
predictor <- data.frame(RID=adps$RID, FEATURE=adps$bl)
predictor <- subset(predictor, !is.na(predictor$FEATURE))
data_set <- merge(predictor, label, by="RID")
i <- length(biomarkers)+2
result[i,"Biomarker"] <- "ADPS"
#bootstrap accuracy
accuracy <- vector(length=num)
for (j in 1:num) {
    index <- sample(1:nrow(data_set), nrow(data_set), replace=T)
    temp <- classify(data_set[index,])
    accuracy[j] <- temp$accuracy
}
accuracy <- sort(accuracy)
result[i,"Mean"] <- mean(accuracy)
result[i,"Std"] <- sd(accuracy)
result[i,"LB"] <- accuracy[as.integer(num*0.025)]
result[i,"UB"] <- accuracy[as.integer(num*0.975)]
#real accuracy
real <- classify(data_set)
result[i,"ACC"] <- real$accuracy
result[i,"PPV"] <- real$ppv
result[i,"NPV"] <- real$npv
result[i,"SEN"] <- real$sen
result[i,"SPE"] <- real$spe
#permutation test
p <- 0
for (j in 1:num)  {
    permute <- classify(permutation(data_set))
    if (permute$accuracy>real$accuracy) {
        p <- p+1
    }
}
result[i,"P"] <- p/num
#statistics
result[i,"Convert"] <- sum(data_set$LABEL==1)/nrow(data_set)
result[i,"NonConvert"] <- sum(data_set$LABEL==0)/nrow(data_set)

###Output###
write.csv(x=result, file=paste0("classification_", time, ".csv"), row.names=F)



