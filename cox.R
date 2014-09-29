#!/usr/bin/Rscript

############################
###Cox Model Significance###
############################

###Zhou Ye###
###08/04/2014###

rm(list=ls())
setwd("/Users/zhouye/Google Drive/ADPS_PPMI/Program/")
suppressMessages(library(survival))

###Process Data###
adni <- read.csv("ADNI_Data.csv", stringsAsFactors=F)
adni <- subset(adni, adni$PHASE==1)
biomarkers <- c("ADAS", "RHIPPO", "MMSE", "RAVLT30", "TMTB", "ATEC", "SVV")
adps <- read.csv("ADPS_bl_7.csv", stringsAsFactors=F)
lr <- read.csv("LR_m24.csv", stringsAsFactors=F)
demg <- read.csv("demog.csv", stringsAsFactors=F)
bl_data <- subset(adni, adni$VISCODE=="bl")

###Process Dignostic###
#ADNI Part
adni_temp <- data.frame(adni$RID, adni$VISCODE, adni$AGE, adni$CURRENTDIAG)
colnames(adni_temp) <- c("RID", "VISCODE", "AGE", "CURRENTDIAG")
bl <- subset(adni_temp, adni_temp$VISCODE=="bl") #Patient at baseline
m06 <- subset(adni_temp, adni_temp$VISCODE=="m06") #Patient at 6 month
m12 <- subset(adni_temp, adni_temp$VISCODE=="m12") #Patient at 12 month
m18 <- subset(adni_temp, adni_temp$VISCODE=="m18") #Patient at 18 month
m24 <- subset(adni_temp, adni_temp$VISCODE=="m24") #Patient at 24 month
m36 <- subset(adni_temp, adni_temp$VISCODE=="m36") #Patient at 36 month
data1 <- data.frame(bl$RID, bl$AGE, bl$CURRENTDIAG, m06$CURRENTDIAG, m12$CURRENTDIAG, m18$CURRENTDIAG, m24$CURRENTDIAG, m36$CURRENTDIAG)
colnames(data1) <- c("RID", "AGE", "Baseline_DIAG", "Month06_DIAG", "Month12_DIAG", "Month18_DIAG", "Month24_DIAG", "Month36_DIAG")
#DEMOG Part
demg_temp <- data.frame(demg$RID, demg$VISCODE, demg$PTGENDER, demg$PTMARRY)
demg_temp <- unique(demg_temp)
colnames(demg_temp) <- c("RID","VISCODE","GENDER","MARRY")
data2 <- subset(demg_temp, demg_temp$VISCODE=="sc" | demg_temp$VISCODE=="f" | demg_temp$VISCODE=="bl")
colnames(data2) <- c("RID", "VISCODE", "GENDER", "MARRY")
data2$VISCODE <- NULL
#Merge
diag <- merge(data1, data2, "RID")

###Construct Survival Table###
survival_table <- function(data) {
    patient <- subset(data, Baseline_DIAG==2 & !is.na(Month06_DIAG) & !is.na(Month12_DIAG) & !is.na(Month18_DIAG) & !is.na(Month24_DIAG) & !is.na(Month36_DIAG) & !is.na(FEATURE) & !is.na(AGE) & GENDER>0 & GENDER!=4 & MARRY>0 & MARRY!=4 & MARRY!=5)
    colnames(patient) <- c("RID", "AGEBL", "DIAGBL", "DIAG06", "DIAG12", "DIAG18", "DIAG24", "DIAG36", "GENDER", "MARRY", "FeatureBL")
    survival_table <- data.frame(matrix(NA,nrow(patient),7))
    colnames(survival_table) <- c("RID","Time","Status","Feature","AGE","GENDER","MARRY")
    survival_table$RID <- patient$RID
    survival_table$Feature <- patient$FeatureBL
    survival_table$AGE <- patient$AGEBL
    survival_table$GENDER <- patient$GENDER
    survival_table$MARRY <- patient$MARRY
    for (i in 1:nrow(patient)) {
        temp <- patient[i,]
        if (temp$DIAG06==3) {
            survival_table[i,]$Time <- 6
            survival_table[i,]$Status <- 1
        }
        else {
            if (temp$DIAG12==3) {
                survival_table[i,]$Time <- 12
                survival_table[i,]$Status <- 1
            }
            else {
                if (temp$DIAG18==3) {
                    survival_table[i,]$Time <- 18
                    survival_table[i,]$Status <- 1
                }
                else {
                    if (temp$DIAG24==3) {
                        survival_table[i,]$Time <- 24
                        survival_table[i,]$Status <- 1
                    }
                    else {
                        if (temp$DIAG36==3) {
                            survival_table[i,]$Time <- 36
                            survival_table[i,]$Status <- 1
                        }
                        else {
                            survival_table[i,]$Time <- 36
                            survival_table[i,]$Status <- 0
                        }
                    }
                }
            }
        }
    }
    return(survival_table)
}

###Result Data Set###
result <- data.frame(matrix(NA, length(biomarkers)+2, 3))
colnames(result) <- c("Biomarker", "P", "LogP")

###Single Biomarkers###
for (i in 1:length(biomarkers)) {
    predictor <- data.frame(RID=bl_data$RID, FEATURE=bl_data[,biomarkers[i]])
    data_set <- merge(diag, predictor, by="RID")
    surv_table <- survival_table(data_set)
    model <- coxph(Surv(Time,Status)~Feature+AGE+GENDER, data=surv_table)
    result[i,"Biomarker"] <- biomarkers[i]
    temp <- summary(model)
    result[i,"P"] <- temp$coefficients["Feature","Pr(>|z|)"]
    result[i,"LogP"] <- log(result[i,"P"])
}

###Logistic Regression###
predictor <- data.frame(RID=lr$RID, FEATURE=lr$bl)
data_set <- merge(diag, predictor, by="RID")
surv_table <- survival_table(data_set)
model <- coxph(Surv(Time,Status)~Feature+AGE+GENDER, data=surv_table)
i <- length(biomarkers)+1
result[i,"Biomarker"] <- "LR"
temp <- summary(model)
result[i,"P"] <- temp$coefficients["Feature","Pr(>|z|)"]
result[i,"LogP"] <- log(result[i,"P"])

###ADPS###
predictor <- data.frame(RID=adps$RID, FEATURE=adps$bl)
data_set <- merge(diag, predictor, by="RID")
surv_table <- survival_table(data_set)
model <- coxph(Surv(Time,Status)~Feature+AGE+GENDER, data=surv_table)
i <- length(biomarkers)+2
result[i,"Biomarker"] <- "ADPS"
temp <- summary(model)
result[i,"P"] <- temp$coefficients["Feature","Pr(>|z|)"]
result[i,"LogP"] <- log(result[i,"P"])

###Output###
write.csv(x=result, file=paste0("cox_siginificance.csv"), row.names=F)




