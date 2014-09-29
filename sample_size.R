#!/usr/bin/Rscript

#############################
###Sample Size Computation###
#############################

###Zhou Ye###
###08/04/2014###

rm(list=ls())
setwd("/Users/zhouye/Google Drive/ADPS_PPMI/Program/")
suppressMessages(library(MBESS))
suppressMessages(library(outliers))

###Process Data###
adni <- read.csv("ADNI_Data.csv", stringsAsFactors=F)
adni <- subset(adni, adni$PHASE==1)
biomarkers <- c("ADAS", "RHIPPO", "MMSE", "RAVLT30", "TMTB", "ATEC", "SVV")
adps <- read.csv("ADPS_m24_7.csv", stringsAsFactors=F)
lr <- read.csv("LR_m24.csv", stringsAsFactors=F)
bl_data <- subset(adni, adni$VISCODE=="bl")
m24_data <- subset(adni, adni$VISCODE=="m24")
diag <- data.frame(RID=bl_data$RID, DIAG=bl_data$CURRENTDIAG)

###Compute Sample Size###
#Model 1
model1 <- function(bl_diag, bl_value, m24_value) {
	data <- data.frame(bl_diag, bl_value, m24_value)
	colnames(data) <- c("diag", "bl", "m24")
	data_mci <- subset(data, data$diag==2) 
	diff <- data_mci$m24-data_mci$bl
	diff <- rm.outlier(diff[complete.cases(diff)])
	n <- length(diff)
	sample_size <- ((0.842+1.96)^2)*2*(sd(diff)^2)/((0.25*mean(diff))^2)
	std_diff <- 0.25*mean(diff)/sd(diff)*sqrt(n/2)
	ci <- ci.smd(smd=std_diff, n.1=n, n.2=n, conf.level=0.95)
	upper <- ci$Upper.Conf.Limit.smd/sqrt(n/2)
	lower <- ci$Lower.Conf.Limit.smd/sqrt(n/2)
	sample_size_upper <- ((0.842+1.96)^2)*2/(lower^2)
	sample_size_lower <- ((0.842+1.96)^2)*2/(upper^2)
	return(c(round(sample_size), round(sample_size_lower), round(sample_size_upper)))
}
#Model 2
model2 <- function(bl_diag, bl_value, m24_value) {
	data <- data.frame(bl_diag, bl_value, m24_value)
	colnames(data) <- c("diag", "bl", "m24")
	data_mci <- subset(data, data$diag==2)
	data_normal <- subset(data, data$diag==1) 
	diff <- data_mci$m24-data_mci$bl
	diff <- rm.outlier(diff[complete.cases(diff)])
	n <- length(diff)
	nature <- data_normal$m24-data_normal$bl
	nature <- rm.outlier(nature[complete.cases(nature)])
	sample_size <- ((0.842+1.96)^2)*2*(sd(diff)^2)/((0.25*(mean(diff)-mean(nature)))^2)
	std_diff <- 0.25*(mean(diff)-mean(nature))/sd(diff)*sqrt(n/2)
	ci <- ci.smd(smd=std_diff, n.1=n, n.2=n, conf.level=0.95)
	upper <- ci$Upper.Conf.Limit.smd/sqrt(n/2)
	lower <- ci$Lower.Conf.Limit.smd/sqrt(n/2)
	sample_size_upper <- ((0.842+1.96)^2)*2/(lower^2)
	sample_size_lower <- ((0.842+1.96)^2)*2/(upper^2)
	return(c(round(sample_size), round(sample_size_lower), round(sample_size_upper)))
}

###Result Data Set###
result <- data.frame(matrix(NA, length(biomarkers)+2, 7))
colnames(result) <- c("Biomarker", "Model1_Size", "Model1_LB", "Model1_UB", "Model2_Size", "Model2_LB", "Model2_UB")

###Single Biomarkers###
for (i in 1:length(biomarkers)) {
	bl_value <- subset(adni, adni$VISCODE=="bl")[,biomarkers[i]]
	m24_value <- subset(adni, adni$VISCODE=="m24")[,biomarkers[i]]
	bl_diag <- subset(adni, adni$VISCODE=="bl")$CURRENTDIAG
    result[i,"Biomarker"] <- biomarkers[i]
    result[i,c("Model1_Size", "Model1_LB", "Model1_UB")] <- model1(bl_diag, bl_value, m24_value)
    result[i,c("Model2_Size", "Model2_LB", "Model2_UB")] <- model2(bl_diag, bl_value, m24_value)
}

###Logistic Regression###
data_set <- merge(lr, diag, by="RID")
bl_value <- data_set$bl
m24_value <- data_set$m24
bl_diag <- data_set$DIAG
i <- length(biomarkers)+1
result[i,"Biomarker"] <- "LR"
result[i,c("Model1_Size", "Model1_LB", "Model1_UB")] <- model1(bl_diag, bl_value, m24_value)
result[i,c("Model2_Size", "Model2_LB", "Model2_UB")] <- model2(bl_diag, bl_value, m24_value)

###ADPS###
data_set <- merge(adps, diag, by="RID")
bl_value <- data_set$bl
m24_value <- data_set$m24
bl_diag <- data_set$DIAG
i <- length(biomarkers)+2
result[i,"Biomarker"] <- "ADPS"
result[i,c("Model1_Size", "Model1_LB", "Model1_UB")] <- model1(bl_diag, bl_value, m24_value)
result[i,c("Model2_Size", "Model2_LB", "Model2_UB")] <- model2(bl_diag, bl_value, m24_value)

###Output###
write.csv(x=result, file=paste0("sample_size_m24.csv"), row.names=F)
