#################
###Data Counts###
#################

###Zhou Ye###
###2013 New Year###

rm(list=ls())
setwd("E:/Research_Bruno/Data") #Change The Route Here!!!
Data <- read.csv("ADNI_ALL_DATA.csv")

###Basic Information###
NameRow1 <- c("Normal_ADNI1_GO","MCI_ADNI1_GO","AD_ADNI1_GO")
N1 <- length(NameRow1)
NameRow2 <- c("Normal_ADNI2","MCI_ADNI2","AD_ADNI2")
N2 <- length(NameRow2)
NameCol <- c("CLASS","VISCODE","AGE","EXAMDATE","DOB","CURRENTDIAG","ADAS11","ADAS13","Hippo","ICV","RHippo","MMSE","CDRSB","SPAREAD","ANARTE","ANARTN","RAVLT1","RAVLT2","RAVLT3","RAVLT4","RAVLT5","RAVLT6","RAVLT30","RAVLTB","RAVLTR","BNT","TMTA","TMTB","TMTE","TAU1","TAU2","TAU3","TAU4","TAU5","TAU","ABETA1","ABETA2","ABETA3","ABETA4","ABETA5","ABETA","PTAU1","PTAU2","PTAU3","PTAU4","PTAU5","PTAU","PIBPET","APGEN1","APGEN2","APOE4","EDU","MRHIPPO","MRICV","MRRHippo","MRPT","MROA","SUVR","VV1","ATEC1","VV2","ATEC2","SVV")
n <- length(NameCol)
Time1 <- c("bl","m03","m06","m12","m18","m24","m30","m36","m42","m48","m54","m60","m66","m72","m78")
T1 <- length(Time1)
Time2 <- c("v01","v02","v03","v04","v05","v06","v07","v11","v12","v21","v22","v31","v32","v41","v42","v51","v52")
T2 <- length(Time2)
##########

###Initialize Result###
#Result#
R1 <- matrix(NA,(N1*T1),n)
Result1 <- data.frame(R1)
colnames(Result1) <- NameCol
Result1$CLASS <- rep(NameRow1,T1)
R2 <- matrix(NA,(N2*T2),n)
Result2 <- data.frame(R2)
colnames(Result2) <- NameCol
Result2$CLASS <- rep(NameRow2,T2)
#Visits#
for (i in 1:T1)
{
	Result1[(N1*(i-1)+1):(N1*i),]$VISCODE <- Time1[i]
}
for (i in 1:T2)
{
	Result2[(N2*(i-1)+1):(N2*i),]$VISCODE <- Time2[i]
}
##########

###ADNI1/ADNIGO/ADNI2###
Data1 <- subset(Data,Data$PHASE==1) #ADNI1/ADNIGO
Data2 <- subset(Data,Data$PHASE==2) #ADNI2
##########

###ADNI1/ADNIGO Counts###
for (i in 1:T1)
{
	time <- Time1[i]
	D1 <- subset(Data1,Data1$VISCODE==time & Data1$CURRENTDIAG==1)
	D2 <- subset(Data1,Data1$VISCODE==time & Data1$CURRENTDIAG==2)
	D3 <- subset(Data1,Data1$VISCODE==time & Data1$CURRENTDIAG==3)
	index1 <- which(Result1$VISCODE==time & Result1$CLASS=="Normal_ADNI1_GO")
	index2 <- which(Result1$VISCODE==time & Result1$CLASS=="MCI_ADNI1_GO")
	index3 <- which(Result1$VISCODE==time & Result1$CLASS=="AD_ADNI1_GO")

	age_1 <- length(which(!is.nan(D1$AGE)))
	age_2 <- length(which(!is.nan(D2$AGE)))
	age_3 <- length(which(!is.nan(D3$AGE)))
	Result1[index1,]$AGE <- age_1
	Result1[index2,]$AGE <- age_2
	Result1[index3,]$AGE <- age_3

	currentdiag_1 <- length(which(!is.nan(D1$CURRENTDIAG)))
	currentdiag_2 <- length(which(!is.nan(D2$CURRENTDIAG)))
	currentdiag_3 <- length(which(!is.nan(D3$CURRENTDIAG)))
	Result1[index1,]$CURRENTDIAG <- currentdiag_1
	Result1[index2,]$CURRENTDIAG <- currentdiag_2
	Result1[index3,]$CURRENTDIAG <- currentdiag_3

	adas11_1 <- length(which(!is.nan(D1$ADAS11)))
	adas11_2 <- length(which(!is.nan(D2$ADAS11)))
	adas11_3 <- length(which(!is.nan(D3$ADAS11)))
	Result1[index1,]$ADAS11 <- adas11_1
	Result1[index2,]$ADAS11 <- adas11_2
	Result1[index3,]$ADAS11 <- adas11_3

	adas13_1 <- length(which(!is.nan(D1$ADAS13)))
	adas13_2 <- length(which(!is.nan(D2$ADAS13)))
	adas13_3 <- length(which(!is.nan(D3$ADAS13)))
	Result1[index1,]$ADAS13 <- adas13_1
	Result1[index2,]$ADAS13 <- adas13_2
	Result1[index3,]$ADAS13 <- adas13_3

	icv_1 <- length(which(!is.nan(D1$ICV)))
	icv_2 <- length(which(!is.nan(D2$ICV)))
	icv_3 <- length(which(!is.nan(D3$ICV)))
	Result1[index1,]$ICV <- icv_1
	Result1[index2,]$ICV <- icv_2
	Result1[index3,]$ICV <- icv_3

	hippo_1 <- length(which(!is.nan(D1$Hippo)))
	hippo_2 <- length(which(!is.nan(D2$Hippo)))
	hippo_3 <- length(which(!is.nan(D3$Hippo)))
	Result1[index1,]$Hippo <- hippo_1
	Result1[index2,]$Hippo <- hippo_2
	Result1[index3,]$Hippo <- hippo_3

	relativehippo_1 <- length(which(!is.nan(D1$RHippo)))
	relativehippo_2 <- length(which(!is.nan(D2$RHippo)))
	relativehippo_3 <- length(which(!is.nan(D3$RHippo)))
	Result1[index1,]$RHippo <- relativehippo_1
	Result1[index2,]$RHippo <- relativehippo_2
	Result1[index3,]$RHippo <- relativehippo_3

	mmse_1 <- length(which(!is.nan(D1$MMSE)))
	mmse_2 <- length(which(!is.nan(D2$MMSE)))
	mmse_3 <- length(which(!is.nan(D3$MMSE)))
	Result1[index1,]$MMSE <- mmse_1
	Result1[index2,]$MMSE <- mmse_2
	Result1[index3,]$MMSE <- mmse_3

	cdrsb_1 <- length(which(!is.nan(D1$CDRSB)))
	cdrsb_2 <- length(which(!is.nan(D2$CDRSB)))
	cdrsb_3 <- length(which(!is.nan(D3$CDRSB)))
	Result1[index1,]$CDRSB <- cdrsb_1
	Result1[index2,]$CDRSB <- cdrsb_2
	Result1[index3,]$CDRSB <- cdrsb_3

	sparead_1 <- length(which(!is.nan(D1$SPAREAD)))
	sparead_2 <- length(which(!is.nan(D2$SPAREAD)))
	sparead_3 <- length(which(!is.nan(D3$SPAREAD)))
	Result1[index1,]$SPAREAD <- sparead_1
	Result1[index2,]$SPAREAD <- sparead_2
	Result1[index3,]$SPAREAD <- sparead_3

	anarterr_1 <- length(which(!is.nan(D1$ANARTERR)))
	anarterr_2 <- length(which(!is.nan(D2$ANARTERR)))
	anarterr_3 <- length(which(!is.nan(D3$ANARTERR)))
	Result1[index1,]$ANARTERR <- anarterr_1
	Result1[index2,]$ANARTERR <- anarterr_2
	Result1[index3,]$ANARTERR <- anarterr_3

anartnd_1 <- length(which(!is.nan(D1$ANARTND)))
anartnd_2 <- length(which(!is.nan(D2$ANARTND)))
anartnd_3 <- length(which(!is.nan(D3$ANARTND)))
Result1[index1,]$ANARTND <- anartnd_1
Result1[index2,]$ANARTND <- anartnd_2
Result1[index3,]$ANARTND <- anartnd_3

ravlt_tot1_1 <- length(which(!is.nan(D1$RAVLT_TOT1)))
ravlt_tot1_2 <- length(which(!is.nan(D2$RAVLT_TOT1)))
ravlt_tot1_3 <- length(which(!is.nan(D3$RAVLT_TOT1)))
Result1[index1,]$RAVLT_TOT1 <- ravlt_tot1_1
Result1[index2,]$RAVLT_TOT1 <- ravlt_tot1_2
Result1[index3,]$RAVLT_TOT1 <- ravlt_tot1_3

ravlt_tot2_1 <- length(which(!is.nan(D1$RAVLT_TOT2)))
ravlt_tot2_2 <- length(which(!is.nan(D2$RAVLT_TOT2)))
ravlt_tot2_3 <- length(which(!is.nan(D3$RAVLT_TOT2)))
Result1[index1,]$RAVLT_TOT2 <- ravlt_tot2_1
Result1[index2,]$RAVLT_TOT2 <- ravlt_tot2_2
Result1[index3,]$RAVLT_TOT2 <- ravlt_tot2_3

ravlt_tot3_1 <- length(which(!is.nan(D1$RAVLT_TOT3)))
ravlt_tot3_2 <- length(which(!is.nan(D2$RAVLT_TOT3)))
ravlt_tot3_3 <- length(which(!is.nan(D3$RAVLT_TOT3)))
Result1[index1,]$RAVLT_TOT3 <- ravlt_tot3_1
Result1[index2,]$RAVLT_TOT3 <- ravlt_tot3_2
Result1[index3,]$RAVLT_TOT3 <- ravlt_tot3_3

ravlt_tot4_1 <- length(which(!is.nan(D1$RAVLT_TOT4)))
ravlt_tot4_2 <- length(which(!is.nan(D2$RAVLT_TOT4)))
ravlt_tot4_3 <- length(which(!is.nan(D3$RAVLT_TOT4)))
Result1[index1,]$RAVLT_TOT4 <- ravlt_tot4_1
Result1[index2,]$RAVLT_TOT4 <- ravlt_tot4_2
Result1[index3,]$RAVLT_TOT4 <- ravlt_tot4_3

ravlt_tot5_1 <- length(which(!is.nan(D1$RAVLT_TOT5)))
ravlt_tot5_2 <- length(which(!is.nan(D2$RAVLT_TOT5)))
ravlt_tot5_3 <- length(which(!is.nan(D3$RAVLT_TOT5)))
Result1[index1,]$RAVLT_TOT5 <- ravlt_tot5_1
Result1[index2,]$RAVLT_TOT5 <- ravlt_tot5_2
Result1[index3,]$RAVLT_TOT5 <- ravlt_tot5_3

ravlt_tot6_1 <- length(which(!is.nan(D1$RAVLT_TOT6)))
ravlt_tot6_2 <- length(which(!is.nan(D2$RAVLT_TOT6)))
ravlt_tot6_3 <- length(which(!is.nan(D3$RAVLT_TOT6)))
Result1[index1,]$RAVLT_TOT6 <- ravlt_tot6_1
Result1[index2,]$RAVLT_TOT6 <- ravlt_tot6_2
Result1[index3,]$RAVLT_TOT6 <- ravlt_tot6_3

ravlt_30min_1 <- length(which(!is.nan(D1$RAVLT_30MIN)))
ravlt_30min_2 <- length(which(!is.nan(D2$RAVLT_30MIN)))
ravlt_30min_3 <- length(which(!is.nan(D3$RAVLT_30MIN)))
Result1[index1,]$RAVLT_30MIN <- ravlt_30min_1
Result1[index2,]$RAVLT_30MIN <- ravlt_30min_2
Result1[index3,]$RAVLT_30MIN <- ravlt_30min_3

ravlt_totb_1 <- length(which(!is.nan(D1$RAVLT_TOTB)))
ravlt_totb_2 <- length(which(!is.nan(D2$RAVLT_TOTB)))
ravlt_totb_3 <- length(which(!is.nan(D3$RAVLT_TOTB)))
Result1[index1,]$RAVLT_TOTB <- ravlt_totb_1
Result1[index2,]$RAVLT_TOTB <- ravlt_totb_2
Result1[index3,]$RAVLT_TOTB <- ravlt_totb_3

ravlt_deltot_1 <- length(which(!is.nan(D1$RAVLT_DELTOT)))
ravlt_deltot_2 <- length(which(!is.nan(D2$RAVLT_DELTOT)))
ravlt_deltot_3 <- length(which(!is.nan(D3$RAVLT_DELTOT)))
Result1[index1,]$RAVLT_DELTOT <- ravlt_deltot_1
Result1[index2,]$RAVLT_DELTOT <- ravlt_deltot_2
Result1[index3,]$RAVLT_DELTOT <- ravlt_deltot_3

bnttotal_1 <- length(which(!is.nan(D1$BNTTOTAL)))
bnttotal_2 <- length(which(!is.nan(D2$BNTTOTAL)))
bnttotal_3 <- length(which(!is.nan(D3$BNTTOTAL)))
Result1[index1,]$BNTTOTAL <- bnttotal_1
Result1[index2,]$BNTTOTAL <- bnttotal_2
Result1[index3,]$BNTTOTAL <- bnttotal_3

traascor_1 <- length(which(!is.nan(D1$TRAASCOR)))
traascor_2 <- length(which(!is.nan(D2$TRAASCOR)))
traascor_3 <- length(which(!is.nan(D3$TRAASCOR)))
Result1[index1,]$TRAASCOR <- traascor_1
Result1[index2,]$TRAASCOR <- traascor_2
Result1[index3,]$TRAASCOR <- traascor_3

trabscor_1 <- length(which(!is.nan(D1$TRABSCOR)))
trabscor_2 <- length(which(!is.nan(D2$TRABSCOR)))
trabscor_3 <- length(which(!is.nan(D3$TRABSCOR)))
Result1[index1,]$TRABSCOR <- trabscor_1
Result1[index2,]$TRABSCOR <- trabscor_2
Result1[index3,]$TRABSCOR <- trabscor_3

traberrcom_1 <- length(which(!is.nan(D1$TRABERRCOM)))
traberrcom_2 <- length(which(!is.nan(D2$TRABERRCOM)))
traberrcom_3 <- length(which(!is.nan(D3$TRABERRCOM)))
Result1[index1,]$TRABERRCOM <- traberrcom_1
Result1[index2,]$TRABERRCOM <- traberrcom_2
Result1[index3,]$TRABERRCOM <- traberrcom_3

tau_1 <- length(which(!is.nan(D1$TAU)))
tau_2 <- length(which(!is.nan(D2$TAU)))
tau_3 <- length(which(!is.nan(D3$TAU)))
Result1[index1,]$TAU <- tau_1
Result1[index2,]$TAU <- tau_2
Result1[index3,]$TAU <- tau_3

abeta_1 <- length(which(!is.nan(D1$ABETA)))
abeta_2 <- length(which(!is.nan(D2$ABETA)))
abeta_3 <- length(which(!is.nan(D3$ABETA)))
Result1[index1,]$ABETA <- abeta_1
Result1[index2,]$ABETA <- abeta_2
Result1[index3,]$ABETA <- abeta_3

ptau_1 <- length(which(!is.nan(D1$PTAU)))
ptau_2 <- length(which(!is.nan(D2$PTAU)))
ptau_3 <- length(which(!is.nan(D3$PTAU)))
Result1[index1,]$PTAU <- ptau_1
Result1[index2,]$PTAU <- ptau_2
Result1[index3,]$PTAU <- ptau_3

pibpet_1 <- length(which(!is.nan(D1$PIBPET)))
pibpet_2 <- length(which(!is.nan(D2$PIBPET)))
pibpet_3 <- length(which(!is.nan(D3$PIBPET)))
Result1[index1,]$PIBPET <- pibpet_1
Result1[index2,]$PIBPET <- pibpet_2
Result1[index3,]$PIBPET <- pibpet_3

apgen1_1 <- length(which(!is.nan(D1$APGEN1)))
apgen1_2 <- length(which(!is.nan(D2$APGEN1)))
apgen1_3 <- length(which(!is.nan(D3$APGEN1)))
Result1[index1,]$APGEN1 <- apgen1_1
Result1[index2,]$APGEN1 <- apgen1_2
Result1[index3,]$APGEN1 <- apgen1_3

apgen2_1 <- length(which(!is.nan(D1$APGEN2)))
apgen2_2 <- length(which(!is.nan(D2$APGEN2)))
apgen2_3 <- length(which(!is.nan(D3$APGEN2)))
Result1[index1,]$APGEN2 <- apgen2_1
Result1[index2,]$APGEN2 <- apgen2_2
Result1[index3,]$APGEN2 <- apgen2_3

apoe4_1 <- length(which(!is.nan(D1$APOE4)))
apoe4_2 <- length(which(!is.nan(D2$APOE4)))
apoe4_3 <- length(which(!is.nan(D3$APOE4)))
Result1[index1,]$APOE4 <- apoe4_1
Result1[index2,]$APOE4 <- apoe4_2
Result1[index3,]$APOE4 <- apoe4_3

edu_1 <- length(which(!is.nan(D1$EDU)))
edu_2 <- length(which(!is.nan(D2$EDU)))
edu_3 <- length(which(!is.nan(D3$EDU)))
Result1[index1,]$EDU <- edu_1
Result1[index2,]$EDU <- edu_2
Result1[index3,]$EDU <- edu_3

mrlhippo_1 <- length(which(!is.nan(D1$MRLHIPPO)))
mrlhippo_2 <- length(which(!is.nan(D2$MRLHIPPO)))
mrlhippo_3 <- length(which(!is.nan(D3$MRLHIPPO)))
Result1[index1,]$MRLHIPPO <- mrlhippo_1
Result1[index2,]$MRLHIPPO <- mrlhippo_2
Result1[index3,]$MRLHIPPO <- mrlhippo_3

mrrhippo_1 <- length(which(!is.nan(D1$MRRHIPPO)))
mrrhippo_2 <- length(which(!is.nan(D2$MRRHIPPO)))
mrrhippo_3 <- length(which(!is.nan(D3$MRRHIPPO)))
Result1[index1,]$MRRHIPPO <- mrrhippo_1
Result1[index2,]$MRRHIPPO <- mrrhippo_2
Result1[index3,]$MRRHIPPO <- mrrhippo_3

mricv_1 <- length(which(!is.nan(D1$MRICV)))
mricv_2 <- length(which(!is.nan(D2$MRICV)))
mricv_3 <- length(which(!is.nan(D3$MRICV)))
Result1[index1,]$MRICV <- mricv_1
Result1[index2,]$MRICV <- mricv_2
Result1[index3,]$MRICV <- mricv_3

mrlent_1 <- length(which(!is.nan(D1$MRLENT)))
mrlent_2 <- length(which(!is.nan(D2$MRLENT)))
mrlent_3 <- length(which(!is.nan(D3$MRLENT)))
Result1[index1,]$MRLENT <- mrlent_1
Result1[index2,]$MRLENT <- mrlent_2
Result1[index3,]$MRLENT <- mrlent_3

mrrent_1 <- length(which(!is.nan(D1$MRRENT)))
mrrent_2 <- length(which(!is.nan(D2$MRRENT)))
mrrent_3 <- length(which(!is.nan(D3$MRRENT)))
Result1[index1,]$MRRENT <- mrrent_1
Result1[index2,]$MRRENT <- mrrent_2
Result1[index3,]$MRRENT <- mrrent_3

mrlpc_1 <- length(which(!is.nan(D1$MRLPC)))
mrlpc_2 <- length(which(!is.nan(D2$MRLPC)))
mrlpc_3 <- length(which(!is.nan(D3$MRLPC)))
Result1[index1,]$MRLPC <- mrlpc_1
Result1[index2,]$MRLPC <- mrlpc_2
Result1[index3,]$MRLPC <- mrlpc_3

mrrpc_1 <- length(which(!is.nan(D1$MRRPC)))
mrrpc_2 <- length(which(!is.nan(D2$MRRPC)))
mrrpc_3 <- length(which(!is.nan(D3$MRRPC)))
Result1[index1,]$MRRPC <- mrrpc_1
Result1[index2,]$MRRPC <- mrrpc_2
Result1[index3,]$MRRPC <- mrrpc_3

mrloa_1 <- length(which(!is.nan(D1$MRLOA)))
mrloa_2 <- length(which(!is.nan(D2$MRLOA)))
mrloa_3 <- length(which(!is.nan(D3$MRLOA)))
Result1[index1,]$MRLOA <- mrloa_1
Result1[index2,]$MRLOA <- mrloa_2
Result1[index3,]$MRLOA <- mrloa_3

mrroa_1 <- length(which(!is.nan(D1$MRROA)))
mrroa_2 <- length(which(!is.nan(D2$MRROA)))
mrroa_3 <- length(which(!is.nan(D3$MRROA)))
Result1[index1,]$MRROA <- mrroa_1
Result1[index2,]$MRROA <- mrroa_2
Result1[index3,]$MRROA <- mrroa_3

mrhippo_1 <- length(which(!is.nan(D1$MRHIPPO)))
mrhippo_2 <- length(which(!is.nan(D2$MRHIPPO)))
mrhippo_3 <- length(which(!is.nan(D3$MRHIPPO)))
Result1[index1,]$MRHIPPO <- mrhippo_1
Result1[index2,]$MRHIPPO <- mrhippo_2
Result1[index3,]$MRHIPPO <- mrhippo_3

mrrelativehippo_1 <- length(which(!is.nan(D1$MRRelativeHippo)))
mrrelativehippo_2 <- length(which(!is.nan(D2$MRRelativeHippo)))
mrrelativehippo_3 <- length(which(!is.nan(D3$MRRelativeHippo)))
Result1[index1,]$MRRelativeHippo <- mrrelativehippo_1
Result1[index2,]$MRRelativeHippo <- mrrelativehippo_2
Result1[index3,]$MRRelativeHippo <- mrrelativehippo_3

mrent_1 <- length(which(!is.nan(D1$MRENT)))
mrent_2 <- length(which(!is.nan(D2$MRENT)))
mrent_3 <- length(which(!is.nan(D3$MRENT)))
Result1[index1,]$MRENT <- mrent_1
Result1[index2,]$MRENT <- mrent_2
Result1[index3,]$MRENT <- mrent_3

mrpc_1 <- length(which(!is.nan(D1$MRPC)))
mrpc_2 <- length(which(!is.nan(D2$MRPC)))
mrpc_3 <- length(which(!is.nan(D3$MRPC)))
Result1[index1,]$MRPC <- mrpc_1
Result1[index2,]$MRPC <- mrpc_2
Result1[index3,]$MRPC <- mrpc_3

mroa_1 <- length(which(!is.nan(D1$MROA)))
mroa_2 <- length(which(!is.nan(D2$MROA)))
mroa_3 <- length(which(!is.nan(D3$MROA)))
Result1[index1,]$MROA <- mroa_1
Result1[index2,]$MROA <- mroa_2
Result1[index3,]$MROA <- mroa_3

suvr_1 <- length(which(!is.nan(D1$SUVR)))
suvr_2 <- length(which(!is.nan(D2$SUVR)))
suvr_3 <- length(which(!is.nan(D3$SUVR)))
Result1[index1,]$SUVR <- suvr_1
Result1[index2,]$SUVR <- suvr_2
Result1[index3,]$SUVR <- suvr_3

ventricles1_1 <- length(which(!is.nan(D1$VENTRICLES1)))
ventricles1_2 <- length(which(!is.nan(D2$VENTRICLES1)))
ventricles1_3 <- length(which(!is.nan(D3$VENTRICLES1)))
Result1[index1,]$VENTRICLES1 <- ventricles1_1
Result1[index2,]$VENTRICLES1 <- ventricles1_2
Result1[index3,]$VENTRICLES1 <- ventricles1_3

entorhin_1 <- length(which(!is.nan(D1$ENTORHIN)))
entorhin_2 <- length(which(!is.nan(D2$ENTORHIN)))
entorhin_3 <- length(which(!is.nan(D3$ENTORHIN)))
Result1[index1,]$ENTORHIN <- entorhin_1
Result1[index2,]$ENTORHIN <- entorhin_2
Result1[index3,]$ENTORHIN <- entorhin_3

ventricles2_1 <- length(which(!is.nan(D1$VENTRICLES2)))
ventricles2_2 <- length(which(!is.nan(D2$VENTRICLES2)))
ventricles2_3 <- length(which(!is.nan(D3$VENTRICLES2)))
Result1[index1,]$VENTRICLES2 <- ventricles2_1
Result1[index2,]$VENTRICLES2 <- ventricles2_2
Result1[index3,]$VENTRICLES2 <- ventricles2_3

erc_1 <- length(which(!is.nan(D1$ERC)))
erc_2 <- length(which(!is.nan(D2$ERC)))
erc_3 <- length(which(!is.nan(D3$ERC)))
Result1[index1,]$ERC <- erc_1
Result1[index2,]$ERC <- erc_2
Result1[index3,]$ERC <- erc_3

ventvol_1 <- length(which(!is.nan(D1$VENTVOL)))
ventvol_2 <- length(which(!is.nan(D2$VENTVOL)))
ventvol_3 <- length(which(!is.nan(D3$VENTVOL)))
Result1[index1,]$VENTVOL <- ventvol_1
Result1[index2,]$VENTVOL <- ventvol_2
Result1[index3,]$VENTVOL <- ventvol_3
}
##########

###ADNI2 Counts###
for (i in 1:T2)
{
time <- Time2[i]
D1 <- subset(Data2,Data2$VISCODE==time & Data2$CURRENTDIAG==1)
D2 <- subset(Data2,Data2$VISCODE==time & Data2$CURRENTDIAG==2)
D3 <- subset(Data2,Data2$VISCODE==time & Data2$CURRENTDIAG==3)
index1 <- which(Result2$VISCODE==time & Result2$CLASS=="Normal_ADNI2")
index2 <- which(Result2$VISCODE==time & Result2$CLASS=="MCI_ADNI2")
index3 <- which(Result2$VISCODE==time & Result2$CLASS=="AD_ADNI2")

age_1 <- length(which(!is.nan(D1$AGE)))
age_2 <- length(which(!is.nan(D2$AGE)))
age_3 <- length(which(!is.nan(D3$AGE)))
Result2[index1,]$AGE <- age_1
Result2[index2,]$AGE <- age_2
Result2[index3,]$AGE <- age_3

currentdiag_1 <- length(which(!is.nan(D1$CURRENTDIAG)))
currentdiag_2 <- length(which(!is.nan(D2$CURRENTDIAG)))
currentdiag_3 <- length(which(!is.nan(D3$CURRENTDIAG)))
Result2[index1,]$CURRENTDIAG <- currentdiag_1
Result2[index2,]$CURRENTDIAG <- currentdiag_2
Result2[index3,]$CURRENTDIAG <- currentdiag_3

adas11_1 <- length(which(!is.nan(D1$ADAS11)))
adas11_2 <- length(which(!is.nan(D2$ADAS11)))
adas11_3 <- length(which(!is.nan(D3$ADAS11)))
Result2[index1,]$ADAS11 <- adas11_1
Result2[index2,]$ADAS11 <- adas11_2
Result2[index3,]$ADAS11 <- adas11_3

adas13_1 <- length(which(!is.nan(D1$ADAS13)))
adas13_2 <- length(which(!is.nan(D2$ADAS13)))
adas13_3 <- length(which(!is.nan(D3$ADAS13)))
Result2[index1,]$ADAS13 <- adas13_1
Result2[index2,]$ADAS13 <- adas13_2
Result2[index3,]$ADAS13 <- adas13_3

lhippo_1 <- length(which(!is.nan(D1$LHippo)))
lhippo_2 <- length(which(!is.nan(D2$LHippo)))
lhippo_3 <- length(which(!is.nan(D3$LHippo)))
Result2[index1,]$LHippo <- lhippo_1
Result2[index2,]$LHippo <- lhippo_2
Result2[index3,]$LHippo <- lhippo_3

rhippo_1 <- length(which(!is.nan(D1$RHippo)))
rhippo_2 <- length(which(!is.nan(D2$RHippo)))
rhippo_3 <- length(which(!is.nan(D3$RHippo)))
Result2[index1,]$RHippo <- rhippo_1
Result2[index2,]$RHippo <- rhippo_2
Result2[index3,]$RHippo <- rhippo_3

icv_1 <- length(which(!is.nan(D1$ICV)))
icv_2 <- length(which(!is.nan(D2$ICV)))
icv_3 <- length(which(!is.nan(D3$ICV)))
Result2[index1,]$ICV <- icv_1
Result2[index2,]$ICV <- icv_2
Result2[index3,]$ICV <- icv_3

hippo_1 <- length(which(!is.nan(D1$Hippo)))
hippo_2 <- length(which(!is.nan(D2$Hippo)))
hippo_3 <- length(which(!is.nan(D3$Hippo)))
Result2[index1,]$Hippo <- hippo_1
Result2[index2,]$Hippo <- hippo_2
Result2[index3,]$Hippo <- hippo_3

relativehippo_1 <- length(which(!is.nan(D1$RelativeHippo)))
relativehippo_2 <- length(which(!is.nan(D2$RelativeHippo)))
relativehippo_3 <- length(which(!is.nan(D3$RelativeHippo)))
Result2[index1,]$RelativeHippo <- relativehippo_1
Result2[index2,]$RelativeHippo <- relativehippo_2
Result2[index3,]$RelativeHippo <- relativehippo_3

mmse_1 <- length(which(!is.nan(D1$MMSE)))
mmse_2 <- length(which(!is.nan(D2$MMSE)))
mmse_3 <- length(which(!is.nan(D3$MMSE)))
Result2[index1,]$MMSE <- mmse_1
Result2[index2,]$MMSE <- mmse_2
Result2[index3,]$MMSE <- mmse_3

cdrsb_1 <- length(which(!is.nan(D1$CDRSB)))
cdrsb_2 <- length(which(!is.nan(D2$CDRSB)))
cdrsb_3 <- length(which(!is.nan(D3$CDRSB)))
Result2[index1,]$CDRSB <- cdrsb_1
Result2[index2,]$CDRSB <- cdrsb_2
Result2[index3,]$CDRSB <- cdrsb_3

sparead_1 <- length(which(!is.nan(D1$SPAREAD)))
sparead_2 <- length(which(!is.nan(D2$SPAREAD)))
sparead_3 <- length(which(!is.nan(D3$SPAREAD)))
Result2[index1,]$SPAREAD <- sparead_1
Result2[index2,]$SPAREAD <- sparead_2
Result2[index3,]$SPAREAD <- sparead_3

anarterr_1 <- length(which(!is.nan(D1$ANARTERR)))
anarterr_2 <- length(which(!is.nan(D2$ANARTERR)))
anarterr_3 <- length(which(!is.nan(D3$ANARTERR)))
Result2[index1,]$ANARTERR <- anarterr_1
Result2[index2,]$ANARTERR <- anarterr_2
Result2[index3,]$ANARTERR <- anarterr_3

anartnd_1 <- length(which(!is.nan(D1$ANARTND)))
anartnd_2 <- length(which(!is.nan(D2$ANARTND)))
anartnd_3 <- length(which(!is.nan(D3$ANARTND)))
Result2[index1,]$ANARTND <- anartnd_1
Result2[index2,]$ANARTND <- anartnd_2
Result2[index3,]$ANARTND <- anartnd_3

ravlt_tot1_1 <- length(which(!is.nan(D1$RAVLT_TOT1)))
ravlt_tot1_2 <- length(which(!is.nan(D2$RAVLT_TOT1)))
ravlt_tot1_3 <- length(which(!is.nan(D3$RAVLT_TOT1)))
Result2[index1,]$RAVLT_TOT1 <- ravlt_tot1_1
Result2[index2,]$RAVLT_TOT1 <- ravlt_tot1_2
Result2[index3,]$RAVLT_TOT1 <- ravlt_tot1_3

ravlt_tot2_1 <- length(which(!is.nan(D1$RAVLT_TOT2)))
ravlt_tot2_2 <- length(which(!is.nan(D2$RAVLT_TOT2)))
ravlt_tot2_3 <- length(which(!is.nan(D3$RAVLT_TOT2)))
Result2[index1,]$RAVLT_TOT2 <- ravlt_tot2_1
Result2[index2,]$RAVLT_TOT2 <- ravlt_tot2_2
Result2[index3,]$RAVLT_TOT2 <- ravlt_tot2_3

ravlt_tot3_1 <- length(which(!is.nan(D1$RAVLT_TOT3)))
ravlt_tot3_2 <- length(which(!is.nan(D2$RAVLT_TOT3)))
ravlt_tot3_3 <- length(which(!is.nan(D3$RAVLT_TOT3)))
Result2[index1,]$RAVLT_TOT3 <- ravlt_tot3_1
Result2[index2,]$RAVLT_TOT3 <- ravlt_tot3_2
Result2[index3,]$RAVLT_TOT3 <- ravlt_tot3_3

ravlt_tot4_1 <- length(which(!is.nan(D1$RAVLT_TOT4)))
ravlt_tot4_2 <- length(which(!is.nan(D2$RAVLT_TOT4)))
ravlt_tot4_3 <- length(which(!is.nan(D3$RAVLT_TOT4)))
Result2[index1,]$RAVLT_TOT4 <- ravlt_tot4_1
Result2[index2,]$RAVLT_TOT4 <- ravlt_tot4_2
Result2[index3,]$RAVLT_TOT4 <- ravlt_tot4_3

ravlt_tot5_1 <- length(which(!is.nan(D1$RAVLT_TOT5)))
ravlt_tot5_2 <- length(which(!is.nan(D2$RAVLT_TOT5)))
ravlt_tot5_3 <- length(which(!is.nan(D3$RAVLT_TOT5)))
Result2[index1,]$RAVLT_TOT5 <- ravlt_tot5_1
Result2[index2,]$RAVLT_TOT5 <- ravlt_tot5_2
Result2[index3,]$RAVLT_TOT5 <- ravlt_tot5_3

ravlt_tot6_1 <- length(which(!is.nan(D1$RAVLT_TOT6)))
ravlt_tot6_2 <- length(which(!is.nan(D2$RAVLT_TOT6)))
ravlt_tot6_3 <- length(which(!is.nan(D3$RAVLT_TOT6)))
Result2[index1,]$RAVLT_TOT6 <- ravlt_tot6_1
Result2[index2,]$RAVLT_TOT6 <- ravlt_tot6_2
Result2[index3,]$RAVLT_TOT6 <- ravlt_tot6_3

ravlt_30min_1 <- length(which(!is.nan(D1$RAVLT_30MIN)))
ravlt_30min_2 <- length(which(!is.nan(D2$RAVLT_30MIN)))
ravlt_30min_3 <- length(which(!is.nan(D3$RAVLT_30MIN)))
Result2[index1,]$RAVLT_30MIN <- ravlt_30min_1
Result2[index2,]$RAVLT_30MIN <- ravlt_30min_2
Result2[index3,]$RAVLT_30MIN <- ravlt_30min_3

ravlt_totb_1 <- length(which(!is.nan(D1$RAVLT_TOTB)))
ravlt_totb_2 <- length(which(!is.nan(D2$RAVLT_TOTB)))
ravlt_totb_3 <- length(which(!is.nan(D3$RAVLT_TOTB)))
Result2[index1,]$RAVLT_TOTB <- ravlt_totb_1
Result2[index2,]$RAVLT_TOTB <- ravlt_totb_2
Result2[index3,]$RAVLT_TOTB <- ravlt_totb_3

ravlt_deltot_1 <- length(which(!is.nan(D1$RAVLT_DELTOT)))
ravlt_deltot_2 <- length(which(!is.nan(D2$RAVLT_DELTOT)))
ravlt_deltot_3 <- length(which(!is.nan(D3$RAVLT_DELTOT)))
Result2[index1,]$RAVLT_DELTOT <- ravlt_deltot_1
Result2[index2,]$RAVLT_DELTOT <- ravlt_deltot_2
Result2[index3,]$RAVLT_DELTOT <- ravlt_deltot_3

bnttotal_1 <- length(which(!is.nan(D1$BNTTOTAL)))
bnttotal_2 <- length(which(!is.nan(D2$BNTTOTAL)))
bnttotal_3 <- length(which(!is.nan(D3$BNTTOTAL)))
Result2[index1,]$BNTTOTAL <- bnttotal_1
Result2[index2,]$BNTTOTAL <- bnttotal_2
Result2[index3,]$BNTTOTAL <- bnttotal_3

traascor_1 <- length(which(!is.nan(D1$TRAASCOR)))
traascor_2 <- length(which(!is.nan(D2$TRAASCOR)))
traascor_3 <- length(which(!is.nan(D3$TRAASCOR)))
Result2[index1,]$TRAASCOR <- traascor_1
Result2[index2,]$TRAASCOR <- traascor_2
Result2[index3,]$TRAASCOR <- traascor_3

trabscor_1 <- length(which(!is.nan(D1$TRABSCOR)))
trabscor_2 <- length(which(!is.nan(D2$TRABSCOR)))
trabscor_3 <- length(which(!is.nan(D3$TRABSCOR)))
Result2[index1,]$TRABSCOR <- trabscor_1
Result2[index2,]$TRABSCOR <- trabscor_2
Result2[index3,]$TRABSCOR <- trabscor_3

traberrcom_1 <- length(which(!is.nan(D1$TRABERRCOM)))
traberrcom_2 <- length(which(!is.nan(D2$TRABERRCOM)))
traberrcom_3 <- length(which(!is.nan(D3$TRABERRCOM)))
Result2[index1,]$TRABERRCOM <- traberrcom_1
Result2[index2,]$TRABERRCOM <- traberrcom_2
Result2[index3,]$TRABERRCOM <- traberrcom_3

tau_1 <- length(which(!is.nan(D1$TAU)))
tau_2 <- length(which(!is.nan(D2$TAU)))
tau_3 <- length(which(!is.nan(D3$TAU)))
Result2[index1,]$TAU <- tau_1
Result2[index2,]$TAU <- tau_2
Result2[index3,]$TAU <- tau_3

abeta_1 <- length(which(!is.nan(D1$ABETA)))
abeta_2 <- length(which(!is.nan(D2$ABETA)))
abeta_3 <- length(which(!is.nan(D3$ABETA)))
Result2[index1,]$ABETA <- abeta_1
Result2[index2,]$ABETA <- abeta_2
Result2[index3,]$ABETA <- abeta_3

ptau_1 <- length(which(!is.nan(D1$PTAU)))
ptau_2 <- length(which(!is.nan(D2$PTAU)))
ptau_3 <- length(which(!is.nan(D3$PTAU)))
Result2[index1,]$PTAU <- ptau_1
Result2[index2,]$PTAU <- ptau_2
Result2[index3,]$PTAU <- ptau_3

pibpet_1 <- length(which(!is.nan(D1$PIBPET)))
pibpet_2 <- length(which(!is.nan(D2$PIBPET)))
pibpet_3 <- length(which(!is.nan(D3$PIBPET)))
Result2[index1,]$PIBPET <- pibpet_1
Result2[index2,]$PIBPET <- pibpet_2
Result2[index3,]$PIBPET <- pibpet_3

apgen1_1 <- length(which(!is.nan(D1$APGEN1)))
apgen1_2 <- length(which(!is.nan(D2$APGEN1)))
apgen1_3 <- length(which(!is.nan(D3$APGEN1)))
Result2[index1,]$APGEN1 <- apgen1_1
Result2[index2,]$APGEN1 <- apgen1_2
Result2[index3,]$APGEN1 <- apgen1_3

apgen2_1 <- length(which(!is.nan(D1$APGEN2)))
apgen2_2 <- length(which(!is.nan(D2$APGEN2)))
apgen2_3 <- length(which(!is.nan(D3$APGEN2)))
Result2[index1,]$APGEN2 <- apgen2_1
Result2[index2,]$APGEN2 <- apgen2_2
Result2[index3,]$APGEN2 <- apgen2_3

apoe4_1 <- length(which(!is.nan(D1$APOE4)))
apoe4_2 <- length(which(!is.nan(D2$APOE4)))
apoe4_3 <- length(which(!is.nan(D3$APOE4)))
Result2[index1,]$APOE4 <- apoe4_1
Result2[index2,]$APOE4 <- apoe4_2
Result2[index3,]$APOE4 <- apoe4_3

edu_1 <- length(which(!is.nan(D1$EDU)))
edu_2 <- length(which(!is.nan(D2$EDU)))
edu_3 <- length(which(!is.nan(D3$EDU)))
Result2[index1,]$EDU <- edu_1
Result2[index2,]$EDU <- edu_2
Result2[index3,]$EDU <- edu_3

mrlhippo_1 <- length(which(!is.nan(D1$MRLHIPPO)))
mrlhippo_2 <- length(which(!is.nan(D2$MRLHIPPO)))
mrlhippo_3 <- length(which(!is.nan(D3$MRLHIPPO)))
Result2[index1,]$MRLHIPPO <- mrlhippo_1
Result2[index2,]$MRLHIPPO <- mrlhippo_2
Result2[index3,]$MRLHIPPO <- mrlhippo_3

mrrhippo_1 <- length(which(!is.nan(D1$MRRHIPPO)))
mrrhippo_2 <- length(which(!is.nan(D2$MRRHIPPO)))
mrrhippo_3 <- length(which(!is.nan(D3$MRRHIPPO)))
Result2[index1,]$MRRHIPPO <- mrrhippo_1
Result2[index2,]$MRRHIPPO <- mrrhippo_2
Result2[index3,]$MRRHIPPO <- mrrhippo_3

mricv_1 <- length(which(!is.nan(D1$MRICV)))
mricv_2 <- length(which(!is.nan(D2$MRICV)))
mricv_3 <- length(which(!is.nan(D3$MRICV)))
Result2[index1,]$MRICV <- mricv_1
Result2[index2,]$MRICV <- mricv_2
Result2[index3,]$MRICV <- mricv_3

mrlent_1 <- length(which(!is.nan(D1$MRLENT)))
mrlent_2 <- length(which(!is.nan(D2$MRLENT)))
mrlent_3 <- length(which(!is.nan(D3$MRLENT)))
Result2[index1,]$MRLENT <- mrlent_1
Result2[index2,]$MRLENT <- mrlent_2
Result2[index3,]$MRLENT <- mrlent_3

mrrent_1 <- length(which(!is.nan(D1$MRRENT)))
mrrent_2 <- length(which(!is.nan(D2$MRRENT)))
mrrent_3 <- length(which(!is.nan(D3$MRRENT)))
Result2[index1,]$MRRENT <- mrrent_1
Result2[index2,]$MRRENT <- mrrent_2
Result2[index3,]$MRRENT <- mrrent_3

mrlpc_1 <- length(which(!is.nan(D1$MRLPC)))
mrlpc_2 <- length(which(!is.nan(D2$MRLPC)))
mrlpc_3 <- length(which(!is.nan(D3$MRLPC)))
Result2[index1,]$MRLPC <- mrlpc_1
Result2[index2,]$MRLPC <- mrlpc_2
Result2[index3,]$MRLPC <- mrlpc_3

mrrpc_1 <- length(which(!is.nan(D1$MRRPC)))
mrrpc_2 <- length(which(!is.nan(D2$MRRPC)))
mrrpc_3 <- length(which(!is.nan(D3$MRRPC)))
Result2[index1,]$MRRPC <- mrrpc_1
Result2[index2,]$MRRPC <- mrrpc_2
Result2[index3,]$MRRPC <- mrrpc_3

mrloa_1 <- length(which(!is.nan(D1$MRLOA)))
mrloa_2 <- length(which(!is.nan(D2$MRLOA)))
mrloa_3 <- length(which(!is.nan(D3$MRLOA)))
Result2[index1,]$MRLOA <- mrloa_1
Result2[index2,]$MRLOA <- mrloa_2
Result2[index3,]$MRLOA <- mrloa_3

mrroa_1 <- length(which(!is.nan(D1$MRROA)))
mrroa_2 <- length(which(!is.nan(D2$MRROA)))
mrroa_3 <- length(which(!is.nan(D3$MRROA)))
Result2[index1,]$MRROA <- mrroa_1
Result2[index2,]$MRROA <- mrroa_2
Result2[index3,]$MRROA <- mrroa_3

mrhippo_1 <- length(which(!is.nan(D1$MRHIPPO)))
mrhippo_2 <- length(which(!is.nan(D2$MRHIPPO)))
mrhippo_3 <- length(which(!is.nan(D3$MRHIPPO)))
Result2[index1,]$MRHIPPO <- mrhippo_1
Result2[index2,]$MRHIPPO <- mrhippo_2
Result2[index3,]$MRHIPPO <- mrhippo_3

mrrelativehippo_1 <- length(which(!is.nan(D1$MRRelativeHippo)))
mrrelativehippo_2 <- length(which(!is.nan(D2$MRRelativeHippo)))
mrrelativehippo_3 <- length(which(!is.nan(D3$MRRelativeHippo)))
Result2[index1,]$MRRelativeHippo <- mrrelativehippo_1
Result2[index2,]$MRRelativeHippo <- mrrelativehippo_2
Result2[index3,]$MRRelativeHippo <- mrrelativehippo_3

mrent_1 <- length(which(!is.nan(D1$MRENT)))
mrent_2 <- length(which(!is.nan(D2$MRENT)))
mrent_3 <- length(which(!is.nan(D3$MRENT)))
Result2[index1,]$MRENT <- mrent_1
Result2[index2,]$MRENT <- mrent_2
Result2[index3,]$MRENT <- mrent_3

mrpc_1 <- length(which(!is.nan(D1$MRPC)))
mrpc_2 <- length(which(!is.nan(D2$MRPC)))
mrpc_3 <- length(which(!is.nan(D3$MRPC)))
Result2[index1,]$MRPC <- mrpc_1
Result2[index2,]$MRPC <- mrpc_2
Result2[index3,]$MRPC <- mrpc_3

mroa_1 <- length(which(!is.nan(D1$MROA)))
mroa_2 <- length(which(!is.nan(D2$MROA)))
mroa_3 <- length(which(!is.nan(D3$MROA)))
Result2[index1,]$MROA <- mroa_1
Result2[index2,]$MROA <- mroa_2
Result2[index3,]$MROA <- mroa_3

suvr_1 <- length(which(!is.nan(D1$SUVR)))
suvr_2 <- length(which(!is.nan(D2$SUVR)))
suvr_3 <- length(which(!is.nan(D3$SUVR)))
Result2[index1,]$SUVR <- suvr_1
Result2[index2,]$SUVR <- suvr_2
Result2[index3,]$SUVR <- suvr_3

ventricles1_1 <- length(which(!is.nan(D1$VENTRICLES1)))
ventricles1_2 <- length(which(!is.nan(D2$VENTRICLES1)))
ventricles1_3 <- length(which(!is.nan(D3$VENTRICLES1)))
Result2[index1,]$VENTRICLES1 <- ventricles1_1
Result2[index2,]$VENTRICLES1 <- ventricles1_2
Result2[index3,]$VENTRICLES1 <- ventricles1_3

entorhin_1 <- length(which(!is.nan(D1$ENTORHIN)))
entorhin_2 <- length(which(!is.nan(D2$ENTORHIN)))
entorhin_3 <- length(which(!is.nan(D3$ENTORHIN)))
Result2[index1,]$ENTORHIN <- entorhin_1
Result2[index2,]$ENTORHIN <- entorhin_2
Result2[index3,]$ENTORHIN <- entorhin_3

ventricles2_1 <- length(which(!is.nan(D1$VENTRICLES2)))
ventricles2_2 <- length(which(!is.nan(D2$VENTRICLES2)))
ventricles2_3 <- length(which(!is.nan(D3$VENTRICLES2)))
Result2[index1,]$VENTRICLES2 <- ventricles2_1
Result2[index2,]$VENTRICLES2 <- ventricles2_2
Result2[index3,]$VENTRICLES2 <- ventricles2_3

erc_1 <- length(which(!is.nan(D1$ERC)))
erc_2 <- length(which(!is.nan(D2$ERC)))
erc_3 <- length(which(!is.nan(D3$ERC)))
Result2[index1,]$ERC <- erc_1
Result2[index2,]$ERC <- erc_2
Result2[index3,]$ERC <- erc_3

ventvol_1 <- length(which(!is.nan(D1$VENTVOL)))
ventvol_2 <- length(which(!is.nan(D2$VENTVOL)))
ventvol_3 <- length(which(!is.nan(D3$VENTVOL)))
Result2[index1,]$VENTVOL <- ventvol_1
Result2[index2,]$VENTVOL <- ventvol_2
Result2[index3,]$VENTVOL <- ventvol_3
}
##########

###Final Result###
R <- matrix(NA,(N1*T1+N2*T2+1),n)
Result <- data.frame(R)
names(Result) <- NameCol
Result[1:(N1*T1+N2*T2),] <- rbind(Result1,Result2)
Result[(N1*T1+N2*T2+1),1] <- "Summary"
Result[(N1*T1+N2*T2+1),2] <- "All Visits"
Result[(N1*T1+N2*T2+1),3:n] <- colSums(Result[1:(N1*T1+N2*T2),3:n])
##########

write.table(Result,file="ADNI_DATA_COUNTS.csv",sep=",",row.names=FALSE)
