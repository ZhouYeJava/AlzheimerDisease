############################
###ADNI DATA Construction###
############################

###Zhou Ye###
###01/21/2014###

#The code to merge ADNI data also with Martin Reuter's data
#This code is fine for ADNI data before 04/05/2013

rm(list=ls())
setwd("E:/Research_Bruno/Data/Data_2013_04_05") #Change the route for your own computer or server

#####ALL USEFUL FILES#####
NEURO <- read.csv("NEUROBAT.csv") #1
ADAS1 <- read.csv("ADASSCORES.csv") #2
UCSF <- read.csv("UCSFFSL.csv") #3
UCSD <- read.csv("UCSDVOL.csv") #4
MMSE <- read.csv("MMSE.csv") #5
CDR <- read.csv("CDR.csv") #6
DXSUM <- read.csv("DXSUM_PDXCONV_ADNIALL.csv") #7
ADAS2 <- read.csv("ADAS_ADNIGO2.csv") #8
UPENN <- read.csv("UPENNSPARE.csv") #9
DEMOG <- read.csv("PTDEMOG.csv") #10
UPEB1 <- read.csv("UPENNBIOMK.csv") #11
UPEB2 <- read.csv("UPENNBIOMK2.csv") #12
UPEB3 <- read.csv("UPENNBIOMK3.csv") #13
UPEB4 <- read.csv("UPENNBIOMK4.csv") #14
UPEB5 <- read.csv("UPENNBIOMK5.csv") #15
PIB <- read.csv("PIBPETSUVR.csv") #16
APO <- read.csv("APOERES.csv") #17
MR <- read.csv("Martin_Reuter_Modified.csv") #18
BER <- read.csv("UCBERKELEYAV45.csv") #19
BSI <- read.csv("BSI.csv") #20
FOX <- read.csv("FOXLABBSI.csv") #21
##########

#####RID#####
#This way is very safe since it will obtain all RIDs. 
D <- list()
D[[1]] <- unique(NEURO$RID)
D[[2]] <- unique(ADAS1$RID)
D[[3]] <- unique(UCSF$RID)
D[[4]] <- unique(UCSD$RID)
D[[5]] <- unique(MMSE$RID)
D[[6]] <- unique(CDR$RID)
D[[7]] <- unique(DXSUM$RID)
D[[8]] <- unique(ADAS2$RID)
D[[9]] <- unique(UPENN$RID)
D[[10]] <- unique(DEMOG$RID)
D[[11]] <- unique(UPEB1$RID)
D[[12]] <- unique(UPEB2$RID)
D[[13]] <- unique(UPEB3$RID)
D[[14]] <- unique(UPEB4$RID)
D[[15]] <- unique(UPEB5$RID)
D[[16]] <- unique(PIB$RID)
D[[17]] <- unique(APO$RID)
D[[18]] <- unique(MR$RID)
D[[19]] <- unique(BER$RID)
D[[20]] <- unique(BSI$RID)
D[[21]] <- unique(FOX$RID)
id <- c()
for (i in 1:length(D))
{
	id <- union(id,D[[i]]) #you might add new packages to do it without the loop
}
ID <- unique(id)
ID <- sort(ID) 
N <- length(ID)
##########

#####VISCODE_ADNI1/GO#####
Time1 <- c("f","bl","m03","m06","m12","m18","m24","m30","m36","m42","m48","m54","m60","m66","m72","m78","uns1","nv")
T1 <- length(Time1) #number of stages
##########

#####VISCODE_ADNI2#####
Time2 <- c("v01","v02","v03","v04","v05","v06","v07","v11","v12","v21","v22","v31","v32","v41","v42","v51","v52")
T2 <- length(Time2) #number of stages
##########

#####INITIALIZE RESULT#####
#I have modified the names and they should be updated in a timely fashion
Name <- c("PHASE","RID","VISCODE","AGE","EXAMDATE","DOB","CURRENTDIAG","ADAS11","ADAS13","Hippo","ICV","RHippo","MMSE","CDRSB","SPAREAD","ANARTE","ANARTN","RAVLT1","RAVLT2","RAVLT3","RAVLT4","RAVLT5","RAVLT6","RAVLT30","RAVLTB","RAVLTR","BNT","TMTA","TMTB","TMTE","TAU1","TAU2","TAU3","TAU4","TAU5","TAU","ABETA1","ABETA2","ABETA3","ABETA4","ABETA5","ABETA","PTAU1","PTAU2","PTAU3","PTAU4","PTAU5","PTAU","PIBPET","APGEN1","APGEN2","APOE4","EDU","MRHIPPO","MRICV","MRRHippo","MRPT","MROA","SUVR","VV1","ATEC1","VV2","ATEC2","SVV")
V <- length(Name) #number of variables in columns
R <- matrix(NA,(N*T1+N*T2),V)
Result <- data.frame(R)
colnames(Result) <- Name
Result[(1:(N*T1)),]$RID <- rep(ID,T1)
Result[((N*T1+1):(N*T1+N*T2)),]$RID <- rep(ID,T2)
for (i in 1:T1)
{
	Result[((N*(i-1)+1):(N*i)),]$VISCODE <- Time1[i]
}
for (i in 1:T2)
{
	Result[((N*T1+N*(i-1)+1):(N*T1+N*i)),]$VISCODE <- Time2[i]
}
##########

#####PHASE#####
Result[(1:(N*T1)),]$PHASE <- 1 #ADNI1/ADNIGO
Result[((N*T1+1):(N*T1+N*T2)),]$PHASE <- 2 #ADNI2
##########

#####ADAS#####
num <- nrow(ADAS1)
for (i in 1:num)
{
	rid <- ADAS1[i,]$RID
	viscode <- as.character(ADAS1[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$ADAS11 <- ADAS1[i,]$TOTAL11 #source 1
}

num <- nrow(ADAS2)
for (i in 1:num)
{
	rid <- ADAS2[i,]$RID
	viscode <- as.character(ADAS2[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$ADAS13 <- ADAS2[i,]$TOTAL13 #source 2
}

#set ADAS<0 as NA
Result$ADAS11[Result$ADAS11<0] <- NA
Result$ADAS13[Result$ADAS13<0] <- NA
##########

#####RelativeHippo#####
num <- nrow(UCSF)
for (i in 1:num)
{
	rid <- UCSF[i,]$RID
	viscode <- as.character(UCSF[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	ICV <- UCSF[i,]$ST10CV
	HIPPOL <- UCSF[i,]$ST29SV
	HIPPOR <- UCSF[i,]$ST88SV
	Result[index,]$RHippo <- (HIPPOL+HIPPOR)/ICV
	Result[index,]$Hippo <- HIPPOL+HIPPOR
	Result[index,]$ICV <- ICV
}

#set RelativeHippo<0 as NA
Result$RHippo[Result$RHippo<0] <- NA
Result$ICV[Result$ICV<0] <- NA
Result$Hippo[Result$Hippo<0] <- NA
##########

#####MMSE#####
num <- nrow(MMSE)
for (i in 1:num)
{
	rid <- MMSE[i,]$RID
	viscode <- as.character(MMSE[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$MMSE <- MMSE[i,]$MMSCORE
}

#set MMSE<0 as NA
Result$MMSE[Result$MMSE<0] <- NA
##########

#####CDRSB#####
num <- nrow(CDR)
for (i in 1:num)
{
	rid <- CDR[i,]$RID
	viscode <- as.character(CDR[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	#set negative value to NA
	if (CDR[i,]$CDCARE<0) {CDR[i,]$CDCARE <- NA}
	if (CDR[i,]$CDCOMMUN<0) {CDR[i,]$CDCOMMUN <- NA}
	if (CDR[i,]$CDHOME<0) {CDR[i,]$CDHOME <- NA}
	if (CDR[i,]$CDJUDGE<0) {CDR[i,]$CDJUDGE <- NA}
	if (CDR[i,]$CDMEMORY<0) {CDR[i,]$CDMEMORY <- NA}
	if (CDR[i,]$CDORIENT<0) {CDR[i,]$CDORIENT <- NA}
	Result[index,]$CDRSB <- CDR[i,]$CDCARE+CDR[i,]$CDCOMMUN+CDR[i,]$CDHOME+CDR[i,]$CDJUDGE+CDR[i,]$CDMEMORY+CDR[i,]$CDORIENT
}
##########

#####CURRENTDIAG/EXAMDATE#####
num <- nrow(DXSUM)
for (i in 1:num)
{
	rid <- DXSUM[i,]$RID
	viscode <- as.character(DXSUM[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	phase <- as.character(DXSUM[i,]$Phase)
	if (phase == "ADNI1")
	{
		Result[index,]$CURRENTDIAG <- DXSUM[i,]$DXCURREN}
	else
	{
		diagn <- DXSUM[i,]$DXCHANGE
		if (diagn==7 | diagn==9) {diagn <- 1}
		if (diagn==4 | diagn==8) {diagn <- 2}
		if (diagn==5 | diagn==6) {diagn <- 3}
		Result[index,]$CURRENTDIAG <- diagn
	}
	exam_date <- as.Date(DXSUM[i,]$EXAMDATE,format='%Y-%m-%d')
	Result[index,]$EXAMDATE <- as.character(exam_date)
}

#set CURRENTDIAG<0 as NA
Result$CURRENTDIAG[Result$CURRENTDIAG<0] <- NA
##########

#####SPAREAD#####
num <- nrow(UPENN)
for (i in 1:num)
{
	rid <- UPENN[i,]$RID
	viscode <- as.character(UPENN[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$SPAREAD <- UPENN[i,]$SPAREAD
}

#set SPAREAD<0 as NA
Result$SPAREAD[Result$SPAREAD<0] <- NA
##########

#####ANART/RAVLT/ETC#####
num <- nrow(NEURO)
for (i in 1:num)
{
	rid <- NEURO[i,]$RID
	viscode <- as.character(NEURO[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$ANARTE <- NEURO[i,]$ANARTERR
	Result[index,]$ANARTN <- NEURO[i,]$ANARTND
	Result[index,]$RAVLT1 <- NEURO[i,]$AVTOT1
	Result[index,]$RAVLT30 <- NEURO[i,]$AVDEL30MIN
	Result[index,]$RAVLT2 <- NEURO[i,]$AVTOT2
	Result[index,]$RAVLT3 <- NEURO[i,]$AVTOT3
	Result[index,]$RAVLT4 <- NEURO[i,]$AVTOT4
	Result[index,]$RAVLT5 <- NEURO[i,]$AVTOT5
	Result[index,]$RAVLT6 <- NEURO[i,]$AVTOT6
	Result[index,]$RAVLTB <- NEURO[i,]$AVTOTB
	Result[index,]$RAVLTR <- NEURO[i,]$AVDELTOT
	Result[index,]$BNT <- NEURO[i,]$BNTTOTAL
	Result[index,]$TMTA <- NEURO[i,]$TRAASCOR
	Result[index,]$TMTB <- NEURO[i,]$TRABSCOR
	Result[index,]$TMTE <- NEURO[i,]$TRABERRCOM
}

#set ANART/RAVLT/ETC<0 as NA
Result$ANARTE[Result$ANARTE<0] <- NA
Result$ANARTN[Result$ANARTN<0] <- NA
Result$RAVLT1[Result$RAVLT1<0] <- NA
Result$RAVLT30[Result$RAVLT30<0] <- NA
Result$RAVLT2[Result$RAVLT2<0] <- NA
Result$RAVLT3[Result$RAVLT3<0] <- NA
Result$RAVLT4[Result$RAVLT4<0] <- NA
Result$RAVLT5[Result$RAVLT5<0] <- NA
Result$RAVLT6[Result$RAVLT6<0] <- NA
Result$RAVLTB[Result$RAVLTB<0] <- NA
Result$RAVLTR[Result$RAVLTR<0] <- NA
Result$BNT[Result$BNT<0] <- NA
Result$TMTA[Result$TMTA<0] <- NA
Result$TMTB[Result$TMTB<0] <- NA
Result$TMTB[Result$TMTB==0] <- NA
Result$TMTB[Result$TMTB>300] <- NA
Result$TMTE[Result$TMTE<0] <- NA
##########

#####DOB/AGE#####
DEMOG$VISCODE <- as.character(DEMOG$VISCODE)
DEMOG <- subset(DEMOG,DEMOG$VISCODE=="sc" | DEMOG$VISCODE=="f" | DEMOG$VISCODE=="scmri" | DEMOG$VISCODE=="bl")
num <- nrow(DEMOG)
for (i in 1:num)
{
	rid <- DEMOG[i,]$RID
	index <- which(Result$RID==rid)
	year <- DEMOG[i,]$PTDOBYY
	month <- DEMOG[i,]$PTDOBMM
	day <- 15 #All birthday is treated as 15
	dob <- paste(month,day,year,sep="/",collapse=NULL)
	dob <- as.Date(dob,format='%m/%d/%Y')
	dob <- as.character(dob)
	Result[index,]$DOB <- dob
}

Difference <- difftime(Result$EXAMDATE,Result$DOB,units='days')
Result$AGE <- as.numeric(Difference) #The unit of age is "day"
Result$AGE <- (Result$AGE)/365.25 #Convert to "year"
##########

#####TAU/ABETA/PTAU#####
num <- nrow(UPEB1)
for (i in 1:num)
{
	rid <- UPEB1[i,]$RID
	viscode <- as.character(UPEB1[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$TAU1 <- UPEB1[i,]$TAU
	Result[index,]$ABETA1 <- UPEB1[i,]$ABETA142 #The column name is ABETA142 in UPENNBIOMK!
	Result[index,]$PTAU1 <- UPEB1[i,]$PTAU181P #The column name is PTAU181P in UPENNBIOMK!
}

num <- nrow(UPEB2)
for (i in 1:num)
{
	rid <- UPEB2[i,]$RID
	viscode <- as.character(UPEB2[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$TAU2 <- UPEB2[i,]$TTAU #The column name is TTAU in UPENNBIOMK2!
	Result[index,]$ABETA2 <- UPEB2[i,]$ABETA142 #The column name is ABETA142 in UPENNBIOMK2!
	#No information about PTAU in UPENNBIOMK2
}

num <- nrow(UPEB3)
for (i in 1:num)
{
	rid <- UPEB3[i,]$RID
	viscode <- as.character(UPEB3[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$TAU3 <- UPEB3[i,]$TAU
	Result[index,]$ABETA3 <- UPEB3[i,]$ABETA
	Result[index,]$PTAU3 <- UPEB3[i,]$PTAU
}

num <- nrow(UPEB4)
for (i in 1:num)
{
	rid <- UPEB4[i,]$RID
	viscode <- as.character(UPEB4[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$TAU4 <- UPEB4[i,]$TAU
	Result[index,]$ABETA4 <- UPEB4[i,]$ABETA
	Result[index,]$PTAU4 <- UPEB4[i,]$PTAU
}

num <- nrow(UPEB5)
for (i in 1:num)
{
	rid <- UPEB5[i,]$RID
	viscode <- as.character(UPEB1[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$TAU5 <- UPEB5[i,]$TAU
	Result[index,]$ABETA5 <- UPEB5[i,]$ABETA
	Result[index,]$PTAU5 <- UPEB5[i,]$PTAU
}

#set TAU/ABETA/PTAU<0 as NA
Result$TAU1[Result$TAU1<0] <- NA
Result$TAU2[Result$TAU2<0] <- NA
Result$TAU3[Result$TAU3<0] <- NA
Result$TAU4[Result$TAU4<0] <- NA
Result$TAU5[Result$TAU5<0] <- NA
Result$ABETA1[Result$ABETA1<0] <- NA
Result$ABETA2[Result$ABETA2<0] <- NA
Result$ABETA3[Result$ABETA3<0] <- NA
Result$ABETA4[Result$ABETA4<0] <- NA
Result$ABETA5[Result$ABETA5<0] <- NA
Result$PTAU1[Result$PTAU1<0] <- NA
Result$PTAU2[Result$PTAU2<0] <- NA
Result$PTAU3[Result$PTAU3<0] <- NA
Result$PTAU4[Result$PTAU4<0] <- NA
Result$PTAU5[Result$PTAU5<0] <- NA

TA <- matrix(NA,N*(T1+T2),5)
TA[,1] <- Result$TAU1
TA[,2] <- Result$TAU2
TA[,3] <- Result$TAU3
TA[,4] <- Result$TAU4
TA[,5] <- Result$TAU5
Result$TAU <- rowMeans(TA,na.rm=TRUE)

AB <- matrix(NA,N*(T1+T2),5)
AB[,1] <- Result$ABETA1
AB[,2] <- Result$ABETA2
AB[,3] <- Result$ABETA3
AB[,4] <- Result$ABETA4
AB[,5] <- Result$ABETA5
Result$ABETA <- rowMeans(AB,na.rm=TRUE)

PT <- matrix(NA,N*(T1+T2),5)
PT[,1] <- Result$PTAU1
PT[,2] <- Result$PTAU2
PT[,3] <- Result$PTAU3
PT[,4] <- Result$PTAU4
PT[,5] <- Result$PTAU5
Result$PTAU <- rowMeans(PT,na.rm=TRUE)

Result$TAU1 <- NULL
Result$TAU2 <- NULL
Result$TAU3 <- NULL
Result$TAU4 <- NULL
Result$TAU5 <- NULL
Result$ABETA1 <- NULL
Result$ABETA2 <- NULL
Result$ABETA3 <- NULL
Result$ABETA4 <- NULL
Result$ABETA5 <- NULL
Result$PTAU1 <- NULL
Result$PTAU2 <- NULL
Result$PTAU3 <- NULL
Result$PTAU4 <- NULL
Result$PTAU5 <- NULL
##########

#####PIBPET#####
num <- nrow(PIB)
for (i in 1:num)
{
	rid <- PIB[i,]$RID
	viscode <- as.character(PIB[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	#set negative value to NA
	if (PIB[i,]$ACG<0) {PIB[i,]$ACG <- NA}
	if (PIB[i,]$FRC<0) {PIB[i,]$FRC <- NA}
	if (PIB[i,]$PAR<0) {PIB[i,]$PAR <- NA}
	if (PIB[i,]$PRC<0) {PIB[i,]$PRC <- NA}
	Result[index,]$PIBPET <- (PIB[i,]$ACG+PIB[i,]$FRC+PIB[i,]$PAR+PIB[i,]$PRC)/4
}
##########

#####APGEN/APOE4#####
num <- nrow(APO)
for (i in 1:num)
{
	rid <- APO[i,]$RID
	index <- which(Result$RID==rid)
	Result[index,]$APGEN1 <- APO[i,]$APGEN1
	Result[index,]$APGEN2 <- APO[i,]$APGEN2
	a1 <- APO[i,]$APGEN1
	a2 <- APO[i,]$APGEN2
	if (a1==4 | a2==4)
	{
		Result[index,]$APOE4 <- 1
	}
	else
	{
		Result[index,]$APOE4 <- 0
	}
}

#set APGEN<0 as NA
Result$APGEN1[Result$APGEN1<0] <- NA
Result$APGEN2[Result$APGEN2<0] <- NA
##########

#####EDU#####
num <- nrow(DEMOG)
for (i in 1:num)
{
	rid <- DEMOG[i,]$RID
	index <- which(Result$RID==rid)
	Result[index,]$EDU <- DEMOG[i,]$PTEDUCAT
}

#set EDU<0 as NA
Result$EDU[Result$EDU<0] <- NA
##########

#####MARTIN REUTER#####
#I have modified the Martin Reuter file by R
num <- nrow(MR)
for (i in 1:num)
{
	rid <- MR[i,]$RID
	viscode <- as.character(MR[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$MRHIPPO <- MR[i,]$Left_Hippocampus+MR[i,]$Right_Hippocampus
	Result[index,]$MRICV <- MR[i,]$IntraCranialVol
	Result[index,]$MRET <- MR[i,]$lh_entorhinal_thickness+MR[i,]$rh_entorhinal_thickness
	Result[index,]$MRPT <- MR[i,]$lh_pc_thickness+MR[i,]$rh_pc_thickness
	Result[index,]$MROA <- MR[i,]$lh_oasis.chubs.retrosplenial_thickness+MR[i,]$rh_oasis.chubs.retrosplenial_thickness
}

Result$MRRHippo <- Result$MRHIPPO/Result$MRICV
##########

#####SUVR#####
num <- nrow(BER)
for (i in 1:num)
{
	rid <- BER[i,]$RID
	viscode <- as.character(BER[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$SUVR <- BER[i,]$SUMMARYSUVR_WHOLECEREBNORM
}

#set SUVR<0 as NA
Result$SUVR[Result$SUVR<0] <- NA
##########

#####VENTRICLES1/ENTORHIN1#####
num <- nrow(UCSD)
for (i in 1:num)
{
	rid <- UCSD[i,]$RID
	viscode <- as.character(UCSD[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$VV1 <- UCSD[i,]$VENTRICLES
	Result[index,]$ATEC1 <- (UCSD[i,]$LENTORHIN+UCSD[i,]$RENTORHIN)/2
}

#set VENTRICLES1/ENTORHIN1<0 as NA
Result$VV1[Result$VV1<0] <- NA
Result$ATEC1[Result$ATEC1<0] <- NA
##########

#####VENTRICLES2/ENTORHIN2#####
num <- nrow(UCSF)
for (i in 1:num)
{
	rid <- UCSF[i,]$RID
	viscode <- as.character(UCSF[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$VV2 <- (UCSF[i,]$ST30SV+UCSF[i,]$ST89SV+UCSF[i,]$ST37SV+UCSF[i,]$ST96SV)/2
	Result[index,]$ATEC2 <- (UCSF[i,]$ST24TA+UCSF[i,]$ST83TA)/2
}

#set VENTRICLES2/ENTORHIN2<0 as NA
Result$VV2[Result$VV2<0] <- NA
Result$ATEC2[Result$ATEC2<0] <- NA
##########

#####SVV#####
num <- nrow(BSI)
for (i in 1:num)
{
	rid <- BSI[i,]$RID
	viscode <- as.character(BSI[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$SVV <- BSI[i,]$SVV
}

num <- nrow(FOX)
for (i in 1:num)
{
	rid <- FOX[i,]$RID
	viscode <- as.character(FOX[i,]$VISCODE)
	if (viscode=="sc" | viscode=="scmri") {viscode <- "bl"}
	index <- which(Result$RID==rid & Result$VISCODE==viscode)
	Result[index,]$SVV <- FOX[i,]$VENTVOL
}

#set SVV<0 as NA
Result$SVV[Result$SVV<0] <- NA
##########

#####Final Output#####
Result <- subset(Result,Result$VISCODE!="uns1" & Result$VISCODE!="nv" & Result$VISCODE!="f")
write.table(Result,"ADNI_ALL_DATA.csv",sep=",",na="NaN",row.names=FALSE)
Result1 <- subset(Result,Result$PHASE==1)
Result2 <- subset(Result,Result$PHASE==2)
write.table(Result1,"ADNI1_ADNIGO_ALL_DATA.csv",sep=",",na="NaN",row.names=FALSE)
write.table(Result2,"ADNI2_ALL_DATA.csv",sep=",",na="NaN",row.names=FALSE)
##########
