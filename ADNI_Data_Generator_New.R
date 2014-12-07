###########################
####ADNI Data Generator####
###########################

####Zhou Ye####
####11/08/2014####

#New ADNI data generator code using a different construction approach
#60 times faster than the old code
#Start from 10/09/2014

#Attention:
#DXSUM_PDXCONV_ADNIALL.csv is damaged. There are embedded nulls in the file.
#Use "tr -d '\0' < DXSUM_PDXCONV_ADNIALL.csv > DXSUM_PDXCONV_ADNIALL_CLEAN.csv" in command line to remove nulls
#If you find there is no problem to read DXSUM_PDXCONV_ADNIALL.csv, then do nothing.

rm(list=ls())
setwd("/Users/yezhou/Documents/data-10-9-2014/") #ADNI data folder

####Read Data Files####
visit <- read.csv("VISITS.csv", stringsAsFactors=F) #Visit
adas1 <- read.csv("ADASSCORES.csv", stringsAsFactors=F) #ADNI1 ADAS
adas2 <- read.csv("ADAS_ADNIGO2.csv", stringsAsFactors=F) #ADNIGO/ADNI2 ADAS
ucsf <- read.csv("UCSFFSL_11_01_13.csv", stringsAsFactors=F) #UCSF Longitudinal FreeSurfer (use most recent data)
mmse <- read.csv("MMSE.csv", stringsAsFactors=F) #MMSE Score
cdr <- read.csv("CDR.csv", stringsAsFactors=F) #CDR Sum Score
dxsum <- read.csv("DXSUM_PDXCONV_ADNIALL_CLEAN.csv", stringsAsFactors=F) #Diagnostic 
demog <- read.csv("PTDEMOG.csv", stringsAsFactors=F) #Demographic Information
neuro <- read.csv("NEUROBAT.csv", stringsAsFactors=F) #Neuropsychological Test Information
upen1 <- read.csv("UPENNBIOMK.csv", stringsAsFactors=F) #Biospecimen Results 1
upen2 <- read.csv("UPENNBIOMK2.csv", stringsAsFactors=F) #Biospecimen Results 2
upen3 <- read.csv("UPENNBIOMK3.csv", stringsAsFactors=F) #Biospecimen Results 3
upen4 <- read.csv("UPENNBIOMK4_09_06_12.csv", stringsAsFactors=F) #Biospecimen Results 4
upen5 <- read.csv("UPENNBIOMK5_10_31_13.csv", stringsAsFactors=F) #Biospecimen Results 5
upen6 <- read.csv("UPENNBIOMK6_07_02_13.csv", stringsAsFactors=F) #Biospecimen Results 6
upen7 <- read.csv("UPENNBIOMK7.csv", stringsAsFactors=F) #Biospecimen Results 7
apoe <- read.csv("APOERES.csv", stringsAsFactors=F) #Genetype Information
bsi <- read.csv("BSI.csv", stringsAsFactors=F) #Boundary Shift Integral Summary
fox1 <- read.csv("FOXLABBSI_11_01_13.csv", stringsAsFactors=F) #Fox Lab 1
fox2 <- read.csv("FOXLABBSI_01_31_14.csv", stringsAsFactors=F) #Fox Lab 2
fox3 <- read.csv("FOXLABBSI_06_23_14.csv", stringsAsFactors=F) #Fox Lab 3
fox4 <- read.csv("FOXLABBSI_11_01_13.csv", stringsAsFactors=F) #Fox Lab 4

####Roster ID####
id_list <- list(adas1$RID, adas2$RID, ucsf$RID, mmse$RID, cdr$RID, dxsum$RID,
                demog$RID, neuro$RID, upen1$RID, upen2$RID, upen3$RID, 
                upen4$RID, upen5$RID, upen6$RID, upen7$RID, apoe$RID, 
                bsi$RID, fox1$RID, fox2$RID, fox3$RID, fox4$RID)
id <- sort(Reduce("union", id_list)) #merge RID from all the files
n <- length(id)

####Visit Code####
adni1_visit <- unique(subset(visit, Phase!="ADNI2")$VISCODE) #ADNI1/ADNIGO
adni2_visit <- unique(subset(visit, Phase=="ADNI2")$VISCODE) #ADNI2
trash_visit <- c("sc", "f", "scmri", "uns1", "nv") #useless visits
#sc: screen which is included in bl
#f: screen fail which is included in bl
#scmri: screen MRI which is included in bl
#uns1: unscheduled visit
#nv: no visit information
adni1_visit <- adni1_visit[!adni1_visit%in%trash_visit]
adni2_visit <- adni2_visit[!adni2_visit%in%trash_visit]
t1 <- length(adni1_visit)
t2 <- length(adni2_visit)

####Initialize Final Data Frame####
name <- c("PHASE", "RID", "VISCODE")
d <- length(name)
adni <- data.frame(matrix(NA, n*t1+n*t2, d))
colnames(adni) <- name
adni$RID <- id #RID
adni[1:(n*t1),]$PHASE <- 1 #ADNI1/ADNIGO identifier
adni[(n*t1+1):(n*t1+n*t2),]$PHASE <- 2 #ADNI2 identifier
adni[1:(n*t1),]$VISCODE <- rep(adni1_visit, each=n) #ADNI1/ADNIGO visit
adni[(n*t1+1):(n*t1+n*t2),]$VISCODE <- rep(adni2_visit, each=n) #ADNI2 visit
#define time order to help sort
visit_order <- c(adni1_visit, adni2_visit)
adni$VISCODE <- ordered(adni$VISCODE, visit_order)

####Exam Date/Diagnostic####
dxsum <- with(dxsum, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    exam_date <- as.Date(EXAMDATE,format="%Y-%m-%d")
    temp <- data.frame(RID=RID, VISCODE=VISCODE, DIAG=NA, ED=as.character(exam_date))
    index0 <- which(Phase=="ADNI1") #ADNI1 patient
    index1 <- which(Phase!="ADNI1" & (DXCHANGE==7 | DXCHANGE==9)) #diag = 1
    index2 <- which(Phase!="ADNI1" & (DXCHANGE==4 | DXCHANGE==8)) #diag = 2
    index3 <- which(Phase!="ADNI1" & (DXCHANGE==5 | DXCHANGE==6)) #diag = 3
    if (length(index0)!=0) {
        temp[index0,]$DIAG <- DXCURREN[index0]
    }
    if (length(index1)!=0) {
        temp[index1,]$DIAG <- 1
    }
    if (length(index2)!=0) {
        temp[index2,]$DIAG <- 2
    }
    if (length(index3)!=0) {
        temp[index3,]$DIAG <- 3
    }
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    temp[!duplicated(temp$RID, temp$VISCODE),] #eliminate duplicates
})
adni <- merge(adni, dxsum, by=c("VISCODE", "RID"), all.x=T)

####Birthday/Age####
demog <- subset(demog, demog$VISCODE=="sc" | demog$VISCODE=="f" | demog$VISCODE=="scmri" | demog$VISCODE=="bl")
demog <- with(demog, {
    dob <- as.Date(paste(PTDOBMM, 15, PTDOBYY, sep="-"), format="%m-%d-%Y") #all birthday is treated on 15
    temp <- data.frame(RID=RID, DOB=as.character(dob))
    temp[!duplicated(temp$RID),]
})
adni <- merge(adni, demog, by="RID", all.x=T)
difference <- difftime(adni$ED, adni$DOB, units='days')
adni$AGE <- as.numeric(difference)/365.25 #convert to "year"
adni <- adni[order(adni$VISCODE),] #order will change after the merge

####APOE4####
apoe <- with(apoe, {
    APGEN1[APGEN1<0] <- NA
    APGEN2[APGEN2<0] <- NA
    temp <- data.frame(RID=RID, APGEN1=APGEN1, APGEN2=APGEN2, APOE4=0)
    temp <- within(temp, {
        APOE4[APGEN1==4 | APGEN2==4] <- 1
    })
})
adni <- merge(adni, apoe, by="RID", all.x=T)
adni <- adni[order(adni$VISCODE),] #order will change after the merge

####ADAS####
adas1 <- with(adas1, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    TOTAL11[TOTAL11<0] <- NA 
    TOTALMOD[TOTALMOD<0] <- NA 
    temp <- data.frame(RID=RID, VISCODE=VISCODE, ADAS11=TOTAL11, ADAS13=TOTALMOD)
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=cbind(ADAS11, ADAS13)~RID+VISCODE, data=temp, FUN=mean, na.rm=T) #remove duplicates
})
adas2 <- with(adas2, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    TOTSCORE[TOTSCORE<0] <- NA 
    TOTAL13[TOTAL13<0] <- NA 
    temp <- data.frame(RID=RID, VISCODE=VISCODE, ADAS11=TOTSCORE, ADAS13=TOTAL13)
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=cbind(ADAS11, ADAS13)~RID+VISCODE, data=temp, FUN=mean, na.rm=T)
})
adas <- rbind(adas1, adas2)
adni <- merge(adni, adas, by=c("VISCODE", "RID"), all.x=T)

####Hippocampus/Ventricle/Entorhin####
ucsf <- with(ucsf, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    ST10CV[ST10CV<0] <- NA #ICV
    ST29SV[ST29SV<0] <- NA #Left Hippocampus
    ST88SV[ST88SV<0] <- NA #Right Hippocampus
    ST30SV[ST30SV<0] <- NA #Ventrocle Part 1
    ST89SV[ST89SV<0] <- NA #Ventrocle Part 2
    ST37SV[ST37SV<0] <- NA #Ventrocle Part 3
    ST96SV[ST96SV<0] <- NA #Ventrocle Part 4
    ST24TA[ST24TA<0] <- NA #Left Entorhin
    ST83TA[ST83TA<0] <- NA #Right Entorhin
    temp <- data.frame(RID=RID, VISCODE=VISCODE, HC=ST29SV+ST88SV, VV=ST30SV+ST89SV+ST37SV+ST96SV, ICV=ST10CV, EC=ST24TA+ST83TA)
    temp <- within(temp, {
        RHC <- HC/ICV #Relative Hippocampus Volume
        RVV <- VV/ICV #Relative Ventrocle Volume
    })
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=cbind(HC, VV, ICV, EC, RHC, RVV)~RID+VISCODE, data=temp, FUN=mean, na.rm=T)
})
adni <- merge(adni, ucsf, by=c("VISCODE", "RID"), all.x=T)

####MMSE####
mmse <- with(mmse, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    MMSCORE[MMSCORE<0] <- NA
    temp <- data.frame(RID=RID, VISCODE=VISCODE, MMSE=MMSCORE)
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=MMSE~RID+VISCODE, data=temp, FUN=mean, na.rm=T)
})
adni <- merge(adni, mmse, by=c("VISCODE", "RID"), all.x=T)

####CDRSB####
cdr <- with(cdr, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    CDMEMORY[CDMEMORY<0] <- NA #memory
    CDORIENT[CDORIENT<0] <- NA #orientation
    CDJUDGE[CDJUDGE<0] <- NA #judgement & problem solving 
    CDCOMMUN[CDCOMMUN<0] <- NA #community affairs
    CDHOME[CDHOME<0] <- NA #home and hobbies
    CDCARE[CDCARE<0] <- NA #personal care
    temp <- data.frame(RID=RID, VISCODE=VISCODE, 
                       CDRSS=CDMEMORY+CDORIENT+CDJUDGE+CDCOMMUN+CDHOME+CDCARE)
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=CDRSS~RID+VISCODE, data=temp, FUN=mean, na.rm=T)
})
adni <- merge(adni, cdr, by=c("VISCODE", "RID"), all.x=T)

####RAVLT/TMT####
ravlt <- with(neuro, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    AVDEL30MIN[AVDEL30MIN<0] <- NA
    TRAASCOR[TRAASCOR<0] <- NA
    TRABSCOR[TRABSCOR<0] <- NA
    TRABERRCOM[TRABERRCOM<0] <- NA
    temp <- data.frame(RID=RID, VISCODE=VISCODE, RAVLT30=AVDEL30MIN, TMTA=TRAASCOR, TMTB=TRABSCOR, TMTE=TRABERRCOM)
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=cbind(RAVLT30, TMTA, TMTB, TMTE)~RID+VISCODE, data=temp, FUN=mean, na.rm=T)
})
adni <- merge(adni, ravlt, by=c("VISCODE", "RID"), all.x=T)

####TAU/ABETA/PTAU####
#use outer join to include all id
upen <- data.frame(RID=upen1$RID, VISCODE=upen1$VISCODE, TAU1=upen1$TAU, ABETA1=upen1$ABETA142, PTAU1=upen1$PTAU181P)
upen <- merge(upen, data.frame(RID=upen2$RID, VISCODE=upen2$VISCODE, TAU2=upen2$TTAU, ABETA2=upen2$ABETA142, PTAU2=NA), by=c("VISCODE", "RID"), all=T)
upen <- merge(upen, data.frame(RID=upen3$RID, VISCODE=upen3$VISCODE, TAU3=upen3$TAU, ABETA3=upen3$ABETA, PTAU3=upen3$PTAU), by=c("VISCODE", "RID"), all=T)
upen <- merge(upen, data.frame(RID=upen4$RID, VISCODE=upen4$VISCODE, TAU4=upen4$TAU, ABETA4=upen4$ABETA, PTAU4=upen4$PTAU), by=c("VISCODE", "RID"), all=T)
upen <- merge(upen, data.frame(RID=upen5$RID, VISCODE=upen5$VISCODE, TAU5=upen5$TAU, ABETA5=upen5$ABETA, PTAU5=upen5$PTAU), by=c("VISCODE", "RID"), all=T)
upen <- merge(upen, data.frame(RID=upen6$RID, VISCODE=upen6$VISCODE, TAU6=upen6$TAU, ABETA6=upen6$ABETA, PTAU6=upen6$PTAU), by=c("VISCODE", "RID"), all=T)
upen <- merge(upen, data.frame(RID=upen7$RID, VISCODE=upen7$VISCODE, TAU7=upen7$TAU, ABETA7=upen7$ABETA, PTAU7=upen7$PTAU), by=c("VISCODE", "RID"), all=T)
upen <- with(upen, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    #NA need not be considered here
    temp <- data.frame(
        RID=RID, VISCODE=VISCODE, 
        TAU=rowMeans(data.frame(TAU1, TAU2, TAU3, TAU4, TAU5, TAU6, TAU7), na.rm=T),
        ABETA=rowMeans(data.frame(ABETA1, ABETA2, ABETA3, ABETA4, ABETA5, ABETA6, ABETA7), na.rm=T),
        PTAU=rowMeans(data.frame(PTAU1, PTAU2, PTAU3, PTAU4, PTAU5, PTAU6, PTAU7), na.rm=T)
    )
    temp <- subset(temp, !VISCODE%in%c("uns1", "nv"))
    aggregate(formula=cbind(TAU, ABETA, PTAU)~RID+VISCODE, data=temp, FUN=mean, na.rm=T)
})
adni <- merge(adni, upen, by=c("VISCODE", "RID"), all.x=T)

####Segmented Ventricle Volume####
svv <- rbind(data.frame(RID=bsi$RID, VISCODE=bsi$VISCODE, SVV=bsi$VENTVOL),
             data.frame(RID=fox1$RID, VISCODE=fox1$VISCODE, SVV=fox1$VENTVOL),
             data.frame(RID=fox2$RID, VISCODE=fox2$VISCODE, SVV=fox2$VENTVOL),
             data.frame(RID=fox3$RID, VISCODE=fox3$VISCODE, SVV=fox3$VENTVOL),
             data.frame(RID=fox4$RID, VISCODE=fox4$VISCODE, SVV=fox4$VENTVOL))
svv <- within(svv, {
    VISCODE[VISCODE=="sc" | VISCODE=="scmri" | VISCODE=="f"] <- "bl"
    SVV[SVV<0] <- NA
})
svv <- subset(svv, !VISCODE%in%c("uns1", "nv"))
svv <- aggregate(formula=SVV~RID+VISCODE, data=svv, FUN=mean, na.rm=T)
adni <- merge(adni, svv, by=c("VISCODE", "RID"), all.x=T)

####Final Output####
write.csv(adni, "ADNI_Data_10_9_2014.csv", na="NaN", row.names=F)

