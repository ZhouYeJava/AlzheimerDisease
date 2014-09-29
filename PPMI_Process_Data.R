##################
###PPMI Process###
##################

###Zhou Ye###
###02/27/2014###

rm(list=ls())
#change the working directory based on your machine
setwd("/Users/zhouye/Documents/Bruno/PPMI/data-10-04-2013/")

###Data Sets###
random <- read.csv("Consent_and_Enrollment.csv", stringsAsFactors=F)
screen <- read.csv("Screening_Demographics.csv", stringsAsFactors=F)
family <- read.csv("Family_History__PD_.csv", stringsAsFactors=F)
feature <- read.csv("PD_Features.csv", stringsAsFactors=F)
nupdrs2p <- read.csv("MDS_UPDRS_Part_II__Patient_Questionnaire.csv", stringsAsFactors=F)
nupdrs3 <- read.csv("MDS_UPDRS_Part_III__Post_Dose_.csv", stringsAsFactors=F)
nupdrs1 <- read.csv("MDS_UPDRS_Part_I.csv", stringsAsFactors=F)
nupdrs1p <- read.csv("MDS_UPDRS_Part_I__Patient_Questionnaire.csv", stringsAsFactors=F)
bjlos <- read.csv("Benton_Judgment_of_Line_Orientation.csv", stringsAsFactors=F)
ess <- read.csv("Epworth_Sleepiness_Scale.csv", stringsAsFactors=F)
gdsrs <- read.csv("Geriatric_Depression_Scale__Short_.csv", stringsAsFactors=F)
hvlt <- read.csv("Hopkins_Verbal_Learning_Test.csv", stringsAsFactors=F)
lns <- read.csv("Letter-Number_Sequencing__PD_.csv", stringsAsFactors=F)
moca <- read.csv("Montreal_Cognitive_Assessment__MoCA_.csv", stringsAsFactors=F)
edu <- read.csv("Socio-Economics.csv", stringsAsFactors=F)
quip <- read.csv("QUIP_Current_Short.csv", stringsAsFactors=F)
rbd <- read.csv("REM_Sleep_Disorder_Questionnaire.csv", stringsAsFactors=F)
scopa <- read.csv("SCOPA-AUT.csv", stringsAsFactors=F)
sft <- read.csv("Semantic_Fluency.csv", stringsAsFactors=F)
stai <- read.csv("State-Trait_Anxiety_Inventory.csv", stringsAsFactors=F)
upsit <- read.csv("Univ._of_Pennsylvania_Smell_ID_Test.csv", stringsAsFactors=F)
cogcatg <- read.csv("Clinical_Cognitive_Categorization.csv", stringsAsFactors=F)
sbr <- read.csv("Striatal_Binding_Ratio_Results.csv", stringsAsFactors=F)
bio <- read.csv("Biospecimen_Analysis_Results.csv", stringsAsFactors=F)
cat("Finish importing data sets!")

###Initialization###
#column names
name <- c("PATNO", "TIME", "DIAG", "AGE", "RACE", "HISTORY", "DURATION",
          "TD_PIGD", "MDS_UPDRS_I", "MDS_UPDRS_II", "MDS_UPDRS_III",
          "MDS_UPDRS", "BJLOS", "ESS", "GDSRS", "HVLTIR",
          "HVLTDG", "LNS", "MOCA", "QUIP", "REMSBD", "SCOPA_AUT",
          "SFT", "STAT", "STASS", "STATS", "UPSITRS", "MCI",
          "CONTRALATERAL", "IPSILATERAL", "CAUDATE", "PUTAMEN",
          "STRIATUM", "CDR", "AI", "ABETA", "TAU", "PTAU",
          "ALPHA", "TAU_ABETA", "PTAU_ABETA", "PTAU_TAU",
          "URATE")
#The followings are well-defined time visits
time <- c("SC", "BL", "V01", "V02", "V03", "V04", "V05", "V06", 
          "V07", "V08", "V09", "V10", "V11", "PW") #PW=V12
randomData <- cbind(random$PATNO, random$ENROLLDT)
colnames(randomData) <- c("PATNO", "ENROLLDT")
screenData <- cbind(screen$PATNO, screen$APPRDX)
colnames(screenData) <- c("PATNO", "APPRDX")
temp <- merge(randomData, screenData)
temp <- subset(temp, !is.na(temp$ENROLLDT)) #choose patients with enrollment dates
id <- temp$PATNO #enrolled patients
diag <- temp$APPRDX #dignostic for each patient
numItem <- length(name) #number of columns
numTime <- length(time) #number of stages
numPatient <- length(id) #number of patients
ppmi <- data.frame(matrix(NA, numPatient*numTime, numItem)) #final data set
colnames(ppmi) <- name
ppmi$PATNO <- id #update patient ID
ppmi$TIME <- rep(time, each=numPatient) #update time visits
ppmi$DIAG <- diag #update dignostic
cat("Finish Initialization!")

###Age(RANDOM)###
#The data set does not provide the "day" of birth. I assume
#They are all at the first day of the year since only
#the difference matters. 
enroll <- paste(random$ENROLLDT, "/01", sep="")
birth <- paste(random$BIRTHDT, "/01", sep="")
enroll <- strftime(as.Date(enroll, "%m/%Y/%d"), "%Y-%m-%d")
birth <- strftime(as.Date(birth, "%m/%Y/%d"), "%Y-%m-%d")
age <- as.numeric(difftime(enroll, birth))/365 #age in units year
random$AGE <- age
for (i in 1:nrow(random)) {
  id <- random[i,]$PATNO
  index <- which(ppmi$PATNO==id & ppmi$TIME%in%c("SC", "BL"))
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$AGE <- random[i,]$AGE
}
timeDiff <- c(3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 54, 60)
#age at V01 to V12
ppmi[(2*numPatient+1):nrow(ppmi),]$AGE <- age+rep(timeDiff/12, each=numPatient)
cat("Age has been updated!")

###Race(SCREEN)###
#1 -> RAINDALS; 2 -> RAASIAN; 3 -> RABLACK; 4 -> RAHAWOPI;
#5 -> RAWHITE; 6 -> RANOS; 7 -> more than one races; 8 -> no races
indexRace <- which(names(screen)%in%c("RAINDALS", "RAASIAN", 
             "RABLACK", "RAHAWOPI", "RAWHITE", "RANOS"))
for (i in 1:nrow(screen)) {
  id <- screen[i,]$PATNO
  index <- which(ppmi$PATNO==id)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  race <- which(screen[i,indexRace]==1)
  if (length(race)==1) {
    ppmi[index,]$RACE <- race
  }
  else if (length(race)>1){
    ppmi[index,]$RACE <- 7
  }
  else {
    ppmi[index,]$RACE <- 8
  }
}
cat("Race has been updated!")

###Family History(FAMHXPD)###
#1 -> yes; 0 -> no
indexHistory <- which(names(family)%in%c("BIOMOMPD", "BIODADPD",
               "FULSIBPD", "HALFSIBPD", "MAGPARPD", "PAGPARPD",
               "MATAUPD", "PATAUPD", "KIDSPD"))
for (i in 1:nrow(family)) {
  id <- family[i,]$PATNO
  index <- which(ppmi$PATNO==id)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  fh <- which(family[i,indexHistory]==1)
  if (length(fh)!=0) {
    ppmi[index,]$HISTORY <- 1 #has history
  }
  else {
    ppmi[index,]$HISTORY <- 0 #no history
  }
}
cat("Family History has been updated!")

###Duration of Disease(PDFEAT, RANDOM)###
#The data set does not provide the "day" of birth. I assume
#They are all at the first day of the year since only
#the difference matters. 
randomSub <- random$PATNO
enroll <- paste(random$ENROLLDT, "/01", sep="")
enroll <- strftime(as.Date(enroll, "%m/%Y/%d"), "%Y-%m-%d")
pdfeatSub <- feature$PATNO
pddxdt <- paste(feature$PDDXDT, "/01", sep="")
pddxdt <- strftime(as.Date(pddxdt, "%m/%Y/%d"), "%Y-%m-%d")
randomData <- data.frame(cbind(randomSub, enroll))
colnames(randomData) <- c("PATNO", "ENROLL")
featData <- data.frame(cbind(pdfeatSub, pddxdt))
colnames(featData) <- c("PATNO", "PDDATE")
temp <- merge(randomData, featData)
temp$DURATION <- as.numeric(difftime(temp$ENROLL, temp$PDDATE), units="weeks")
#Duration = number of months between two dates
for (i in 1:nrow(temp)) {
  id <- temp[i,]$PATNO
  index <- which(ppmi$PATNO==id)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  if (!is.na(temp[i,]$DURATION)) {
    ppmi[index,]$DURATION <- temp[i,]$DURATION
  }
}
cat("Disease Duration has been updated!")

###TD/PIGD Classification(NUPDRS2P, NUPDRS3)###
#1 -> TD; 2 -> PIGD; 3 -> Indeterminate; 
part2Name <- c("PATNO", "EVENT_ID", "NP2TRMR", "NP2WALK", "NP2FREZ")
part2Index <- which(names(nupdrs2p)%in%part2Name)
part2 <- nupdrs2p[,part2Index]
part3Name <- c("PATNO", "EVENT_ID", "NP3PTRMR", "NP3PTRML",
               "NP3KTRMR", "NP3KTRML", "NP3RTARU",
               "NP3RTALU", "NP3RTARL", "NP3RTALL",
               "NP3RTALJ", "NP3RTCON", "NP3GAIT",
               "NP3FRZGT", "NP3PSTBL")
part3Index <- which(names(nupdrs3)%in%part3Name)
part3 <- nupdrs3[,part3Index]
temp <- merge(part2, part3)
temp$TREMOR <- rowMeans(cbind(temp$NP3PTRMR, temp$NP3PTRML,
                    temp$NP3KTRMR, temp$NP3KTRML, temp$NP3RTARU,
                    temp$NP3RTALU, temp$NP3RTARL, temp$NP3RTALL,
                    temp$NP3RTALJ, temp$NP3RTCON))
temp$PIGD <- rowMeans(cbind(temp$NP2WALK, temp$NP2FREZ, temp$NP3GAIT,
                     temp$NP3FRZGT, temp$NP3PSTBL))
for (i in 1:nrow(temp)) {
  id <- temp[i,]$PATNO
  time <- temp[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  if (is.na(temp[i,]$TREMOR) | is.na(temp[i,]$PIGD)) {
    next #missing value is useless
  }
  if (is.nan(temp[i,]$TREMOR/temp[i,]$PIGD)) {
    next #NaN is useless
  }
  if (temp[i,]$TREMOR/temp[i,]$PIGD>=1.15) {
    ppmi[index,]$TD_PIGD <- 1
  }
  else if (temp[i,]$TREMOR/temp[i,]$PIGD<=0.9) {
    ppmi[index,]$TD_PIGD <- 2
  }
  else if ((temp[i,]$TREMOR/temp[i,]$PIGD > 0.9 
           & temp[i,]$TREMOR/temp[i,]$PIGD<1.15) 
           | (temp[i,]$PIGD==0 & temp[i,]$TREMOR==0)) {
    ppmi[index,]$TD_PIGD <- 3
  }
  else if (temp[i,]$PIGD==0 & temp[i,]$TREMOR>0) {
    ppmi[index,]$TD_PIGD <- 1
  }
  else {
    next
  }
}
cat("TD/PIGD Classification has been updated!")

###MDS-UPDRS Part I(NUPDRS1, NUPDRS1P)###
part1Name <- c("PATNO", "EVENT_ID", "NP1COG", "NP1HALL",
               "NP1DPRS", "NP1ANXS", "NP1APAT", "NP1DDS")
part1Index <- which(names(nupdrs1)%in%part1Name)
part1 <- nupdrs1[,part1Index]
part1qName <- c("PATNO", "EVENT_ID", "NP1SLPN", "NP1SLPD",
                "NP1PAIN", "NP1URIN", "NP1CNST", "NP1LTHD", "NP1FATG")
part1qIndex <- which(names(nupdrs1p)%in%part1qName)
part1q <- nupdrs1p[,part1qIndex]
temp <- merge(part1, part1q)
temp$SUM <- rowSums(temp[,-which(names(temp)%in%c("PATNO", "EVENT_ID"))])
for (i in 1:nrow(temp)) {
  id <- temp[i,]$PATNO
  time <- temp[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$MDS_UPDRS_I <- temp[i,]$SUM
}
cat("MDS-UPDRS Part I has been updated!")

###MDS-UPDRS Part II(NUPDRS2P)###
part2qName <- c("PATNO", "EVENT_ID", "NP2SPCH", "NP2SALV",
                "NP2SWAL", "NP2EAT", "NP2DRES", "NP2HYGN",
                "NP2HWRT", "NP2HOBB", "NP2TURN", "NP2TRMR",
                "NP2RISE", "NP2WALK", "NP2FREZ")
part2qIndex <- which(names(nupdrs2p)%in%part2qName)
part2q <- nupdrs2p[,part2Index]
part2q$SUM <- rowSums(part2q[,-which(names(part2q)%in%c("PATNO", "EVENT_ID"))])
for (i in 1:nrow(part2q)) {
  id <- part2q[i,]$PATNO
  time <- part2q[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$MDS_UPDRS_II <- part2q[i,]$SUM
}
cat("MDS-UPDRS Part II has been updated!")

###MDS-UPDRS Part III(NUPDRS3)###
part3Name <- c("PATNO", "EVENT_ID", "NP3SPCH", "NP3FACXP",
               "NP3RIGN", "NP3RIGRU", "NP3RIGLU", "PN3RIGRL",
               "NP3RIGLL", "NP3FTAPR", "NP3FTAPL", "NP3HMOVR",
               "NP3HMOVL", "NP3PRSPR", "NP3PRSPL", "NP3TTAPR",
               "NP3TTAPL", "NP3LGAGR", "NP3LGAGL", "NP3RISNG",
               "NP3GAIT", "NP3FRZGT", "NP3PSTBL", "NP3POSTR",
               "NP3BRADY", "NP3PTRMR", "NP3PTRML", "NP3KTRMR",
               "NP3KTRML", "NP3RTARU", "NP3RTALU", "NP3RTARL",
               "NP3RTALL", "NP3RTALJ", "NP3RTCON")
part3Index <- which(names(nupdrs3)%in%part3Name)
part3 <- nupdrs3[,part3Index]
part3$SUM <- rowSums(part3[,-which(names(part3)%in%c("PATNO", "EVENT_ID"))])
for (i in 1:nrow(part3)) {
  id <- part3[i,]$PATNO
  time <- part3[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$MDS_UPDRS_III <- part3[i,]$SUM
}
cat("MDS-UPDRS Part III has been updated!")

###MDS-UPDRS Total Score(NUPDRS1, NUPDRS1P, NUPDRS2P, NUPDRS3)###
ppmi$MDS_UPDRS <- rowSums(cbind(ppmi$MDS_UPDRS_I, ppmi$MDS_UPDRS_II, ppmi$MDS_UPDRS_III))
cat("MDS-UPDRS Total Score has been updated!")

###Benton Judgement of Line Orientation Score(LINEORNT)###
bjlosName <- paste("BJLOT", 1:30, sep="")
bjlosIndex <- which(names(bjlos)%in%bjlosName)
bjlos$SUM <- rowSums(bjlos[,bjlosIndex])
for (i in 1:nrow(bjlos)) {
  id <- bjlos[i,]$PATNO
  time <- bjlos[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$BJLOS <- bjlos[i,]$SUM
}
cat("Benton Judgement of Line Orientation Score has been updated!")

###Epworth Sleepiness Scale##(EPWORTH)###
#1 -> sleepy; 0 -> not sleepy
essName <- paste("ESS", 1:8, sep="")
essIndex <- which(names(ess)%in%essName)
ess$SUM <- rowSums(ess[,essIndex])
for (i in 1:nrow(ess)) {
  id <- ess[i,]$PATNO
  time <- ess[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$ESS <- ess[i,]$SUM
}
cat("Epworth Sleepiness Scale has been updated!")

###GDS Raw Score(GDSSHORT)###
#1 -> depressed; 0 -> not depressed
gdsName1 <- c("GDSSATIS", "GDSGSPIR", "GDSHAPPY", "GDSALIVE", "GDSENRGY")
gdsIndex1 <- which(names(gdsrs)%in%gdsName1)
gdsName2 <- c("GDSDROPD", "GDSEMPTY", "GDSBORED", "GDSAFRAD",
              "GDSHLPLS", "GDSHOME", "GDSMEMRY", "GDSWRTLS",
              "GDSHOPLS", "GDSBETER")
gdsIndex2 <- which(names(gdsrs)%in%gdsName2)
for (i in 1:nrow(gdsrs)) {
  id <- gdsrs[i,]$PATNO
  time <- gdsrs[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  score1 <- length(which(gdsrs[i,gdsIndex1]==0))
  score2 <- length(which(gdsrs[i,gdsIndex2]==1))
  ppmi[index,]$GDSRS <- score1+score2
}
cat("GDS Raw Score has been updated!")

###HVLT Immediate Recall/HVLT Discrimination Recognition(HVLT)###
hvlt$RECALL <- hvlt$HVLTRT1+hvlt$HVLTRT2+hvlt$HVLTRT3
hvlt$RECOG <- hvlt$HVLTREC-(hvlt$HVLTFPRL+hvlt$HVLTFPUN)
for (i in 1:nrow(hvlt)) {
  id <- hvlt[i,]$PATNO
  time <- hvlt[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$HVLTIR <- hvlt[i,]$RECALL
  ppmi[index,]$HVLTDG <- hvlt[i,]$RECOG
}
cat("HVLT Recall and Recognition have been updated!")

###Letter Number Sequencing(LNSPD)###
lnsName <- c(paste("LNS", 1:7, "A", sep=""), 
             paste("LNS", 1:7, "B", sep=""),
             paste("LNS", 1:7, "C", sep=""))
lnsIndex <- which(names(lns)%in%lnsName)
lns$SUM <- rowSums(lns[,lnsIndex])
for (i in 1:nrow(lns)) {
  id <- lns[i,]$PATNO
  time <- lns[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$LNS <- lns[i,]$SUM
}
cat("Letter Number Sequencing has been updated!")

###MOCA Total Score(MOCA, SOCIOECO)###
mocaName <- c("PATNO", "EVENT_ID", "MCAALTTM", "MCACUBE", "MCACLCKC", "MCACLCKN", 
              "MCACLCKH", "MCALION", "MCARHINO", "MCACAMEL", 
              "MCAFDS", "MCABDS", "MCAVIGIL", "MCASER7", "MCASNTNC", 
              "MCAVF", "MCAABSTR", "MCAREC1", "MCAREC2", "MCAREC3", 
              "MCAREC4", "MCAREC5", "MCADATE", "MCAMONTH", "MCAYR", 
              "MCADAY", "MCAPLACE", "MCACITY")
mocaIndex <- which(names(moca)%in%mocaName)
mocaData <- moca[,mocaIndex]
mocaData$SUM <- rowSums(mocaData[,-which(names(mocaData)%in%c("PATNO", "EVENT_ID"))])
eduName <- c("PATNO", "EVENT_ID", "EDUCYRS")
eduIndex <- which(names(edu)%in%eduName)
eduData <- edu[,eduIndex]
temp <- merge(mocaData, eduData)
for (i in 1:nrow(temp)) {
  id <- temp[i,]$PATNO
  time <- temp[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  if (is.na(temp[i,]$SUM)) {
    next #missing value is useless
  }
  if (!is.na(temp[i,]$EDUCYRS)) {
    if (temp[i,]$EDUCYRS<=12 & temp[i,]$SUM<30) {
      ppmi[index,]$MOCA <- temp[i,]$SUM+1 #adjustment
    }
    else {
      ppmi[index,]$MOCA <- temp[i,]$SUM
    }
  }
  else {
    ppmi[index,]$MOCA <- temp[i,]$SUM
  }
}
cat("MOCA Total Score has been updated!")

###QUIP(QUIPCS)###
secA <- c("CNTRLGMB", "TMGAMBLE")
secAIndex <- which(names(quip)%in%secA)
secB <- c("CNTRLSEX", "TMSEX")
secBIndex <- which(names(quip)%in%secB)
secC <- c("CNTRLBUY", "TMBUY")
secCIndex <- which(names(quip)%in%secC)
secD <- c("CNTRLEAT", "TMEAT")
secDIndex <- which(names(quip)%in%secD)
secE <- c("TMTORACT", "TMTMTACT", "TMTRWD")
secEIndex <- which(names(quip)%in%secE)
for (i in 1:nrow(quip)) {
  id <- quip[i,]$PATNO
  time <- quip[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  score <- 0
  if (length(which(quip[i,secAIndex]==1))!=0) {
    score <- score+1
  }
  if (length(which(quip[i,secBIndex]==1))!=0) {
    score <- score+1
  }
  if (length(which(quip[i,secCIndex]==1))!=0) {
    score <- score+1
  }
  if (length(which(quip[i,secDIndex]==1))!=0) {
    score <- score+1
  }
  score <- score+length(which(quip[i,secEIndex]==1))
  ppmi[index,]$QUIP <- score
}
cat("QUIP has been updated!")

###REM Sleep Behavior Disorder(REMSLEEP)###
group1 <- c("DRMVIVID", "DRMAGRAC", "DRMNOCTB", "SLPLMBMV", "SLPINJUR", 
            "DRMVERBL", "DRMFIGHT", "DRMUMV", "DRMOBJFL", "MVAWAKEN", 
            "DRMREMEM", "SLPDSTRB")
group1Index <- which(names(rbd)%in%group1)
group2 <- c("STROKE", "HETRA", "PARKISM", "RLS", "NARCLPSY", "DEPRS", 
            "EPILEPSY", "BRNINFM", "CNSOTH")
group2Index <- which(names(rbd)%in%group2)
for (i in 1:nrow(rbd)) {
  id <- rbd[i,]$PATNO
  time <- rbd[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  if (is.na(sum(rbd[i,group1Index], rbd[i,group2Index]))) {
    next #any missing value leads to RBD score missing
  }
  score <- length(which(rbd[i,group1Index]==1))
  if (length(which(rbd[i,group2Index]==1))!=0) {
    score <- score+1
  }
  ppmi[index,]$REMSBD <- score
}
cat("RBD has been updated!")

###SCOPA-AUT Total Score(SCOPA)###
group1 <- paste("SCAU", 1:21, sep="")
group1Index <- which(names(scopa)%in%group1)
group2 <- paste("SCAU", 22:25, sep="")
group2Index <- which(names(scopa)%in%group2)
for (i in 1:nrow(scopa)) {
  id <- scopa[i,]$PATNO
  time <- scopa[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  score <- 0
  group1Score <- scopa[i,group1Index]
  group1_9 <- which(group1Score==9)
  if (length(group1_9)!=0) {
    score <- score+sum(group1Score[-group1_9])+3*length(group1_9)
  }
  else {
    score <- score+sum(group1Score)
  }
  group2Score <- scopa[i,group2Index]
  group2_9 <- which(group2Score==9)
  if (length(group2_9)!=0) {
    score <- score+sum(group2Score[-group2_9])
  }
  else {
    score <- score+sum(group2Score)
  }
  ppmi[index,]$SCOPA_AUT <- score
}
cat("SCOPA-AUT Total Score hae been updated!")

###Semantic Fluency(SFT)###
sft$SUM <- sft$VLTANIM+sft$VLTVEG+sft$VLTFRUIT
for (i in 1:nrow(sft)) {
  id <- sft[i,]$PATNO
  time <- sft[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[i,]$SFT <- sft[i,]$SUM
}
cat("Semantic Fluency has been updated!")

###State Trait Anxiety Score and Subscore(STAI)###
group1AName <- paste("STAIAD", c(3, 4, 6, 7, 9, 12, 13, 14, 17, 18), sep="")
group1AIndex <- which(names(stai)%in%group1AName)
group1RName <- paste("STAIAD", c(1:20)[-c(3, 4, 6, 7, 9, 12, 13, 14, 17, 18)], sep="")
group1RIndex <- which(names(stai)%in%group1RName)
group2AName <- paste("STAIAD", c(22, 24, 25, 28, 29, 31, 32, 35, 37, 38, 40), sep="")
group2AIndex <- which(names(stai)%in%group2AName)
group2RName <- paste("STAIAD", c(21:40)[-c(22, 24, 25, 28, 29, 31, 32, 35, 37, 38, 40)], sep="")
group2RIndex <- which(names(stai)%in%group2RName)
for (i in 1:nrow(stai)) {
  id <- stai[i,]$PATNO
  time <- stai[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  score1 <- sum(stai[i,group1AIndex])+sum(5-stai[i,group1RIndex])
  score2 <- sum(stai[i,group1AIndex])+sum(5-stai[i,group2RIndex])
  score <- score1+score2
  ppmi[index,]$STAT <- score
  ppmi[index,]$STASS <- score1
  ppmi[index,]$STATS <- score2
}
cat("State Trait Anxiety Score and Sub-score have been updated!")

###UPSIT Raw Score(UPSIT)###
upsit$SUM <- upsit$UPSITBK1+upsit$UPSITBK2+upsit$UPSITBK3+upsit$UPSITBK4
for (i in 1:nrow(upsit)) {
  id <- upsit[i,]$PATNO
  time <- upsit[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[index,]$UPSITRS <- upsit[i,]$SUM
}
cat("UPSIT Raw Score hae been updated!")

###Mild Coginitive Impairment(COGCATG)###
#1 -> MCI, 0 -> no MCI
for (i in 1:nrow(cogcatg)) {
  id <- cogcatg[i,]$PATNO
  time <- cogcatg[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  rule1 <- cogcatg[i,]$COGDECLN
  #missing rule 2
  rule3 <- cogcatg[i,]$FNCDTCOG
  if (rule1==1 & rule3==0) {
    ppmi[index,]$MCI <- 1
  }
  else {
    ppmi[index,]$MCI <- 0
  }
}
cat("MCI has been updated!")

###Contralateral, Lpsilateral(SBR, PDFET)###
sbrName <- c("PATNO", "EVENT_ID", "CAUDATE_R", "CAUDATE_L", "PUTAMEN_R", "PUTAMEN_L")
sbrIndex <- which(names(sbr)%in%sbrName)
sbrData <- sbr[,sbrIndex]
featName <- c("PATNO", "EVENT_ID", "DOMSIDE")
featIndex <- which(names(feature)%in%featName)
featData <- feature[,featIndex]
temp <- merge(sbrData, featData)
for (i in 1:nrow(temp)) {
  id <- temp[i,]$PATNO
  time <- temp[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  health <- ppmi[index,]$DIAG
  left <- temp[i,]$CAUDATE_L+temp[i,]$PUTAMEN_L
  right <- temp[i,]$CAUDATE_R+temp[i,]$PUTAMEN_R
  domside <- temp[i,]$DOMSIDE
  if (is.na(health) | is.na(left) | is.na(right) | is.na(domside)) {
    next #cannot decide value without all information
  }
  if (health==1) {
    ppmi[index,]$CONTRALATERAL <- max(left, right)
    ppmi[index,]$IPSILATERAL <- min(left, right)
  }
  else {
    if (domside==2) {
      ppmi[index,]$CONTRALATERAL <- left
      ppmi[index,]$IPSILATERAL <- right
    }
    else if (domside==1) {
      ppmi[index,]$CONTRALATERAL <- right
      ppmi[index,]$IPSILATERAL <- left
    }
    else if (domside==3) {
      ppmi[index,]$CONTRALATERAL <- max(left, right)
      ppmi[index,]$IPSILATERAL <- min(left, right)
    }
    else {
      next
    }
  }
}
cat("Contralateral and Ipsilateral have been updated!")

###Other DaTSCAN values(SBR)###
for (i in 1:nrow(sbr)) {
  id <- sbr[i,]$PATNO
  time <- sbr[i,]$EVENT_ID
  index <- which(ppmi$PATNO==id & ppmi$TIME==time)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  ppmi[i,]$CAUDATE <- (sbr[i,]$CAUDATE_L+sbr[i,]$CAUDATE_R)/2
  ppmi[i,]$PUTAMEN <- (sbr[i,]$PUTAMEN_L+sbr[i,]$PUTAMEN_R)/2
  ppmi[i,]$STRIATUM <- (sbr[i,]$CAUDATE_L+sbr[i,]$CAUDATE_R+sbr[i,]$PUTAMEN_L+sbr[i,]$PUTAMEN_R)/4
  ppmi[i,]$CDR <- (sbr[i,]$CAUDATE_L+sbr[i,]$CAUDATE_R)/(sbr[i,]$PUTAMEN_L+sbr[i,]$PUTAMEN_R)
  ppmi[i,]$AI <- abs((sbr[i,]$CAUDATE_L+sbr[i,]$PUTAMEN_L-sbr[i,]$CAUDATE_R-sbr[i,]$PUTAMEN_R)/mean(sbr[i,]$CAUDATE_L+sbr[i,]$CAUDATE_R+sbr[i,]$PUTAMEN_L+sbr[i,]$PUTAMEN_R)*100)
}
cat("Caudate, Putamen, Striatum, CDR, and AI have been updated!")

###CSF(BIOSPECAN)###
for (i in 1:nrow(bio)) {
  id <- bio[i,]$PATNO
  index <- which(ppmi$PATNO==id)
  if (length(index)==0) {
    next #this is not an enrolled patient
  }
  if (is.na(bio[i,]$TESTNAME)) {
    next #cannot decide what test
  }
  if (bio[i,]$TESTNAME=="Abeta 42") {
    ppmi[index,]$ABETA <- as.numeric(bio[i,]$TESTVALUE)
  }
  else if (bio[i,]$TESTNAME=="Total tau") {
    ppmi[index,]$TAU <- as.numeric(bio[i,]$TESTVALUE)
  }
  else if (bio[i,]$TESTNAME=="p-Tau 181P") {
    ppmi[index,]$PTAU <- as.numeric(bio[i,]$TESTVALUE)
  }
  else if (bio[i,]$TESTNAME=="CSF Alpha-synuclein") {
    ppmi[index,]$ALPHA <- as.numeric(bio[i,]$TESTVALUE)
  }
  else {
    next
  }
}
ppmi$TAU_ABETA <- ppmi$TAU/ppmi$ABETA
ppmi$PTAU_ABETA <- ppmi$PTAU/ppmi$ABETA
ppmi$PTAU_TAU <- ppmi$PTAU/ppmi$TAU
cat("CSFs have been updated!")

###Urate(COVANCE)###
#Cannot find COVANCE

###Output###
write.csv(x=ppmi, file="../PPMI_Data.csv", row.names=F)
