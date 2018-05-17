##########################################################################

#                             Clean up Mort Data  


##########################################################################

setwd("C:/Users/allison/Documents/Project/WolfData")

mort1 <- read.csv("C:/Users/allison/Documents/Project/WolfData/Complete MT_Wolf DATA/MT_MORTALITY.csv")
mort2 <- read.csv("C:/Users/allison/Documents/Project/WolfData/Complete MT_Wolf DATA/MT_MTWOLFMORTALITY.csv")
mort3 <- read.csv("C:/Users/allison/Documents/Project/WolfData/Complete MT_Wolf DATA/WesternMT_MTWOLFMORTALITY.csv")
mort4 <- read.csv("C:/Users/allison/Documents/Project/WolfData/Complete MT_Wolf DATA/WesternMT_MTWOLFMORTALITY2009.csv")
mort5 <- read.csv("C:/Users/allison/Documents/Project/WolfData/Complete MT_Wolf DATA/WesternMT_WolfMortality.csv")
mort6 <- read.csv("C:/Users/allison/Documents/Project/WolfData/Complete MT_Wolf DATA/Wolf MRRE data 8.19.2016 to UM.csv")


colnames(mort1)
colnames(mort2)
colnames(mort3)
colnames(mort4)
colnames(mort5)


mort1 <-mort1[,-c(8,15,16,21,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42,43,44,48,49)]
colnames(mort1)

mort2 <-mort2[,-c(6,16,17,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,44,45)]
colnames(mort2)

mort3 <-mort3[,-c(6,16,17,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43)]
colnames(mort3)

mort4 <-mort4[,-c(7,16,17,18,26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44)]
colnames(mort4)

mort5 <-mort5[,-c(5,14,15,17,23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41)]
colnames(mort5)


colnames(mort1)
colnames(mort2)
colnames(mort3)
colnames(mort4)
colnames(mort5)

colnames(mort2)[4] <- "WOLFNO" 

colnames(mort3)[2] <- "RECOVERY_AREA" 
colnames(mort3)[4] <- "WOLFNO" 
colnames(mort3)[6] <- "PACK_AT_TIME_OF_DEATH" 
colnames(mort3)[7] <- "AGE_CLASS"
colnames(mort3)[9] <- "SOCIAL_RANK" 
colnames(mort3)[10] <- "BREEDING_STATUS" 
colnames(mort3)[14] <- "DATE_OF_DISCOVERY" 
colnames(mort3)[16] <- "UTM_E" 
colnames(mort3)[17] <- "UTM_N" 

colnames(mort4)[3] <- "RECOVERY_AREA" 
colnames(mort4)[5] <- "WOLFNO" 
colnames(mort4)[7] <- "PACK_AT_TIME_OF_DEATH" 
colnames(mort4)[8] <- "AGE_CLASS"
colnames(mort4)[10] <- "SOCIAL_RANK" 
colnames(mort4)[11] <- "BREEDING_STATUS" 
colnames(mort4)[16] <- "UTM_E" 
colnames(mort4)[17] <- "UTM_N"

colnames(mort5)[1] <- "RECOVERY_AREA" 
colnames(mort5)[2] <- "FREQ"
colnames(mort5)[3] <- "WOLFNO"
colnames(mort5)[4] <- "SEX"
colnames(mort5)[5] <- "PACK_AT_TIME_OF_DEATH" 
colnames(mort5)[6] <- "AGE_CLASS"
colnames(mort5)[7] <- "AGE"
colnames(mort5)[8] <- "SOCIAL_RANK"
colnames(mort5)[9] <- "YEAR"
colnames(mort5)[10] <- "DATE_LAST_ALIVE" 
colnames(mort5)[12] <- "DATE_OF_DEATH" 
colnames(mort5)[14] <- "UTM_E" 
colnames(mort5)[15] <- "UTM_N"
colnames(mort5)[16] <- "DATUM"
colnames(mort5)[20] <- "COMMENTS"


str(mort1)
mort1$LATITUDE<-as.factor(mort1$LATITUDE)
mort1$LONGITUDE<-as.factor(mort1$LONGITUDE)
str(mort2)
str(mort3)
str(mort4)
mort4$LATITUDE<-as.factor(mort4$LATITUDE)
mort4$LONGITUDE<-as.factor(mort4$LONGITUDE)
str(mort5)
mort5$LATITUDE<-as.factor(mort5$LATITUDE)
mort5$LONGITUDE<-as.factor(mort5$LONGITUDE)

library(dplyr)
library(tidyr)

totalmort <- bind_rows(mort1, mort2, mort3, mort4, mort5)
write.csv(totalmort, "totalwolfmort1.csv")


totalmort <- read.csv("C:/Users/allison/Documents/Project/WolfData/totalwolfmort1.csv")
colnames(totalmort)


table(totalmort$RECOVERY_AREA)
totalmort$RECOVERY_AREA[totalmort$RECOVERY_AREA=="Canada"] <- "CANADA"
totalmort$RECOVERY_AREA[totalmort$RECOVERY_AREA=="GYE/YNP"] <- "GYE"
totalmort$RECOVERY_AREA[totalmort$RECOVERY_AREA=="CID/ID?"] <- "CID"
totalmort$RECOVERY_AREA[totalmort$RECOVERY_AREA=="NWMT IDAHO"] <- "NWMT"
totalmort$RECOVERY_AREA <- factor(totalmort$RECOVERY_AREA)
table(totalmort$RECOVERY_AREA)


table(totalmort$YEAR)
table(totalmort$DATE_OF_DEATH) # FORTMAT IN EXCEL
table(totalmort$DATE_LAST_ALIVE) # FORTMAT IN EXCEL
table(totalmort$FREQ)
totalmort$FREQ[totalmort$FREQ == "Freq?"] <- "UNKNOWN"
totalmort$FREQ[totalmort$FREQ == "YES"] <- "UNKNOWN"
totalmort$FREQ[totalmort$FREQ == "GPS"] <- "UNKNOWN"
totalmort$FREQ[totalmort$FREQ == "WDFWCollar"] <- "UNKNOWN"
totalmort$FREQ[totalmort$FREQ == "S/N 646105"] <- "UNKNOWN"
totalmort$FREQ[totalmort$FREQ == "??"] <- "UNKNOWN"
totalmort$FREQ[totalmort$FREQ == "NO COLLAR"] <- ""
totalmort$FREQ[totalmort$FREQ == "Uncollared"] <- ""
totalmort$FREQ <- factor(totalmort$FREQ)

table(totalmort$SEX)
totalmort$SEX[totalmort$SEX == "?"] <- ""
totalmort$SEX[totalmort$SEX == "u"] <- ""
totalmort$SEX[totalmort$SEX == "U"] <- ""
totalmort$SEX[totalmort$SEX == "unk"] <- ""
totalmort$SEX[totalmort$SEX == "UNK"] <- ""
totalmort$SEX[totalmort$SEX == "UNKNOWN"] <- ""
totalmort$SEX[totalmort$SEX == "F"] <- "FEMALE"
totalmort$SEX[totalmort$SEX == "M"] <- "MALE"
totalmort$SEX <- factor(totalmort$SEX)


table(totalmort$PACK_AT_TIME_OF_DEATH) #Clean up in Excel#
table(totalmort$AGE_CLASS)
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "pup"] <- "PUP"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Pup"] <- "PUP"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Pup?"] <- "PUP"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Adult"] <- "ADULT"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Adult?"] <- "ADULT"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Yearling?"] <- "YEARLING"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Yearling"] <- "YEARLING"
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "u"] <- ""
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "unk"] <- ""
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Unk"] <- ""
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "Unknown"] <- ""
totalmort$AGE_CLASS[totalmort$AGE_CLASS == "UNKNOWN"] <- ""
totalmort$AGE_CLASS <- factor(totalmort$AGE_CLASS)

table(totalmort$AGE) ### FIX THIS MORE ; STANDARDIZE WORDS IN EXCEL##
totalmort$AGE[totalmort$AGE == "6 months"] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "6 Months"] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "6 mos"] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "6 MOS"] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "6 mos."] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "6mos"] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "6MOS"] <- "6 MONTHS"
totalmort$AGE[totalmort$AGE == "?"] <- ""
totalmort$AGE[totalmort$AGE == "u"] <- ""
totalmort$AGE[totalmort$AGE == "U"] <- ""
totalmort$AGE[totalmort$AGE == "unk"] <- ""
totalmort$AGE[totalmort$AGE == "Unk"] <- ""
totalmort$AGE[totalmort$AGE == "unknown"] <- ""
totalmort$AGE[totalmort$AGE == "Unknown"] <- ""
totalmort$AGE[totalmort$AGE == "UNKNOWN"] <- ""
totalmort$AGE[totalmort$AGE == "10 months"] <- "10 MONTHS"
totalmort$AGE[totalmort$AGE == "10 mos"] <- "10 MONTHS"
totalmort$AGE[totalmort$AGE == "10 MOS"] <- "10 MONTHS"
totalmort$AGE[totalmort$AGE == "8 months"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "8 MOSNTHS"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "8 mos"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "8 MOS"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "8mos"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "8MONTHS"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "11 months"] <- "11 MOS"
totalmort$AGE[totalmort$AGE == "11 mos"] <- "11 MOS"
totalmort$AGE[totalmort$AGE == "9 months"] <- "9 MONTHS"
totalmort$AGE[totalmort$AGE == "9 mos"] <- "9 MONTHS"
totalmort$AGE[totalmort$AGE == "9 MOS"] <- "9 MONTHS"
totalmort$AGE[totalmort$AGE == "9 mos."] <- "9 MONTHS"
totalmort$AGE[totalmort$AGE == "9mos"] <- "9 MONTHS"
totalmort$AGE[totalmort$AGE == "7 months"] <- "7 MONTHS"
totalmort$AGE[totalmort$AGE == "7 mos"] <- "7 MONTHS"
totalmort$AGE[totalmort$AGE == "7 MOS"] <- "7 MONTHS"
totalmort$AGE[totalmort$AGE == "7mos"] <- "7 MONTHS"
totalmort$AGE[totalmort$AGE == "4 months"] <- "4 MONTHS"
totalmort$AGE[totalmort$AGE == "4 mos"] <- "4 MONTHS"
totalmort$AGE[totalmort$AGE == "4 MOS"] <- "4 MONTHS"
totalmort$AGE[totalmort$AGE == "4mos"] <- "4 MONTHS"
totalmort$AGE[totalmort$AGE == "3 months"] <- "3 MONTHS"
totalmort$AGE[totalmort$AGE == "3 mos"] <- "3 MONTHS"
totalmort$AGE[totalmort$AGE == "3 MOS"] <- "3 MONTHS"
totalmort$AGE[totalmort$AGE == "3.5 mos"] <- "3 MONTHS"
totalmort$AGE[totalmort$AGE == "5 months"] <- "5 MONTHS"
totalmort$AGE[totalmort$AGE == "5 mos"] <- "5 MONTHS"
totalmort$AGE[totalmort$AGE == "5 MOS"] <- "5 MONTHS"
totalmort$AGE[totalmort$AGE == "5 MTHS"] <- "5 MONTHS"
totalmort$AGE[totalmort$AGE == "1 year"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "1 yr"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "1 Yr"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "1yr"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "1YR"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "1"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "2 years"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2 yr"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2 yrs"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2 Yrs"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2 YRS"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2yr"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2yrs"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "2YR"] <- "2 YR"
totalmort$AGE[totalmort$AGE == "3 years"] <- "3 YR"
totalmort$AGE[totalmort$AGE == "3 yrs"] <- "3 YR"
totalmort$AGE[totalmort$AGE == "3 YRS"] <- "3 YR"
totalmort$AGE[totalmort$AGE == "3"] <- "3 YR"
totalmort$AGE[totalmort$AGE == "5 yrs"] <- "5 YR"
totalmort$AGE[totalmort$AGE == "5 Yrs"] <- "5 YR"
totalmort$AGE[totalmort$AGE == "5 YRS"] <- "5 YR"
totalmort$AGE[totalmort$AGE == "5"] <- "5 YR"
totalmort$AGE[totalmort$AGE == "5yrs"] <- "5 YR"
totalmort$AGE[totalmort$AGE == "5 YR"] <- "5 YR"
totalmort$AGE[totalmort$AGE == "6 yrs"] <- "6 YR"
totalmort$AGE[totalmort$AGE == "6 YRS"] <- "6 YR"
totalmort$AGE[totalmort$AGE == "6"] <- "6 YR"
totalmort$AGE[totalmort$AGE == "6yr"] <- "6 YR"
totalmort$AGE[totalmort$AGE == "6 YR"] <- "6 YR"
totalmort$AGE[totalmort$AGE == "7 YRS"] <- "7"
totalmort$AGE[totalmort$AGE == "8 YRS"] <- "8"
totalmort$AGE[totalmort$AGE == "8 yrs"] <- "8"
totalmort$AGE[totalmort$AGE == "3 +"] <- "3+"
totalmort$AGE[totalmort$AGE == "3+ yrs"] <- "3+"
totalmort$AGE[totalmort$AGE == "4-5 yrs"] <- "4-5 YR"
totalmort$AGE[totalmort$AGE == "4-5 YRS"] <- "4-5 YR"
totalmort$AGE[totalmort$AGE == "4 yr"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "4"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "4 yrs"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "4 YRS"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "4yrs"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "5 +"] <- "5+"
totalmort$AGE[totalmort$AGE == "5+ years"] <- "5+"
totalmort$AGE[totalmort$AGE == "5 + years"] <- "5+"
totalmort$AGE[totalmort$AGE == "5+ YR"] <- "5+"
totalmort$AGE[totalmort$AGE == "5+ yrs"] <- "5+"
totalmort$AGE[totalmort$AGE == "6 +"] <- "6+"
totalmort$AGE[totalmort$AGE == "6+YR"] <- "6+"
totalmort$AGE[totalmort$AGE == "6+ YRS"] <- "6+"
totalmort$AGE[totalmort$AGE == "8+ yrs"] <- "8+"
totalmort$AGE[totalmort$AGE == "8+ YRS"] <- "8+"
totalmort$AGE[totalmort$AGE == "9"] <- "9 YR"
totalmort$AGE[totalmort$AGE == "5-6 yrs"] <- "5-6 YR"
totalmort$AGE[totalmort$AGE == "5-6 YRS"] <- "5-6 YR"
totalmort$AGE[totalmort$AGE == "30 MOS"] <- "30 months"
totalmort$AGE[totalmort$AGE == "3-4 yrs"] <- "3-4 YR"
totalmort$AGE[totalmort$AGE == "3/4 yrs"] <- "3-4 YR"
totalmort$AGE[totalmort$AGE == "3/4 YRS EST"] <- "3-4 YR"
totalmort$AGE[totalmort$AGE == "23 MOS"] <- "23 Months"
totalmort$AGE[totalmort$AGE == "23 months"] <- "23 Months"
totalmort$AGE[totalmort$AGE == "22 months?"] <- "22 months"
totalmort$AGE[totalmort$AGE == "2 mos"] <- "2 MOS"
totalmort$AGE[totalmort$AGE == "2-3 years"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2-3 yr"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2-3 yrs"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2-3 YRS"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2-3?"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "16 MOS"] <- "16 MONTHS"
totalmort$AGE[totalmort$AGE == "16 months"] <- "16 MONTHS"
totalmort$AGE[totalmort$AGE == "15 MOS"] <- "15 months"
totalmort$AGE[totalmort$AGE == "14 mos"] <- "14 months"
totalmort$AGE[totalmort$AGE == "13 mos"] <- "13 MONTHS"
totalmort$AGE[totalmort$AGE == "13 months"] <- "13 MONTHS"
totalmort$AGE[totalmort$AGE == "12 MOS"] <- "12 months"
totalmort$AGE[totalmort$AGE == "11 yrs"] <- "11 YR"
totalmort$AGE[totalmort$AGE == "1-2 yrs"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "1-2 YRS"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "~ 3 yrs"] <- "~3 yrs"
totalmort$AGE[totalmort$AGE == "~3 yrs"] <- "3 YR"
totalmort$AGE[totalmort$AGE == "~10 yrs"] <- "10+ yrs"
totalmort$AGE[totalmort$AGE == "~ 4 YRS"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "-4"] <- "4 YR"
totalmort$AGE[totalmort$AGE == "<2Yrs"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "1.1"] <- "1.1 yrs"
totalmort$AGE[totalmort$AGE == "1.76"] <- "1.7 yrs"
totalmort$AGE[totalmort$AGE == "1.83"] <- "1.8"
totalmort$AGE[totalmort$AGE == "1.84"] <- "1.8"
totalmort$AGE[totalmort$AGE == "1.9 yrs?"] <- "1.9"
totalmort$AGE[totalmort$AGE == "1.92"] <- "1.9"
totalmort$AGE[totalmort$AGE == "1.95"] <- "1.9"
totalmort$AGE[totalmort$AGE == "0.44"] <- "5 MONTHS"
totalmort$AGE[totalmort$AGE == "0.64"] <- "8 MONTHS"
totalmort$AGE[totalmort$AGE == "0.8"] <- "10 MONTHS"
totalmort$AGE[totalmort$AGE == "0.83"] <- "10 MONTHS"
totalmort$AGE[totalmort$AGE == "0.9"] <- "11 MOS"
totalmort$AGE[totalmort$AGE == "12 months"] <- "1 YR"
totalmort$AGE[totalmort$AGE == "1.1 yrs"] <- "13 MONTHS"
totalmort$AGE[totalmort$AGE == "1.3 yrs"] <- "16 MONTHS"
totalmort$AGE[totalmort$AGE == "1.5"] <- "18 MONTHS"
totalmort$AGE[totalmort$AGE == "1.8"] <- "22 months"
totalmort$AGE[totalmort$AGE == "1.9 yrs"] <- "23 Months"
totalmort$AGE[totalmort$AGE == "2.46"] <- "30 months"
totalmort$AGE[totalmort$AGE == "2.5"] <- "30 months"
totalmort$AGE[totalmort$AGE == "2.87"] <- "34 months"
totalmort$AGE[totalmort$AGE == "2.91"] <- "35 month est"
totalmort$AGE[totalmort$AGE == "Adult"] <- ""
totalmort$AGE[totalmort$AGE == "MOS"] <- ""

totalmort$AGE[totalmort$AGE == "1.4 yrs"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "1.6 yrs"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "1.7 yrs"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "1.9"] <- "1-2 YR"
totalmort$AGE[totalmort$AGE == "2.17"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2.3"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2.34"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2.57"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2.6"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "2.76"] <- "2-3 YR"
totalmort$AGE[totalmort$AGE == "3.47"] <- "3-4 YR"
totalmort$AGE[totalmort$AGE == "3.5"] <- "3-4 YR"
totalmort$AGE[totalmort$AGE == "4.61"] <- "4-5 YR"
totalmort$AGE[totalmort$AGE == "4.86"] <- "4-5 YR"
totalmort$AGE[totalmort$AGE == "5.1"] <- "5-6 YR"
totalmort$AGE[totalmort$AGE == "5.23"] <- "5-6 YR"
totalmort$AGE[totalmort$AGE == "6.5"] <- "6-7 yrs"
totalmort$AGE <-factor(totalmort$AGE)

table(totalmort$SOCIAL_RANK)
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "Alpha"] <- "ALPHA"
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "Alpha?"] <- "ALPHA"
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "Breeder?"] <- "Breeder"
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "Poss breeder"] <- "PAST BREEDER"
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "Subordinate"] <- "SUBORDINATE"
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "UKNOWN"] <- ""
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "UNKNOWN"] <- ""
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "unknown"] <- ""
totalmort$SOCIAL_RANK[totalmort$SOCIAL_RANK == "Unknown"] <- ""
totalmort$SOCIAL_RANK <- factor(totalmort$SOCIAL_RANK)

table(totalmort$BREEDING_STATUS)
totalmort$BREEDING_STATUS[totalmort$BREEDING_STATUS == "UNKNOWN"] <- ""
totalmort$BREEDING_STATUS[totalmort$BREEDING_STATUS == "Unknown"] <- ""
totalmort$BREEDING_STATUS[totalmort$BREEDING_STATUS == "Breeder"] <- "BREEDER"
totalmort$BREEDING_STATUS[totalmort$BREEDING_STATUS == "breeder?"] <- "BREEDER"
totalmort$BREEDING_STATUS[totalmort$BREEDING_STATUS == "nonbreeder"] <- "NON-BREEDER"
totalmort$BREEDING_STATUS[totalmort$BREEDING_STATUS == "Nonbreeder"] <- "NON-BREEDER"
totalmort$BREEDING_STATUS <- factor(totalmort$BREEDING_STATUS)


table(totalmort$DATE_OF_DISCOVERY) ## FORTMAT LATER ##
table(totalmort$ZONE)
table(totalmort$UTM_E)
table(totalmort$UTM_N)
table(totalmort$DATUM)
totalmort$DATUM[totalmort$DATUM == "NAD83"] <- "NAD 83"
totalmort$DATUM <- factor(totalmort$DATUM)

table(totalmort$LATITUDE)
table(totalmort$LONGITUDE)
table(totalmort$COD)
totalmort$COD[totalmort$COD == "u"] <- ""
totalmort$COD[totalmort$COD == "unk"] <- ""
totalmort$COD[totalmort$COD == "Unkn"] <- ""
totalmort$COD[totalmort$COD == "Unknow"] <- ""
totalmort$COD[totalmort$COD == "unknown"] <- ""
totalmort$COD[totalmort$COD == "Unknown"] <- ""
totalmort$COD[totalmort$COD == "UNKNOWN"] <- ""
totalmort$COD[totalmort$COD == "VEHICLE"] <- "Vehicle"
totalmort$COD[totalmort$COD == "TRAIN"] <- "Train"
totalmort$COD[totalmort$COD == "suspect illegally shot"] <- "Illegal Shot"
totalmort$COD[totalmort$COD == "NATURAL"] <- "Natural"
totalmort$COD[totalmort$COD == "MT HARVEST"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "Legal Hunter Take"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "LEGAL HARVEST"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "LEGAL"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "INCIDENTAL"] <- "Incidental"
totalmort$COD[totalmort$COD == "incidental"] <- "Incidental"
totalmort$COD[totalmort$COD == "ILLEGAL"] <- "Illegal"
totalmort$COD[totalmort$COD == "illegal"] <- "Illegal"
totalmort$COD[totalmort$COD == "SELF DEFENSE"] <- "Self-defense"
totalmort$COD[totalmort$COD == "Shot-self defense claim"] <- "Self-defense"
totalmort$COD[totalmort$COD == "Claim self defense"] <- "Self-defense"
totalmort$COD[totalmort$COD == "CLAIM SELF DEFENSE"] <- "Self-defense"
totalmort$COD[totalmort$COD == "CONTROL"] <- "Control"
totalmort$COD[totalmort$COD == "HARVEST"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "harvest"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "Hunter"] <- "Legal Harvest"
totalmort$COD[totalmort$COD == "ELECTROCUTED"] <- "Electrocuted"
totalmort$COD[totalmort$COD == "Legal Trapper Take"] <- "Legally Trapped"
totalmort$COD <- factor(totalmort$COD)

table(totalmort$METHOD_OF_CONTROL)
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "AERIAL GUINNED"] <- "AERIAL GUNNED"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "AERIAL GUNNED PLANE"] <- "AERIAL GUNNED"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "ARIAL GUNNED"] <- "AERIAL GUNNED"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "EARIAL GUINNED"] <- "AERIAL GUNNED"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "call and shoot"] <- "Call and shoot"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "Call and Shoot"] <- "Call and shoot"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "SNARED"] <- "Snared"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "SNARE"] <- "Snared"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "trapped and shot"] <- "Trapped and shot"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "Trapped and Shot"] <- "Trapped and shot"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "trapped/shot"] <- "Trapped and shot"
totalmort$METHOD_OF_CONTROL[totalmort$METHOD_OF_CONTROL == "Shot from ground"] <- "GROUND SHOT"
totalmort$METHOD_OF_CONTROL <- factor(totalmort$METHOD_OF_CONTROL)

write.csv(totalmort, "totalwolfmort1.csv")




# Clean up MRRE database

colnames(mort6)
table(mort6$REFERENCENO)
table(mort6$LICYR)
table(mort6$SEX)
table(mort6$WMU)
table(mort6$SUBUNIT)
table(mort6$DEERELKDISTRICT)
table(mort6$COUNTY)
table(mort6$LOCATION)
table(mort6$HARVESTDATE)
table(mort6$HARVESTTIME)
table(mort6$REPORTEDDATE)
table(mort6$TYPEOFHARVEST)
table(mort6$AGECLASS)
table(mort6$TAGCOLLARTYPE)
table(mort6$TAGCOLLARDESC)
table(mort6$DISEASEPARASITESINJURIES)
table(mort6$STATUS)
table(mort6$WEAPONTYPE)
table(mort6$ESTDISTANCEYDS)
table(mort6$LANDOWNERTYPE)
table(mort6$CEM_AGE)
table(mort6$TOWNSHIP)
table(mort6$WOLFPGM2009PACK)
table(mort6$WOLFPGM2009ANIMALID)
table(mort6$WOLFID)
table(mort6$RESIDENCYSTATUS)



##########################################################################

#                             Clean up Capture Data  


##########################################################################

setwd("C:/Users/allison/Documents/Project/WolfData")

cap1 <- read.csv("C:/Users/allison/Documents/Project/WolfData/MT CAPTURES.csv")
cap2 <- read.csv("C:/Users/allison/Documents/Project/WolfData/MT_CAPTURE.csv")
cap3 <- read.csv("C:/Users/allison/Documents/Project/WolfData/WesternMT_Captures.csv")


colnames(cap1)
colnames(cap2)
colnames(cap3)

colnames(cap1)[2] <- "AREA"
colnames(cap1)[5] <- "PACK"
colnames(cap1)[7] <- "WOLFNO"
colnames(cap1)[14] <- "AGE_EoA"
colnames(cap1)[16] <- "WEIGHT"
colnames(cap1)[17] <- "WEIGHT_EoA"

colnames(cap2)[58] <- "EARTAG"
colnames(cap2)[18] <- "WEIGHT_EoA"
colnames(cap2)[24] <- "UTM_EASTING"
colnames(cap2)[25] <- "UTM_NORTHING"
colnames(cap2)[20] <- "TRAP_TYPE"
colnames(cap2)[36] <- "PIT_TAG_ID"
colnames(cap2)[47] <- "LAB_SAMPLE_NO"

colnames(cap3) <- c("AREA", "RECAPTURE", "PACK", "YEAR", "DATE", "WOLFNO",
                    "EARTAG", "OBSERVER", "MANGE", "SEX", "COLOR", "AGE",
                    "AGE_CLASS", "WEIGHT", "BREEDER", "LOCATION_OF_CAPTURE",
                    "UTM_EASTING", "UTM_NORTHING", "DATUM", "LATITUDE", 
                    "LONGITUDE", "METHOD_OF_CAPTURE", "TRAP_TYPE", 
                    "OTHER_TRAP", "AMBIENT_TEMP", "WEATHER_CONDITIONS", 
                    "METHOD_OF_DRUG_ADMINISTRATION", "COLLAR_FREQUENCY",
                    "SERIAL_NO", "PIT_TAG_ID", "CANINES", "VULVA", 
                    "INGUINAL_TEATS", "TESTES", "SET_TYPE", "LURE", 
                    "TIME_TRANS_SET", "TIME_IN_TRAP","FOOT_PIC",
                    "TRAP_INJURIES")

library(dplyr)
library(tidyr)

str(cap1)
str(cap2)
str(cap3)

cap2$YEAR <- as.factor(cap2$YEAR)
cap3$YEAR <- as.factor(cap3$YEAR)
cap2$ZONE <- as.factor(cap2$ZONE)
cap2$LAB_SAMPLE_NO <- as.factor(cap2$LAB_SAMPLE_NO)
cap3$MANGE <- as.factor(cap3$MANGE)
cap3$WEIGHT <- as.factor(cap3$WEIGHT)

totalcap <- bind_rows(cap1, cap2, cap3)

colnames(totalcap)

table(totalcap$AREA)
table(totalcap$DATE) # FORMAT DATE IN EXCEL
table(totalcap$RECAPTURE)
totalcap$RECAPTURE[totalcap$RECAPTURE=="Yes"] <- "YES"
totalcap$RECAPTURE[totalcap$RECAPTURE=="No"] <- "NO"
totalcap$RECAPTURE <- factor(totalcap$RECAPTURE)

table(totalcap$PACK) # FIX IN EXCEL
table(totalcap$YEAR)
table(totalcap$WOLFNO)
table(totalcap$AGENCY)
totalcap$AGENCY[totalcap$AGENCY == "NEZPERCE"] <- "NEZ PERCE"
totalcap$AGENCY[totalcap$AGENCY == "NPT"] <- "NEZ PERCE"
totalcap$AGENCY[totalcap$AGENCY == "UNKN"] <- ""
totalcap$AGENCY[totalcap$AGENCY == "UNKNOWN"] <- ""
totalcap$AGENCY[totalcap$AGENCY == "Wildife Services"] <- "WILDLIFE SERVICES"
totalcap$AGENCY[totalcap$AGENCY == "Wildlife Services"] <- "WILDLIFE SERVICES"
totalcap$AGENCY[totalcap$AGENCY == "WS"] <- "WILDLIFE SERVICES"
totalcap$AGENCY <- factor(totalcap$AGENCY)

table(totalcap$SEX)
totalcap$SEX[totalcap$SEX == "F"] <- "FEMALE"
totalcap$SEX[totalcap$SEX == "M"] <- "MALE"
totalcap$SEX[totalcap$SEX == "UNKNOWN"] <- ""
totalcap$SEX <- factor(totalcap$SEX)

table(totalcap$COLOR)
totalcap$COLOR[totalcap$COLOR == "Black"] <- "BLACK"
totalcap$COLOR[totalcap$COLOR == "Drk Gray"] <- "DARK GRAY"
totalcap$COLOR[totalcap$COLOR == "Gray"] <- "GRAY"
totalcap$COLOR[totalcap$COLOR == "UNKNOWN"] <- ""
totalcap$COLOR <- factor(totalcap$COLOR)

table(totalcap$AGE) ## COME BACK TO THIS ONE
table(totalcap$AGE_EoA)
table(totalcap$AGE_CLASS)
totalcap$AGE_CLASS[totalcap$AGE_CLASS == "Adult"] <- "ADULT"
totalcap$AGE_CLASS[totalcap$AGE_CLASS == "Pup"] <- "PUP"
totalcap$AGE_CLASS[totalcap$AGE_CLASS == "Yearling"] <- "YEARLING"
totalcap$AGE_CLASS[totalcap$AGE_CLASS == "UNKNOWN"] <- ""
totalcap$AGE_CLASS <- factor(totalcap$AGE_CLASS)

table(totalcap$WEIGHT)
table(totalcap$WEIGHT_EoA)
totalcap$WEIGHT_EoA[totalcap$WEIGHT_EoA == "A"] <- "ACTUAL"
totalcap$WEIGHT_EoA[totalcap$WEIGHT_EoA == "E"] <- "ESTIMATED"
totalcap$WEIGHT_EoA[totalcap$WEIGHT_EoA == "unknown"] <- ""
totalcap$WEIGHT_EoA <- factor(totalcap$WEIGHT_EoA)

table(totalcap$BREEDER)
totalcap$BREEDER[totalcap$BREEDER == "No"] <- "NO"
totalcap$BREEDER[totalcap$BREEDER == "NO COLLAR"] <- "NO"
totalcap$BREEDER[totalcap$BREEDER == "UNK"] <- ""
totalcap$BREEDER[totalcap$BREEDER == "UNKN"] <- ""
totalcap$BREEDER[totalcap$BREEDER == "UNKNOWN"] <- ""
totalcap$BREEDER <- factor(totalcap$BREEDER)

table(totalcap$COLLAR_FREQUENCY)
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "no collar"] <- "NO COLLAR"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "n/a"] <- "NO COLLAR"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "N/A"] <- "NO COLLAR"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "No Collar"] <- "NO COLLAR"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "none"] <- "NO COLLAR"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "NONE"] <- "NO COLLAR"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "YES"] <- "COLLARED"
totalcap$COLLAR_FREQUENCY[totalcap$COLLAR_FREQUENCY == "UNKN"] <- "UNKNOWN"
totalcap$COLLAR_FREQUENCY <- factor(totalcap$COLLAR_FREQUENCY)

table(totalcap$LOCATION_OF_CAPTURE)
table(totalcap$ZONE)
table(totalcap$UTM_EASTING)
table(totalcap$UTM_NORTHING)
table(totalcap$DATUM)
totalcap$DATUM[totalcap$DATUM == "NAD27"] <- "NAD 27"
totalcap$DATUM[totalcap$DATUM == "NAD83"] <- "NAD 83"
totalcap$DATUM[totalcap$DATUM == "UNKNOWN"] <- ""
totalcap$DATUM <- factor(totalcap$DATUM)

table(totalcap$LATITUDE)
table(totalcap$LONGITUDE)
table(totalcap$METHOD_OF_CAPTURE)
totalcap$METHOD_OF_CAPTURE[totalcap$METHOD_OF_CAPTURE == "Aerial darted"] <- "AERIAL DARTED"
totalcap$METHOD_OF_CAPTURE[totalcap$METHOD_OF_CAPTURE == "GRIZZLY SNARE"] <- "BEAR SNARE"
totalcap$METHOD_OF_CAPTURE[totalcap$METHOD_OF_CAPTURE == "Net gunned"] <- "NET GUN"
totalcap$METHOD_OF_CAPTURE[totalcap$METHOD_OF_CAPTURE == "Trapped"] <- "TRAPPED"
totalcap$METHOD_OF_CAPTURE[totalcap$METHOD_OF_CAPTURE == "UNKN"] <- "UNKNOWN"
totalcap$METHOD_OF_CAPTURE <- factor(totalcap$METHOD_OF_CAPTURE)


table(totalcap$TRAP_TYPE)
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "McBride button"] <- "MCBRIDE BUTTON"
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "McBride rubber"] <- "MCBRIDE RUBBER"
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "Montana Special #3"] <- "MONTANA SPECIAL #3"
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "NEWHOUSE 114"] <- "NEWHOUSE #114"
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "NEWHOUSE 14"] <- "NEWHOUSE #14"
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "Other"] <- "OTHER"
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "u"] <- ""
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "unk"] <- ""
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "UNKN"] <- ""
totalcap$TRAP_TYPE[totalcap$TRAP_TYPE == "UNKNOWN"] <- ""
totalcap$TRAP_TYPE <- factor(totalcap$TRAP_TYPE)

write.csv(totalcap, "totalcapture1.csv")






##########################################################################

#                      Clean more data


##########################################################################


setwd("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData")
totalcap <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/2007to2013_CaptureData.csv")
totalmort <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/2007to2015_MortalityData.csv")
newcap <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/2014to2017_UMCaptureData.csv")



colnames(totalcap)

totalcap <- totalcap[,-c(9,13,16,30,31,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,52,53,54,55,57,58,59,61,62,63)]
colnames(totalcap)

table(totalcap$AREA)
table(totalcap$DATE)
table(totalcap$RECAPTURE)
table(totalcap$PACK)
table(totalcap$YEAR)
table(totalcap$WOLFNO)
table(totalcap$EARTAG)
totalcap$EARTAG[totalcap$EARTAG=="N/A"] <- "NONE"
totalcap$EARTAG[totalcap$EARTAG==""] <- "NONE"
totalcap$EARTAG <- factor(totalcap$EARTAG)

table(totalcap$AGENCY)
table(totalcap$SEX)
table(totalcap$COLOR)
table(totalcap$AGE)
totalcap$AGE[totalcap$AGE == "~6 yrs"] <- "~ 6 yrs"
totalcap$AGE[totalcap$AGE == "~5 yrs"] <- "~ 5 yrs"
totalcap$AGE[totalcap$AGE == "1 yr"] <- "1 YR"
totalcap$AGE[totalcap$AGE == "1yr"] <- "1 YR"
totalcap$AGE[totalcap$AGE == "1"] <- "1 YR"
totalcap$AGE[totalcap$AGE == "13 months"] <- "13 MOS"
totalcap$AGE[totalcap$AGE == "13 mos"] <- "13 MOS"
totalcap$AGE[totalcap$AGE == "14 months"] <- "14 MOS"
totalcap$AGE[totalcap$AGE == "15 months"] <- "15 MOS"
totalcap$AGE[totalcap$AGE == "2 yrs"] <- "2 YR"
totalcap$AGE[totalcap$AGE == "2 years"] <- "2 YR"
totalcap$AGE[totalcap$AGE == "2"] <- "2 YR"
totalcap$AGE[totalcap$AGE == "2YR"] <- "2 YR"
totalcap$AGE[totalcap$AGE == "2-3 yrs"] <- "2-3 YR"
totalcap$AGE[totalcap$AGE == "2 months"] <- "2 MOS"
totalcap$AGE[totalcap$AGE == "3"] <- "3 yrs"
totalcap$AGE[totalcap$AGE == "3 years"] <- "3 yrs"
totalcap$AGE[totalcap$AGE == "3-4 years"] <- "3-4 yrs"
totalcap$AGE[totalcap$AGE == "3 months"] <- "3 MOS"
totalcap$AGE[totalcap$AGE == "3 mos"] <- "3 MOS"
totalcap$AGE[totalcap$AGE == "4-5 yrs"] <- "4-5 YR"
totalcap$AGE[totalcap$AGE == "4 months"] <- "4 MOS"
totalcap$AGE[totalcap$AGE == "4 mos"] <- "4 MOS"
totalcap$AGE[totalcap$AGE == "4"] <- "4 YR"
totalcap$AGE[totalcap$AGE == "4 yrs"] <- "4 YR"
totalcap$AGE[totalcap$AGE == "5 months"] <- "5 MOS"
totalcap$AGE[totalcap$AGE == "5 mos"] <- "5 MOS"
totalcap$AGE[totalcap$AGE == "5"] <- "5 yrs"
totalcap$AGE[totalcap$AGE == "5-6 years"] <- "5-6 yr"
totalcap$AGE[totalcap$AGE == "5+"] <- "5+ YR"
totalcap$AGE[totalcap$AGE == "5+ yrs"] <- "5+ YR"
totalcap$AGE[totalcap$AGE == "6 months"] <- "6 MOS"
totalcap$AGE[totalcap$AGE == "6"] <- "6yr"
totalcap$AGE[totalcap$AGE == "7 mos"] <- "7 MOS"
totalcap$AGE[totalcap$AGE == "8 mos"] <- "8 MOS"
totalcap$AGE[totalcap$AGE == "unk"] <- ""
totalcap$AGE[totalcap$AGE == "UNKN"] <- ""
totalcap$AGE[totalcap$AGE == "UNKNOWN"] <- ""
totalcap$AGE <- factor(totalcap$AGE)

table(totalcap$AGE_CLASS)
table(totalcap$WEIGHT)
table(totalcap$BREEDER)
table(totalcap$COLLAR_FREQUENCY)
table(totalcap$LOCATION_OF_CAPTURE)
table(totalcap$ZONE)
table(totalcap$UTM_EASTING)
table(totalcap$UTM_NORTHING)
table(totalcap$DATUM)
table(totalcap$LATITUDE)
table(totalcap$LONGITUDE)
table(totalcap$METHOD_OF_CAPTURE)
table(totalcap$TRAP_TYPE)
table(totalcap$METHOD_OF_DRUG_ADMINISTRATION)
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "Dart"] <- "DART"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "Hand Inject"] <- "HAND INJECT"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "Jab Stick"] <- "JAB STICK"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "JAB STICK, HAND INJECT"] <- "JAB STICK"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "JABSTICK"] <- "JAB STICK"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "n/a"] <- "NONE"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "N/A"] <- "NONE"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "N/Z"] <- "NONE"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "none"] <- "NONE"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "No Drugs"] <- "NONE"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "NO DRUGS"] <- "NONE"
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "u"] <- ""
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "unk"] <- ""
totalcap$METHOD_OF_DRUG_ADMINISTRATION[totalcap$METHOD_OF_DRUG_ADMINISTRATION == "UNKNOWN"] <- ""
totalcap$METHOD_OF_DRUG_ADMINISTRATION <- factor(totalcap$METHOD_OF_DRUG_ADMINISTRATION)

table(totalcap$SERIAL_NO)
totalcap$SERIAL_NO[totalcap$SERIAL_NO == "u"] <- ""
totalcap$SERIAL_NO[totalcap$SERIAL_NO == "unk"] <- ""
totalcap$SERIAL_NO[totalcap$SERIAL_NO == "UNKNOWN"] <- ""
totalcap$SERIAL_NO[totalcap$SERIAL_NO == "n/a"] <- "N/A"
totalcap$SERIAL_NO[totalcap$SERIAL_NO == "none"] <- "N/A"
totalcap$SERIAL_NO <- factor(totalcap$SERIAL_NO)

table(totalcap$PIT_TAG_ID)
totalcap$PIT_TAG_ID[totalcap$PIT_TAG_ID == "n/a"] <- "NONE"
totalcap$PIT_TAG_ID[totalcap$PIT_TAG_ID == "N/A"] <- "NONE"
totalcap$PIT_TAG_ID[totalcap$PIT_TAG_ID == "NO"] <- "NONE"
totalcap$PIT_TAG_ID[totalcap$PIT_TAG_ID == "u"] <- ""
totalcap$PIT_TAG_ID[totalcap$PIT_TAG_ID == "UNKNOWN"] <- ""
totalcap$PIT_TAG_ID <- factor(totalcap$PIT_TAG_ID)

table(totalcap$COMMENTS)
table(totalcap$CAPTUREID)
table(totalcap$latlongdatum)
table(totalcap$TRAP_INJURIES)


colnames(newcap) <- c("ID", "PACK", "WOLFNO", "COLLARID", "FREQ", "DATE", "REMOVAL_DATE", "CAUSE_OF_REMOVAL", "AoE", "MOS_DEPLOYED", "SEX", "AGE_CLASS", "COMMENTS")

library(dplyr)
library(tidyr)

allcap <- bind_rows(newcap, totalcap)
write.csv(allcap, "allcaptures.csv")


colnames(totalmort)

table(totalmort$REFNO)
table(totalmort$WOLFNO)
table(totalmort$RECOVERY_AREA)
table(totalmort$YEAR)
table(totalmort$DATE_OF_DEATH)
table(totalmort$DATE_LAST_ALIVE)
table(totalmort$FREQ)
table(totalmort$SEX)
table(totalmort$PACK_AT_TIME_OF_DEATH)
table(totalmort$AGE_CLASS)
table(totalmort$AGE)
table(totalmort$SOCIAL_RANK)
table(totalmort$BREEDING_STATUS)
table(totalmort$DATE_OF_DISCOVERY)
table(totalmort$ZONE)
table(totalmort$UTM_E)
table(totalmort$UTM_N)
table(totalmort$DATUM)
table(totalmort$LATITUDE)
table(totalmort$LONGITUDE)
table(totalmort$COD)
table(totalmort$METHOD_OF_CONTROL)
table(totalmort$COMMENTS)

write.csv(totalmort, "totalmortality.csv")



##########################################################################

#                      Merge data for survival stuff


##########################################################################

setwd("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData")
cap <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/allcaptures.csv")
mort <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/totalmortality.csv")
mrre <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/MRRE.csv")


mrre2<-mrre[,-c(3,15,16,19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,53,57,58,59,60,61,62,64,65,66,67,68,69,70,71)]
trial<- merge(cap, mrre2, by.x="WOLFNO", by.y="WOLFID")
View(trial)
trial<-trial[-(1:1516),]
View(cap)
View(mort)
trial1<- merge(cap, mort, by="WOLFNO")
View(trial1)
trial1<-trial1[-c(1:806),]

length(which(cap$WOLFNO == ""))
unique(cap$WOLFNO)
unique(trial$WOLFNO)
unique(trial1$WOLFNO)


# merge(cap[,1:6], trial[,1:6], by="WOLFNO")
trial$LATITUDE.x<-as.factor(trial$LATITUDE.x)
trial$LONGITUDE.x<-as.factor(trial$LONGITUDE.x)
trial1$LATITUDE.y<-as.factor(trial1$LATITUDE.y)
trial1$LONGITUDE.y<-as.factor(trial1$LONGITUDE.y)

matchedcaps <- bind_rows(trial, trial1)

write.csv(matchedcaps, "matchedCaptures.csv")

matched <- read.csv("matchedCaptures2.csv")

matched <- matched[,order(names(matched))]
write.csv(matched, "matchedCaptures2.csv")



matched <- read.csv("matchedCaptures2.csv")
matched <- matched[,order(names(matched))]

dput(colnames(matched))

##arrange df vars by position
##'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}

matched <- arrange.vars(matched, c("WOLFNO"=1, "MRRENO"=2, "SEX"=3, 
                      "BREEDING_STATUS"=4, 
                      "CAPTURE_DATE"=5, "CAPTURE_YEAR"=6, 
                      "CAPTURE_AGE_CLASS"=7, "CAPTURE_AGE"=8, 
                      "DATE_OF_DEATH"=9, "DATE_LAST_ALIVE"=10, 
                      "MORT_YEAR"=11, "DEATH_AGE_CLASS"=12, 
                      "DEATH_AGE"=13, "PACK"=14, 
                      "PACK_AT_TIME_OF_DEATH"=15))

matched <- arrange.vars(matched, c("RECOVERY_AREA"=3, "Comments"=62, 
                                   "TRAP_INJURIES"=61, 
                                   "HARVEST_LICYR"=16, "HARVESTDATE"=17, 
                                   "HARVESTTIME"=18, 
                                   "HARVEST_LOCATION"=19, 
                                   "CAPTURE_DATUM"=20, 
                                   "CAPTURE_LATITUDE"=21, 
                                   "CAPTURE_LONGITUDE"=22, 
                                   "CAPTURE_ZONE"=23, 
                                   "CAPTURE_UTM_EASTING"=24, 
                                   "CAPTURE_UTM_NORTHING"=25, 
                                   "CAPTURE_LOCATION"=26, 
                                   "CAPTURE_METHOD"=27, "CAPTUREID"=28, 
                                   "COLLARID"=29, "FREQ"=30, 
                                   "RECAPTURE"=31))


write.csv(matched, "matchedCaptures2.csv")


matched <- read.csv("matchedCaptures2.csv")
captures <- read.csv("allcaptures.csv")

notmatched <- data.frame("WOLFNO" = c(setdiff(captures$WOLFNO, matched$WOLFNO)))
nomatchcap <- merge(notmatched, captures, by="WOLFNO")

morts <- read.csv("totalmortality.csv")

nomatchmort <- data.frame("WOLFNO" = c(setdiff(morts$WOLFNO, matched$WOLFNO)))
nomatchmort <- merge(nomatchmort, morts, by="WOLFNO")

mrre <- read.csv("MRRE.csv")

nomatchmrre <- data.frame("WOLFNO" = c(setdiff(mrre$WOLFID, matched$WOLFNO)))
nomatchmrre <- merge(nomatchmrre, morts, by="WOLFNO")

write.csv(nomatchcap, "nomatchcap.csv")
write.csv(nomatchmort, "nomatchmort.csv")
write.csv(nomatchmrre, "nomatchmrre.csv")







###########################################################

# Getting rid of duplicates in total mort data

###########################################################

totalmort <- read.csv(file.choose())
  

nondup_totalmort<-totalmort[duplicated(totalmort),]
write.csv(nondup_totalmort, "duplicatedtotalmort.csv")
