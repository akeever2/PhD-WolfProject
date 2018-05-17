
####Prepping data for survival####



cap<-read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/matchedCaptures2.csv")
str(cap)



y.surv <- cap[,c(1,4,5,7,10,12,19,20,21)]

y.surv$start <- as.Date(cap$CAPTURE_DATE, "%m/%d/%Y")
y.surv$stop <- as.Date(y.surv$stop, "%m/%d/%Y")

str(y.surv)
table(y.surv$PACK)


# ytemp <- data.frame(intstartmo= as.integer(substr(y.surv$start, 6,7)))
# ytemp$intstopmo <- as.integer(substr(y.surv$stop, 6,7))
# ytemp$startYr <- as.integer(substr(y.surv$start, 1,4))
# ytemp$stopYr <- as.integer(substr(y.surv$stop, 1,4))
# ytemp$srP <- ifelse(ytemp$intstartmo > 3 & ytemp$intstartmo < 6, 1,
#                     ifelse(ytemp$intstartmo > 5 & ytemp$intstartmo < 9, 2,
#                            ifelse(ytemp$intstartmo > 8 & ytemp$intstartmo < 12,
#                                   3, 4)))
# ytemp$spP <- ifelse(ytemp$intstopmo > 3 & ytemp$intstopmo < 6, 1,
#                     ifelse(ytemp$intstopmo > 5 & ytemp$intstopmo < 9, 2,
#                            ifelse(ytemp$intstopmo > 8 & ytemp$intstopmo < 12,
#                                   3, 4)))
# ytemp$freq <- ytemp$spP - ytemp$srP

# y.surv$freq <- ytemp$freq
# write.csv(y.surv, "ysurv.csv")
y.surv2 <- read.csv("ysurv.csv")
y.surv$freq <- y.surv2$freq
y.surv$startint <- ytemp$srP
y.surv <- y.surv %>% drop_na(freq)

y.surv <- expandRows(y.surv, "freq", drop = FALSE)

# y.surv <- y.surv %>% group_by(WOLFNO) %>%
#   mutate(period = for(i in 1:n()){
#     ifelse(abs(first(freq) - 4 + first(startint) - 1) < first(freq),
#            c(first(startint):4, 
#              rep_len(1:4, (first(freq) - 4 + first(startint) - 1))),
#            first(startint):(first(startint) + first(freq) - 1))
#     })
# 
# y.surv$period <- NA
# for(i in levels(y.surv$WOLFNO)){
#   for(j in 1:n()){
#   period[i[j]] <- ifelse(abs(y.surv$freq - 4 + y.surv$startint - 1) < y.surv$freq,
#                    c(y.surv$startint:4, 
#                      rep_len(1:4, (y.surv$freq - 4 + y.surv$startint - 1))),
#                    y.surv$startint:(y.surv$startint + y.surv$freq - 1))
#   }
# }



