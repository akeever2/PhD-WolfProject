
library(R2jags)
library(mcmcplots)
library(survival)
library(dplyr)
library(tidyr)



sink("survival.txt")
cat("
    model{
    
    # Priors

      for(i in 1:ints){
        b.int[i] ~ dnorm(0,0.001)
      }

      for(i in 1:10){
        b[i] ~ dnorm(0, 0.001)
      }
    
    
    
    # Likelihood

      for(i in 1:nobs){
        event[i] ~ dbern(mu[i])
        cloglog(mu[i]) <- b.int[interval[i]] + b[1] * workprg[i] + 
                          b[2] * priors[i] + b[3] * tserved[i] + b[4] * 
                          felon[i] + b[5] * alcohol[i] + b[6] * 
                          drugs[i] + b[7] * black[i] + b[8] * 
                          married[i] + b[9] * educ[i] + b[10] * age[i]
      }
    


    # Predicted values
      for(i in 1:nobs){
        hazard[i] <- -log(1 - mu[i])
      }
    
    # Cumulative hazard and survival 
      
    
    }
    ", fill=TRUE)

sink()


win.data <- list("ints" = length(unique(recidx$interval)), 
                 "nobs" = nrow(recidx), "event" = recidx$fail, 
                 "interval" = recidx$interval, "workprg" = 
                   recidx$workprg, "priors" = recidx$priors, 
                 "tserved" = recidx$tserved, "felon" = recidx$felon, 
                 "alcohol" = recidx$alcohol, "drugs" = recidx$drugs, 
                 "black" = recidx$black, "married" = recidx$married, 
                 "educ" = recidx$educ, "age" = recidx$age)


#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("hazard", "b", "b.int") 


# MCMC Settings 
ni <- 10000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)











sink("survival.txt")
cat("
    model{
    
      # Priors
        #b0.surv ~ dnorm(0,0.001)

        for(i in 1:nyrs){
          bYR.surv[i] ~ dnorm(0,0.001)
        }

        for(i in 1:npds){
          bPD.surv[i] ~ dnorm(0,0.001)
        }

        

      # Likelihood
        for(i in 1:n){
          death[i] ~ dpois(mu.surv[i])
          log(mu.surv[i]) <- logEX[i] + bYR.surv[Year[i]] + 
                            bPD.surv[Period[i]]
        }
    
      # Predicted values
        for(i in 1:n){
          preds[i] <- logEX[i] + bYR.surv[Year[i]] + bPD.surv[Period[i]]
          hazard[i] <- exp(preds[i]-logEX[i])
        }

       # Cumulative hazard and survival 
        for(yr in 1:nyrs){
          for(i in (yr * npds - (npds - 1))){
            H[i] <- hazard[i] * w[Period[i]]
          }
            
          for(i in (yr * npds - (npds - 1) + 1):(npds * yr)){
            H[i] <- H[i-1] + hazard[i] * w[Period[i]]
          }
        }

        for(i in 1:n){
          S[i] <- exp(-H[i])
        }

    }
    ", fill=TRUE)

sink()


win.data <- list("logEX"=log(datum$exposure), "n"=length(datum$exposure), "Year"=datum$Year, 
                 "Period"=datum$Period, "death"=datum$deaths, "npds"=length(unique(datum$Period)),
                 "nyrs"=length(unique(datum$Year)), "w"=w)


#  Initial Values	
inits <- function(){list(bYR.surv=runif(3, -1, 1), 
                         bPD.surv=runif(5, -5, -1))}


# Parameters to keep track of and report
params <- c("H", "hazard", "bYR.surv", "bPD.surv", "S") 


# MCMC Settings 
ni <- 10000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)







data.fn <- function(bYR.surv.range = c(-1, 0.5), bPD.surv.range = c(-5,-2), 
                    bHarv.surv = -2.5, w = c(2, 3, 3, 3, 1), nyears = 3, 
                    nperiods = 5, Nsurv = 100){
  
  # Make empty vectors to hold values
  hazard <-0
  survival<-0
  deaths<-0
  exposure<-0
  
  # Make Year and Period covariate
  Year <- rep(1:nyears, nperiods)
  Year <- Year[order(Year)]
  Period <- rep(1:nperiods, nyears)
  
  # Make coefficients based on ranges
  bYR.surv <- runif(nyears, bYR.surv.range[1], bYR.surv.range[2])
  bPD.surv <- runif(nperiods, bPD.surv.range[1], bPD.surv.range[2])
  bPD.surv <- bPD.surv[order(bPD.surv)]
  
  # Baseline hazard
  for(i in 1:(nyears * nperiods)){
    hazard[i] <- exp(bYR.surv[Year[i]] + bPD.surv[Period[i]])
  }
  
  # Associated cumulative hazard and survival
  Year <- as.factor(Year)
  
  for(yr in levels(Year)) {
    i <- which(Year == yr)
    Hazard <- cumsum(hazard[i] * w)
    survival[i] <- exp(-Hazard)
  }
  
  
  # Expected survival time for individuals??
  
  for(yr in 1:nyears){
    for(i in (yr * nperiods - (nperiods - 1))){
      exposure[i] <- Nsurv * w[Period[i]] / 2
      deaths[i] <- rpois(1, hazard[i] * exposure[i])
    }
    
    for(i in (yr * nperiods - (nperiods - 1) + 1):(nperiods * yr)){
      exposure[i] <- (Nsurv - sum(deaths[(yr * nperiods - (nperiods -1)):(i-1)])) * 
        w[Period[i]] / 2
      deaths[i] <- rpois(1, hazard[i] * exposure[i])
    }
  }
  
  return(list(deaths=deaths, exposure=exposure, hazard=hazard, survival=survival,
              bYR.surv=bYR.surv, bPD.surv=bPD.surv, Year=Year, Period=Period)) 
  
}


datum<-data.fn()
















data.fn <- function(bYR.surv.range = c(-1, 0.5), bPD.surv.range = c(-5,-2), 
                    bHarv.surv = -2.5, w = c(2, 3, 3, 3, 1), nyears = 3, 
                    nperiods = 5, Nsurv = 100){
  
  # Make empty vectors to hold values
  hazard <-0
  survival<-0
  deaths<-0
  exposure<-0
  
  # Make Year and Period covariate
  Year <- rep(1:nyears, nperiods)
  Year <- Year[order(Year)]
  Period <- rep(1:nperiods, nyears)
  
  # Make coefficients based on ranges
  bYR.surv <- runif(nyears, bYR.surv.range[1], bYR.surv.range[2])
  bPD.surv <- runif(nperiods, bPD.surv.range[1], bPD.surv.range[2])
  bPD.surv <- bPD.surv[order(bPD.surv)]
  
  # Baseline hazard
  for(i in 1:(nyears * nperiods)){
    hazard[i] <- exp(bYR.surv[Year[i]] + bPD.surv[Period[i]])
  }
  
  # Associated cumulative hazard and survival
  Year <- as.factor(Year)
  
  for(yr in levels(Year)) {
    i <- which(Year == yr)
    Hazard <- cumsum(hazard[i] * w)
    survival[i] <- exp(-Hazard)
  }
  
  
  # Expected survival time for individuals??
  
  survtime.1 <- rexp(Nsurv, hazard[1:5])
  survtime.2 <- rexp(Nsurv, hazard[6:10])
  survtime.3 <- rexp(Nsurv, hazard[11:15])
  censort.1 = 12 * runif(Nsurv,0,1)
  censort.2 = 12 * runif(Nsurv,0,1)
  censort.3 = 12 * runif(Nsurv,0,1)
  
  event.1 = as.numeric(censort.1 > survtime.1)
  event.2 = as.numeric(censort.2 > survtime.2)
  event.3 = as.numeric(censort.3 > survtime.3)
  
  y.1 = survtime.1; 
  y.1[event.1==0] = censort.1[event.1==0]
  y.2 = survtime.2; 
  y.2[event.2==0] = censort.2[event.2==0]
  y.3 = survtime.3; 
  y.3[event.3==0] = censort.3[event.3==0]
  
  
  breaks=c(2,5,8,11,12)
  
  dat.1 <- data.frame("duration"=round(y.1 +.5), "event"=event.1, 
                      "id"=c(1:Nsurv))
  dat.2 <- data.frame("duration"=round(y.2 +.5), "event"=event.2, 
                      "id"=c(1:Nsurv))
  dat.3 <- data.frame("duration"=round(y.3 +.5), "event"=event.3, 
                      "id"=c(1:Nsurv))
  
  
  dat.1 <- survSplit(Surv(duration, event) ~ ., data=dat.1, cut=breaks,
                         episode="interval", start="start")
  dat.2 <- survSplit(Surv(duration, event) ~ ., data=dat.2, cut=breaks,
                         episode="interval", start="start")
  dat.3 <- survSplit(Surv(duration, event) ~ ., data=dat.3, cut=breaks,
                         episode="interval", start="start")

  dat.1 <- mutate(dat.1, 
                 interval=factor(as.numeric(interval), 
                                 labels=paste("(", c(0,breaks[1:4]), ",", 
                                              c(breaks[1:4],breaks[5]), 
                                              "]", sep="")), 
                 num.months=factor(as.numeric(interval), 
                                   labels=paste(c(2,3,3,3,1))), 
                 exposure=duration - start)
  dat.2 <- mutate(dat.2, 
                  interval=factor(as.numeric(interval), 
                                  labels=paste("(", c(0,breaks[1:4]), ",", 
                                               c(breaks[1:4],breaks[5]), 
                                               "]", sep="")), 
                  num.months=factor(as.numeric(interval), 
                                    labels=paste(c(2,3,3,3,1))), 
                  exposure=duration - start)
  dat.3 <- mutate(dat.3, 
                  interval=factor(as.numeric(interval), 
                                  labels=paste("(", c(0,breaks[1:4]), ",", 
                                               c(breaks[1:4],breaks[5]), 
                                               "]", sep="")), 
                  num.months=factor(as.numeric(interval), 
                                    labels=paste(c(2,3,3,3,1))), 
                  exposure=duration - start)
  
  
  
  dat.1 <- dat.1[order(dat.1$id, as.numeric(dat.1$interval)),]
  dat.2 <- dat.2[order(dat.2$id, as.numeric(dat.2$interval)),]
  dat.3 <- dat.3[order(dat.3$id, as.numeric(dat.3$interval)),]
  
  agr.1 <- aggregate(exposure ~ interval, data=dat.1, FUN=sum) %>% 
    mutate(year=1)
  agr.2 <- aggregate(exposure ~ interval, data=dat.2, FUN=sum) %>% 
    mutate(year=2) 
  agr.3 <- aggregate(exposure ~ interval, data=dat.3, FUN=sum) %>% 
    mutate(year=3) 
  
  total.dat <- bind_rows(agr.1,agr.2,agr.3)

  for(yr in 1:nyears){
    for(i in (yr * nperiods - (nperiods - 1))){
      deaths[i] <- rpois(1, hazard[i] * total.dat$exposure[i])
    }
    
    for(i in (yr * nperiods - (nperiods - 1) + 1):(nperiods * yr)){
      deaths[i] <- rpois(1, hazard[i] * total.dat$exposure[i])
    }
  }
  
  
  
  # # Determine number of deaths using GLM with offset
  # lm.1 <- exp(log(agr.1$exposure) + bYR.surv[1] + sum(bPD.surv))
  # lm.2 <- exp(log(agr.2$exposure) + bYR.surv[2] + sum(bPD.surv))
  # lm.3 <- exp(log(agr.3$exposure) + bYR.surv[3] + sum(bPD.surv))
  # deaths.1 <- rpois(nperiods, lm.1)
  
  return(list(deaths=deaths, exposure=total.dat$exposure, hazard=hazard, survival=survival,
              bYR.surv=bYR.surv, bPD.surv=bPD.surv, Year=Year, Period=Period)) 
  
}


datum<-data.fn()



















n.years <- 4
periods <- c("Apr-May", "Jun-Aug", "Sep-Nov", "Dec-Feb", "Mar")
w <- c(2, 3, 3, 3, 1)/12




fakewolf<-data.frame(wolfID=1:100)
fakewolf$start<-round(runif(100,0,47))
fakewolf$stop<-round(runif(100, fakewolf$start, 68))
fakewolf$fail<-ifelse(fakewolf$stop<49, rbinom(100,1,0.5), 0)
fakewolf$durat<-fakewolf$stop-fakewolf$start+1
fakewolf$durat[fakewolf$durat>48]=48
fakewolf$year1<-NA
fakewolf$year2<-NA
fakewolf$year3<-NA
fakewolf$year4<-NA

fakewolf$year1<-ifelse(fakewolf$start<12, 1, 0)
fakewolf$year2<-ifelse( (fakewolf$start<24 & fakewolf$stop>11), 1, 0)
fakewolf$year3<-ifelse( (fakewolf$start<36 & fakewolf$stop>23), 1, 0)
fakewolf$year4<-ifelse( (fakewolf$start<48 & fakewolf$stop>35), 1, 0)


fakewolf$int1<-ifelse(fakewolf$start<breaks[1], 1, 0)
fakewolf$int2<-ifelse((fakewolf$start<breaks[2] & fakewolf$stop>breaks[2]), 1, 0)
fakewolf$int3<-ifelse((fakewolf$start<breaks[3] & fakewolf$stop>breaks[3]), 1, 0)
fakewolf$int4<-ifelse((fakewolf$start<breaks[4] & fakewolf$stop>breaks[4]), 1, 0)
fakewolf$int5<-ifelse((fakewolf$start<breaks[5] & fakewolf$stop>breaks[5]), 1, 0)
fakewolf$int6<-ifelse((fakewolf$start<breaks[6] & fakewolf$stop>breaks[6]), 1, 0)
fakewolf$int7<-ifelse((fakewolf$start<breaks[7] & fakewolf$stop>breaks[7]), 1, 0)
fakewolf$int8<-ifelse((fakewolf$start<breaks[8] & fakewolf$stop>breaks[8]), 1, 0)
fakewolf$int9<-ifelse((fakewolf$start<breaks[9] & fakewolf$stop>breaks[9]), 1, 0)
fakewolf$int10<-ifelse((fakewolf$start<breaks[10] & fakewolf$stop>breaks[10]), 1, 0)
fakewolf$int11<-ifelse((fakewolf$start<breaks[11] & fakewolf$stop>breaks[11]), 1, 0)
fakewolf$int12<-ifelse((fakewolf$start<breaks[12] & fakewolf$stop>breaks[12]), 1, 0)
fakewolf$int13<-ifelse((fakewolf$start<breaks[13] & fakewolf$stop>breaks[13]), 1, 0)
fakewolf$int14<-ifelse((fakewolf$start<breaks[14] & fakewolf$stop>breaks[14]), 1, 0)
fakewolf$int15<-ifelse((fakewolf$start<breaks[15] & fakewolf$stop>breaks[15]), 1, 0)
fakewolf$int16<-ifelse((fakewolf$start<breaks[16] & fakewolf$stop>breaks[16]), 1, 0)
fakewolf$int17<-ifelse((fakewolf$start<breaks[17] & fakewolf$stop>breaks[17]), 1, 0)
fakewolf$int18<-ifelse((fakewolf$start<breaks[18] & fakewolf$stop>breaks[18]), 1, 0)
fakewolf$int19<-ifelse((fakewolf$start<breaks[19] & fakewolf$stop>breaks[19]), 1, 0)
fakewolf$int20<-ifelse((fakewolf$start<breaks[20] & fakewolf$stop>breaks[20]), 1, 0)

wolfdat2<-gather(fakewolf, "Interval", "Value", 10:29 )
wolfdat2 <- wolfdat2[order(wolfdat2$wolfID, wolfdat2$Interval),]
wolfdat2<-wolfdat2[!wolfdat2$Value==0,]

wolfdat2$Interval[wolfdat2$Interval=="int1"]<-1
wolfdat2$Interval[wolfdat2$Interval=="int2"]<-2
wolfdat2$Interval[wolfdat2$Interval=="int3"]<-3
wolfdat2$Interval[wolfdat2$Interval=="int4"]<-4
wolfdat2$Interval[wolfdat2$Interval=="int5"]<-5
wolfdat2$Interval[wolfdat2$Interval=="int6"]<-6
wolfdat2$Interval[wolfdat2$Interval=="int7"]<-7
wolfdat2$Interval[wolfdat2$Interval=="int8"]<-8
wolfdat2$Interval[wolfdat2$Interval=="int9"]<-9
wolfdat2$Interval[wolfdat2$Interval=="int10"]<-10
wolfdat2$Interval[wolfdat2$Interval=="int11"]<-11
wolfdat2$Interval[wolfdat2$Interval=="int12"]<-12
wolfdat2$Interval[wolfdat2$Interval=="int13"]<-13
wolfdat2$Interval[wolfdat2$Interval=="int14"]<-14
wolfdat2$Interval[wolfdat2$Interval=="int15"]<-15
wolfdat2$Interval[wolfdat2$Interval=="int16"]<-16
wolfdat2$Interval[wolfdat2$Interval=="int17"]<-17
wolfdat2$Interval[wolfdat2$Interval=="int18"]<-18
wolfdat2$Interval[wolfdat2$Interval=="int19"]<-19
wolfdat2$Interval[wolfdat2$Interval=="int20"]<-20
wolfdat2<-wolfdat2[,-11]


wolfdat2 <- mutate(wolfdat2, 
                   interval=factor(as.numeric(Interval), labels=paste("(", c(0,breaks[1:4]), 
                                                          ",", 
                                                          c(breaks[1:4],breaks[5]), 
                                                          "]", sep="")))
wolfdat2<-mutate(wolfdat2, num.months=factor(as.numeric(Interval), 
                                             labels=paste(c2,3,3,3,1)))

wolfdat2 <- wolfdat2[order(wolfdat2$wolfID, as.numeric(wolfdat2$Interval)),]




breaks<-c(c(2,5,8,11,12), 12+c(2,5,8,11,12), 24+c(2,5,8,11,12),36+c(2,5,8,11,12))
breaks2<-c(c(3,6,9,12,13), 12+c(3,6,9,12,13), 24+c(3,6,9,12,13),36+c(3,6,9,12,13))

wolfdat<-survSplit(Surv(durat, fail)~., data=fakewolf, cut=breaks, 
                   episode="interval", start="start.int")

wolfdat <- mutate(wolfdat, exposure = durat - start.int, 
                  interval = factor(interval,  
                                    labels = paste("(", c(0,breaks[1:19]), ",", 
                                                   c(breaks[1:19],breaks[20]), "]", sep=""))) %>%
  rename(events = fail)
wolfdat$Year<-ifelse(wolfdat$start.int<12, 1, 
                     ifelse(wolfdat$start.int<24, 2, 
                            ifelse(wolfdat$start.int<36, 3, 4)))
wolfdat3<-wolfdat[!((wolfdat$start.int+3)<wolfdat$start),]












wolfdat2<-data.frame(Period=rep(periods, 4), Year=c(rep(1,5), rep(2,5), rep(3,5),
                                                    rep(4,5)))








N.surv = 100
round=2
x1 = rbinom(N.surv,size=1,prob=0.5)
x = t(rbind(x1))
censortime = runif(N.surv,0,1)
survtime= rexp(N.surv,rate=exp(-1.5*x1))
survtime = round(survtime,digits=round)
event = as.numeric(censortime>survtime)
y = survtime; 
y[event==0] = censortime[event==0]
t=sort(unique(y[event==1]))
t=c(t,max(censortime))
bigt=length(t)-1


fit <- coxph(Surv(y, event)~x)
survs <- survfit(fit, newdata=data.frame(x=1))
survs2 <- survfit(fit, newdata=data.frame(x=0))
plot(survs)
lines(survs2)





##########################################################

# Cox proportional hazards model code with a counting process. 
# Adapted from Luek Winbugs example

##########################################################





library(R2jags)
library(mcmcplots)
################################################################################
#  Specify model in BUGS language

sink("survival.txt")
cat("
    
    data {
    
    # Set up the data
    
    for(i in 1:N.surv) {
    for(j in 1:T) {
    # Create risk set; Yij = 1 if obs.t >= t and start.t <= t
    Y[i,j] <- step(obs.t[i] - t[j] + eps) #* step(t[j] - start.t[i] + eps) 
    
    # Counting process jump; dN = 1 if obs.t is in interval [t[j], t[j+1) 
    # and fail = 1
    dN[i,j] <- Y[i,j] * step(t[j+1] - obs.t[i] - eps) * fail[i] 
    }
    }
    
    
    }    
    
    
    model {
    
    # Now that the data are set up, run the model
    
    
    # Specify priors
    
    # for(j in 1:T){
    # Prior for baseline hazard if using the Poisson trick
    # b0.surv[j] ~ dnorm(0, 0.001)
    # dL0[j] <- exp(b0.surv[j])
    # }
    
    # Prior for beta coefficient
    beta ~ dnorm(0, 0.0001)
    
    
    # Likelihood and Intensity process...What is the intensity process?
    
    for(j in 1:T) {
    for(i in 1:N.surv) {
    
    dN[i,j] ~ dpois(Idt[i,j])
    Idt[i,j] <- Y[i,j] * exp(beta * x[i]) * dL0[j]
    
    # Intensity process if using the Poisson trick - independent log-normal
    # increments which means drop dL0, c, r, and mu from model
    #Idt[i,j] <- Y[i,j] * exp(b0.surv[j] + beta * x[i])
    }
    
    # Prior mean hazard
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c
    
    # Survival function for both groups of fake data; multiply beta coefficient
    # by covariate value (1 or 0 in this example)
    S.treat[j] <- pow(exp(-sum(dL0[1:j])), exp(beta * 1))
    S.no[j] <- pow(exp(-sum(dL0[1:j])), exp(beta * 0))
    
    }
    
    # Confidence in guess for dL0
    c <- 0.001
    
    # Prior guess at failure rate
    r <- 0.1
    
    for(j in 1:T){
    dL0.star[j] <- r * (t[j+1] - t[j])
    }
    
    # Annual survival in this example because observations end after 1 year
    S.annual.treat <- S.treat[T]
    S.annual.no <- S.no[T]
    
    
    }
    ", fill=TRUE)
sink()


win.data <- list(x=x1, obs.t=y, t=t, T=bigt, N.surv=N.surv,
                 fail=event, eps=1E-10)



#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("S.annual.treat", "S.annual.no","S.no", "S.treat", "dL0", "beta") 


# MCMC Settings 
ni <- 5000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)



