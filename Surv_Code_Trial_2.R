

library(R2jags)
library(mcmcplots)


sink("survival.txt")
cat("
    model{
    
      # Priors
        b0.surv ~ dnorm(0,0.001)

        for(i in 1:nyrs){
          bYR.surv[i] ~ dnorm(0,0.001)
        }

        for(i in 1:npds){
          bPD.surv[i] ~ dnorm(0,0.001)
        }

        

      # Likelihood
        #for(t in 1:nyears){
          for(i in 1:n){
            S[i] ~ dpois(mu.surv[i])
            log(mu.surv[i]) <- logEX[i] + b0.surv + bYR.surv[Year[i]] + 
                              bPD.surv[Period[i]]
          }
        #}
    
      # Predicted values
        for(i in 1:n){
          preds[i] <- logEX[i] +b0.surv + bYR.surv[Year[i]] + bPD.surv[Period[i]]
          hazard[i] <- exp(preds[i]-logEX[i])
        }

      # # Cumulative hazard and survival 
      #   for(yr in Year){
      #     H[Year[yr[1]]] <- hazard[Year[yr[1]]]
      # 
      #     for(i in 2:npds){
      #       H[Year[yr[i]]] <- H[Year[yr[i-1]]] + hazard[Year[yr[i]]] * w[Period[i]]
      #     }
      #   }
    


    
    }
    ", fill=TRUE)

sink()


win.data <- list("logEX"=log(co$exposure), "n"=nrow(co), "Year"=co$cohort, 
                 "Period"=co$age, "S"=co$deaths, "npds"=length(unique(co$age)),
                 "nyrs"=length(unique(co$cohort)), "w"=w)


#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("H", "hazard", "b0.surv", "bYR.surv", "bPD.surv") 


# MCMC Settings 
ni <- 1000
nt <- 2
nb <- 250
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)






for(cohort in levels(co$cohort)) {
  i <- which(co$cohort == cohort)
  Hazard <- cumsum(co$hazard[i] * w)
  co$survival[i] <- exp(-Hazard)
}






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
                   interval=factor(as.numeric(Interval), labels=paste("(", c(0,breaks[1:19]), 
                                                          ",", 
                                                          c(breaks[1:19],breaks[20]), 
                                                          "]", sep="")))
wolfdat2<-mutate(wolfdat2, num.months=factor(as.numeric(Interval), 
                                             labels=paste(rep(c(2,3,3,3,1),4))))

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
survtime= rexp(N.surv,rate=exp(-.5*x1))
survtime = round(survtime,digits=round)
event = as.numeric(censortime>survtime)
y = survtime; 
y[event==0] = censortime[event==0]
t=sort(unique(y[event==1]))
t=c(t,max(censortime))
bigt=length(t)-1






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



