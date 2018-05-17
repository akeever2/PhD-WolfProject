##############################################################################################################

#                         Survival submodel for IPM

##############################################################################################################



# Bring in data

y.surv <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/ysurv_subset.csv")
y.surv$STAGE=ifelse(y.surv$stage == "ADULT", 2, ifelse(y.surv$stage == "YEARLING", 2, ifelse(y.surv$stage == "PUP", 1, 2)))

y.surv<-y.surv[y.surv$CAPTURE_AGE_CLASS != "ADULT" & y.surv$CAPTURE_AGE_CLASS != "YEARLING" | y.surv$freq > 6,]


# Call for appropriate packages


library(R2jags)
library(mcmcplots)


################################################################################I
#  Specify model in BUGS language

sink("Survival.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################


    # Random effect for year
    for(k in 1:nyears){
      eps.surv[k] ~ dnorm (0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    # Beta coefficients
    b0.surv ~ dnorm(0,0.001)
    
    for(p in 1:nperiods){
      b.period.surv[p] ~ dnorm(0,0.001)
    }

    for(i in 1:2){
      b.stage[i] ~ dnorm(0,0.001)
    }

    

    ############################################################
    
    #             2. Likelihoods
    
    ############################################################


    #####################
    
    # 2.3. Survival likelihood 
    
    # Output is survival indexed by year (k)
    
    ####################
    
    # Estimate the harzard. 
    # This part transforms the linear predictor (mu.surv)
    # using the cloglog link and relates it to the data (event) for each 
    # observation
    
    for(i in 1:nobs){
      event[i] ~ dbern(mu.surv[i])
      cloglog(mu.surv[i]) <- b0.surv + b.period.surv[Period[i]] + b.stage[Stage[i]] + eps.surv[Year[i]]
      # cloglog(mu.surv[i]) <- b0.surv + b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard for pups (stage 1), yearlings (stage 2) and adults (stage 3)
    
     for(k in 1:nyears){
       for(p in 1:nperiods){
          # cloglog(mu.pred.adult[p,k]) <- b0.surv + b.period.surv[p] + b.stage[3] + eps.surv[k]
          # hazard.adult[p,k] <- -log(1-mu.pred.adult[p,k])

          cloglog(mu.pred.yrl[p,k]) <- b0.surv + b.period.surv[p] + b.stage[2] + eps.surv[k]
          hazard.yrl[p,k] <- -log(1-mu.pred.yrl[p,k])

          cloglog(mu.pred.pup[p,k]) <- b0.surv + b.period.surv[p] + b.stage[1] + eps.surv[k]
          hazard.pup[p,k] <- -log(1-mu.pred.pup[p,k])
       }#p
     }#k
    
    
    # # Cumulative hazard and survival with yearlings and adults grouped
    #   for(k in 1:nyears){
    #     for(p in 1:nperiods){
    #       cloglog(mu.pred[p,k]) <- b0.surv + b.period.surv[p] + eps.surv[k]
    #       hazard[p,k] <- -log(1-mu.pred[p,k])
    #     }#p
    #   }#k
    # 
    # 
    # for(k in 1:nyears){
    #   base.H[1,k] <- hazard[1,k] * width.interval[1]
    #   for(p in 2:nperiods){
    #     base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    #   }#p
    # }#k
    # 
    # for(k in 1:nyears){
    #   for(p in 1:nperiods){
    #     base.s[p,k] <- exp(-base.H[p,k])
    #   }#p
    # 
    #   annual.s[k] <- base.s[length(width.interval), k]
    # }#k


    # Cumulative hazard and survival by stage
    
    # for(k in 1:nyears){
    #   base.H.adult[1,k] <- hazard.adult[1,k] * width.interval[1]
    #   for(p in 2:nperiods){
    #     base.H.adult[p,k] <- base.H.adult[p-1, k] + hazard.adult[p,k] * width.interval[p]
    #   }#p
    # }#k
    # 
    # for(k in 1:nyears){
    #   for(p in 1:nperiods){
    #     base.s.adult[p,k] <- exp(-base.H.adult[p,k])
    #   }#p
    # 
    #   annual.s.adult[k] <- base.s.adult[length(width.interval), k]
    # }#k


    for(k in 1:nyears){
      base.H.yrl[1,k] <- hazard.yrl[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H.yrl[p,k] <- base.H.yrl[p-1, k] + hazard.yrl[p,k] * width.interval[p]
      }#p
    }#k

    for(k in 1:nyears){
      for(p in 1:nperiods){
        base.s.yrl[p,k] <- exp(-base.H.yrl[p,k])
      }#p

      annual.s.yrl[k] <- base.s.yrl[length(width.interval), k]
    }#k


    for(k in 1:nyears){
      base.H.pup[1,k] <- hazard.pup[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H.pup[p,k] <- base.H.pup[p-1, k] + hazard.pup[p,k] * width.interval[p]
      }#p
    }#k

    for(k in 1:nyears){
      for(p in 1:nperiods){
        base.s.pup[p,k] <- exp(-base.H.pup[p,k])
      }#p

      annual.s.pup[k] <- base.s.pup[length(width.interval), k]
    }#k
    



  #   Calculate log likelihood for each iteration (for WAIC computations)
  #   Calculates log likelihood for each site, sampling occasion, and year
  #   based on the data, detection probability, and truth (occupied/unoccupied)

    for(i in 1:nobs){
      loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
      #loglik2[i] <- logdensity.bin(event[i])
    }#i

    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nyears"=length(unique(y.surv$year)), "event" = y.surv$Event, 
                 "nperiods" = length(unique(y.surv$period)),  "Period" = y.surv$period, 
                 "Year" = as.factor(y.surv$year), "width.interval" = c(2, 3, 3, 4), 
                 "nobs" =  nrow(y.surv), "Stage"=as.factor(y.surv$STAGE))


#  Initial Values	
inits <- function(){list(b0.surv=runif(1,-1,1))}


# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "hazard.adult", "base.H.adult", 
            "base.s.adult", "annual.s.adult", "hazard.yrl", "base.H.yrl", 
            "base.s.yrl", "annual.s.yrl", "hazard.pup", "base.H.pup", 
            "base.s.pup", "annual.s.pup") 


# MCMC Settings 
ni <- 50000
nt <- 2
nb <- 5000
nc <- 3


# Call JAGS 
out3 <- jags(win.data, inits, params, "Survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out3, dig=2)

mcmcplot(out3)

datum <- data.frame(s=c(out3$BUGSoutput$mean$annual.s.pup, 
                        out3$BUGSoutput$mean$annual.s.yrl), 
                        #out3$BUGSoutput$mean$annual.s.adult), 
                    sd=c(out3$BUGSoutput$sd$annual.s.pup, 
                         out3$BUGSoutput$sd$annual.s.yrl), 
                         #out3$BUGSoutput$sd$annual.s.adult), 
                    year=rep(2007:2017,2), 
                    stage=c(rep("pup", 11), rep("yearling", 11)))#, rep("adult", 11)))

ggplot(datum, aes(x=year, y=s))+
  geom_point()+
  geom_errorbar(aes(ymin=s-sd, ymax=s+sd), colour="black", width=.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(breaks=c(2007:2017))
