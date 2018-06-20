#####################################################################/

#                         Survival submodel for IPM

#####################################################################/


#### Call Packages and Bring in Data ####
# Call for appropriate packages

library(R2jags)
library(mcmcplots)
library(loo)
library(dplyr)
library(tidyr)


# Set working directory
setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results")


# Bring in data

y.surv <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/ysurv_subset.csv")
y.surv$STAGE=ifelse(y.surv$stage == "ADULT", 3, ifelse(y.surv$stage == "YEARLING", 2, ifelse(y.surv$stage == "PUP", 1, 3)))

# y.surv.ad<-y.surv[y.surv$CAPTURE_AGE_CLASS != "ADULT" & y.surv$CAPTURE_AGE_CLASS != "YEARLING" | y.surv$freq > 6,]
# y.surv.ad<-y.surv[y.surv$CAPTURE_AGE_CLASS == "ADULT" | y.surv$CAPTURE_AGE_CLASS == "YEARLING",]
# y.surv.ad<-y.surv.ad[y.surv.ad$freq > 4, ]
# y.surv.pup<-y.surv[y.surv$CAPTURE_AGE_CLASS == "PUP",]
# y.surv2 <- bind_rows(y.surv.ad, y.surv.pup)

y.surv2 <- y.surv[y.surv$STAGE == 3 | y.surv$STAGE == 2,]
y.surv2 <- y.surv2[y.surv2$freq > 4, ]

# Bring in covariate data
packlist <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/PackList.csv")
harv <- read.csv("C:/Users/allison/Documents/Project/WolfData/HarvestCovariates.csv")

# Merge the covariate data
y.surv3 <- merge(y.surv2, packlist, sort=FALSE, by.x="PACK", by.y="Revised.Pack.Name")
y.surv3 <- merge(y.surv3, harv, sort=FALSE, by.x="year", by.y="X")
y.surv3$HarvMethod <- ifelse(y.surv3$Methods == "None", 1, ifelse(y.surv3$Methods == "Harvest", 2, 3))
y.surv3$Region <- as.factor(y.surv3$FWP.Region)

write.csv(y.surv3, "ysurv_subset2.csv")

################################################################################I
#  Specify model in BUGS language



####   M1; S ~ Random Year ####

sink("Survival_m1.txt")
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
    # b0.surv ~ dnorm(0,0.001)
    
    for(p in 1:nperiods){
      b.period.surv[p] ~ dnorm(0,0.001)
    }

    # for(i in 1:2){
    #   b.stage[i] ~ dnorm(0,0.001)
    # }

    

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
      # cloglog(mu.surv[i]) <- b0.surv + b.period.surv[Period[i]] + b.stage[Stage[i]] + eps.surv[Year[i]]
      cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard for pups (stage 1), yearlings (stage 2) and adults (stage 3)
    
     for(k in 1:nyears){
       for(p in 1:nperiods){
          cloglog(mu.pred.adult[p,k]) <- b.period.surv[p] + eps.surv[k]
          hazard.adult[p,k] <- -log(1-mu.pred.adult[p,k])

          # cloglog(mu.pred.yrl[p,k]) <- b0.surv + b.period.surv[p] + b.stage[2] + eps.surv[k]
          # hazard.yrl[p,k] <- -log(1-mu.pred.yrl[p,k])

          # cloglog(mu.pred.pup[p,k]) <- b0.surv + b.period.surv[p] + b.stage[1] + eps.surv[k]
          # hazard.pup[p,k] <- -log(1-mu.pred.pup[p,k])
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
    
    for(k in 1:nyears){
      base.H.adult[1,k] <- hazard.adult[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H.adult[p,k] <- base.H.adult[p-1, k] + hazard.adult[p,k] * width.interval[p]
      }#p
    }#k

    for(k in 1:nyears){
      for(p in 1:nperiods){
        base.s.adult[p,k] <- exp(-base.H.adult[p,k])
      }#p

      annual.s.adult[k] <- base.s.adult[length(width.interval), k]
    }#k


    # for(k in 1:nyears){
    #   base.H.yrl[1,k] <- hazard.yrl[1,k] * width.interval[1]
    #   for(p in 2:nperiods){
    #     base.H.yrl[p,k] <- base.H.yrl[p-1, k] + hazard.yrl[p,k] * width.interval[p]
    #   }#p
    # }#k
    # 
    # for(k in 1:nyears){
    #   for(p in 1:nperiods){
    #     base.s.yrl[p,k] <- exp(-base.H.yrl[p,k])
    #   }#p
    # 
    #   annual.s.yrl[k] <- base.s.yrl[length(width.interval), k]
    # }#k
    # 
    # 
    # for(k in 1:nyears){
    #   base.H.pup[1,k] <- hazard.pup[1,k] * width.interval[1]
    #   for(p in 2:nperiods){
    #     base.H.pup[p,k] <- base.H.pup[p-1, k] + hazard.pup[p,k] * width.interval[p]
    #   }#p
    # }#k
    # 
    # for(k in 1:nyears){
    #   for(p in 1:nperiods){
    #     base.s.pup[p,k] <- exp(-base.H.pup[p,k])
    #   }#p
    # 
    #   annual.s.pup[k] <- base.s.pup[length(width.interval), k]
    # }#k
    



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
win.data <- list("nyears"=length(unique(y.surv3$year)), "event" = y.surv3$Event, 
                 "nperiods" = length(unique(y.surv3$period)),  "Period" = y.surv3$period, 
                 "Year" = as.factor(y.surv3$year), "width.interval" = c(2, 3, 3, 4), 
                 "nobs" =  nrow(y.surv3), "Stage"=as.factor(y.surv3$STAGE), 
                 "Region"=y.surv3$Region, "Method"=as.factor(y.surv3$HarvMethod), 
                 "nregions"=length(unique(y.surv3$FWP.Region)))


#  Initial Values	
inits <- function(){list(b0.surv=runif(1,-1,1))}


# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "hazard.adult", "base.H.adult", 
            "base.s.adult", "annual.s.adult", "hazard.yrl", "base.H.yrl", 
            "base.s.yrl", "annual.s.yrl", "hazard.pup", "base.H.pup", 
            "base.s.pup", "annual.s.pup", "b0.surv", "b.period.surv", 
            "b.stage") 


# MCMC Settings 
ni <- 1000000
nt <- 3
nb <- 100000
nc <- 3


# Call JAGS 
out_m1 <- jags(win.data, inits, params, "Survival_m1.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out_m1, dig=2)

#mcmcplot(out3)

jag.sum <- out_m1$BUGSoutput$summary
write.table(x=jag.sum, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/m1_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m1 <- as.mcmc(out_m1)
mcmc_all_m1 <- rbind(mcmc_out_m1[[1]], mcmc_out_m1[[2]], mcmc_out_m1[[3]])
dim(mcmc_all_m1)

#  Computes WAIC value for model comparison with loo() package
m1_Survival_WAIC <- waic(mcmc_all_m1)

print(m1_Survival_WAIC)
#save(SS_Wolf_basep_nwolf_pHARV_WAIC, file ="C:/Users/sarah.bassing/POM & data/Outputs/SS_Wolf_basep_nwolf_pHARV_WAIC.RData")



####   M2; S ~ Random Year + Region ####

sink("Survival_m2.txt")
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
    # b0.surv ~ dnorm(0,0.001)
    
    for(p in 1:nperiods){
      b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    # for(i in 1:2){
    #   b.stage[i] ~ dnorm(0,0.001)
    # }

    for(r in 1:nregions){
      b.region[r] ~ dnorm(0,0.001)
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
      cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]] + b.region[Region[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard for pups (stage 1), yearlings (stage 2) and adults (stage 3)
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred.adult[p,k]) <- b.period.surv[p] + eps.surv[k] + b.region[1]
    hazard.adult[p,k] <- -log(1-mu.pred.adult[p,k])
    
    # cloglog(mu.pred.yrl[p,k]) <- b0.surv + b.period.surv[p] + b.stage[2] + eps.surv[k] + b.region[1]
    # hazard.yrl[p,k] <- -log(1-mu.pred.yrl[p,k])
    # 
    # cloglog(mu.pred.pup[p,k]) <- b0.surv + b.period.surv[p] + b.stage[1] + eps.surv[k] + b.region[1]
    # hazard.pup[p,k] <- -log(1-mu.pred.pup[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival by stage
    
    for(k in 1:nyears){
      base.H.adult[1,k] <- hazard.adult[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H.adult[p,k] <- base.H.adult[p-1, k] + hazard.adult[p,k] * width.interval[p]
      }#p
    }#k

    for(k in 1:nyears){
      for(p in 1:nperiods){
        base.s.adult[p,k] <- exp(-base.H.adult[p,k])
      }#p

      annual.s.adult[k] <- base.s.adult[length(width.interval), k]
    }#k
    
    
    # for(k in 1:nyears){
    # base.H.yrl[1,k] <- hazard.yrl[1,k] * width.interval[1]
    # for(p in 2:nperiods){
    # base.H.yrl[p,k] <- base.H.yrl[p-1, k] + hazard.yrl[p,k] * width.interval[p]
    # }#p
    # }#k
    # 
    # for(k in 1:nyears){
    # for(p in 1:nperiods){
    # base.s.yrl[p,k] <- exp(-base.H.yrl[p,k])
    # }#p
    # 
    # annual.s.yrl[k] <- base.s.yrl[length(width.interval), k]
    # }#k
    # 
    # 
    # for(k in 1:nyears){
    # base.H.pup[1,k] <- hazard.pup[1,k] * width.interval[1]
    # for(p in 2:nperiods){
    # base.H.pup[p,k] <- base.H.pup[p-1, k] + hazard.pup[p,k] * width.interval[p]
    # }#p
    # }#k
    # 
    # for(k in 1:nyears){
    # for(p in 1:nperiods){
    # base.s.pup[p,k] <- exp(-base.H.pup[p,k])
    # }#p
    # 
    # annual.s.pup[k] <- base.s.pup[length(width.interval), k]
    # }#k
    
    
    
    
    #   Calculate log likelihood for each iteration (for WAIC computations)
    #   Calculates log likelihood for each site, sampling occasion, and year
    #   based on the data, detection probability, and truth (occupied/unoccupied)
    
    for(i in 1:nobs){
    loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
    #loglik2[i] <- logdensity.bin(event[i])
    }#i
    
    
    }", fill=TRUE)
sink()



# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "hazard.adult", "base.H.adult", 
            "base.s.adult", "annual.s.adult", "hazard.yrl", "base.H.yrl", 
            "base.s.yrl", "annual.s.yrl", "hazard.pup", "base.H.pup", 
            "base.s.pup", "annual.s.pup", "b0.surv", "b.period.surv", 
            "b.stage", "b.region") 


# Call JAGS 
out_m2 <- jags(win.data, inits, params, "Survival_m2.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m2, dig=2)

#mcmcplot(out3)

jag.sum <- out_m2$BUGSoutput$summary
write.table(x=jag.sum, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/m2_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m2 <- as.mcmc(out_m2)
mcmc_all_m2 <- rbind(mcmc_out_m2[[1]], mcmc_out_m2[[2]], mcmc_out_m2[[3]])
dim(mcmc_all_m2)

#  Computes WAIC value for model comparison with loo() package
m2_Survival_WAIC <- waic(mcmc_all_m2)

print(m2_Survival_WAIC)
#save(SS_Wolf_basep_nwolf_pHARV_WAIC, file ="C:/Users/sarah.bassing/POM & data/Outputs/SS_Wolf_basep_nwolf_pHARV_WAIC.RData")



####   M3; S ~ Random Year + Method ####

sink("Survival_m3.txt")
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
    # b0.surv ~ dnorm(0,0.001)
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    # for(i in 1:2){
    # b.stage[i] ~ dnorm(0,0.001)
    # }
    
    for(r in 1:3){
    b.method[r] ~ dnorm(0,0.001)
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]] + b.method[Method[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard for pups (stage 1), yearlings (stage 2) and adults (stage 3)
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred.adult[p,k]) <- b.period.surv[p] + eps.surv[k] + b.method[3]
    hazard.adult[p,k] <- -log(1-mu.pred.adult[p,k])
    
    # cloglog(mu.pred.yrl[p,k]) <- b0.surv + b.period.surv[p] + b.stage[2] + eps.surv[k] + b.method[3]
    # hazard.yrl[p,k] <- -log(1-mu.pred.yrl[p,k])
    # 
    # cloglog(mu.pred.pup[p,k]) <- b0.surv + b.period.surv[p] + b.stage[1] + eps.surv[k] + b.method[3]
    # hazard.pup[p,k] <- -log(1-mu.pred.pup[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival by stage
    
    for(k in 1:nyears){
      base.H.adult[1,k] <- hazard.adult[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H.adult[p,k] <- base.H.adult[p-1, k] + hazard.adult[p,k] * width.interval[p]
      }#p
    }#k

    for(k in 1:nyears){
      for(p in 1:nperiods){
        base.s.adult[p,k] <- exp(-base.H.adult[p,k])
      }#p

      annual.s.adult[k] <- base.s.adult[length(width.interval), k]
    }#k
    
    # 
    # for(k in 1:nyears){
    # base.H.yrl[1,k] <- hazard.yrl[1,k] * width.interval[1]
    # for(p in 2:nperiods){
    # base.H.yrl[p,k] <- base.H.yrl[p-1, k] + hazard.yrl[p,k] * width.interval[p]
    # }#p
    # }#k
    # 
    # for(k in 1:nyears){
    # for(p in 1:nperiods){
    # base.s.yrl[p,k] <- exp(-base.H.yrl[p,k])
    # }#p
    # 
    # annual.s.yrl[k] <- base.s.yrl[length(width.interval), k]
    # }#k
    # 
    # 
    # for(k in 1:nyears){
    # base.H.pup[1,k] <- hazard.pup[1,k] * width.interval[1]
    # for(p in 2:nperiods){
    # base.H.pup[p,k] <- base.H.pup[p-1, k] + hazard.pup[p,k] * width.interval[p]
    # }#p
    # }#k
    # 
    # for(k in 1:nyears){
    # for(p in 1:nperiods){
    # base.s.pup[p,k] <- exp(-base.H.pup[p,k])
    # }#p
    # 
    # annual.s.pup[k] <- base.s.pup[length(width.interval), k]
    # }#k
    # 
    
    
    
    #   Calculate log likelihood for each iteration (for WAIC computations)
    #   Calculates log likelihood for each site, sampling occasion, and year
    #   based on the data, detection probability, and truth (occupied/unoccupied)
    
    for(i in 1:nobs){
    loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
    #loglik2[i] <- logdensity.bin(event[i])
    }#i
    
    
    }", fill=TRUE)
sink()



# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "hazard.adult", "base.H.adult", 
            "base.s.adult", "annual.s.adult", "hazard.yrl", "base.H.yrl", 
            "base.s.yrl", "annual.s.yrl", "hazard.pup", "base.H.pup", 
            "base.s.pup", "annual.s.pup", "b0.surv", "b.period.surv", 
            "b.stage", "b.method") 


# Call JAGS 
out_m3 <- jags(win.data, inits, params, "Survival_m3.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m2, dig=2)

#mcmcplot(out3)

jag.sum <- out_m3$BUGSoutput$summary
write.table(x=jag.sum, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/m3_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m3 <- as.mcmc(out_m3)
mcmc_all_m3 <- rbind(mcmc_out_m3[[1]], mcmc_out_m3[[2]], mcmc_out_m3[[3]])
dim(mcmc_all_m3)

#  Computes WAIC value for model comparison with loo() package
m3_Survival_WAIC <- waic(mcmc_all_m3)

print(m3_Survival_WAIC)
#save(SS_Wolf_basep_nwolf_pHARV_WAIC, file ="C:/Users/sarah.bassing/POM & data/Outputs/SS_Wolf_basep_nwolf_pHARV_WAIC.RData")





#### Extras ####











datum <- data.frame(s=c(out3$BUGSoutput$mean$annual.s.pup, 
                        out3$BUGSoutput$mean$annual.s.yrl, 
                        out3$BUGSoutput$mean$annual.s.adult), 
                    sd=c(out3$BUGSoutput$sd$annual.s.pup, 
                         out3$BUGSoutput$sd$annual.s.yrl, 
                         out3$BUGSoutput$sd$annual.s.adult), 
                    year=rep(2007:2017,3), 
                    stage=c(rep("pup", 11), rep("yearling", 11), rep("adult", 11)))


datum2 <- data.frame(s=c(out3$BUGSoutput$mean$annual.s.pup, 
                        out3$BUGSoutput$mean$annual.s.yrl), 
                    sd=c(out3$BUGSoutput$sd$annual.s.pup, 
                         out3$BUGSoutput$sd$annual.s.yrl), 
                    year=rep(2007:2017,2), 
                    stage=c(rep("pup", 11), rep("yearling", 11)))


ggplot(datum, aes(x=year, y=s))+
  geom_point()+
  geom_errorbar(aes(ymin=s-sd, ymax=s+sd), colour="black", width=.1)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(breaks=c(2007:2017))




