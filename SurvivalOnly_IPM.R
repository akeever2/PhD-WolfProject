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

y.surv <- read.csv("ysurv_subset2.csv")
# y.surv$STAGE=ifelse(y.surv$stage == "ADULT", 3, ifelse(y.surv$stage == "YEARLING", 2, ifelse(y.surv$stage == "PUP", 1, 3)))
# 
# # y.surv.ad<-y.surv[y.surv$CAPTURE_AGE_CLASS != "ADULT" & y.surv$CAPTURE_AGE_CLASS != "YEARLING" | y.surv$freq > 6,]
# # y.surv.ad<-y.surv[y.surv$CAPTURE_AGE_CLASS == "ADULT" | y.surv$CAPTURE_AGE_CLASS == "YEARLING",]
# # y.surv.ad<-y.surv.ad[y.surv.ad$freq > 4, ]
# # y.surv.pup<-y.surv[y.surv$CAPTURE_AGE_CLASS == "PUP",]
# # y.surv2 <- bind_rows(y.surv.ad, y.surv.pup)
# 
# y.surv2 <- y.surv[y.surv$STAGE == 3 | y.surv$STAGE == 2,]
# y.surv2 <- y.surv2[y.surv2$freq > 4, ]
# 
# # Bring in covariate data
# packlist <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/PackList.csv")
# harv <- read.csv("C:/Users/allison/Documents/Project/WolfData/HarvestCovariates.csv")
# 
# # Merge the covariate data
# y.surv3 <- merge(y.surv2, packlist, sort=FALSE, by.x="PACK", by.y="Revised.Pack.Name")
# y.surv3 <- merge(y.surv3, harv, sort=FALSE, by.x="year", by.y="X")
# y.surv3$HarvMethod <- ifelse(y.surv3$Methods == "None", 1, ifelse(y.surv3$Methods == "Harvest", 2, 3))
# y.surv3$Region <- as.factor(y.surv3$FWP.Region)
# 
# write.csv(y.surv3, "ysurv_subset2.csv")

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
    loglik[i] <- logdensity.bern(event[i], mu.surv[i])
    #loglik2[i] <- logdensity.bin(event[i])
    }#i
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nyears"=length(unique(y.surv$year)), "event" = y.surv$Event, 
                 "nperiods" = length(unique(y.surv$period)),  "Period" = y.surv$period, 
                 "Year" = as.factor(y.surv$year), "width.interval" = c(2, 3, 3, 4), 
                 "nobs" =  nrow(y.surv), "Stage"=as.factor(y.surv$STAGE), 
                 "Region"=y.surv$Region, "Method"=as.factor(y.surv$HarvMethod), 
                 "nregions"=length(unique(y.surv$FWP.Region)), 
                 "Sex"=as.factor(Sex))


#  Initial Values	
inits <- function(){list(b0.surv=runif(1,-1,1))}


# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "hazard.adult", "base.H.adult", 
            "base.s.adult", "annual.s.adult", "hazard.yrl", "base.H.yrl", 
            "base.s.yrl", "annual.s.yrl", "hazard.pup", "base.H.pup", 
            "base.s.pup", "annual.s.pup", "b0.surv", "b.period.surv", 
            "b.stage") 


# MCMC Settings 
ni <- 100000
nt <- 3
nb <- 10000
nc <- 3


# Call JAGS 
out_m1 <- jags(win.data, inits, params, "Survival_m1.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m1, dig=2)

#mcmcplot(out3)

jag.sum <- out_m1$BUGSoutput$summary
write.table(x=jag.sum, file="m1_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m1 <- as.mcmc(out_m1)
mcmc_all_m1 <- rbind(mcmc_out_m1[[1]], mcmc_out_m1[[2]], mcmc_out_m1[[3]])
dim(mcmc_all_m1)

#  Computes WAIC value for model comparison with loo() package
m1_Survival_WAIC <- waic(mcmc_all_m1)

print(m1_Survival_WAIC)
save(m1_Survival_WAIC, file ="m1_WAIC.RData")



####   M2; S ~ Random Year + Random Region ####

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
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    for(r in 1:nregions){
    b.region[r] ~ dnorm(0,tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    cloglog(mu.pred.r1[p,k]) <- b.period.surv[p] + eps.surv[k] + b.region[1]
    hazard.r1[p,k] <- -log(1-mu.pred.r1[p,k])
    
    cloglog(mu.pred.r4[p,k]) <- b.period.surv[p] + eps.surv[k] + b.region[4]
    hazard.r4[p,k] <- -log(1-mu.pred.r4[p,k])
    
    }
    }
    
    
    
    # Cumulative hazard and survival by stage
    
    for(k in 1:nyears){
    base.H.r1[1,k] <- hazard.r1[1,k] * width.interval[1]
    for(p in 2:nperiods){
    base.H.r1[p,k] <- base.H.r1[p-1, k] + hazard.r1[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s.r1[p,k] <- exp(-base.H.r1[p,k])
    }#p
    
    annual.s.r1[k] <- base.s.r1[length(width.interval), k]
    }#k
    
    
    for(k in 1:nyears){
    base.H.r4[1,k] <- hazard.r4[1,k] * width.interval[1]
    for(p in 2:nperiods){
    base.H.r4[p,k] <- base.H.r4[p-1, k] + hazard.r4[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s.r4[p,k] <- exp(-base.H.r4[p,k])
    }#p
    
    annual.s.r4[k] <- base.s.r4[length(width.interval), k]
    }#k
    
    
    
    
    
    #   Calculate log likelihood for each iteration (for WAIC computations)
    #   Calculates log likelihood for each site, sampling occasion, and year
    #   based on the data, detection probability, and truth (occupied/unoccupied)
    
    for(i in 1:nobs){
    loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
    }#i
    
    
    }", fill=TRUE)
sink()



# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "annual.s.r1", "annual.s.r4", "b.period.surv", 
            "var.reg", "b.region") 


# Call JAGS 
out_m2 <- jags(win.data, inits, params, "Survival_m2.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m2, dig=2)

#mcmcplot(out3)

jag.sum <- out_m2$BUGSoutput$summary
write.table(x=jag.sum, file="m2_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m2 <- as.mcmc(out_m2)
mcmc_all_m2 <- rbind(mcmc_out_m2[[1]], mcmc_out_m2[[2]], mcmc_out_m2[[3]])
dim(mcmc_all_m2)

#  Computes WAIC value for model comparison with loo() package
m2_Survival_WAIC <- waic(mcmc_all_m2)

print(m2_Survival_WAIC)
save(m2_Survival_WAIC, file ="m2_WAIC.RData")


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
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
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
    
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred.no[p,k]) <- b.period.surv[p] + eps.surv[k] + b.method[1]
    hazard.no[p,k] <- -log(1-mu.pred.no[p,k])
    
    cloglog(mu.pred.hat[p,k]) <- b.period.surv[p] +eps.surv[k] + b.method[3]
    hazard.hat[p,k] <- -log(1-mu.pred.hat[p,k])
    
    }#p
    }#k
    
    
    # Cumulative hazard and survival by stage
    
    for(k in 1:nyears){
    base.H.no[1,k] <- hazard.no[1,k] * width.interval[1]
    for(p in 2:nperiods){
    base.H.no[p,k] <- base.H.no[p-1, k] + hazard.no[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s.no[p,k] <- exp(-base.H.no[p,k])
    }#p
    
    annual.s.no[k] <- base.s.no[length(width.interval), k]
    }#k
    
    
    for(k in 1:nyears){
    base.H.hat[1,k] <- hazard.hat[1,k] * width.interval[1]
    for(p in 2:nperiods){
    base.H.hat[p,k] <- base.H.hat[p-1, k] + hazard.hat[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s.hat[p,k] <- exp(-base.H.hat[p,k])
    }#p
    
    annual.s.hat[k] <- base.s.hat[length(width.interval), k]
    }#k
    
    
    
    #   Calculate log likelihood for each iteration (for WAIC computations)
    #   Calculates log likelihood for each site, sampling occasion, and year
    #   based on the data, detection probability, and truth (occupied/unoccupied)
    
    for(i in 1:nobs){
    loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
    }#i
    
    
    }", fill=TRUE)
sink()



# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "annual.s.no", "annual.s.hat", 
            "b.period.surv", "b.method") 


# Call JAGS 
out_m3 <- jags(win.data, inits, params, "Survival_m3.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m2, dig=2)

#mcmcplot(out3)

jag.sum <- out_m3$BUGSoutput$summary
write.table(x=jag.sum, file="m3_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m3 <- as.mcmc(out_m3)
mcmc_all_m3 <- rbind(mcmc_out_m3[[1]], mcmc_out_m3[[2]], mcmc_out_m3[[3]])
dim(mcmc_all_m3)

#  Computes WAIC value for model comparison with loo() package
m3_Survival_WAIC <- waic(mcmc_all_m3)

print(m3_Survival_WAIC)
save(m3_Survival_WAIC, file ="m3_WAIC.RData")


####   M4; S ~ Random Year + Stage ####

sink("Survival_m4.txt")
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
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    for(r in 1:2){
    b.stage[r] ~ dnorm(0,0.001)
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]] + b.stage[Stage[i]]
    }#i
    
    
    # Predicted values
    
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred.yrl[p,k]) <- b.period.surv[p] + eps.surv[k] + b.stage[1]
    hazard.yrl[p,k] <- -log(1-mu.pred.yrl[p,k])
    
    cloglog(mu.pred.adult[p,k]) <- b.period.surv[p] +eps.surv[k] + b.stage[2]
    hazard.adult[p,k] <- -log(1-mu.pred.adult[p,k])
    
    }#p
    }#k
    
    
    # Cumulative hazard and survival by stage
    
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
    
    
    
    #   Calculate log likelihood for each iteration (for WAIC computations)
    #   Calculates log likelihood for each site, sampling occasion, and year
    #   based on the data, detection probability, and truth (occupied/unoccupied)
    
    for(i in 1:nobs){
    loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
    }#i
    
    
    }", fill=TRUE)
sink()



# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "annual.s.adult", "annual.s.yrl", 
            "b.period.surv", "b.stage") 


# Call JAGS 
out_m4 <- jags(win.data, inits, params, "Survival_m4.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m2, dig=2)

#mcmcplot(out3)

jag.sum <- out_m4$BUGSoutput$summary
write.table(x=jag.sum, file="m4_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m4 <- as.mcmc(out_m4)
mcmc_all_m4 <- rbind(mcmc_out_m4[[1]], mcmc_out_m4[[2]], mcmc_out_m4[[3]])
dim(mcmc_all_m4)

#  Computes WAIC value for model comparison with loo() package
m4_Survival_WAIC <- waic(mcmc_all_m4)

print(m4_Survival_WAIC)
save(m4_Survival_WAIC, file ="m4_WAIC.RData")



####   M5; S ~ Random Year + Sex ####

sink("Survival_m5.txt")
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
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    for(r in 1:2){
    b.sex[r] ~ dnorm(0,0.001)
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]] + b.sex[Sex[i]]
    }#i
    
    
    # Predicted values
    
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred.f[p,k]) <- b.period.surv[p] + eps.surv[k] + b.sex[1]
    hazard.f[p,k] <- -log(1-mu.pred.f[p,k])
    
    cloglog(mu.pred.m[p,k]) <- b.period.surv[p] +eps.surv[k] + b.sex[2]
    hazard.m[p,k] <- -log(1-mu.pred.m[p,k])
    
    }#p
    }#k
    
    
    # Cumulative hazard and survival by stage
    
    for(k in 1:nyears){
    base.H.f[1,k] <- hazard.f[1,k] * width.interval[1]
    for(p in 2:nperiods){
    base.H.f[p,k] <- base.H.f[p-1, k] + hazard.f[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s.f[p,k] <- exp(-base.H.f[p,k])
    }#p
    
    annual.s.f[k] <- base.s.f[length(width.interval), k]
    }#k
    
    
    for(k in 1:nyears){
    base.H.m[1,k] <- hazard.m[1,k] * width.interval[1]
    for(p in 2:nperiods){
    base.H.m[p,k] <- base.H.m[p-1, k] + hazard.m[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s.m[p,k] <- exp(-base.H.m[p,k])
    }#p
    
    annual.s.m[k] <- base.s.m[length(width.interval), k]
    }#k
    
    
    
    #   Calculate log likelihood for each iteration (for WAIC computations)
    #   Calculates log likelihood for each site, sampling occasion, and year
    #   based on the data, detection probability, and truth (occupied/unoccupied)
    
    for(i in 1:nobs){
    loglik[i] <- logdensity.bin(event[i], mu.surv[i],1)
    }#i
    
    
    }", fill=TRUE)
sink()



# Parameters to keep track of and report
params <- c("eps.surv", "var.surv", "annual.s.f", "annual.s.m", 
            "b.period.surv", "b.sex") 


# Call JAGS 
out_m5 <- jags(win.data, inits, params, "Survival_m5.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

*#print(out_m2, dig=2)
  
  #mcmcplot(out3)
  
  jag.sum <- out_m5$BUGSoutput$summary
write.table(x=jag.sum, file="m5_SurvivalSummary.txt", sep="\t")

#  Pull MCMC chains together to compute WAIC
mcmc_out_m5 <- as.mcmc(out_m5)
mcmc_all_m5 <- rbind(mcmc_out_m5[[1]], mcmc_out_m5[[2]], mcmc_out_m5[[3]])
dim(mcmc_all_m5)

#  Computes WAIC value for model comparison with loo() package
m5_Survival_WAIC <- waic(mcmc_all_m5)

print(m5_Survival_WAIC)
save(m5_Survival_WAIC, file ="m5_WAIC.RData")


##### PLOTTING  #####

library(ggplot2)

ggplot(data=datum, aes(x=year, y=Surv, colour=Reg))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)


  scale_x_continuous(name="Year", breaks=c(1:10))+
  scale_y_continuous(name="Survival probability", breaks=c(0.2, 0.4, 0.6, 0.8, 1))+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14))

