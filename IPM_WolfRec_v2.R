#####################################################################################I

#             Integrated Population Model to Estimate Recruitment of Wolves


# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2018

#             This model uses group count, collar, and occupancy data from MT
#             to estiamte recruitment from 2007-2014 using an IPM. The IPM model 
#             is outlined in word document ?????


#####################################################################################I

# Set the working directory I want to save files to

setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results")

#### Bring in data ####

# site covariates for occupancy
sitecovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/sitecovariates.csv")
# survey covariates for occupancy
survcovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/surveycovariates.csv")
# encounter histories for occupancy
encounter <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/encounterhistory.csv")

# group count data for the group level model
group <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/GroupCountData.csv")

# collar data for survival model
y.surv <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/ysurv_subset.csv")




#### Organize/format data for use in model ####

# Load packages for manipulating data
library(tidyr)
library(dplyr)

# y.occ is an array by 695 sites (rows), 5 survey occasions (columns), and 8 years (3D)
# unlist allows me to put it back together properly, I add 1 at the end because
# encounter data can't contain 0s for the model
y.occ <- array(unlist(encounter[,2:41]), dim=c(695, 5, 8))+1
nsites <- nrow(y.occ)
noccs <- ncol(y.occ)
nyears <- 8


# group.g is a new dataframe with only the good (G) and moderate (M) counts used
group.g <- group[group$Count.Quality=="G" | group$Count.Quality=="M", ]

# g2 is a new dataframe with only the year, pack, recreation area, and count of groups
# from the good and moderate counts
g2 <- data.frame(year=group.g$YEAR, pack=group.g$Revised.Pack.Name, area=group.g$RECAREA,  
                 count=group.g$CountPlusRemovals)

# the data are currently in long format so I change it to wide format
g2 <-g2 %>% spread(year, count)

# I only need years from 2007-2014, so I only use columns corresponding to those years 
y.group <- g2[,22:29]

# Some packs were only around prior to 2007, so I get rid of packs/group counts that 
# do not have any counts in years 2007-2014
g3 <- y.group[as.logical((rowSums(is.na(y.group))-8)),]

# My final group count data frame where the number of packs is the number of rows
y.group <- g3
ngroups <- nrow(y.group)


# The survival data currently includes 2007-2017, however I only need until 2014
# (7 years {nyears - 1} worth of data so I need to remove the data that I do not
# need
y.surv <- y.surv[y.surv$year < 2014]

# The event is whether the animal died (1) or not (0). If the animal was censored
# then the last full time period of observation is used and the event is 0. If
# the animal only lives through part of another period, then the previous period
# is used as the last observation
event <- y.surv$Event

# Period is which period each observation is in
Period = y.surv$period

# Year is which year each observation is in
Year = as.factor(y.surv$year)

# The width interval is the number of months in each period. Because they are 
# not all the same, the hazard is adjusted by multiplying by the number of 
# months in each period to account for differences when calculating the 
# cumulative harzard (H) and survival
width.interval = c(2, 3, 3, 4)
 
# The number of periods (which is currently 4) and the number of observations
nperiods = length(unique(y.surv$period))
nobs =  nrow(y.surv)
 

# Set my constants: em.group is a 10% dispersal rate from a pack, 
# territory size is set to 600 sq km, and territory overlap is set to the 
# values determined by FWP
em.group <- 0.1
mu.T <- 599.83
sd.T <- 368.21
T.overlap <- c(1.12, 1.08, 1.13, 1.16, 1.26, 1.27, 1.33, 1.24, 1.26, 1.32)





#### Model code ####


# Call for appropriate packages


library(R2jags)
library(mcmcplots)


################################################################################I
#  Specify model in BUGS language

sink("IPMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Occupancy priors
    
    # psi1 coefficient (occupancy in year 1)
        B0.psi1 ~ dnorm(0,0.001)	  
    
    # Priors for transition probabilities (survival and colonization)
        for(k in 1:(nyears-1)){
          B0.phi[k] ~ dnorm(0,0.001)
          B0.colo[k] ~ dnorm(0,0.001)
        }#k
    
    # Priors for detection probabilities (only varies by year)
        for(k in 1:nyears){
          B0.p11[k] ~ dnorm(0,0.001)
          B0.p10[k] ~ dnorm(0,0.001)
        }#k
    
        B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
        b.pc1.psi ~ dnorm(0,0.001) 
        b.pcks5yrs.psi ~ dnorm(0,0.001) 
        b.pc1.colo ~ dnorm(0,0.001) 
        b.pcks5yrs.colo ~ dnorm(0,0.001) 
        b.area ~ dnorm(0,0.001) 
        b.huntdays.p11 ~ dnorm(0,0.001) 
        b.acv ~ dnorm(0,0.001) 
        b.huntdays.p10 ~ dnorm(0,0.001) 
        b.nonfrrds.p10 ~ dnorm(0,0.001)
        b.nonfrrds.p11 ~ dnorm(0,0.001)
        b.frrds.p11 ~ dnorm(0,0.001) 
        b.frrds.p10 ~ dnorm(0,0.001) 
    

    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors

    # Random effect for year
      for(k in 1:(nyears-1)){
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



    ## 1.4 Group priors
    
    # Initial group sizes
        for(i in 1:ngroups){
          G[i,1] ~ dpois(5)
        }
    
    # Process error
        tauy.group <- pow(sigma.group, -2)
        sigma.group ~ dunif(0,100)
        var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
        B0.gam ~ dunif(-10,10)
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Occupancy likelihood 
    
    # Adapted from Sarah B Bassing
    # Montana Cooperative Wildlife Research Unit
    # August 2016. 
    # This is a DYNAMIC FALSE-POSITVE MULTI-SEASON occupancy model 
    # Encounter histories include:
    #     1 = no detection 
    #     2 = uncertain detection
    #     3 = certain detection
    
    ####################
    
    # Ecological process/submodel
    # Define State (z) conditional on parameters- Nmbr sites occupied
    
    for(i in 1:nsites){
      logit.psi1[i] <- B0.psi1 + b.pc1.psi * PC1[i] + b.pcks5yrs.psi * pcks5yrs[i,1]       
      logit(psi1[i]) <- logit.psi1[i]                                     
      z[i,1] ~ dbern(psi1[i])
      
      for(k in 1:(nyears-1)){
        logit.phi[i,k] <- B0.phi[k] 
        logit.colo[i,k] <- B0.colo[k] + b.pc1.colo * PC1[i] + b.pcks5yrs.colo * pcks5yrs[i,k] 
        logit(phi[i,k]) <- logit.phi[i,k]
        logit(colo[i,k]) <- logit.colo[i,k]
      }#k
      
      for(k in 2:nyears){
        muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*colo[i,k-1]
        z[i,k] ~ dbern(muZ[i,k])
      }#k
    }#i
    
    
    #  Observation process/submodel
    #  z is either 0 or 1 (unoccupied or occupied)
    #  y (observation dependent on the state z) can be 0,1,2 (no obs, uncertain obs,
    #  certain obs) but JAGS's multinomial link function, dcat(), needs y to be 1,2,3
    
    #  Observation process: define observations [y,i,j,k,z]
    #  y|z has a probability of...
    #  Detection probabilities are site, occasion, and year specific
    
    for(i in 1:nsites){
      for (j in 1:noccs){
        for(k in 1:nyears){			                                  
          p[1,i,j,k,1] <- (1-p10[i,j,k])
          p[1,i,j,k,2] <- (1-p11[i,j,k])
          p[2,i,j,k,1] <- p10[i,j,k]
          p[2,i,j,k,2] <- (1-b[i,j,k])*p11[i,j,k]
          p[3,i,j,k,1] <- 0
          p[3,i,j,k,2] <- b[i,j,k]*p11[i,j,k]
        }#k
      }#j
    }#i
    
    # Need mulitnomial link function for false positive detections: dcat function in JAGS
    # p11 is normal detction, p10 is false-positive dection (i.e., detected but wrong), 
    # b is certain detection
    

    # Observation model
    
    for(i in 1:nsites){
      for(j in 1:noccs){
        for(k in 1:nyears){
          multilogit.p11[i,j,k] <- B0.p11[k] + b.area * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i]
          logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
          multilogit.p10[i,j,k] <- B0.p10[k] + b.acv * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
          logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
          multilogit.b[i,j,k] <- B0.b
          logit(b[i,j,k]) <- multilogit.b[i,j,k]
          
          y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
        }#k
      }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
      psi[i,1] <- psi1[i]
      growthr[i,1] <- 1  
      
      for (k in 2:nyears){                                          
        psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*colo[i,k-1]
        growthr[i,k] <- psi[i,k]/psi[i,k-1]
        turnover[i,k-1] <- (1 - psi[i,k-1]) * colo[i,k-1]/psi[i,k]
      }#k
    }#i
    
    # Area occpupied indexed by year and region
    
    for(k in 1:nyears){
      A[k] <- sum(psi[,k] * area[])
    }
    
    
    
    
    
    #####################
    
    # 2.2. Territory model 
    
    # Input includes area occupied (A) indexed by year (k) from occupancy model (2.1.)
    # Area occupied is divided by territory size (T), which is currently set to 600
    # km squared, but will be based on data from Rich et al. 2012 for territory size

    # Output is number of packs (P) indexed by year (k)
    
    ####################
    
    # Pull in data for the mean for territory size
        T ~ dnorm(mu.T, 1/(sd.T * sd.T))
    
    # Estimate number of packs from area occupied (A) and territory size (T)
        for(k in 1:nyears){
          P[k] <- (A[k] / T)*T.overlap[k]
        }
    
    
    
    

    #####################
    
    # 2.3. Survival likelihood 
    
    # Current model is 
    # Output is survival indexed by year (k)
    
    ####################
    
    # Estimate the harzard. 
    # This part transforms the linear predictor (mu.surv)
    # using the cloglog link and relates it to the data (event) for each 
    # observation
    
    for(i in 1:nobs){
      event[i] ~ dbern(mu.surv[i])
      cloglog(mu.surv[i]) <- b0.surv + b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:(nyears-1)){
      for(p in 1:nperiods){
        cloglog(mu.pred[p,k]) <- b0.surv + b.period.surv[p] + eps.surv[k]
        hazard[p,k] <- -log(1-mu.pred[p,k])
      }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:(nyears-1)){
      base.H[1,k] <- hazard[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
      }#p
    }#k
    
    for(k in 1:(nyears-1)){
      for(p in 1:nperiods){
        base.s[p,k] <- exp(-base.H[p,k])
        }#p
      annual.s[k] <- base.s[length(width.interval), k]
    }#k
    

    
    
    
    #####################
    
    # 2.4. Group level counts likelihood 
    
    # Input data are group counts (y.group)
    # Input estimates are survival (s) from survival model indexed by year (k) and
    #   recruitment (number of pups per pack, gamma) indexed by year and group (i)
    
    # Output is mean estimate of group size (G) which are indexed by year and group
    
    ####################
    
    # Ecological model/ system process
        for(i in 1:ngroups){
          for(k in 2:nyears){
            g.mu[i,k] <- G[i,k-1] * s[k-1] * (1 - em.group) + gamma[i,k-1]
            G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+.01))T(0,)
          }
        }
    
    # Observation proccess
        for(i in 1:ngroups){
          for(k in 1:nyears){
            y.group[i,k] ~ dnorm(G[i,k], tauy.group)T(0,)
          }
        }
    
    # Derived parameters
        for(k in 1:nyears){
          G.mean[k] <- mean(G[,k])
          gamma.mean[k] <- mean(gamma[,k])
          n.est[k] <- P[k] * G.mean[k]
        }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
        for(i in 1:ngroups){
          for(k in 1:nyears){
            log(mu.gamma[i,k]) <- B0.gam
            gamma[i,k] ~ dpois(mu.gamma[i,k])
          }
        }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "ngroups"=ngroups, 
                 "area"=sitecovs$AREASAMP, "noccs"=noccs, "y.occ"=y.occ, 
                 "PC1"=sitecovs$PC1, "pcks5yrs"=sitecovs[,24:31], 
                 "huntdays"=array(unlist(survcovs$HUNTDAYS), dim=c(695, 5, 8)),
                 "nonforrds"=sitecovs$LOWUSENONFORESTRDS, 
                 "forrds"=sitecovs$LOWUSEFORESTRDS, 
                 "acv"=array(unlist(survcovs$ACV), dim=c(695, 5, 8)),  
                 "y.surv"=y.surv, "em.group"=em.group, "y.group"=y.group, 
                 "T.overlap"=T.overlap, "event"=event, "Period"=Period, 
                 "Year"=Year, "width.interval"=width.interval, 
                 "nperiods"=nperiods, "nobs"=nobs, "mu.T" =mu.T,
                 "sd.T"=sd.T)


#  Initial Values	
zst <- apply(y.occ,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 1
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(B0.gam=runif(1,-1,1), sd.proc=runif(1,0,10), 
                         sigma.group=runif(1,0,10), z=zst)}
#### add initial values ####

# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", "s", "phi", 
            "colo", "psi", "B0.phi", "B0.colo", "b.pc1.colo", "b.pcks5yrs.colo",
            "B0.gam", "A", "gamma", "G", "p11", "p10", "b", "eps.surv", 
            "hazard", "var.surv", "base.H", "base.s", "annual.s") 


# MCMC Settings 
ni <- 100000
nt <- 2
nb <- 10000
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)


#### format output for population level code ####

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 2))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean
phi2 <- array(NA, dim=c(nsites, nyears-1, 2))
phi2[,,1] <- out$BUGSoutput$mean$phi
phi2[,,2] <- out$BUGSoutput$sd$phi
colo2 <- array(NA, dim=c(nsites, nyears-1, 2))
colo2[,,1] <- out$BUGSoutput$mean$colo
colo2[,,2] <- out$BUGSoutput$sd$colo
b2 <- array(NA, dim=c(nsites, noccs, nyears, 2))
b2[,,,1] <- out$BUGSoutput$mean$b
b2[,,,2] <- out$BUGSoutput$sd$b
p102 <- array(NA, dim=c(nyears,1, 2))
p102[,,1] <- out$BUGSoutput$mean$p10
p102[,,2] <- out$BUGSoutput$sd$p10
p112 <- array(NA, dim=c(nyears,1, 2))
p112[,,1] <- out$BUGSoutput$mean$p11
p112[,,2] <- out$BUGSoutput$sd$p11
psi2 <- array(NA, dim=c(nsites, nyears, 2))
psi2[,,1] <- out$BUGSoutput$mean$psi
psi2[,,2] <- out$BUGSoutput$sd$psi
betas <- data.frame("B0.phi" = out$BUGSoutput$mean$B0.phi, 
                    "B0.phi.sd" = out$BUGSoutput$sd$B0.phi,
                    "B0.colo" = out$BUGSoutput$mean$B0.colo, 
                    "B0.colo.sd" = out$BUGSoutput$sd$B0.colo,
                    "b.pc1.colo" = out$BUGSoutput$mean$b.pc1.colo, 
                    "b.pc1.colo.sd" = out$BUGSoutput$sd$b.pc1.colo,
                    "b.pcks5yrs.colo" = out$BUGSoutput$mean$b.pcks5yrs.colo, 
                    "b.pcks5yrs.colo.sd" = out$BUGSoutput$sd$b.pcks5yrs.colo,
                    "B0.gam" = out$BUGSoutput$mean$B0.gam, 
                    "B0.gam.sd" = out$BUGSoutput$sd$B0.gam)

#### population level code ####

sink("PopLevel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## Population priors
    
    # Initial population size
    
    N.tot[1] ~ dnorm(800, 0.0001)I(0,)
    
    
    ## Bring in data s, G.mean, gamma.mean, P, colo, and phi
    
    for(k in 1:nyears){
    s[k] <- 0.72
    P[k] ~ dnorm(P2[k,1], 1 / (P2[k,2] * P2[k,2]))
    G.mean[k] ~ dnorm(G.mean2[k,1], 1 / (G.mean2[k,2] * G.mean2[k,2]))
    gamma.mean[k] ~ dnorm(gamma2[k,1], 1 / (gamma2[k,2] * gamma2[k,2]))
    }
    
    for(k in 1:(nyears-1)){
    B0.phi[k] ~ dnorm(betas[k,1], 1 / (betas[k,2] * betas[k,2]))
    B0.colo[k] ~ dnorm(betas[k,3], 1 / (betas[k,4] * betas[k,4]))
    }
    
    b.pc1.colo ~ dnorm(betas[1,5], 1 / (betas[1,6] * betas[1,6]))
    b.pcks5yrs.colo ~ dnorm(betas[1,7], 1 / (betas[1,8] * betas[1,8]))
    
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    
    # Ecological model/ system process
    
    # First determine colonization and extinction 
    
    for(i in 1:nsites){
    for(k in 1:(nyears-1)){
    colo[i,k] <- B0.colo[k] + b.pc1.colo * PC1[i] + b.pcks5yrs.colo * pcks5yrs[i,k]
    }
    }
    
    for(k in 1:(nyears-1)){
    phi[k] <- B0.phi[k]      
    }
    
    for(k in 2:nyears){
    N.rec[k] ~ dpois(P[k-1] * gamma.mean[k-1])
    N.ps[k] ~ dpois((P[k-1] + sum(colo[,k-1]*(area[]/T - P[k-1])) - P[k-1] * (1 - phi[k-1])) * G.mean[k-1])
    N.ad[k] ~ dbin(s[k-1], N.ps[k])
    N.tot[k] <- N.ad[k] + N.rec[k]
    #n.mu[k] <- N.tot[k-1]*(1+colo[k]-(1-phi[k]))*s[k-1]+P[k-1]*gamma.mean[k-1]
    #N.tot[k] ~ dpois(n.mu[k])
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    n.est[k,1] ~ dnorm(N.tot[k], (1 / (n.est[k,2]*n.est[k,2]+.001)))
    }
    
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data2 <- list("nyears"=nyears, "phi2"=phi2, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "s2"=s2, "gamma2"=gamma2, "P2"=P2, 
                  "colo2"=colo2, "T"=Tr, "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "pcks5yrs"=sitecovs[,24:31], "betas"=betas)


#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("P", "s", "N.tot", "gamma.mean") 


# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

print(out2, dig=2)
mcmcplot(out2)

