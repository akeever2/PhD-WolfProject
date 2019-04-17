#########################################################################################$

###         MAKE FAKE DATA

#########################################################################################$

#### DATA GENERATION ####
library(snowfall)
memory.limit(size = 7500000) #just shy of 8 tb limit i think

# Load package for categorical distribution to create false-positive encounter history
library(LaplacesDemon)

data.fn <- function(nsites=20, nyears=4, ngroups=12, s.prob.range=c(0.62, 0.67), Tr=600, 
                    em.range=c(0.08, 0.18), b0.gam.range=c(1.4, 1.7), sd.proc=4, 
                    sigma.group=1, nobs=20, psi1=0.4, p11.range=c(0.4,0.6), p10.range=c(0,0.1), 
                    b.range=c(0.2,0.4), patch.surv.range=c(0.7,0.9), patch.col.range=c(0,0.1),
                    noccs=3){
  
  # Function to generate data for IPM model to estimate recruitment using dynamic, false-positive occupancy
  # model and group counts with fixed survival. Annual varaition in parameters specified by uniform distr.
  
  # nsites = number of sites
  # nyears = number of years
  # ngroups = number packs monitored
  # nobs = number of observations for survival
  # s.prob.range = range for survival probability
  # Tr = territory size
  # em.range = bounds of uniform distribution to draw from for emigration out of group
  # b0.gam.range = range for intercept/ mean recruitment of pups
  # sd.proc = sd of process error for group counts
  # sigma.group = sd of observation error for group counts
  # psi1 = first year occupancy probability 
  # p11.range = bounds of uniform distribution for detection probability 
  # p10.range = bounds for false-positive detection probability 
  # b.range = bounds for certain detection probability 
  # patch.surv.range = bounds of uniform distribution for patch survival probability 
  # patch.col.range = bounds of uniform distribution for patch colinization probability
  # noccs = number of occasions for occupancy surveys
  
  
  # Arrays
  y.group <- array(dim=c(ngroups, nyears))
  G <- array(dim=c(ngroups, nyears))
  sigma.proc.group <- array(dim=c(nyears))
  y.temp <- array(dim=c(ngroups, nyears))
  G.mean <- array(dim=c(nyears))
  psi <- array(dim=c(nsites, nyears))
  area <- array(dim=c(nsites, nyears))
  A <- array(dim=c(nyears))
  P <- array(dim=c(nyears))
  N.tot <- array(dim=c(nyears))
  N.rec <- array(dim=c(nyears))
  y.occ <- array(dim=c(nsites, noccs, nyears))
  z <- muZ <- array(dim=c(nsites, nyears))
  p.0 <- p.1 <- array(dim=c(nyears, 3))
  prob <- array(dim=c(nsites, 3, nyears))
  y.surv <- array(dim=c(nobs,nyears))

  # Recruitment and group size
  # Recruitment
  b0.gam <- runif(n=nyears-1, min=b0.gam.range[1], max=b0.gam.range[2])
  
  gamma <- exp(b0.gam)
  
  # Random noise for first year group size
  for(i in 1:ngroups){
    G[i,1] <- max(2, rpois(1, 7))
  }
  
  # Process noise for group counts
  # for(k in 1:nyears){
  #   for(r in 1:nregions){
  #     sigma.proc.group[k,r] <- rnorm(1, 0, sd.proc)
  #   }
  # }
  
  # Determine annual group immigration and emigration based on bounds
  
  s.prob <- runif(n=nyears-1, min=s.prob.range[1], max=s.prob.range[2])
  em.group <- runif(n=nyears-1, min=em.range[1], max=em.range[2])
  
  # Group size in successive years
  for(i in 1:ngroups){
    for(k in 2:nyears){
      G[i,k] <- ifelse(G[i,k-1]==0, 0, rpois(1, (G[i,k-1]*s.prob[k-1]*(1-em.group[k-1])+gamma[k-1]))) #+sigma.proc.group[k-1,r])))
    }
  }
  
  # Group counts
  for(i in 1:ngroups){
    for(k in 1:nyears){
      y.temp[i,k] <- round(rnorm(1, G[i,k], sigma.group))
      y.temp[i,k] <- ifelse(y.temp[i,k] < 0, 0, y.temp[i,k])
      y.group[i,k] <- ifelse(y.temp[i,k] > G[i,k], G[i,k], y.temp[i,k])
    }
  }
  
  # Observation error
  var.group <- sigma.group*sigma.group
  
  # Mean group size
  for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
  }

  # survival data
  for(k in 1:(nyears-1)){
    for(i in 1:nobs){
      y.surv[i,k] <- rbinom(1,1,(1-s.prob[k]))
    }
  }
  
  # Occupancy and territory model to determine number of packs
  # Initial occupancy, psi
  for(i in 1:nsites){
    psi[i,1] <- psi1
  }
  
  # Determine detection and patch colinization and survivial
  p11 <- runif(n=nyears, min=p11.range[1], max=p11.range[2])
  p10 <- runif(n=nyears, min=p10.range[1], max=p10.range[2])
  b <- runif(n=nyears, min=b.range[1], max=b.range[2])
  patch.col <- runif(n=nyears-1, min=patch.col.range[1], max=patch.col.range[2])
  patch.surv <- runif(n=nyears-1, min=patch.surv.range[1], max=patch.surv.range[2])
  
  # Generate latent states of occurance
  # Year 1
  z[,1] <- rbinom(nsites, 1, psi[1])
  
  # Later years
  for(i in 1:nsites){
    for(k in 2:nyears){
      muZ[k] <- z[i,k-1]*patch.surv[k-1]+(1-z[i,k-1])*patch.col[k-1]
      z[i,k] <- rbinom(1,1,muZ[k])
    }
  }
  
  # Generate encounter histories using detection and true state of occurance
  # Set up detection vectors to be used in categorical distribution for when z=0 and z=1
  # p.0 is z=0, p.1 is z=1. The order is no detection (y=0), uncertain det (y=1) and 
  # certain det (y=2) for the probabilites
  for(k in 1:nyears){
    p.0[k,] <- c(1-p10[k], p10[k], 0)
    p.1[k,] <- c(1-p11[k], p11[k]*(1-b[k]), b[k]*p11[k])
  }
  
  # Create occupancy data
  for(i in 1:nsites){
    for(k in 1:nyears){
      if(z[i,k]==0){
        prob[i,,k] <- p.0[k,]
      } else {
        prob[i,,k] <- p.1[k,]
      }
      for(j in 1:noccs){
        y.occ[i,j,k] <- rcat(1, prob[i,,k])
      }
    }
  }
  
  # Compute annual occupancy so I can get area occupied
  for(i in 1:nsites){
    for(k in 2:nyears){
      psi[i,k] <- psi[i,k-1]*patch.surv[k-1]+(1-psi[i,k-1])*patch.col[k-1] 
    }
  }
  
  # Area
  for(i in 1:nsites){
    area[i,] <- 600
  }
  
  # Area occupied
  for(k in 1:nyears){
    A[k] <- sum(psi[,k]*area[,k])
  }
  
  # Number of packs
  for(k in 1:nyears){
    P[k] <- A[k]/Tr
  }
  
  # Population level
  # Total population in the first year
  N.tot <- P*G.mean


  
  
  return(list(nsites=nsites, nyears=nyears, ngroups=ngroups, nobs=nobs, noccs=noccs,
              P=P, Tr=Tr, patch.col=patch.col, patch.surv=patch.surv, A=A, 
              s.prob=s.prob, em.group=em.group, y.surv=y.surv, area=area,
              y.group=y.group, b0.gam=b0.gam, G.mean=G.mean, 
              gamma.mean=gamma, N.tot=N.tot, var.group=var.group,
              psi=psi, y.occ=y.occ, p10=p10, p11=p11, b=b, z=z))
}



datum<-data.fn(nsites=500, nyears=15, ngroups=100, s.prob.range=c(.55, .75),  
               em.range=c(0.08, 0.18), b0.gam.range=c(1.4, 1.7), sd.proc=4, 
               sigma.group=1, nobs=50, Tr=600, psi1=0.3, p11.range=c(0.4,0.6), 
               p10.range=c(0,0.1), b.range=c(0.1,0.2), 
               patch.surv.range=c(0.7,0.9), patch.col.range=c(0.1,0.6),
               noccs=5)

#detach()
attach(datum)
str(datum)

save(datum, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/datum_truth.RData")


# Make the datasets by sampling from the full data
library(dplyr)
library(tidyverse)

y.group.50 <- sample_n(as.data.frame(y.group), 50)
y.group.25 <- sample_n(as.data.frame(y.group), 25)
y.group.15 <- sample_n(as.data.frame(y.group), 15)

y.surv.20 <- sample_n(as.data.frame(y.surv), 20)
y.surv.10 <- sample_n(as.data.frame(y.surv), 10)

second <- seq(1, ncol(y.surv.10.2), 2)
y.surv.20.2 <- sample_n(as.data.frame(y.surv), 20)
y.surv.20.2[,second] <- NA
y.surv.10.2 <- sample_n(as.data.frame(y.surv), 10)
y.surv.10.2[,second] <- NA

fifth <- c(2,3,4,5,7,8,9,10,12,13,14,15)
y.surv.20.5 <- sample_n(as.data.frame(y.surv), 20)
y.surv.20.5[,fifth] <- NA
y.surv.10.5 <- sample_n(as.data.frame(y.surv), 10)
y.surv.10.5[,fifth] <- NA

# Write data as output to have for later
write.csv(y.group.50, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_group_50.csv")
write.csv(y.group.25, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_group_25.csv")
write.csv(y.group.15, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_group_15.csv")

write.csv(y.surv.20, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_surv_20.csv")
write.csv(y.surv.20.5, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_surv_20_5.csv")
write.csv(y.surv.20.2, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_surv_20_2.csv")
write.csv(y.surv.10, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_surv_10.csv")
write.csv(y.surv.10.5, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_surv_10_5.csv")
write.csv(y.surv.10.2, file="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/y_surv_10_2.csv")


#### IPM CODE ####
##############################################################################$

#             Integrated population model to estimate recruitment

# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2017

##############################################################################$

# Pull in data and call for appropriate packages


library(R2jags)
library(mcmcplots)

setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults")
################################################################################$
#  Specify model in BUGS language

sink("IPMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Occupancy priors

    #  psi1 coefficient (occupancy in year 1)
    
    psi1 ~ dunif(0,1)
    
    
    #  Priors for transition probabilities (survival and colonization)
    
    for(k in 1:(nyears-1)){
      phi[k] ~ dunif(0,1)
      colo[k] ~ dunif(0,1)
    }#k
    
    
    # Priors for detection probabilities (only varies by year)
    
    for(k in 1:nyears){
      p11[k] ~ dunif(0,1)
      p10[k] ~ dunif(0,1)
      b[k] ~ dunif(0,1)
    }#k

    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
      
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
      
      
   ## 1.4 Group priors
    
      # Initial group sizes
    
      for(i in 1:ngroups){
        G[i,1] ~ dpois(7)T(2,)
      }
    
    
    # Process and observation error

      tauy.group <- pow(sigma.group, -2)
      sigma.group ~ dunif(0,100)
      var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
      # Priors for beta coefficients
    
      for(k in 1:nyears){
        b0.gam[k] ~ dunif(-10,10)
      }
    
    
    
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
      z[i,1] ~ dbern(psi1)
    
      for(k in 2:nyears){
        muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*colo[k-1]
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
          p[1,i,j,k,1] <- (1-p10[k])
          p[1,i,j,k,2] <- (1-p11[k])
          p[2,i,j,k,1] <- p10[k]
          p[2,i,j,k,2] <- (1-b[k])*p11[k]
          p[3,i,j,k,1] <- 0
          p[3,i,j,k,2] <- b[k]*p11[k]
        }#k
      }#j
    }#i
    
    # Need mulitnomial link function for false positive detections: dcat function in JAGS
    # p11 is normal detction, p10 is false-positive dection (i.e., detected but wrong), 
    # b is certain detection
    
    #  Observation model
    
    for(i in 1:nsites){
      for(j in 1:noccs){
        for(k in 1:nyears){
          y.occ[i,j,k] ~ dcat(p[,i,j,k,z[i,k]+1])
        }#k
      }#j
    }#i
    
    
    # Derived parameters
    
    for(i in 1:nsites){
      psi[i,1] <- psi1
      growthr[i,1] <- 1  
    
      for (k in 2:nyears){                                          
        psi[i,k] <- psi[i,k-1]*phi[k-1] + (1-psi[i,k-1])*colo[k-1]
      }#k
    }#i
    
    #  Area occpupied indexed by year and region
    
    for(k in 1:nyears){
      A[k] <- sum(psi[,k] * area[,k])
    }
    #####################
    
    # 2.2. Territory model 
    
    ####################

    # Estimate number of packs from area occupied (A) and territory size (T)
    
    for(k in 1:nyears){
      P[k] <- A[k] / T
    }

    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:nyears){
      for(i in 1:nobs){
        event[i,k] ~ dbern(mu.surv[i,k])
        cloglog(mu.surv[i,k]) <- b0.surv + eps.surv[k]
      }
    }
    
    
    # Predicted values
    
    for(k in 1:(nyears-1)){
      cloglog(mu.pred[k]) <- b0.surv + eps.surv[k]
      hazard[k] <- -log(1-mu.pred[k])
    }
    
    # Cumulative hazard and survival 
    
    for(k in 1:(nyears-1)){
      H[k] <- hazard[k]
      annual.s[k] <- exp(-H[k])
    }
    
    
    #####################
    
    # 2.4. Group level counts likelihood 
    
    # Input data are group counts (y.group)
    # Input estimates are survival (s) from survival model indexed by year (k) and region (r) and
    #   recruitment (number of pups per pack, gamma) indexed by year, region, and group (i)
    
    # Output is mean estimate of group size (G) which are indexed by year, region, and group
    
    ####################
    
    # Ecological model/ system process
    
    for(i in 1:ngroups){
      for(k in 2:nyears){
        g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group[k-1]) + gamma[i,k-1] 
        G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))
        # G[i,k] ~ dpois(g.mu[i,k])T(2,)
      }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        y.group[i,k] ~ dnorm(G[i,k], tauy.group)T(2,)
      }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
      G.mean[k] <- mean(G[,k])
      n.est[k] <- P[k] * G.mean[k]
    }#k
    
    for(k in 1:(nyears-1)){
      gamma.mean[k] <- mean(gamma[,k])
    }
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(i in 1:ngroups){
      for(k in 1:(nyears-1)){
        mu.gamma[i,k] <- exp(b0.gam[k])
        gamma[i,k] ~ dpois(mu.gamma[i,k])
      }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.20, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


#  Initial Values	
zst <- apply(y.occ,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 0
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(b0.gam=runif(nyears, 1.4, 1.7), z=zst, colo=runif(nyears-1, 0.05, 0.15), phi=runif(nyears-1, 0.7,0.8),
                         b=runif(nyears, 0.1, 0.3), p10=runif(nyears, 0.05, 0.15),
                         p11=runif(nyears, 0.4, 0.6))} #, b0.surv=runif(nyears, 0.71, 0.74))}


# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", "annual.s",
            "phi", "colo", "psi", "p11", "p10", "b", "b0.surv", "eps.surv") 


# MCMC Settings 
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=3, n.thin=3, n.iter=50000, 
            n.burnin=20000, jags.module = c("glm", "dic"))

#print(out, dig=2)

#mcmcplot(out)

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
  }

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv, 14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv, 14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_full.RData")

#### POPULATION LEVEL CODE ####

sink("PopLevel_Sim.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## Population priors
    
    # Initial population size
    
    N.tot[1] ~ dnorm(1000, 0.0001)I(0,)
    
    
    ## Bring in data s, G.mean, gamma.mean, P, colo, and phi
    
    for(k in 1:nyears){
      P[k] ~ dnorm(P2[k,1], 1 / (P2[k,3] * P2[k,3]+ 0.0000001))
      # G.mean[k] ~ dnorm(G.mean2[k,1], 1 / (G.mean2[k,2] * G.mean2[k,2] + 0.0000001))
    }
    
    for(k in 1:(nyears-1)){
      # colo2[k] ~ dbeta(betas[k,1], betas[k,2])
      eps.surv[k] ~ dnorm(betas[k,7], 1 / (betas[k,8] * betas[k,8] + 0.0000001))
      # phi2[k] ~ dbeta(betas[k,3], betas[k,4])
      gamma.mean[k] ~ dnorm(gamma2[k,1], 1 / (gamma2[k,2] * gamma2[k,2] + 0.0000001))
    }

    b0.surv ~ dnorm(betas[1,5], 1 / (betas[1,6] * betas[1,6] + 0.0000001))
    
    
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    
    # Ecological model/ system process
    
    # First determine colonization and extinction 
    
    # for(i in 1:nsites){
    #   for(k in 1:(nyears-1)){
    #     colo[i,k] <- colo2[k]
    #     phi[i,k] <- phi2[k]
    #   }
    # }
    
    # Then determine survival
    # Baseline hazard
    
    for(k in 1:(nyears-1)){
      cloglog(mu.pred[k]) <- b0.surv + eps.surv[k]
      hazard[k] <- -log(1-mu.pred[k])
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:(nyears-1)){
      H[k] <- hazard[k]
      annual.s[k] <- exp(-H[k])
    }
    
    
    for(k in 2:nyears){
      N.rec[k] ~ dpois(P[k-1] * gamma.mean[k-1])
      # P.col[k] <- sum(colo[,k-1]*(area[,1]/T - P[k-1]))
      # P.ext[k] <- P[k-1] * (1 - sum(phi[,k-1]))
      # Ps[k] <- P[k-1] + P.col[k] - P.ext[k]
      # N.ps[k] ~ dpois(Ps[k] * G.mean[k-1])
      N.ad[k] ~ dbin(annual.s[k-1], round(N.tot[k-1]))
      N.tot[k] <- N.ad[k] + N.rec[k]
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
      n.est[k,1] ~ dnorm(N.tot[k], (1 / (n.est[k,2]*n.est[k,2]+0.00001)))
    }
    
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)


#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("P", "annual.s", "N.tot", "gamma.mean", 
             "N.rec", "eps.surv", 
             "N.ad") 


# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                n.burnin=nb, jags.module = c("glm", "dic"))


# truth <- data.frame("G.mean"=G.mean,
#                     "n.est"=N.tot, 
#                     "s"=c(s.prob, NA),
#                     "gamma"=c(gamma.mean,NA), 
#                     "dataset"=rep("Truth", nyears), "year"=c(1:15))


full <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                   "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                  "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                  "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                  "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                  "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                  "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                  "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(full, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_full.csv")
# write.csv(truth, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_truth.csv")

#### A: 50 group; 10 surv ####
y.surv.10<-read.csv(file.choose(), row.names=1)
y.surv.10[,15]<-as.integer(NA)

win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.10, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_A.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

A <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                   "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                   "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                   "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                   "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                   "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                   "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                   "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(A, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_A.csv")

#### B: 50 group; 20.2 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.20.2, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2,out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_B.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

B <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(B, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_B.csv")

#### C: 50 group; 10.2 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.10.2, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_C.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

C <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(C, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_C.csv")


#### D: 50 group; 20.5 surv ####
win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.20.5, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_D.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

D <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(D, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_D.csv")

#### E: 50 group; 10.5 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.10.5, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_E.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

E <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(E, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_E.csv")

#### F: 25 group; 20 surv ####
win.data <- list("nyears"=nyears, 
                 "ngroups"=25, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.25, "event"=y.surv.20, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_F.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

Ff <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(Ff, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_F.csv")
#### G: 25 group; 10 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=25, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.25, "event"=y.surv.10, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_G.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

G <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                 "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                 "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                 "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                 "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                 "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                 "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                 "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(G, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_G.csv")

#### H: 25 group; 20.2 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=25, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.25, "event"=y.surv.20.2, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_H.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

H <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(H, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_H.csv")

#### I: 25 group; 10.2 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=25, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.25, "event"=y.surv.10.2, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_I.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

I <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(I, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_I.csv")

#### J: 25 group; 20.5 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=25, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.25, "event"=y.surv.20.5, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_J.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

J <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(J, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_J.csv")
#### K: 25 group; 10.5 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=25, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.25, "event"=y.surv.10.5, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_K.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

K <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(K, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_K.csv")
#### L: 15 group; 20 surv ####
win.data <- list("nyears"=nyears, 
                 "ngroups"=15, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.15, "event"=y.surv.20, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_L.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

L <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(L, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_L.csv")

#### M: 15 group; 10 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=15, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.15, "event"=y.surv.10, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_M.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

M <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(M, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_M.csv")

#### N: 15 group; 20.2 surv ####


win.data <- list("nyears"=nyears, 
                 "ngroups"=15, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.15, "event"=y.surv.20.2, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_N.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

N <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(N, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_N.csv")

#### O: 15 group; 10.2 surv ####


win.data <- list("nyears"=nyears, 
                 "ngroups"=15, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.15, "event"=y.surv.10.2, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_O.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

O <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(O, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_O.csv")
#### P: 15 group; 20.5 surv ####


win.data <- list("nyears"=nyears, 
                 "ngroups"=15, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.15, "event"=y.surv.20.5, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_P.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

P <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(P, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_P.csv")
#### Q: 15 group; 10.5 surv ####

win.data <- list("nyears"=nyears, 
                 "ngroups"=15, "nobs"=10, 
                 "em.group"=em.group, 
                 "y.group"=y.group.15, "event"=y.surv.10.5, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))


n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out$BUGSoutput$mean$G.mean
G.mean2[,2] <- out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c(nyears-1, 2))
s2[,1] <- out$BUGSoutput$mean$annual.s
s2[,2] <- out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears-1, 2))
gamma2[,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out$BUGSoutput$sd$gamma.mean

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta,
                    "b0.surv" = rep(out$BUGSoutput$mean$b0.surv,14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv,14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

datums <-list(betas, n.est, G.mean2, s2, P2, gamma2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/out1_Q.RData")

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

Q <- data.frame("G.mean"=out2$BUGSoutput$mean$G.mean, "G.mean.lci"=out2$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out2$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[30:44,3],
                "n.est.uci"=out2$BUGSoutput$summary[30:44,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[60:73,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[60:73,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[7089:7102,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[7089:7102,7],NA),"dataset"=rep("full", nyears))

write.csv(Q, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_Q.csv")



#### Extras ####

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_I.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]

I <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(I, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_I.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_J.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
J <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(J, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_J.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_K.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
K <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(K, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_K.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_L.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
L <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(L, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_L.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)


load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_M.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
M <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(M, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_M.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_N.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
N <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                 "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                 "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                 "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                 "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                 "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                 "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                 "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(N, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_N.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_O.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
O <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(O, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_O.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)

load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_P.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
P <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(P, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_P.csv")

rm(betas, G.mean2, gamma2, n.est, P2, s2, datums, out2, out1)


load("~/Project/Dissertation/Recruitment/Results/SimulationResults/out1_Q.RData")
betas<-datums[[1]]
n.est<-datums[[2]]
G.mean2<-datums[[3]]
s2<-datums[[4]]
P2<-datums[[5]]
gamma2<-datums[[6]]

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))
out1 <- datums[[7]]
Q <- data.frame("G.mean"=out1$BUGSoutput$mean$G.mean, "G.mean.lci"=out1$BUGSoutput$summary[1:15,3],
                "G.mean.uci"=out1$BUGSoutput$summary[1:15,7],
                "n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[29:43,3],
                "n.est.uci"=out2$BUGSoutput$summary[29:43,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[59:72,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[59:72,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma.mean,NA), "gamma.lci"=c(out2$BUGSoutput$summary[88:101,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[88:101,7],NA),"dataset"=rep("full", nyears))

write.csv(Q, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_Q.csv")














A <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                   "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                   "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                   "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                   "dataset"=rep("A", nyears), "year"=c(1:15))

B <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                 "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                 "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                 "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                 "dataset"=rep("B", nyears), "year"=c(1:15))

C <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                 "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                 "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                 "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                 "dataset"=rep("C", nyears), "year"=c(1:15))

D <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                 "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                 "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                 "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                 "dataset"=rep("D", nyears), "year"=c(1:15))


output <- bind_rows(truth, full, A, B, C, D)
output$year <- rep(1:10, 6)
output[1:10, 7:10] <- 0
# output$G.mean.true <- rep(G.mean, 5)
# output$n.est.true <- rep(N.tot, 5)
# output$gamma.true <- rep(gamma.mean, 5)
# output$s.true <- rep(s.prob, 5)





ggplot(data=output2, aes(x=year, y=gamma, colour=dataset))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=(gamma-gamma.sd), ymax=(gamma+gamma.sd)), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)+
  scale_x_continuous(name="Year", breaks=c(1:10))+
  scale_y_continuous(name="Recruitment", breaks=c(1,2,3,4,5,6,7,8))+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14))


ggplot(data=output2, aes(x=year, y=s, colour=dataset))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=(s-s.sd), ymax=(s+s.sd)), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)+
  scale_x_continuous(name="Year", breaks=c(1:10))+
  scale_y_continuous(name="Survival probability", breaks=c(0.2, 0.4, 0.6, 0.8, 1))+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14))


ggplot(data=output2, aes(x=year, y=n.est, colour=dataset))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=(n.est-n.est.sd), ymax=(n.est+n.est.sd)), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)+
  scale_y_continuous(name="Abundance")+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16))

write.csv(output, "simoutput.csv")
output2 <- read.csv("simoutput.csv")



