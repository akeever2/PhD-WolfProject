#####################################################################I

#             Simulation of recruitment IPM

#####################################################################I

#### Function to create datasets ####

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


#### Load packages and set wd ####

library(dplyr)
library(tidyverse)
library(R2jags)
library(snowfall)
memory.limit(size = 7500000) #just shy of 8 tb limit i think


setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults")

#### Model code ####
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
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)T(0,)
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
    eps.surv[k] ~ dnorm(betas[k,3], 1 / (betas[k,4] * betas[k,4] + 0.0000001))
    # phi2[k] ~ dbeta(betas[k,3], betas[k,4])
    gamma.mean[k] ~ dnorm(gamma2[k,1], 1 / (gamma2[k,2] * gamma2[k,2] + 0.0000001))
    }
    
    b0.surv ~ dnorm(betas[1,1], 1 / (betas[1,2] * betas[1,2] + 0.0000001))
    
    
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



#### Constants ####

# Parameters to keep track of and report
params <- c("P", "G.mean", "gamma.mean", "n.est", "annual.s",
            "b0.surv", "eps.surv") 

#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("annual.s", "N.tot", "gamma.mean") 


#### Function to run models ####

modelrun_fn <- function(ngroups2=12, nobs2=20, nobsyr=1){

# Create and attach data
datum <- data.fn(nsites=500, nyears=15, ngroups=ngroups2, s.prob.range=c(.55, .75),  
               em.range=c(0.08, 0.18), b0.gam.range=c(1.4, 1.7), sd.proc=4, 
               sigma.group=1, nobs=nobs2, Tr=600, psi1=0.3, p11.range=c(0.4,0.6), 
               p10.range=c(0,0.1), b.range=c(0.1,0.2), 
               patch.surv.range=c(0.7,0.9), patch.col.range=c(0.1,0.6),
               noccs=5)
attach(datum)

thin <- ifelse(nobsyr == 1, 15, ifelse(nobsyr == 2, c(1,3,5,7,9,11,13,15), c(2,3,4,5,7,8,9,10,12,13,14,15)))
y.surv[,thin] <- as.integer(NA)

# Set up data for models

win.data <- list("nyears"=nyears, 
                 "ngroups"=ngroups2, "nobs"=nobs2, 
                 "em.group"=em.group, 
                 "y.group"=y.group, "event"=y.surv, 
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
                         p11=runif(nyears, 0.4, 0.6))}

tryCatch({
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=3, n.thin=3, n.iter=30000, 
            n.burnin=10000, jags.module = c("glm", "dic"))

out <- autojags(out, n.iter=10000, n.thin=3, n.update=5)

# Output for the pop level model
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
betas <- data.frame("b0.surv" = rep(out$BUGSoutput$mean$b0.surv, 14),
                    "surv.sd" = rep(out$BUGSoutput$sd$b0.surv, 14), 
                    "eps.surv" = out$BUGSoutput$mean$eps.surv[1:14], 
                    "eps.sd" = out$BUGSoutput$sd$eps.surv[1:14])

win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=area,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "betas"=betas, "T"=Tr)

out2 <- jags(win.data2, inits2, params2, "PopLevel_Sim.txt", n.chains=3, n.thin=3, n.iter=30000, 
             n.burnin=10000, jags.module = c("glm", "dic"))

out2 <- autojags(out2, n.iter=10000, n.thin=3, n.update=5)

results <- c(out2$BUGSoutput$summary[,c(1,2,3,7,8)], out$BUGSoutput$summary[1:15,c(1,2,3,7,8)], N.tot, G.mean, gamma.mean, s.prob)
names(results) <- c(rep("mean_Ntot", 15), rep("mean_s", 14), "mean_deviance", rep("mean_gamma", 14),
                    rep("sd_Ntot", 15), rep("sd_s", 14), "sd_deviance", rep("sd_gamma", 14),
                    rep("LCI_Ntot", 15), rep("LCI_s", 14), "LCI_deviance", rep("LCI_gamma", 14),
                    rep("UCI_Ntot", 15), rep("UCI_s", 14), "UCI_deviance", rep("UCI_gamma", 14),
                    rep("Rhat_Ntot", 15), rep("Rhat_s", 14), "Rhat_deviance", rep("Rhat_gamma", 14),
                    rep("mean_G", 15),rep("sd_G", 15),rep("LCI_G", 15),rep("UCI_G", 15),
                    rep("Rhat_G", 15),rep("Ntot", 15), rep("G", 15), rep("gamma", 14), rep("s", 14))

}, error= function(e) NULL)


detach()

return(tryCatch(results, error = function(e) NA))
}


#### Run the models ####


simoutput <- lapply(rep(50, 100), FUN = modelrun_fn, nobs2 = 10, nobsyr=1)
names(simoutput) <- 1:100

output_datum <- bind_rows(simoutput)
colnames(output) <- c(rep("mean_Ntot", 15), rep("mean_s", 14), "mean_deviance", rep("mean_gamma", 14),
                      rep("sd_Ntot", 15), rep("sd_s", 14), "sd_deviance", rep("sd_gamma", 14),
                      rep("LCI_Ntot", 15), rep("LCI_s", 14), "LCI_deviance", rep("LCI_gamma", 14),
                      rep("UCI_Ntot", 15), rep("UCI_s", 14), "UCI_deviance", rep("UCI_gamma", 14),
                      rep("Rhat_Ntot", 15), rep("Rhat_s", 14), "Rhat_deviance", rep("Rhat_gamma", 14),
                      rep("mean_G", 15),rep("sd_G", 15),rep("LCI_G", 15),rep("UCI_G", 15),
                      rep("Rhat_G", 15),rep("Ntot", 15), rep("G", 15), rep("gamma", 14), rep("s", 14))

write.csv(output_datum, "final_A.csv")




