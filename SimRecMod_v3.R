################################################################M

# Simulation code for group/pop level IPM with data generation
# Version 3
# 2/12/19

################################################################M


#### Function to create datasets ####

# Load package for categorical distribution to create false-positive encounter history
library(LaplacesDemon)

data.fn <- function(nyears=4, ngroups=12, s.prob.range=c(0.62, 0.67), 
                    em.range=c(-2.5, -1.5), b0.gam.range=c(1.4, 1.7), sd.proc=4, 
                    sigma.group=1, nobs=20, noccs=3, init.groups=50, group.ext=0.10, 
                    group.form=0.2, sigma.pop=50, ndisp=30, totgroups=400){
  
  # Function to generate data for IPM model to estimate recruitment using dynamic, false-positive occupancy
  # model and group counts with fixed survival. Annual varaition in parameters specified by uniform distr.
  
  # nyears = number of years
  # ngroups = number packs monitored
  # nobs = number of observations for survival
  # s.prob.range = range for survival probability
  # em.range = bounds of uniform distribution to draw from for emigration out of group
  # b0.gam.range = range for intercept/ mean recruitment of pups
  # sd.proc = sd of process error for group counts
  # sigma.group = sd of observation error for group counts
  # noccs = number of occasions for occupancy surveys
  
  
  # Arrays
  y.group <- array(dim=c(ngroups, nyears))
  G <- array(dim=c(totgroups, nyears))
  G.ad <- array(dim=c(totgroups, nyears))
  G.rec <- array(dim=c(totgroups, nyears))
  sigma.proc.group <- array(dim=c(nyears))
  y.temp <- array(dim=c(ngroups, nyears))
  G.mean <- array(dim=c(nyears))
  N.tot <- array(dim=c(nyears))
  N.rec <- array(dim=c(nyears))
  y.surv <- array(dim=c(nobs,nyears))
  y.temp2 <- array(dim=c(nyears))
  y.pop <- array(dim=c(nyears))
  N <- array(dim=c(nyears))
  y.disp <- array(dim=c(ndisp,nyears))
  P <- array(dim=c(nyears))
  em.group <- array(dim=c(nyears))
  
  # Recruitment and group size
  # Recruitment
  b0.gam <- runif(n=nyears-1, min=b0.gam.range[1], max=b0.gam.range[2])
  
  gamma <- exp(b0.gam)
  
  # Random noise for first year group size
  for(i in 1:totgroups){
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
  b0.disp <- runif(n=nyears-1, min=em.range[1], max=em.range[2])
  
  # Dispersal
  for(k in 1:(nyears-1)){
    em.group[k] <- 1/(1 + exp(-b0.disp[k]))
  }
  
  # Group size in successive years
  for(i in 1:totgroups){
    for(k in 2:nyears){
      G.ad[i,k] <- rbinom(1, G[i,k-1], s.prob[k-1]*(1-em.group[k-1]))
      G.rec[i,k] <- rpois(1, gamma[k-1])
      G[i,k] <- G.ad[i,k] + G.rec[i,k]
      # G[i,k] <- ifelse(G[i,k-1]==0, 0, rpois(1, (G[i,k-1]*s.prob[k-1]*(1-em.group[k-1])+gamma[k-1]))) #+sigma.proc.group[k-1,r])))
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
    G.mean[k] <- mean(G[1:ngroups,k])
  }
  
  # survival data
  for(k in 1:(nyears-1)){
    for(i in 1:nobs){
      y.surv[i,k] <- rbinom(1,1,(1-s.prob[k]))
    }
  }
  
  # Dispersal data
  for(k in 1:(nyears-1)){
    for(i in 1:ndisp){
      y.disp[i,k] <- rbinom(1,1,em.group[k])
    }
  }
  
  # Population level
  P[1] <- init.groups
  N[1] <- sum(G[1:P[1],1])
  for(k in 2:nyears){
    P[k] <- round(P[k-1] + P[k-1] * group.form - P[k-1] * group.ext)
    N[k] <- sum(G[1:P[k],k])
  }
  
  # Pop counts
  for(k in 1:nyears){
    y.temp2[k] <- round(rnorm(1, N[k], sigma.pop))
    y.temp2[k] <- ifelse(y.temp2[k] < 0, 0, y.temp2[k])
    y.pop[k] <- y.temp2[k] #ifelse(y.temp2[k] > N[k], N[k], y.temp2[k])
  }
  
  # Observation error
  var.pop <- sigma.pop*sigma.pop
  
  
  
  
  return(list(nyears=nyears, ngroups=ngroups, nobs=nobs, noccs=noccs, ndisp=ndisp,
              s.prob=s.prob, em.group=em.group[1:14], y.surv=y.surv, y.disp=y.disp,
              y.group=y.group, b0.gam=b0.gam, G.mean=G.mean, 
              gamma.mean=gamma, var.group=var.group, N=N, y.pop=y.pop,  
              var.pop=var.pop, group.form=group.form, group.ext=group.ext, 
              init.groups=init.groups))
}


#### Load packages and set wd ####

library(dplyr)
library(tidyverse)
library(R2jags)
library(snowfall)
memory.limit(size = 7500000) #just shy of 8 tb limit i think


setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults")

#### Model code ####
sink("SimIPMmodel_scen1.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:(nyears-1)){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    # Initial group sizes
    
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    
    # Process and observation error
    
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.gam[k] ~ dunif(-10,10)
    }
    
    
    ## 1.5 Dispersal priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.disp[k] ~ dunif(-10,10)
    }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dpois(mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    for(k in 1:(nyears-1)){
    logit(em.group[k]) <- b0.disp[k]
    
    for(i in 1:ndisp){
    y.disp[i,k] ~ dbern(em.group[k])
    }
    }
    
    
    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:(nyears-1)){
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
    G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+0.001))T(0,)
    }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
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
    gamma[i,k] ~ dpois(mu.gamma[i,k])T(0,)
    }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()







sink("SimIPMmodel_scen2.txt")
cat("
    model{
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    # Initial group sizes
    
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    
    # Process and observation error
    
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.gam[k] ~ dunif(-10,10)
    }
    
    
    ## 1.5 Dispersal priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    em.group[k] ~ dbeta(alpha, beta)
    }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dpois(mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    # for(k in 1:(nyears-1)){
    # logit(em.group[k]) <- b0.disp[k]
    
    # for(i in 1:ndisp){
    # y.disp[i,k] ~ dbern(em.group[k])
    # }
    # }
    
    
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
    G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+0.001))T(0,)
    }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
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
    gamma[i,k] ~ dpois(mu.gamma[i,k])T(0,)
    }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    }
    " , fill=TRUE)
sink()




sink("SimIPMmodel_scen3.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    # Initial group sizes
    
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    
    # Process and observation error
    
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.gam[k] ~ dunif(-10,10)
    }
    
    
    ## 1.5 Dispersal priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.disp[k] ~ dunif(-10,10)
    }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dpois(mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    for(k in 1:(nyears-1)){
    logit(em.group[k]) <- b0.disp[k]
    
    for(i in 1:ndisp){
    y.disp[i,k] ~ dbern(em.group[k])
    }
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
    G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+0.001))T(0,)
    }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
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
    gamma[i,k] ~ dpois(mu.gamma[i,k])T(0,)
    }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()






sink("SimIPMmodel_scen4.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    # Initial group sizes
    
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    
    # Process and observation error
    
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.gam[k] ~ dunif(-10,10)
    }
    
    
    ## 1.5 Dispersal priors
  
    for(k in 1:(nyears-1)){
    em.group[k] ~ dbeta(alpha, beta)
    }
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dpois(mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    # for(k in 1:(nyears-1)){
    # logit(em.group[k]) <- b0.disp[k]
    # 
    # for(i in 1:ndisp){
    # y.disp[i,k] ~ dbern(em.group[k])
    # }
    # }
    
    
    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:(nyears-1)){
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
    G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+0.001))T(0,)
    }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
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
    gamma[i,k] ~ dpois(mu.gamma[i,k])T(0,)
    }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()




#### v2 Pop only models ####
sink("SimIPMmodel_scen1pop.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:(nyears-1)){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(k in 1:nyears){
    G.mean[k] <- mean(y.group[,k])
    }

    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    b0.gam ~ dnorm(0,0.001)

    for(k in 1:(nyears-1)){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }

    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    
    ## 1.5 Dispersal priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.disp[k] ~ dunif(-10,10)
    }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] * (1-em.group[k-1]) + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dnorm(mu.N[k], 1/mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    for(k in 1:(nyears-1)){
    logit(em.group[k]) <- b0.disp[k]
    
    for(i in 1:ndisp){
    y.disp[i,k] ~ dbern(em.group[k])
    }
    }
    
    
    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:(nyears-1)){
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
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(k in 1:(nyears-1)){
    mu.gamma[k] <- exp(b0.gam + eps.gam[k])
    gamma.mean[k] ~ dpois(mu.gamma[k])T(0,)
    }#k

    
    

    }", fill=TRUE)
sink()







sink("SimIPMmodel_scen2pop.txt")
cat("
    model{
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    for(k in 1:nyears){
    G.mean[k] <- mean(y.group[,k])
    }
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    b0.gam ~ dnorm(0,0.001)

    for(k in 1:(nyears-1)){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    
    ## 1.5 Dispersal priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    em.group[k] ~ dbeta(alpha, beta)
    }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] * (1-em.group[k-1]) + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dnorm(mu.N[k], 1/mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    # for(k in 1:(nyears-1)){
    # logit(em.group[k]) <- b0.disp[k]
    
    # for(i in 1:ndisp){
    # y.disp[i,k] ~ dbern(em.group[k])
    # }
    # }
    
    
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
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(k in 1:(nyears-1)){
    mu.gamma[k] <- exp(b0.gam + eps.gam[k])
    gamma.mean[k] ~ dpois(mu.gamma[k])T(0,)
    }#k

    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    }
    " , fill=TRUE)
sink()




sink("SimIPMmodel_scen3pop.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    # Initial group sizes
    
    for(k in 1:nyears){
    G.mean[k] <- mean(y.group[,k])
    }
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    b0.gam ~ dnorm(0,0.001)

    for(k in 1:(nyears-1)){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    
    ## 1.5 Dispersal priors
    
    # Priors for beta coefficients
    
    for(k in 1:(nyears-1)){
    b0.disp[k] ~ dunif(-10,10)
    }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] * (1-em.group[k-1]) + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dnorm(mu.N[k], 1/mu.N[k])
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)T(0,10000)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    for(k in 1:(nyears-1)){
    logit(em.group[k]) <- b0.disp[k]
    
    for(i in 1:ndisp){
    y.disp[i,k] ~ dbern(em.group[k])
    }
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
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(k in 1:(nyears-1)){
    mu.gamma[k] <- exp(b0.gam + eps.gam[k])
    gamma.mean[k] ~ dpois(mu.gamma[k])T(0,)
    }#k

    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()






sink("SimIPMmodel_scen4pop.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Pop Level
    
    # Initial population size
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)
    
    # Process and observation error
    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    
    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    # Initial group sizes
    
    for(k in 1:nyears){
    G.mean[k] <- mean(y.group[,k])
    }
    
    
    
    ## 1.4 Recruitment priors
    
    # Priors for beta coefficients
    
    b0.gam ~ dnorm(0,0.001)

    for(k in 1:(nyears-1)){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    
    ## 1.5 Dispersal priors
    
    for(k in 1:(nyears-1)){
    em.group[k] ~ dbeta(alpha, beta)
    }
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Populaiton level 
    
    ####################
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] * (1-em.group[k-1]) + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dnorm(mu.N[k], 1/mu.N[k])T(0,)
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    
    #####################
    
    # 2.2. Dispersal 
    
    ####################    
    
    # for(k in 1:(nyears-1)){
    # logit(em.group[k]) <- b0.disp[k]
    # 
    # for(i in 1:ndisp){
    # y.disp[i,k] ~ dbern(em.group[k])
    # }
    # }
    
    
    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:(nyears-1)){
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
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(k in 1:(nyears-1)){
    mu.gamma[k] <- exp(b0.gam + eps.gam[k])
    gamma.mean[k] ~ dpois(mu.gamma[k])T(0,)
    }#k

    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### Constants ####

# Parameters to keep track of and report
params <- c("P", "G.mean", "gamma.mean", "annual.s",
            "em.group", 
            "N") 


#### Function to run models ####

beta.MoM.fcn <- function (mean, var, param){
  alpha <-mean*((mean*(1-mean)/var)-1)
  beta <-(1-mean)*((mean*(1-mean)/var)-1)
  y <-rbeta(1:1000, alpha, beta)
  beta.plot <-plot(y, dbeta(y, alpha, beta), ylab="Frequency", xlab=param)
  return(c(list(alpha=alpha, beta=beta)))
}


modelrun_fn <- function(gam.range=c(1.3, 1.7), 
                        modelfile="SimIPMmodel_scen4pop.txt", 
                        ngroups2=50, 
                        nobs2=20, nobsyr=1){
  
  # Clear last runs of anything
  # detach()
  # rm(out, results, datum, y.surv)
  
  # Create and attach data
  datum <- data.fn(nyears=15, ngroups=ngroups2, s.prob.range=c(.55, .75),  
                   em.range=c(-2.5, -1.5), b0.gam.range=gam.range, sd.proc=4, 
                   sigma.group=1, nobs=nobs2, ndisp=50, init.groups=50, group.ext=0.10, 
                   group.form=0.2, sigma.pop=50)

  thin <- ifelse(nobsyr == 1, 15, 
                 ifelse(nobsyr == 2, c(1,3,5,7,9,11,13,15), 
                        c(2,3,4,5,7,8,9,10,12,13,14,15)))
  datum$y.surv[,thin] <- as.integer(NA)
  
  # Set up data for models
  shapes <- ifelse(modelfile=="SimIPMmodel_scen4pop.txt" | 
                     modelfile=="SimIPMmodel_scen4.txt", 
                   list(beta.MoM.fcn(mean(datum$em.group, na.rm=T) * 0.8, var(datum$em.group, na.rm=T), "Survival")), 
                   list(beta.MoM.fcn(mean(datum$em.group, na.rm=T), var(datum$em.group, na.rm=T), "Survival")))
  
  win.data <- list("nyears"=datum$nyears, "ngroups"=datum$ngroups, "nobs"=datum$nobs, 
                   "y.group"=datum$y.group, "event"=datum$y.surv, "ndisp"=datum$ndisp, 
                   "init.groups"=datum$init.groups, "y.pop"=datum$y.pop,
                   "group.form"=datum$group.form, "group.ext"=datum$group.ext, 
                   "y.disp"=datum$y.disp, "alpha"=shapes[[1]][1], "beta"=shapes[[1]][2])
  
  #  Initial Values	
  inits <- function(){list(b0.gam=runif(1, gam.range[1], gam.range[2]), 
                           b0.disp=runif(datum$nyears-1, -2.5, -1.5))}
  
  tryCatch({
    out <- jags(win.data, inits, params, modelfile, n.chains=3, n.thin=3, n.iter=30000, 
                n.burnin=10000, jags.module = c("glm", "dic"))
    
    out <- autojags(out, n.iter=20000, n.thin=3, n.update=3)
    
    
    results <- c(out$BUGSoutput$summary[,c(1,2,3,7,8)], datum$N, datum$G.mean, datum$gamma.mean, datum$s.prob, datum$em.group)
    names(results) <- c(rep("mean_G", 15),rep("mean_Ntot", 15), rep("mean_P", 15), rep("mean_s", 14), 
                        "mean_deviance", rep("mean_disp", 14), rep("mean_gamma", 14),
                        rep("sd_G", 15),rep("sd_Ntot", 15), rep("sd_P", 15), rep("sd_s", 14), 
                        "sd_deviance", rep("sd_disp", 14), rep("sd_gamma", 14),
                        rep("LCI_G", 15),rep("LCI_Ntot", 15), rep("LCI_P", 15), rep("LCI_s", 14), 
                        "LCI_deviance", rep("LCI_disp", 14), rep("LCI_gamma", 14),
                        rep("UCI_G", 15),rep("UCI_Ntot", 15), rep("UCI_P", 15), rep("UCI_s", 14), 
                        "UCI_deviance", rep("UCI_disp", 14), rep("UCI_gamma", 14),
                        rep("Rhat_G", 15),rep("Rhat_Ntot", 15), rep("Rhat_P", 15), rep("Rhat_s", 14), 
                        "Rhat_deviance", rep("Rhat_disp", 14), rep("Rhat_gamma", 14),
                        rep("Ntot", 15), rep("G", 15), rep("gamma", 14), rep("s", 14), rep("em", 14))
    
  }, error= function(e) NULL)
  
  
  
  return(tryCatch(results, error = function(e) NA))
}



#### Run the models and save output ####

# Scen 1 is normal variation (1.3-1.7), dispersal data
# Scen 2 is no dispersal data, accurate dispersal mean
# Scen 3 is high variation (1-2.2), dispersal data
# Scen 4 is no dispersal data, 20% biased dispersal mean

start.time <- Sys.time()
simoutput <- lapply(rep("SimIPMmodel_scen1pop.txt", 250), FUN = modelrun_fn, 
                    gam.range=c(1.3, 1.7), ngroups2=50, nobs2=30, nobsyr=1)
end.time <- Sys.time()

end.time-start.time


fin2<-simoutput[-(which(is.na(simoutput), arr.ind=TRUE))]


save(fin, file ="list_scen4pop.RData")



#### Manipulate results and summarize ####

group.datum <- list()
n.datum <- list()
s.datum <- list()
mod.dev <- list()
em.datum <- list()
gam.datum <- list()

for(i in 1:100){
  group.datum[[i]] <- data.frame(Year=1:15, "G.mean" = simoutput[[i]][1:15], "G.mean.sd"=simoutput[[i]][89:103],
                            "LL"=simoutput[[i]][177:191], "UL"=simoutput[[i]][265:279],
                            "Rhat" = simoutput[[i]][353:367], "truth" = simoutput[[i]][456:470])
  n.datum[[i]] <- data.frame(Year=1:15, "mean" = simoutput[[i]][16:30], "sd"=simoutput[[i]][104:118],
                        "LL"=simoutput[[i]][192:206], "UL"=simoutput[[i]][280:294],
                        "Rhat" = simoutput[[i]][368:382], "truth" = simoutput[[i]][441:455])
  s.datum[[i]] <- data.frame(Year=1:14, "mean" = simoutput[[i]][46:59], "sd"=simoutput[[i]][134:147],
                        "LL"=simoutput[[i]][222:235], "UL"=simoutput[[i]][310:323],
                        "Rhat" = simoutput[[i]][398:411], "truth" = simoutput[[i]][485:498])
  mod.dev[[i]] <- data.frame(Mod="scen1", "mean" = simoutput[[i]][60], "sd"=simoutput[[i]][148],
                        "LL"=simoutput[[i]][236], "UL"=simoutput[[i]][324],
                        "Rhat" = simoutput[[i]][412])
  em.datum[[i]] <- data.frame(Year=1:14, "mean" = simoutput[[i]][61:74], "sd"=simoutput[[i]][149:162],
                         "LL"=simoutput[[i]][237:250], "UL"=simoutput[[i]][325:338],
                         "Rhat" = simoutput[[i]][413:426], "truth" = simoutput[[i]][499:512])
  gam.datum[[i]] <- data.frame(Year=1:14, "mean" = simoutput[[i]][75:88], "sd"=simoutput[[i]][163:176],
                          "LL"=simoutput[[i]][251:264], "UL"=simoutput[[i]][339:352],
                          "Rhat" = simoutput[[i]][427:440], "truth" = simoutput[[i]][471:484])

}


### Check Rhat values
# Single dataframe
em.datum[[1]]$Rhat[em.datum[[1]]$Rhat > 1.1]

# Function to run for all dataframes
rhat.extract <- function(x){
  tmp <- x$Rhat[x$Rhat >= 1.1]
  out <- length(tmp) > 0 
  return(out)
}

# Run function for all datafames
sum(unlist(lapply(s.datum, rhat.extract)))
sum(unlist(lapply(n.datum, rhat.extract)))
sum(unlist(lapply(gam.datum, rhat.extract)))
sum(unlist(lapply(em.datum, rhat.extract)))
sum(unlist(lapply(group.datum, rhat.extract)))


### PRECISION: Calculate CV for each estimate for each year and 
### create new column in dataframes
# Single dataframe
em.datum[[1]]$CV <- 100*em.datum[[1]]$sd/em.datum[[1]]$mean


# Function to run for all dataframes
cv.create <- function(x){
  x$CV <- 100 * x$G.mean.sd / x$G.mean
  return(x)
}

n.datum <- lapply(n.datum, cv.create)
s.datum <- lapply(s.datum, cv.create)
gam.datum <- lapply(gam.datum, cv.create)
em.datum <- lapply(em.datum, cv.create)
group.datum <- lapply(group.datum, cv.create)

# Calculate summary statistics for each scenario
CV.results <- data.frame("Var"="EM", 
                         "cv.mean"=mean(unlist(lapply(em.datum, function(z) z$CV))),
                         "cv.sd"=sd(unlist(lapply(em.datum, function(z) z$CV))))
CV.results <- CV.results %>%
  bind_rows(data.frame("Var"="S", 
             "cv.mean"=mean(unlist(lapply(s.datum, function(z) z$CV))),
             "cv.sd"=sd(unlist(lapply(s.datum, function(z) z$CV))))) 
CV.results <- CV.results %>%
  bind_rows(data.frame("Var"="GAM", 
                       "cv.mean"=mean(unlist(lapply(gam.datum, function(z) z$CV))),
                       "cv.sd"=sd(unlist(lapply(gam.datum, function(z) z$CV)))))
CV.results <- CV.results %>%
  bind_rows(data.frame("Var"="N", 
                       "cv.mean"=mean(unlist(lapply(n.datum, function(z) z$CV))),
                       "cv.sd"=sd(unlist(lapply(n.datum, function(z) z$CV)))))
CV.results <- CV.results %>%
  bind_rows(data.frame("Var"="GROUP", 
                       "cv.mean"=mean(unlist(lapply(group.datum, function(z) z$CV))),
                       "cv.sd"=sd(unlist(lapply(group.datum, function(z) z$CV)))))


### ACCURACY: Calculate SMAE for each year and simulation of variables and
### create 1 datafame to hold all values (1 for each variable and scenario)
# Single dataframe
sum(abs(em.datum[[2]]$mean - em.datum[[2]]$truth)/em.datum[[2]]$truth)*1/14

# Function to run for all dataframes
smae.calc <- function(x){
  SMAE <- abs(x$G.mean - x$truth)/x$truth
  return(SMAE)
}

SMAE.results <- data.frame("Var"="EM", "Scen"=3, "Social"=1, 
                           "smae"=sum(unlist(lapply(em.datum, smae.calc)))*1/(14*100))
SMAE.results <- SMAE.results %>% 
  bind_rows(data.frame("Var"="GAM", "Scen"=3, "Social"=1, 
                       "smae"=sum(unlist(lapply(gam.datum, smae.calc)))*1/(14*100)))
SMAE.results <- SMAE.results %>% 
  bind_rows(data.frame("Var"="S", "Scen"=3, "Social"=1, 
                       "smae"=sum(unlist(lapply(s.datum, smae.calc)))*1/(14*100)))
SMAE.results <- SMAE.results %>% 
  bind_rows(data.frame("Var"="N", "Scen"=3, "Social"=1, 
                       "smae"=sum(unlist(lapply(n.datum, smae.calc)))*1/(15*100)))
SMAE.results <- SMAE.results %>% 
  bind_rows(data.frame("Var"="GROUP", "Scen"=3, "Social"=1, 
                       "smae"=sum(unlist(lapply(group.datum, smae.calc)))*1/(15*100)))


### BIAS: Calculate PAR for each year and simulation and 
### create 1 dataframe to hold all values (1 for each variable and scenrio)
# Single dataframe
sum(100 * em.datum[[1]]$mean / em.datum[[1]]$truth) * 1/14
# sum((em.datum[[1]]$mean - em.datum[[1]]$truth) / em.datum[[1]]$truth) * 1/14

# Function to run for all dataframes
par.calc <- function(x){
  PAR <- 100 * x$G.mean / x$truth
  return(PAR)
}

PAR.results <- data.frame("Var"="EM", "Scen"=3, "Social"=1, 
                           "par"=sum(unlist(lapply(em.datum, par.calc)))*1/(14*23))
PAR.results <- PAR.results %>% 
  bind_rows(data.frame("Var"="GAM", "Scen"=3, "Social"=1, 
                       "par"=sum(unlist(lapply(gam.datum, par.calc)))*1/(14*23)))
PAR.results <- PAR.results %>% 
  bind_rows(data.frame("Var"="S", "Scen"=3, "Social"=1, 
                       "par"=sum(unlist(lapply(s.datum, par.calc)))*1/(14*23)))
PAR.results <- PAR.results %>% 
  bind_rows(data.frame("Var"="N", "Scen"=3, "Social"=1, 
                       "par"=sum(unlist(lapply(n.datum, par.calc)))*1/(15*23)))
PAR.results <- PAR.results %>% 
  bind_rows(data.frame("Var"="GROUP", "Scen"=3, "Social"=1, 
                       "par"=sum(unlist(lapply(group.datum, par.calc)))*1/(15*23)))


### Combine into 1 dataframe for all

TRIALscen3.results <- merge(PAR.results, SMAE.results, by=c("Var", "Social", "Scen"))
TRIALscen3.results <- merge(TRIALscen3.results, CV.results, by="Var")


### Save output summaries
write.csv(TRIALscen3.results, file="scen3TRIAL1.csv")


#### Graphing ####

# Bring in data
scen1 <- read.csv("scen1.csv")[,2:8]
scen2 <- read.csv("scen2.csv")[,2:8]
scen3 <- read.csv("scen3.csv")[,2:8]
scen4 <- read.csv("scen4.csv")[,2:8]
scen1pop <- read.csv("scen1pop.csv")[,2:8]
scen2pop <- read.csv("scen2pop.csv")[,2:8]
scen3pop <- read.csv("scen3pop.csv")[,2:8]
scen4pop <- read.csv("scen4pop.csv")[,2:8]


# Load packages
library(ggplot2)
library(RColorBrewer)

# Combine data
datum <- bind_rows(scen1, scen2, scen3, scen4, 
                   scen1pop, scen2pop, scen3pop, scen4pop)
str(datum)
datum$Var <- as.factor(datum$Var)
datum$Social <- as.factor(datum$Social)
datum$Scen <- as.factor(datum$Scen)


# Plot for Bias
labs2 <- c("EM"="Dispersal", 
           "GAM"="Recruitment",
           "N"="Abundance",
           "S"="Survival")
datum2<-droplevels(datum[datum$Var != "GROUP",])
datum2$Scen <- factor(datum2$Scen, level=c(1,3,2,4))

ggplot(data=datum2, aes(x=Scen, y=par))+
  geom_point(aes(color=Social), 
             position=position_dodge(width=0.5),
             size=2.5)+
  facet_wrap(~Var, ncol=2, labeller=labeller(Var=labs2))+
  geom_hline(yintercept=100, linetype="dashed", 
             color = "#666666", size=1.2)+
  scale_x_discrete(name="Dispersal data scenarios", 
                   breaks=c(1, 3, 2, 4),
                   labels=c("Data\naverage\nvariation", "Data\nhigh\nvariation",
                            "No data\naccurate", "No data\n20% biased"))+
  scale_y_continuous(name="Bias")+
  scale_color_brewer(palette="Dark2", 
                     labels=c("Without Social Structure", 
                              "With Social Structure"))+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.key = element_rect(fill= "white"),
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        # legend.key.size = unit(1.5, 'lines'),
        strip.background=element_rect(color="black", fill="#666666"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

# Plot for Accuracy
ggplot(data=datum2, aes(x=Scen, y=smae))+
  geom_point(aes(color=Social), 
             position=position_dodge(width=0.5),
             size=2.5)+
  facet_wrap(~Var, ncol=2, labeller=labeller(Var=labs2))+
  scale_x_discrete(name="Dispersal data scenarios", 
                   breaks=c(1, 3, 2, 4),
                   labels=c("Data\naverage\nvariation", "Data\nhigh\nvariation",
                            "No data\naccurate", "No data\n20% biased"))+
  scale_y_continuous(name="Accuracy")+
  scale_color_brewer(palette="Dark2", 
                     labels=c("Without Social Structure", 
                              "With Social Structure"))+
  theme(legend.position = "bottom",
        # legend.justification = c(1,0),
        legend.key = element_rect(fill= "white"),
        # legend.direction = "vertical",
        # legend.key.size = unit(1.5, 'lines'),
        strip.background=element_rect(color="black", fill="#666666"),
        strip.text = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=11))

# Precision
ggplot(data=datum2, aes(x=Scen, y=cv.mean))+
  geom_point(aes(color=Social), 
             position=position_dodge(width=0.5),
             size=2.5)+
  facet_wrap(~Var, ncol=2, labeller=labeller(Var=labs2))+
  scale_x_discrete(name="Dispersal data scenarios", 
                   breaks=c(1, 3, 2, 4),
                   labels=c("Data\naverage\nvariation", "Data\nhigh\nvariation",
                            "No data\naccurate", "No data\n20% biased"))+
  scale_y_continuous(name="Precision", breaks=c(seq(0,150,by=25)))+
  scale_color_brewer(palette="Dark2", 
                     labels=c("Without Social Structure", 
                              "With Social Structure"))+
  geom_pointrange(aes(ymin=cv.mean-cv.sd, ymax=cv.mean+cv.sd, colour=Social), 
                  position=position_dodge(width=.5))+
  theme(legend.position = "bottom",
        # legend.justification = c(1,0),
        legend.key = element_rect(fill= "white"),
        # legend.direction = "vertical",
        # legend.key.size = unit(1.5, 'lines'),
        legend.text = element_text(size=12),
        strip.background=element_rect(color="black", fill="#666666"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

#### Graphing Trials ####

# Simple plot for Bias
ggplot(data=datum[datum$Var=="EM",], aes(x=Scen, y=par))+
  geom_point(aes(color=Social), 
             position=position_dodge(width=0.5),
             size=2)+
  geom_hline(yintercept=100, linetype="dashed", 
             color = "#666666", size=1.2)+
  scale_x_discrete(name="Dispersal data scenarios", 
                   breaks=c(1, 3, 2, 4),
                   labels=c("Dispersal data\n average variation", "Dispersal data\n high variation",
                            "No dispersal data\n accurate", "No dispersal data\n 20% biased"))+
  scale_y_continuous(name="Bias", breaks=c(80, 85, 90, 95, 100, 105, 110))+
  scale_color_brewer(palette="Dark2", labels=c("With\nSocial Structure", "Without\nSocial Structure"))+
  theme(legend.position = c(.99,.99),
        legend.justification = c(1,1),
        legend.key = element_rect(fill= "white"),
        legend.key.size = unit(1.5, 'lines'),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))


# Social structure facets
labs <- c("0"="Without Social Structure", "1"="With Social Structrue")

ggplot(data=datum, aes(x=Scen, y=par))+
  geom_point(aes(color=Var), 
             position=position_dodge(width=0.5),
             size=2.5)+
  facet_grid(Social~., labeller = labeller(Social=labs))+
  geom_hline(yintercept=100, linetype="dashed", 
             color = "#666666", size=1.2)+
  scale_x_discrete(name="Dispersal data scenarios", 
                   breaks=c(1, 3, 2, 4),
                   labels=c("Dispersal data\n average variation", "Dispersal data\n high variation",
                            "No dispersal data\n accurate", "No dispersal data\n 20% biased"))+
  scale_y_continuous(name="Bias")+
  scale_color_brewer(palette="Dark2", 
                     labels=c("Dispersal", 
                              "Recruitment", 
                              "Group Size", 
                              "Abundance",
                              "Survival"))+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.key = element_rect(fill= "white"),
        legend.direction = "horizontal",
        # legend.key.size = unit(1.5, 'lines'),
        strip.background=element_rect(color="black", fill="white"),
        strip.text = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=11))


# Var facets
labs2 <- c("EM"="Dispersal", 
           "GAM"="Recruitment",
           "N"="Abundance",
           "S"="Survival")
datum2<-droplevels(datum[datum$Var != "GROUP",])
datum2$Scen <- factor(datum2$Scen, level=c(1,3,2,4))

ggplot(data=datum2, aes(x=Scen, y=par))+
  geom_point(aes(color=Social), 
             position=position_dodge(width=0.5),
             size=2.5)+
  facet_wrap(~Var, ncol=2, labeller=labeller(Var=labs2))+
  geom_hline(yintercept=100, linetype="dashed", 
             color = "#666666", size=1.2)+
  scale_x_discrete(name="Dispersal data scenarios", 
                   breaks=c(1, 3, 2, 4),
                   labels=c("Data\naverage\nvariation", "Data\nhigh\nvariation",
                            "No data\naccurate", "No data\n20% biased"))+
  scale_y_continuous(name="Bias")+
  scale_color_brewer(palette="Dark2", 
                     labels=c("Without Social Structure", 
                              "With Social Structure"))+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.key = element_rect(fill= "white"),
        legend.direction = "vertical",
        # legend.key.size = unit(1.5, 'lines'),
        strip.background=element_rect(color="black", fill="#666666"),
        strip.text = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=11))




#### Testing the prior ####

sink("SimIPMmodel_prior.txt")
cat("
    model {
    
    ##########################
    
    #        1. Priors
    
    ##########################
    
    
    ## 1.4 Recruitment priors
    
    for(k in 1:(nyears-1)){
    b0.gam[k] ~ dunif(-10,10)
    }

    ## 1.1 Pop Level
    P[1] ~ dpois(50)
    N[1] ~ dpois(P[1]*7)

    tauy.pop <- pow(sigma.pop, -2)
    sigma.pop ~ dunif(0,100)
    var.pop <- pow(sigma.pop, 2)
    

    ## 1.2 Survival priors
    
    b0.surv ~ dnorm(0,0.001)
    
    for(k in 1:(nyears-1)){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ## 1.3 Group priors
    
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }

    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    ## 1.5 Dispersal priors

    for(k in 1:(nyears-1)){
    b0.disp[k] ~ dunif(-10,10)
    }
    
    
    #############################
    
    #        2. Likelihoods
    
    #############################

    
    ## 2.1. Populaiton level 
    
    for(k in 2:nyears){
    P[k] <- P[k-1] + P[k-1] * group.form - P[k-1] * group.ext
    
    mu.N[k] <- P[k-1] * gamma.mean[k-1] + N[k-1] * annual.s[k-1] + (P[k] - P[k-1]) * G.mean[k]
    N[k] ~ dpois(mu.N[k])T(0,)
    }
    
    for(k in 1:nyears){
    y.pop[k] ~ dnorm(N[k], tauy.pop)
    }
    

    ## 2.2. Dispersal 

    for(k in 1:(nyears-1)){
    logit(em.group[k]) <- b0.disp[k]
    
    for(i in 1:ndisp){
    y.disp[i,k] ~ dbern(em.group[k])
    }
    }
    
    
    ## 2.3. Survival likelihood 

    for(k in 1:(nyears-1)){
    for(i in 1:nobs){
    event[i,k] ~ dbern(mu.surv[i,k])
    cloglog(mu.surv[i,k]) <- b0.surv + eps.surv[k]
    }
    }
    
    for(k in 1:(nyears-1)){
    cloglog(mu.pred[k]) <- b0.surv + eps.surv[k]
    hazard[k] <- -log(1-mu.pred[k])
    }

    for(k in 1:(nyears-1)){
    H[k] <- hazard[k]
    annual.s[k] <- exp(-H[k])
    }

    
    ## 2.4. Group level counts likelihood 

    for(i in 1:ngroups){
    for(k in 2:nyears){
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group[k-1]) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+0.001))T(0,)
    }#k
    }#i

    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }#k
    }#i

    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
    }#k
    
    for(k in 1:(nyears-1)){
    gamma.mean[k] <- mean(gamma[,k])
    }
    

    ## 2.5. Recruitment model 
    
    for(i in 1:ngroups){
    for(k in 1:(nyears-1)){
    mu.gamma[i,k] <- exp(b0.gam[k])
    gamma[i,k] ~ dpois(mu.gamma[i,k])T(0,)
    }#k
    }#i
    

    
    }", fill=TRUE)
sink()


datum <- data.fn(nyears=5, ngroups=50, s.prob.range=c(.55, .75),  
                 em.range=c(-2.5, -1.5), b0.gam.range=c(1.3, 1.7), sd.proc=4, 
                 sigma.group=1, nobs=30, ndisp=50, init.groups=50, group.ext=0.10, 
                 group.form=0.2, sigma.pop=50)

datum$y.surv[,5] <- as.integer(NA)

params <- c("b0.gam") 

win.data <- list("nyears"=datum$nyears, "ngroups"=datum$ngroups, "nobs"=datum$nobs, 
                 "y.group"=datum$y.group, "event"=datum$y.surv, "ndisp"=datum$ndisp, 
                 "init.groups"=datum$init.groups, "y.pop"=datum$y.pop,
                 "group.form"=datum$group.form, "group.ext"=datum$group.ext, 
                 "y.disp"=datum$y.disp)

#  Initial Values	
inits <- function(){list(b0.gam=runif(datum$nyears-1, 1.3, 1.7), 
                         b0.disp=runif(datum$nyears-1, -2.5, -1.5))}

out <- jags(win.data, inits, params, "SimIPMmodel_prior3.txt", n.chains=3, n.thin=3, n.iter=30000, 
              n.burnin=10000, jags.module = c("glm", "dic"))
  
out <- autojags(out, n.iter=20000, n.thin=3, n.update=3)
  

out_mcmc3 <- as.mcmc(out)

PR <- runif(20000, -2, 2)

MCMCtrace(out_mcmc2, params="b0.gam", 
          priors=PR, gvals = truth.gam, 
          pdf=FALSE)

least <- MCMCchains(out_mcmc, params="b0.gam")
mid <- MCMCchains(out_mcmc2, params="b0.gam")
most <- MCMCchains(out_mcmc3, params="b0.gam")

dens.least <- density(least[,3])
dens.mid <- density(mid[,3])
dens.most <- density(most[,3])

dd <- with(dens.least, data.frame(x,y))
dd2 <- with(dens.mid, data.frame(x,y))
dd3 <- with(dens.most, data.frame(x,y))

colnames(dd) <- c("Val", "Dens")
colnames(dd2) <- c("Val", "Dens")
colnames(dd3) <- c("Val", "Dens")

dens.datum <- bind_rows(dd,dd2,dd3) %>% 
  mutate("Type" = c(rep("Least", nrow(dd)),rep("Mid", nrow(dd2)),rep("Most", nrow(dd3))))

P1 <- runif(20000, -10, 10)
densP1 <- density(P1)
P2 <- runif(20000, -2, 2)
densP2 <- density(P2)
P3 <- rnorm(20000, 0, 0.5)
densP3 <- density(P3)

pdat1 <- with(densP1, data.frame(x,y))
pdat2 <- with(densP2, data.frame(x,y))
pdat3 <- with(densP3, data.frame(x,y))
colnames(pdat1) <- c("Val", "Dens")
colnames(pdat2) <- c("Val", "Dens")
colnames(pdat3) <- c("Val", "Dens")

p.datum <- bind_rows(pdat1,pdat2,pdat3) %>% 
  mutate("Type" = c(rep("Least", nrow(dd)),rep("Mid", nrow(dd2)),rep("Most", nrow(dd3))), 
         "Thing"=rep("Prior", nrow(pdat1)*3))
dens.datum <- dens.datum %>% 
  mutate("Thing"=rep("Post", nrow(dens.datum))) %>%
  bind_rows(p.datum)


dat.text <- data.frame(labels=c("6% overlap", "18% overlap", "20% overlap"), 
                       Type=c("Least", "Mid", "Most"), 
                       Thing=rep("Prior", 3))

ggplot(data=dens.datum, aes(x=Val, y=Dens, color=Thing))+
  geom_line()+
  geom_ribbon(data=dens.datum, 
              aes(ymax=Dens, fill=Thing), ymin=0, alpha=0.5)+
  scale_x_continuous(limit=c(-2,2), name="Coefficient value")+
  scale_y_continuous(name="Density")+
  geom_vline(xintercept=datum$b0.gam[3], linetype="dashed", 
             color = "#666666", size=1.2)+
  facet_grid(.~Type)+
  scale_color_brewer(palette="Dark2", 
                     labels=c("Posterior", 
                              "Prior"))+
  scale_fill_brewer(palette="Dark2", 
                     labels=c("Posterior", 
                              "Prior"))+
  theme(legend.position = "bottom",
        # legend.justification = c(0,1),
        legend.key = element_rect(fill= "white"),
        legend.direction = "horizontal",
        legend.text = element_text(size=12),
        # legend.key.size = unit(1.5, 'lines'),
        strip.background=element_rect(color="black", fill="#666666"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))+
  geom_text(data=dat.text, 
            mapping=aes(x=-Inf, y=Inf, label=labels),
            hjust=-0.2, 
            vjust=2, 
            size=4, 
            color="black")
