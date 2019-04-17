####################################################################I

# A priori models using IPM to estimate recruitment of wolves in Montana

# For more details regarding the IPM see the rscript IPM_WolfRec_V2

# For running models see the rscript Analyses_WolfRec_IPM


# Allison C. Keever
# akeever1122@gmail.com
# github.com/akeever2
# Montana Cooperative Wildlife Research Unit
# 2018

# Occupancy code adapted from Sarah B Bassing
# sarah.bassing@gmail.com
# github.com/SarahBassing
# Montana Cooperative Wildlife Research Unit
# August 2016

######################################################################I


#### Set WD ####

# Set the working directory I want to save files to
setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/ModelResultFiles")


#### M1; FIX R ~ pack size  ####

sink("M1_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit <- sum(event.chi[])
    
    for(i in 1:nobs){
    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])
    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    # pois2[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(2,25)
    # pois1[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(y.group[i,k],25)
    # G[i,k] <- (pois1[i,k] * indicator[i,k]) + (pois2[i,k] * (1 - indicator[i,k]))
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### M2; FIX R ~ pack size + ran year ####

sink("M2_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:(nyears-1)){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]
    
    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)

    fit <- sum(event.chi[])
    
    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### M3; FIX R ~ pack size + ran region ####

sink("M3_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])

    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.reg[GroupReg[i]])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()

#### M4; FIX R ~ pack size + ran year + ran region ####

sink("M4_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    # Random effect for year
    for(k in 1:(nyears-1)){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])

    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### M5: FIX R ~ pack size + ran year + ran region + method ####
sink("M5_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    for(i in 1:3){
    b.method[i] ~ dnorm(0, 0.001)
    }
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])

    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }

    for(k in 2:nyears){
      pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.method[Method[k]])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### M6: FIX R ~ pack size + ran year + ran region + 4wd + 2wd ####
sink("M6_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    b.2wd ~ dnorm(0, 0.001)
    b.4wd ~ dnorm(0, 0.001)
    
    # Prior for missing data
    for(i in 1:ngroups){
    FourWD[i] ~ dnorm(0,1)
    TwoWD[i] ~ dnorm(0,1)
    }
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]
    
    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])
    
    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.2wd * TwoWD[i] + b.4wd * FourWD[i])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()

#### M7: FIX R ~ pack size + ran year + ran region + forcov ####
sink("M7_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    b.forest ~ dnorm(0, 0.001)
    
    # Prior for missing data
    for(i in 1:ngroups){
    Forest[i] ~ dnorm(0,1)
    }
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])

    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }

    for(k in 2:nyears){
      pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.forest * Forest[i])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### M8: FIX R ~ pack size + ran year + ran region +   ####

#### M9: R ~ pack size + ran year + ran region + dd ####

sink("M9_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    b.dd ~ dnorm(0, 0.001)
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)

    # Dispersal
    for(k in 1:nyears){
      em.group[1,k] ~ dbeta(51.21888, 394.1627)
      em.group[2,k] ~ dbeta(47.74287, 525.4008)
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit <- sum(event.chi[])

    for(i in 1:nobs){
    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group[Harv[k-1],k-1]) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group[Harv[k],k]))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est2[k] <- P[k] * G.mean[k]
    n.est[k] <- P[k] * G.dat[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }

    for(k in 1:nyears){
      G.dat[k] ~ dnorm(mu.G[k], 1 / (sd.G[k] * sd.G[k]))T(0,)
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.dd * LogN[k])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()

#### M10: R ~ pack size + ran year + ran region + harv ####

sink("M10_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    for(i in 1:2){
    b.harv[i] ~ dnorm(0, 0.001)
    }
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)

    # Dispersal
    for(k in 1:nyears){
    em.group[1,k] ~ dbeta(51.21888, 394.1627)
    em.group[2,k] ~ dbeta(47.74287, 525.4008)
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]

    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])

    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])


    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group[Harv[k-1],k-1]) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group[Harv[k-1],k-1]))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est2[k] <- P[k] * G.mean[k]
    n.est[k] <- P[k] * G.dat[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    for(k in 1:nyears){
      G.dat[k] ~ dnorm(mu.G[k], 1 / (sd.G[k] * sd.G[k]))T(0,)
    }    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.harv[Harv[k]])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()

#### M11: R ~ pack size + ran year + ran region + winter ####

sink("M11_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    b.winter ~ dnorm(0, 0.001)
    
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)

    # Dispersal
    for(k in 1:nyears){
    em.group[1,k] ~ dbeta(51.21888, 394.1627)
    em.group[2,k] ~ dbeta(47.74287, 525.4008)
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]
    
    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])
    
    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])
    
    
    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group[Harv[k-1],k-1]) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est2[k] <- P[k] * G.mean[k]
    n.est[k] <- P[k] * G.dat[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    for(k in 1:nyears){
      G.dat[k] ~ dnorm(mu.G[k], 1 / (sd.G[k] * sd.G[k]))T(0,)
    }    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.winter * Winter[k])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### M12: FIX R~pack size+ran year+ran region+winter+harv+forest+road ####

sink("M12_GroupRecIPM.txt")
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
    B0.colo[k] ~ dnorm(0,0.001)
    }#k
    
    B0.phi ~ dnorm(0,0.001)
    
    # Priors for detection probabilities
    B0.p11 ~ dnorm(0,0.001)
    B0.p10 ~ dnorm(0,0.001)
    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    b.pc1.psi ~ dnorm(0,0.001) 
    b.recPC.psi ~ dnorm(0,0.001) 
    b.pc1.colo ~ dnorm(0,0.001) 
    b.recPC.colo ~ dnorm(0,0.001)
    b.pc1.phi ~ dnorm(0,0.001)
    b.area.p11 ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv.p11 ~ dnorm(0,0.001) 
    b.map.p11 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001)
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.frrds.p10 ~ dnorm(0,0.001)
    b.acv.p10 ~ dnorm(0,0.001)
    
    
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    # Random effect for year
    for(k in 1:nyears){
    eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(0,0.001)
    }
    
    
    
    ## 1.4 Group priors
    
    # Initial group sizes
    for(i in 1:ngroups){
    G[i,1] ~ dpois(7)T(2,)
    }
    
    # Process error
    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    var.group <- pow(sigma.group, 2)
    
    
    
    ## 1.5 Recruitment priors
    
    # Priors for beta coefficients
    B0.gam ~ dnorm(0,0.001)
    B1.gam ~ dnorm(0,0.001)
    
    for(i in 1:2){
    b.harv[i] ~ dnorm(0, 0.001)
    }
    
    # Random effect for year
    for(k in 1:nyears){
    eps.gam[k] ~ dnorm(0, tau.gam)
    }
    
    sigma.gam ~ dunif(0,100)
    tau.gam <- pow(sigma.gam, -2)
    var.gam <- pow(sigma.gam, 2)
    
    # Random effect for region
    for(r in 1:nregions){
    eps.reg[r] ~ dnorm(0, tau.reg)
    }
    
    sigma.reg ~ dunif(0,100)
    tau.reg <- pow(sigma.reg, -2)
    var.reg <- pow(sigma.reg, 2)
    
    
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
    logit(psi1[i]) <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
    z[i,1] ~ dbern(psi1[i])
    
    for(k in 1:(nyears-1)){
    logit(phi[i,k]) <- B0.phi + b.pc1.phi * PC1[i] 
    logit(colo[i,k]) <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
    }#k
    
    for(k in 2:nyears){
    muZ[i,k] <- z[i,k-1] * phi[i,k-1] + (1-z[i,k-1]) * colo[i,k-1]
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
    p[1,i,j,k,1] <- (1 - p10[i,j,k])
    p[1,i,j,k,2] <- (1 - p11[i,j,k])
    p[2,i,j,k,1] <- p10[i,j,k]
    p[2,i,j,k,2] <- (1 - b[i,j,k]) * p11[i,j,k]
    p[3,i,j,k,1] <- 0
    p[3,i,j,k,2] <- b[i,j,k] * p11[i,j,k]
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
    logit(p11[i,j,k]) <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
    logit(p10[i,j,k]) <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
    logit(b[i,j,k]) <- B0.b
    
    y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
    }#k
    }#j
    }#i
    
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    psi[i,1] <- psi1[i]
    
    for (k in 2:nyears){                                          
    psi[i,k] <- psi[i,k-1] * phi[i,k-1] + (1 - psi[i,k-1]) * colo[i,k-1]
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
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
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
    cloglog(mu.surv[i]) <- b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
    hazard[p,k] <- -log(1 - mu.pred[p,k])
    }#p
    }#k
    
    
    # Cumulative hazard and survival 
    
    for(k in 1:nyears){
    base.H[1,k] <- hazard[1,k] * width.interval[1]
    
    for(p in 2:nperiods){
    base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    }#p
    }#k
    
    for(k in 1:nyears){
    for(p in 1:nperiods){
    base.s[p,k] <- exp(-base.H[p,k])
    }#p
    
    annual.s[k] <- base.s[length(width.interval), k]
    }#k
    
    
    # Compute posterior predictive check statistics
    for(i in 1:nobs){
    # Expected values #### NOTE-Maybe multiple by nobs instead of 1? ####
    event.expected[i] <- 1 * mu.surv[i]
    
    # Fit statistic for actual data
    event.chi[i] <- pow((event[i] - event.expected[i],2))/(event.expected[i] + 0.01)
    fit <- sum(event.chi[])
    
    # Replicate data for GOF
    event.rep[i] ~ dbern(mu.surv[i])
    
    # Fit statistics for replicate data
    event.chi.new[i] <- pow((event.rep[i] - event.expected[i]),2)/(event.expected[i] + 0.01)
    }
    fit.new <- sum(event.chi.new[])
    
    
    
    
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
    g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group) + gamma[i,k-1]
    G[i,k] ~ dnorm(g.mu[i,k], 1 / (g.mu[i,k] + 0.00001))T(0,25)
    }
    }
    
    # Observation proccess
    for(i in 1:ngroups){
    for(k in 1:nyears){
    y.group[i,k] ~ dnorm(G[i,k], tauy.group)
    }
    }
    
    # Derived parameters
    for(k in 1:nyears){
    G.mean[k] <- mean(G[,k] * annual.s[k] * (1 - em.group))
    G.mean.high[k] <- mean(G[,k])
    gamma.mean[k] <- mean(gamma[,k])
    n.est[k] <- P[k] * G.mean[k]
    }
    
    for(k in 2:nyears){
    pop.growth[k] <- n.est[k] / n.est[k-1] 
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    for(i in 1:ngroups){
    for(k in 1:nyears){
    mu.gamma[i,k] <- exp(B0.gam + B1.gam * G[i,k] + eps.gam[k] + eps.reg[GroupReg[i]] + b.harv[Harv[i]])
    gamma[i,k] ~ dpois(mu.gamma[i,k])
    }
    }
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()



#### Population level - same for all models ####

sink("PopRecIPM.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## Population priors
    
    # Initial population size
    
    N.tot[1] ~ dnorm(600, 0.0001)I(0,)
    
    
    ## Bring in data s, G.mean, gamma.mean, P, colo, and phi
    
    for(k in 1:nyears){
    P[k] ~ dnorm(P2[k,1], 1 / (P2[k,2] * P2[k,2]+ 0.0000001))
    # G.mean[k] ~ dnorm(G.mean2[k,1], 1 / (G.mean2[k,2] * G.mean2[k,2]+ 0.0000001))
    gamma.mean[k] ~ dnorm(gamma2[k,1], 1 / (gamma2[k,2] * gamma2[k,2]+ 0.0000001))
    }
    
    for(k in 1:(nyears-1)){
    # B0.colo[k] ~ dnorm(betas[k,3], 1 / (betas[k,4] * betas[k,4]+ 0.0000001))
    eps.surv[k] ~ dnorm(betas[k,13], 1 / (betas[k,14] * betas[k,14]+ 0.0000001))
    }
    
    # B0.phi ~ dnorm(betas[1,1], 1 / (betas[1,2] * betas[1,2]+ 0.0000001))
    # b.pc1.colo ~ dnorm(betas[1,5], 1 / (betas[1,6] * betas[1,6]+ 0.0000001))
    # b.recPC.colo ~ dnorm(betas[1,7], 1 / (betas[1,8] * betas[1,8]+ 0.0000001))
    # b.pc1.phi ~ dnorm(betas[1,9], 1 / (betas[1,10] * betas[1,10]+ 0.0000001))
    
    for(p in 1:nperiods){
    b.period.surv[p] ~ dnorm(betas[p,11], 1 / (betas[p,12] * betas[p,12]+ 0.0000001))
    }
    
    T ~ dlnorm(6.22985815, 1/0.58728123)
    
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    
    # Ecological model/ system process
    
    # First determine colonization and extinction 
    
    # for(i in 1:nsites){
    # for(k in 1:(nyears-1)){
    # colo[i,k] <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1]
    # phi[i,k] <- B0.phi + b.pc1.phi * PC1[i]
    # }
    # }
    
    # Then determine survival
    # Baseline hazard
    
    for(k in 1:(nyears-1)){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
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
    
    
    for(k in 2:nyears){
    N.rec[k] ~ dpois(P[k-1] * gamma.mean[k-1])
    # N.ps[k] ~ dpois((P[k-1] + sum(colo[,k-1]*((area[] * T.overlap[k])/T - P[k-1])) - P[k-1] * (1 - sum(phi[,k-1]))) * G.mean[k-1])
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







