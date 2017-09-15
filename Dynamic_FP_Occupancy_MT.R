sink("IPMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Occupancy priors
    
    #  psi1 coefficient (occupancy in year 1)
    
    #B0.psi1 ~ dnorm(0,0.001)	  
    psi1 ~ dunif(0,1)
    
    
    #  Priors for transition probabilities (survival and colonization)
    
    for(k in 1:(nyears-1)){
    #B0.phi[k] ~ dnorm(0,0.001)
    #B0.gamma[k] ~ dnorm(0,0.001)
    phi[k] ~ dunif(0,1)
    colo[k] ~ dunif(0,1)
    }#k
    
    
    # Priors for detection probabilities (only varies by year)
    
    for(k in 1:nyears){
    #B0.p11[k] ~ dnorm(0,0.001)
    #B0.p10[k] ~ dnorm(0,0.001)
    #B0.b[k] ~ dnorm(0,0.001)
    p11[k] ~ dunif(0,1)
    p10[k] ~ dunif(0,1)
    b[k] ~ dunif(0,1)
    }#k
    
    
    
    
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
    for(r in 1:nregions){
     logit.psi1[i] <- B0.psi1 + b.pc1.psi * PC1[i] + b.pcks5yr.psis + pcks5yrs[i]       
     logit(psi1[i]) <- logit.psi1[i]                                     
     z[i,1] ~ dbern(psi1[i])
    z[i,1,r] ~ dbern(psi1)
    
    for(k in 1:nyears-1){
      logit.phi[i,k] <- B0.phi[k]+ b.year[k] * Year[i,k]
      logit.gamma[i,k] <- B0.gamma[k] + b.year[k] * Year[i,k] + b.pc1.gam * PC1[i] + b.pcks5yrs.gam + pcks5yrs[i,k] 
      logit(phi[i,k]) <- logit.phi[i,k]
      logit(gamma[i,k]) <- logit.gamma[i,k]
    }#k
    
    for(k in 2:nyears){
    muZ[i,k,r] <- z[i,k-1,r]*phi[i,k-1] + (1-z[i,k-1,r])*gamma[i,k-1]
    z[i,k,r] ~ dbern(muZ[i,k,r])
    }#k
    }#r
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
    
    #  Observation model
    
    for(i in 1:nsites){
    for(j in 1:noccs){
    for(k in 1:nyears){
    for(r in 1:nregions){
    multilogit.p11[i,j,k] <- B0.p11[k] + b.area * area[i,j,k] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i,j,k] + b.frrds.p11 * forrds[i,j,k]
    logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
    multilogit.p10[i,j,k] <- B0.p10[k] + b.acv * acv[i,j,k] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i,j,k] + b.frrds.p11 * forrds[i,j,k]
    logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
    multilogit.b[i,j,k] <- B0.b
    logit(b[i,j,k]) <- multilogit.b[i,j,k]
    
    y.occ[i,j,k,r] ~ dcat(p[,i,j,k,z[i,k,r]+1])
    }#r
    }#k
    }#j
    }#i
    
    
    # Derived parameters
    
    for(i in 1:nsites){
    for(r in 1:nregions){
    psi[i,1,r] <- psi1
    growthr[i,1,r] <- 1  
    
    for (k in 2:nyears){                                          
    psi[i,k,r] <- psi[i,k-1,r]*phi[k-1] + (1-psi[i,k-1,r])*colo[k-1]
    growthr[i,k,r] <- psi[i,k,r]/psi[i,k-1,r]
    turnover[i,k-1,r] <- (1 - psi[i,k-1,r]) * colo[k-1]/psi[i,k,r]
    }#k
    }#r
    }#i
    
    #  Area occpupied indexed by year and region
    
    for(k in 1:nyears){
    for(r in 1:nregions){
    A[k,r] <- sum(psi[,k,r] * area[,k,r])
    }
    }
    
    
    
    
    
    #####################
    
    # 2.2. Territory model 
    
    # Input includes area occupied (A) indexed by year (k) and region (r) from occupancy model (2.1.)
    # Area occupied is divided by territory size (T)
    #   Territory size is fixed value
    
    # Output is number of packs (P) indexed by year (k) and region (r)
    
    ####################
    
    # Pull in data for the mean for territory size
    
    #T ~ dnorm(mu.T, tau.t)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    
    for(k in 1:nyears){
    for(r in 1:nregions){
    P[k,r] <- area[k,r] / T
    }
    }
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "nregions"=nregions, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


#  Initial Values	
zst <- apply(y.occ,c(1,3,4), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 1
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(z=zst, colo=runif(nyears-1, 0.05, 0.15), phi=runif(nyears-1, 0.7,0.8),
                         b=runif(nyears, 0.1, 0.3), p10=runif(nyears, 0.05, 0.15),
                         p11=runif(nyears, 0.4, 0.6))}


# Parameters to keep track of and report
params <- c("P", "phi", "colo", "psi", 
            "p11", "p10", "b") 


# MCMC Settings 
ni <- 8000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)

