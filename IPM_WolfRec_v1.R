##############################################################################

#             Integrated population model to estimate recruitment

# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2017

##############################################################################

# Pull in data and call for appropriate packages


library(R2jags)
library(mcmcplots)


################################################################################
#  Specify model in BUGS language

sink("IPMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Occupancy priors
    
      #  psi1 coefficient (occupancy in year 1)
      
      b0.psi1 ~ dnorm(0,0.001)	  
      
      
      #  Priors for transition probabilities (survival and colonization)
      
      for(k in 1:(nyears-1)){
        b0.phi[k] ~ dnorm(0,0.001)
        b0.colo[k] ~ dnorm(0,0.001)
      }#k
      
      
      # Priors for detection probabilities (only varies by year)
      
      for(k in 1:nyears){
        b0.p11[k] ~ dnorm(0,0.001)
        b0.p10[k] ~ dnorm(0,0.001)
        b0.b[k] ~ dnorm(0,0.001)
      }#k
      
      
      
    ## 1.2 Territory priors
   

 
    ## 1.3 Survival priors

      # Intercept and covariates
      
      b0.surv ~ dnorm(0,0.001)

      
      # Prior for time periods
      
      for(i in 1:nperiods){
        b.period[i] ~ dnorm(0, 0.001)
      }
      

      # Random effect for year with year 1 fixed at 0

      b.year[1] <- 0
      for(k in 2:nyears){
        b.year[k] ~ dnorm(0, 0.001)
      }
 

    
    ## 1.4 Group priors
    
      # Initial group sizes
      
      for(r in 1:nregions){
        for(i in 1:ngroups){
          G[i,1,r] ~ dpois(5)
        }
      }
      
      
      # Process and observation error
      
      # for(k in 1:nyears){
      #   for(r in 1:nregions){
      #     sigma.proc.group[k,r] ~ dnorm(0, tau.proc)
      #   }
      # }
      # 
      # tau.proc <- 1/(sd.proc * sd.proc)
      # sd.proc ~ dunif(0, 20)
      # var.proc <- sd.proc*sd.proc
      
      tauy.group <- pow(sigma.group, -2)
      sigma.group ~ dunif(0,100)
      var.group <- pow(sigma.group, 2)
      
      
      
    ## 1.5 Recruitment priors
    
      # Priors for beta coefficients
    
      b0.gam ~ dunif(-10,10)
    


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
        logit.psi1[i] <- b0.psi1       
        logit(psi1[i]) <- logit.psi1[i]                                     
        z[i,1] ~ dbern(psi1[i])
    
        for(k in 1:nyears-1){                                                    
          logit.phi[i,k] <- b0.phi[k]                                       
          logit.gamma[i,k] <- b0.gamma[k]
          logit(phi[i,k]) <- logit.phi[i,k]                         
          logit(gamma[i,k]) <- logit.gamma[i,k]
        }#k
    
        for(k in 2:nyears){
          muZ[i,k,r] <- z[i,k-1,r]*phi[k-1] + (1-z[i,k-1,r])*colo[k-1]
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
          for(r in 1:nregions){
            multilogit.p11[i,j,k] <- b0.p11[k]
            logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
            multilogit.p10[i,j,k] <- b0.p10[k]
            logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
            multilogit.b[i,j,k] <- b0.b[k]
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
      }#r
    }#k
    
    
    
    #####################
    
    # 2.2. Territory model 
    
    # Input includes area occupied (A) indexed by year (k) and region (r) from occupancy model (2.1.)
    # Area occupied is divided by territory size (T)
    #   Territory size is fixed value
    
    # Output is number of packs (P) indexed by year (k) and region (r)
    
    ####################
    
    # Pull in data for the mean for territory size
    
    T ~ dnorm(mu.T, tau.t)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        P[k,r] <- A[k,r] / T
      }#r
    }#k
    
    
    
    #####################
    
    # 2.3. Survival likelihood 
    
    # Current model is simple, based on informative prior and is a simple binomial model. 
    # Output is survival indexed by year (k) and region (r)
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(i in 1:nobs){
        event[i] ~ dbern(mu.surv[i])
        cloglog(mu.surv[i]) <- b0.surv + b.year[Year[i]] + b.period[Period[i]]
    }#i
    

    # Predicted values

    # Baseline hazard    
    for(k in 1:nyears){
      for(p in 1:nperiods){
        cloglog(mu.pred[p,k]) <- b0.surv +b.year[k] + b.period[p]
        hazard[p,k] <- -log(1-mu.pred[p,k])
      }#p
    }#k
    
    # Cumulative hazard and survival 

    for(k in 1:nyears){
      base.H[1,k] <- hazard[1,k] * width.interval[1]
      for(p in 2:nperiods){
        base.H[p,k] <- H[p-1, k] + hazard[p,k] * width.interval[p]
      }#p
    }#k
    
    for(k in 1:nyears){
      for(p in 1:nperiods){
        base.s[p,k] <- exp(-base.H[p,k])
        annual.s[k] <- base.s[length(p), k]
      }#p
    }#k
    
    
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
        for(r in 1:nregions){
          g.mu[i,k,r] <- G[i,k-1,r] * annual.s[k-1,r] * (1 - em.group[k-1]) + gamma[i,k-1,r] # + sigma.proc.group[k-1,r])
          G[i,k,r] ~ dnorm(g.mu[i,k,r], 1/(g.mu[i,k,r]+.01))
        }#r
      }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        for(r in 1:nregions){
          y.group[i,k,r] ~ dnorm(G[i,k,r], tauy.group)T(0,)
        }#r
      }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        G.mean[k,r] <- mean(G[,k,r])
        gamma.mean[k,r] <- mean(gamma[,k,r])
        n.est[k,r] <- P[k,r] * G.mean[k,r]
      }#r
    }#k
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        for(r in 1:nregions){
          log(mu.gamma[i,k,r]) <- B0.gam
          gamma[i,k,r] ~ dpois(mu.gamma[i,k,r])
        }#r
      }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "nregions"=nregions, "area"=area,
                 "ngroups"=ngroups, "noccs"=noccs, "T"=Tr, "s.prob"=s.prob, "imm.group"=imm.group, 
                 "em.group"=em.group, "y.group"=y.group, "y.occ"=y.occ)


#  Initial Values	
zst <- apply(y.occ,c(1,3,4), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 1
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(B0.gam=runif(1,-1,1), sd.proc=runif(1,0,10), 
                         sigma.group=runif(1,0,10), z=zst, colo=runif(nyears-1, 0.05, 0.15), 
                         phi=runif(nyears-1, 0.7,0.8), b=runif(nyears, 0.1, 0.3), 
                         p10=runif(nyears, 0.05, 0.15), p11=runif(nyears, 0.4, 0.6))}


# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", "annual.s", "phi", "colo", "psi", 
            "p11", "p10", "b") 


# MCMC Settings 
ni <- 10000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)




n.est <- array(NA, dim=c(nyears, nregions, 2))
n.est[,,1] <- out$BUGSoutput$mean$n.est
n.est[,,2] <- out$BUGSoutput$sd$n.est
s2 <- array(NA, dim=c(nyears, nregions, 2))
s2[,,1] <- out$BUGSoutput$mean$s
s2[,,2] <- out$BUGSoutput$sd$s
P2 <- array(NA, dim=c(nyears, nregions, 2))
P2[,,1] <- out$BUGSoutput$mean$P
P2[,,2] <- out$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears, nregions, 2))
gamma2[,,1] <- out$BUGSoutput$mean$gamma.mean
gamma2[,,2] <- out$BUGSoutput$sd$gamma.mean
phi2 <- array(NA, dim=c(nyears-1,1, 2))
phi2[,,1] <- out$BUGSoutput$mean$phi
phi2[,,2] <- out$BUGSoutput$sd$phi
colo2 <- array(NA, dim=c(nyears-1,1, 2))
colo2[,,1] <- out$BUGSoutput$mean$colo
colo2[,,2] <- out$BUGSoutput$sd$colo
b2 <- array(NA, dim=c(nyears,1, 2))
b2[,,1] <- out$BUGSoutput$mean$b
b2[,,2] <- out$BUGSoutput$sd$b
p102 <- array(NA, dim=c(nyears,1, 2))
p102[,,1] <- out$BUGSoutput$mean$p10
p102[,,2] <- out$BUGSoutput$sd$p10
p112 <- array(NA, dim=c(nyears,1, 2))
p112[,,1] <- out$BUGSoutput$mean$p11
p112[,,2] <- out$BUGSoutput$sd$p11
psi2 <- array(NA, dim=c(nsites, nyears, 2))
psi2[,,1] <- out$BUGSoutput$mean$psi[,,1]
psi2[,,2] <- out$BUGSoutput$sd$psi[,,1]



sink("PopLevel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## Population priors
    
    # Initial population size
    
    for(r in 1:nregions){
      N.tot[1,r] ~ dnorm(800, 0.0001)I(0,)
    }#r
    
    
    ## Bring in data s, G.mean, gamma.mean, P, colo, and phi
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        s[k,r] ~ dnorm(s2[k,r,1], 1 / (s2[k,r,2] * s2[k,r,2])) 
        P[k,r] ~ dnorm(P2[k,r,1], 1 / (P2[k,r,2] * P2[k,r,2]))
        gamma.mean[k,r] ~ dnorm(gamma2[k,r,1], 1 / (gamma2[k,r,2] * gamma2[k,r,2]))
      }#r
    }#k
    
    for(k in 1:(nyears-1)){
      #Do I need to use binomial approximation? np=mean, np(1-p)=variance
      phi[k] ~ dnorm(phi2[k,1,1], 1 / (phi2[k,1,2] * phi2[k,1,2]))
      colo[k] ~ dnorm(colo2[k,1,1], 1 / (colo2[k,1,2] * colo2[k,1,2]))
    }#k
    
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    
    # Ecological model/ system process
    
    for(k in 2:nyears){
      for(r in 1:nregions){
        N.rec[k,r] ~ dpois(P[k-1,r] * gamma.mean[k-1,r])
        N.ad[k,r] ~ dbin(s[k-1,r], round(N.tot[k-1,r] * (1+colo[k-1]-(1-phi[k-1]))))
        N.tot[k,r] <- N.ad[k,r] + N.rec[k,r]
        #n.mu[k,r] <- N.tot[k-1,r]*(1+colo[k]-(1-phi[k]))*s[k-1,r]+P[k-1,r]*gamma.mean[k-1,r]
        #N.tot[k,r] ~ dpois(n.mu[k,r])
      }#r
    }#k
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        n.est[k,r,1] ~ dnorm(N.tot[k,r], (1 / (n.est[k,r,2]*n.est[k,r,2]+.001)))
      }#r
    }#k
    
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data2 <- list("nyears"=nyears, "nregions"=nregions, "phi2"=phi2,
                  "n.est"=n.est, "s2"=s2, "gamma2"=gamma2, "P2"=P2, 
                  "colo2"=colo2)


#  Initial Values	
inits2 <- function(){list(colo=runif(nyears-1, 0.05, 0.15), 
                          phi=runif(nyears-1, 0.7,0.8))}


# Parameters to keep track of and report
params2 <- c("P", "s", "N.tot", "gamma.mean") 


# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

print(out2, dig=2)
mcmcplot(out2)



