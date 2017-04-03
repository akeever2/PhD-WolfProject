##############################################################################

#             Integrated population model code to estiamte recruitment

# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2016

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
    
    ## 1.2 Territory priors
    
    ## 1.3 Survival priors
    
    ## 1.4 Group priors
    
    # Initial group sizes
    
    for(r in 1:nregions){
      for(i in 1:ngroups){
        G[i,1,r] ~ dpois(5)
      }
    }
    
    # Process and observation error
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        sigma.proc.group[k,r] ~ dnorm(0, tau.proc)
      }
    }
    
    tau.proc <- 1/(sd.proc*sd.proc)
    sd.proc ~ dunif(0, 20)
    var.proc <- sd.proc*sd.proc
    
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
    
    ####################
    
    # Simplified model to get area occupied based on fixed values of psi (occupancy probability)
    
    
    #  Area occpupied indexed by year and region
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        A[k,r] <- sum(psi[,k,r]*area[,k,r])
      }
    }
    
    
    
    #####################
    
    # 2.2. Territory model 
    
    # Input includes area occupied (A) indexed by year (k) and region (r) from occupancy model (2.1.)
    # Area occupied is divided by territory size (T)
    #   Simple model of the mean based on results from Rich et al. 2012
    
    # Output is number of packs (P) indexed by year (k) and region (r)
    
    ####################
    
    # Pull in data for the mean for territory size
    
    #T ~ dnorm(mu.T, tau.t)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        P[k,r] <- A[k,r] / T
      }
    }
    
    
    
    
    #####################
    
    # 2.3. Survival likelihood 
    
    # Current model is simple, based on informative prior and is a simple binomial model. 
    # Output is survival indexed by year (k) and region (r)
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        s[k,r] <- s.prob
      }
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
        for(r in 1:nregions){
          g.mu[i,k,r] <- G[i,k-1,r] * s[k-1,r]*(1 + imm.group - em.group) + gamma[i,k-1,r] + sigma.proc.group[k-1,r]
          G[i,k,r] ~ dpois(g.mu[i,k,r])
        }
      }
    }
    
    # Observation proccess
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        for(r in 1:nregions){
          y.group[i,k,r] ~ dnorm(G[i,k,r], tauy.group)
        }
      }
    }
    
    # Derived parameters
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        G.mean[k,r] <- mean(G[,k,r])
        gamma.mean[k,r] <- mean(gamma[,k,r])
      }
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(i in 1:ngroups){
      for(k in 2:nyears){
        for(r in 1:nregions){
          log(mu.gamma[i,k,r]) <- b0.gam
          gamma[i,k,r] ~ dpois(mu.gamma[i,k,r])
        }
      }
    }
    
    
    
    
    
    
    
    ############################################################
    
    #             3. Derived parameters
    
    ############################################################
    
    
    
    
    ############################################################
    
    #             4. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "nregions"=nregions, "area"=area,
                 "ngroups"=ngroups, "psi"=psi, "T"=T, "s.prob"=s.prob, "imm.group"=imm.group, 
                 "em.group"=em.group, "y.group"=y.group)


#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("P", "s", "var.proc", "var.group", "G.mean", "gamma.mean") 


# MCMC Settings 
ni <- 500
nt <- 2
nb <- 250
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)




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
    }



    #####################
    
    # Population level model 

    ####################

    # Ecological model/ system process

    for(k in 2:nyears){
      for(r in 1:nregions){
        N.rec[k,r] ~ dpois(P[k-1,r]*gamma.mean[k-1,r])
        N.ad[k,r] ~ dbin(s[k-1,r], round(N.tot[k-1,r]*(1+imm.p-em.p)))
        N.tot[k,r] <- N.ad[k,r]+N.rec[k,r]
        #n.mu[k,r] <- N.tot[k-1,r]*(1+imm.p-em.p)*s[k-1,r]+P[k-1,r]*gamma.mean[k-1,r]
        #N.tot[k,r] ~ dpois(n.mu[k,r])
      }
    }

    # Linking pack size (P) and mean group size (G.mean) to abundance (N.tot)

    for(k in 1:nyears){
      for(r in 1:nregions){
        n.est[k,r] <- P[k,r]*G.mean[k,r]
        n.est[k,r] ~ dpois(N.tot[k,r])
      }
    }


    ############################################################
    
    #             3. Derived parameters
    
    ############################################################
    
    
    
    ############################################################
    
    #             4. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "nregions"=nregions, "area"=area,
                 "ngroups"=ngroups, "psi"=psi, "T"=T, "s.prob"=s.prob, "imm.group"=imm.group, 
                 "em.group"=em.group, "y.group"=y.group)


#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("P", "s", "var.proc", "var.group", "G.mean", "gamma.mean") 


# MCMC Settings 
ni <- 500
nt <- 2
nb <- 250
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)














#########################################################################################

###         MAKE FAKE DATA

#########################################################################################



data.fn <- function(nsites=20, nyears=4, nregions=3, ngroups=12, s.prob=0.8, T=600, 
                    imm.group=0.1, em.group=0.15, imm.p=0, em.p=0, b0.gam=1.5, sd.proc=4, 
                    sigma.group=1){
  
  # Arrays
  y.group <- array(dim=c(ngroups, nyears, nregions))
  G <- array(dim=c(ngroups, nyears, nregions))
  sigma.proc.group <- array(dim=c(nyears, nregions))
  y.temp <- array(dim=c(ngroups, nyears, nregions))
  G.mean <- array(dim=c(nyears, nregions))
  psi <- array(dim=c(nsites, nyears, nregions))
  area <- array(dim=c(nsites, nyears, nregions))
  A <- array(dim=c(nyears, nregions))
  P <- array(dim=c(nyears, nregions))
  N.tot <- array(dim=c(nyears, nregions))
  N.rec <- array(dim=c(nyears, nregions))
  N.ad <- array(dim=c(nyears, nregions))
  N.alt <- array(dim=c(nyears, nregions))
  
  # Recruitment and group size
  # Recruitment
  gamma <- exp(b0.gam)
  
  # Random noise for first year group size
  for(i in 1:ngroups){
    for(r in 1:nregions){
      G[i,1,r] <- max(2, rpois(1, 4))
    }
  }
  
  # Process noise for group counts
  for(k in 1:nyears){
    for(r in 1:nregions){
      sigma.proc.group[k,r] <- rnorm(1, 0, sd.proc)
    }
  }
  
  # Group size in successive years
  for(i in 1:ngroups){
    for(k in 2:nyears){
      for(r in 1:nregions){
        G[i,k,r] <- ifelse(G[i,k-1,r]==0, 0, rpois(1, G[i,k-1,r]*s.prob*(1+imm.group-em.group)+gamma+sigma.proc.group[k-1,r]))
      }
    }
  }
  
  # Group counts
  for(i in 1:ngroups){
    for(k in 1:nyears){
      for(r in 1:nregions){
        y.temp [i,k,r] <- round(rnorm(1, G[i,k,r], sigma.group))
        if(y.temp[i,k,r]<0)y.temp[i,k,r]=0
        y.group[i,k,r] <- ifelse(y.temp[i,k,r] > G[i,k,r], G[i,k,r], y.temp[i,k,r])
      }
    }
  }
  
  # Observation error
  var.group <- sigma.group*sigma.group
  
  # Mean group size
  for(k in 1:nyears){
    for(r in 1:nregions){
      G.mean[k,r] <- mean(G[,k,r])
    }
  }
  
  
  # Occupancy and territory model to determine number of packs
  # Psi
  for(i in 1:nsites){
    for(r in 1:nregions){
      psi[i,,r] <- runif(1, 0.35, 0.85)
    }
  }
  
  # Area
  for(i in 1:nsites){
    for(r in 1:nregions){
      area[i,,r] <- 600
    }
  }
  
  # Area occupied
  for(k in 1:nyears){
    for(r in 1:nregions){
      A[k,r] <- sum(psi[,k,r]*area[,k,r])
    }
  }
  
  # Number of packs
  for(k in 1:nyears){
    for(r in 1:nregions){
      P[k,r] <- A[k,r]/T
    }
  }
  
  
  # Population level
  # Total population in the first year
  for(r in 1:nregions){
    N.tot[1,r] <- P[1,r]*G.mean[1,r]
  }
  
  # Population size in successive years
  for(k in 2:nyears){
    for(r in 1:nregions){
      N.rec[k,r] <- rpois(1, P[k-1,r]*gamma)
      N.ad[k,r] <- rbinom(1, round(N.tot[k-1,r]*(1+imm.p-em.p)), s.prob)
      N.tot[k,r] <- N.rec[k,r]+N.ad[k,r]
    }
  }
  
  # Alt pop size
  for(k in 1:nyears){
    for(r in 1:nregions){
      N.alt[k,r] <- P[k,r]*G.mean[k,r]
    }
  }
  
  
  
  
  
  return(list(nsites=nsites, nyears=nyears, nregions=nregions, ngroups=ngroups, 
              s.prob=s.prob, T=T, imm.group=imm.group, em.group=em.group, imm.p=imm.p, 
              em.p=em.p, y.group=y.group, b0.gam=b0.gam, P=P, G.mean=G.mean, 
              gamma.mean=gamma, N.tot=N.tot, psi=psi, area=area, var.group=var.group, 
              N.alt=N.alt))
}



datum<-data.fn(nsites=20, nyears=4, nregions=3, ngroups=12, s.prob=0.8, T=600, 
               imm.group=0.1, em.group=0.15, imm.p=0, em.p=0, b0.gam=1.5, sd.proc=4, 
               sigma.group=1)
attach(datum)
str(datum)



