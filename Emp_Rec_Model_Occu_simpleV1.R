
#########################################################################################

###         MAKE FAKE DATA

#########################################################################################

# Load package for categorical distribution to create false-positive encounter history
library(LaplacesDemon)

data.fn <- function(nsites=20, nyears=4, nregions=3, ngroups=12, s.prob=0.8, Tr=600, 
                    imm.range=c(0,0.1), em.range=c(0.1, 0.2), b0.gam=1.5, sd.proc=4, 
                    sigma.group=1, psi1=0.4, p11.range=c(0.4,0.6), p10.range=c(0,0.1), 
                    b.range=c(0.2,0.4), patch.surv.range=c(0.7,0.9), patch.col.range=c(0,0.1),
                    noccs=3){
  
  # Function to generate data for IPM model to estimate recruitment using dynamic, false-positive occupancy
  # model and group counts with fixed survival. Annual varaition in parameters specified by uniform distr.
  
    # nsites = number of sites
    # nyears = number of years
    # nregions = number of regions
    # ngroups = number packs monitored 
    # s.prob = survival probability
    # Tr = territory size
    # imm.range = bounds of uniform distribution to draw from for immigration into group
    # em.range = bounds of uniform distribution to draw from for emigration out of group
    # b0.gam = intercept/ mean recruitment of pups
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
  y.occ <- array(dim=c(nsites, noccs, nyears, nregions))
  z <- muZ <- array(dim=c(nsites, nyears, nregions))
  p.0 <- p.1 <- array(dim=c(nyears, 3))
  prob <- array(dim=c(nsites, 3, nyears, nregions))
  
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
  # for(k in 1:nyears){
  #   for(r in 1:nregions){
  #     sigma.proc.group[k,r] <- rnorm(1, 0, sd.proc)
  #   }
  # }
  
  # Determine annual group immigration and emigration based on bounds
  
  imm.group <- runif(n=nyears-1, min=imm.range[1], max=imm.range[2])
  em.group <- runif(n=nyears-1, min=em.range[1], max=em.range[2])
  
  # Group size in successive years
  for(i in 1:ngroups){
    for(k in 2:nyears){
      for(r in 1:nregions){
        G[i,k,r] <- ifelse(G[i,k-1,r]==0, 0, rpois(1, (G[i,k-1,r]*s.prob*(1+imm.group[k-1]-em.group[k-1])+gamma))) #+sigma.proc.group[k-1,r])))
      }
    }
  }
  
  # Group counts
  for(i in 1:ngroups){
    for(k in 1:nyears){
      for(r in 1:nregions){
        y.temp [i,k,r] <- round(rnorm(1, G[i,k,r], sigma.group))
        y.temp[i,k,r] <- ifelse(y.temp[i,k,r] < 0, 0, y.temp[i,k,r])
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
  # Initial occupancy, psi
  for(i in 1:nsites){
    for(r in 1:nregions){
      psi[i,1,r] <- psi1
    }
  }
  
  # Determine detection and patch colinization and survivial
  p11 <- runif(n=nyears, min=p11.range[1], max=p11.range[2])
  p10 <- runif(n=nyears, min=p10.range[1], max=p10.range[2])
  b <- runif(n=nyears, min=b.range[1], max=b.range[2])
  patch.col <- runif(n=nyears-1, min=patch.col.range[1], max=patch.col.range[2])
  patch.surv <- runif(n=nyears-1, min=patch.surv.range[1], max=patch.surv.range[2])
  
  # Generate latent states of occurance
  # Year 1
  for(r in 1:nregions){
    z[,1,r] <- rbinom(nsites, 1, psi[1])
  }
  
  # Later years
  for(i in 1:nsites){
    for(r in 1:nregions){
      for(k in 2:nyears){
        muZ[k] <- z[i,k-1,r]*patch.surv[k-1]+(1-z[i,k-1,r])*patch.col[k-1]
        z[i,k,r] <- rbinom(1,1,muZ[k])
      }
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
    for(r in 1:nregions){
      for(k in 1:nyears){
        if(z[i,k,r]==0){
          prob[i,,k,r] <- p.0[k,]
        } else {
          prob[i,,k,r] <- p.1[k,]
        }
        for(j in 1:noccs){
          y.occ[i,j,k,r] <- rcat(1, prob[i,,k,r])
        }
      }
    }
  }
  
  # Compute annual occupancy so I can get area occupied
  for(i in 1:nsites){
    for(r in 1:nregions){
      for(k in 2:nyears){
       psi[i,k,r] <- psi[i,k-1,r]*patch.surv[k-1]+(1-psi[i,k-1,r])*patch.col[k-1] 
      }
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
      P[k,r] <- A[k,r]/Tr
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
      N.ad[k,r] <- rbinom(1, round(N.tot[k-1,r]*(1+patch.col[k-1]-(1-patch.surv[k-1]))), s.prob)
      N.tot[k,r] <- N.rec[k,r]+N.ad[k,r]
    }
  }
  
  # Alt pop size
  for(k in 1:nyears){
    for(r in 1:nregions){
      N.alt[k,r] <- P[k,r]*G.mean[k,r]
    }
  }
  
  
  
  
  
  return(list(nsites=nsites, nyears=nyears, nregions=nregions, ngroups=ngroups, noccs=noccs,
              s.prob=s.prob, Tr=Tr, imm.group=imm.group, em.group=em.group, patch.col=patch.col, 
              patch.surv=patch.surv, y.group=y.group, b0.gam=b0.gam, P=P, G.mean=G.mean, A=A,
              gamma.mean=gamma, N.tot=N.tot, psi=psi, area=area, var.group=var.group, 
              N.alt=N.alt, psi=psi, y.occ=y.occ, p10=p10, p11=p11, b=b, z=z))
}



datum<-data.fn(nsites=50, nyears=10, nregions=4, ngroups=40, s.prob=0.8, Tr=600, 
               imm.range=c(0,0.1), em.range=c(0.1, 0.2), b0.gam=1.5, sd.proc=4, 
               sigma.group=1, psi1=0.4, p11.range=c(0.4,0.6), p10.range=c(0,0.1), 
               b.range=c(0.2,0.4), patch.surv.range=c(0.7,0.9), patch.col.range=c(0,0.1),
               noccs=10)

#detach()
attach(datum)
str(datum)














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
      for(r in 1:nregions){
        # logit.psi1[i] <- B0.psi1       
        # logit(psi1[i]) <- logit.psi1[i]                                     
        # z[i,1] ~ dbern(psi1[i])
        z[i,1,r] ~ dbern(psi1)
    
        # for(k in 1:nyears-1){                                                    
        #   # logit.phi[i,k] <- B0.phi[k]                                       
        #   # logit.gamma[i,k] <- B0.gamma[k]
        #   # logit(phi[i,k]) <- logit.phi[i,k]                         
        #   # logit(gamma[i,k]) <- logit.gamma[i,k]
        # }#k
    
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
                # multilogit.p11[i,j,k] <- B0.p11[k]
                # logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
                # multilogit.p10[i,j,k] <- B0.p10[k]
                # logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
                # multilogit.b[i,j,k] <- B0.b[k]
                # logit(b[i,j,k]) <- multilogit.b[i,j,k]
          
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
          g.mu[i,k,r] <- G[i,k-1,r] * s[k-1,r] * (1 + imm.group[k-1] - em.group[k-1]) + gamma[i,k-1,r] # + sigma.proc.group[k-1,r])
          G[i,k,r] ~ dnorm(g.mu[i,k,r], 1/(g.mu[i,k,r]+.01))
        }
      }
    }
    
    # Observation proccess
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        for(r in 1:nregions){
          y.group[i,k,r] ~ dnorm(G[i,k,r], tauy.group)T(0,)
        }
      }
    }
    
    # Derived parameters
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        G.mean[k,r] <- mean(G[,k,r])
        gamma.mean[k,r] <- mean(gamma[,k,r])
        n.est[k,r] <- P[k,r] * G.mean[k,r]
      }
    }
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        for(r in 1:nregions){
          log(mu.gamma[i,k,r]) <- B0.gam
          gamma[i,k,r] ~ dpois(mu.gamma[i,k,r])
        }
      }
    }
    
    
    
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
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", "s", "phi", "colo", "psi", 
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
    }
    
    
    ## Bring in data s, G.mean, gamma.mean, P, colo, and phi
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        #s[k,r] ~ dnorm(s2[k,r,1], 1 / (s2[k,r,2] * s2[k,r,2])) 
        s[k,r] <- 0.8
        P[k,r] ~ dnorm(P2[k,r,1], 1 / (P2[k,r,2] * P2[k,r,2]))
        gamma.mean[k,r] ~ dnorm(gamma2[k,r,1], 1 / (gamma2[k,r,2] * gamma2[k,r,2]))
      }
    }

    for(k in 1:(nyears-1)){
      #Do I need to use binomial approximation? np=mean, np(1-p)=variance
      phi[k] ~ dnorm(phi2[k,1,1], 1 / (phi2[k,1,2] * phi2[k,1,2]))
      colo[k] ~ dnorm(colo2[k,1,1], 1 / (colo2[k,1,2] * colo2[k,1,2]))
    }
    
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
      }
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        n.est[k,r,1] ~ dnorm(N.tot[k,r], (1 / (n.est[k,r,2]*n.est[k,r,2]+.001)))
      }
    }
    
    
    
    
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







# DATA THINGS

library(dplyr)
library(tidyr)
library(ggplot2)

# chars<-c("mean$", "sd$")
# varss<-c("G.mean", "P")
# for(i in chars){
#   for(j in varss){
#     output[[paste(varss[j],sep="")]]<- paste(c("out$BUGSoutput$", chars[i], varss[j]), sep="")
#   }
# }

# Estimated results

GM<-as.data.frame(out$BUGSoutput$mean$G.mean)
colnames(GM)<-c(1:nregions)
rownames(GM)<-c(1:nyears)
GM<-GM %>% mutate(Year=c(1:nyears))

GM.sd<-as.data.frame(out$BUGSoutput$sd$G.mean) %>% mutate(Year=c(1:nyears))
colnames(GM.sd)<-c(1:nregions, "Year")
GM.sd<-gather(GM.sd, Region, G.mean.sd, 1:nregions)

Ps<-as.data.frame(out2$BUGSoutput$mean$P) %>% mutate(Year=c(1:nyears))
colnames(Ps)<-c(1:nregions, "Year")
Ps<-gather(Ps, Region, P, 1:nregions) 

Ps.sd<-as.data.frame(out2$BUGSoutput$sd$P) %>% mutate(Year=c(1:nyears))
colnames(Ps.sd)<-c(1:nregions, "Year")
Ps.sd<-gather(Ps.sd, Region, P.sd, 1:nregions)

N<-as.data.frame(out2$BUGSoutput$mean$N.tot) %>% mutate(Year=c(1:nyears))
colnames(N)<-c(1:nregions, "Year")
N<-gather(N, Region, N.tot, 1:nregions) 

N.sd<-as.data.frame(out2$BUGSoutput$sd$N.tot) %>% mutate(Year=c(1:nyears))
colnames(N.sd)<-c(1:nregions, "Year")
N.sd<-gather(N.sd, Region, N.tot.sd, 1:nregions)

gam<-as.data.frame(out$BUGSoutput$mean$gamma.mean) %>% mutate(Year=c(1:nyears))
colnames(gam)<-c(1:nregions, "Year")
gam<-gather(gam, Region, gamma.mean, 1:nregions) 

gam.sd<-as.data.frame(out$BUGSoutput$sd$gamma.mean) %>% mutate(Year=c(1:nyears))
colnames(gam.sd)<-c(1:nregions, "Year")
gam.sd<-gather(gam.sd, Region, gamma.mean.sd, 1:nregions)


output<-gather(GM, Region, G.mean, 1:nregions) 
output<- left_join(output, Ps, by=c("Year","Region"))
output<- left_join(output, N, by=c("Year","Region"))
output<- left_join(output, gam, by=c("Year","Region"))
output<- left_join(output, GM.sd, by=c("Year","Region"))
output<- left_join(output, Ps.sd, by=c("Year","Region"))
output<- left_join(output, N.sd, by=c("Year","Region"))
output<- left_join(output, gam.sd, by=c("Year","Region"))


# Truth

GM<-as.data.frame(datum$G.mean) %>% mutate(Year=c(1:nyears))
colnames(GM)<-c(1:nregions, "Year")
GM<-gather(GM, Region, G.mean.true, 1:nregions) 

Ps<-as.data.frame(datum$P) %>% mutate(Year=c(1:nyears))
colnames(Ps)<-c(1:nregions, "Year")
Ps<-gather(Ps, Region, P.true, 1:nregions) 

N<-as.data.frame(datum$N.tot) %>% mutate(Year=c(1:nyears))
colnames(N)<-c(1:nregions, "Year")
N<-gather(N, Region, N.tot.true, 1:nregions) 

# gam<-as.data.frame(datum$gamma.mean)
# colnames(gam)<-c(1:nregions)
# gam<-gather(gam, Region, gamma.mean.true, 1:nregions) %>% mutate(Site=c(1:(nyears*nregions)))


output<- left_join(output, Ps, by=c("Year", "Region"))
output<- left_join(output, N, by=c("Year", "Region"))
output<- output %>% mutate(gamma.mean.true=rep(datum$gamma.mean, nrow(output)))
output<- left_join(output, GM, by=c("Year", "Region"))


ggplot(data=output, aes(x=Year, y=G.mean, color="blue"))+
  geom_point(data=output, aes(x=Year, y=G.mean.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=G.mean, color="blue"))+
  facet_grid(Region~.)

ggplot(data=output, aes(x=Year, y=N.tot, color="blue"))+
  geom_point(data=output, aes(x=Year, y=N.tot.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=N.tot, color="blue"))+
  facet_grid(Region~.)

ggplot(data=output, aes())+
  geom_point(data=output, aes(x=Year, y=gamma.mean.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=gamma.mean, color="blue"))+
  facet_grid(Region~.)




# Now randomly get rid of some observations and run it again :)

y.group.1<-as.data.frame(y.group[,,1])
y.group.2<-as.data.frame(y.group[,,2])
y.group.3<-as.data.frame(y.group[,,3])
y.group.4<-as.data.frame(y.group[,,4])
y.group.5<-as.data.frame(y.group[,,5])
y.group.25.1<-sample_n(y.group.1, 25, replace=FALSE)
y.group.25.2<-sample_n(y.group.2, 25, replace=FALSE)
y.group.25.3<-sample_n(y.group.3, 25, replace=FALSE)
y.group.25.4<-sample_n(y.group.4, 25, replace=FALSE)
y.group.25.5<-sample_n(y.group.5, 25, replace=FALSE)
y.group.25<-array(c(y.group.25.1,y.group.25.2,y.group.25.3,y.group.25.4,y.group.25.5), dim=c(ngroups, nyears, nregions))
y.group.25[,,1]<-y.group.25.1
y.group.25[,,2]<-y.group.25.2
y.group.25[,,3]<-y.group.25.3
y.group.25[,,4]<-y.group.25.4
y.group.25[,,5]<-y.group.25.5


