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

#  Priors for OCCUPANCY parameters that are indexed by year (k)
  
    #  psi1 coefficient (occupancy in year 1)
    
    B0.psi1 ~ dnorm(0,0.001)	        
    
    #  Priors for transition probabilities (survival and colonization)
    
    for(k in 1:K-1){
      B0.phi[k] ~ dnorm(0,0.001)
      B0.gamma[k] ~ dnorm(0,0.001)
    }#k
    
    
    # Priors for detection probabilities (only varies by year)
    
    for(k in 1:K){
      B0.p11[k] ~ dnorm(0,0.001)
      B0.p10[k] ~ dnorm(0,0.001)
      B0.b[k] ~ dnorm(0,0.001)
    }#k
    
    
    #  Priors for covariates




#  Priors for Group level counts

  #  Prior for initial group size

    for(i in 1:ngroups){
      G[1,,i] ~ 
    }


  # Prior for observation error on group counts

    tauy.group <- pow(sigma.group, -2)
    sigma.group ~ dunif(0,100)
    sigma2.group <- pow(sigma.group, 2)  # keep track of this





#  Priors for recruitment model

    #  Priors for covariates
    
      b0.gam ~ dnorm(0,0.001)
    



############################################################

#             2. Likelihoods

############################################################


#####################

# 2.1. Occupancy likelihood 

# Adapted from Sarah B Bassing
# Montana Cooperative Wildlife Research Unit
# August 2016. 
# This is a DYNAMIC FALSE POSITVE MULTI-SEASON occupancy model 
# Encounter histories include:
#     1 = no detection 
#     2 = uncertain detection
#     3 = certain detection

####################

# Ecological process/submodel
# Define State (z) conditional on parameters- Nmbr sites occupied
    
    for(i in 1:nSite){
      logit.psi1[i] <- B0.psi1       
      logit(psi1[i]) <- logit.psi1[i]                                     
      z[i,1] ~ dbern(psi1[i])                                             
    
      for(k in 1:K-1){                                                    
        logit.phi[i,k] <- B0.phi[k]                                       
        logit.gamma[i,k] <- B0.gamma[k]
        logit(phi[i,k]) <- logit.phi[i,k]                         
        logit(gamma[i,k]) <- logit.gamma[i,k]
      }#k
    
      for(k in 2:K){
        muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*gamma[i,k-1]
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
    
    for(i in 1:nSite){
      for (j in 1:nOcc){
        for(k in 1:K){			                                  
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

    #  Observation model
    
    for(i in 1:nSite){
      for(j in 1:nOcc){
        for(k in 1:K){
          multilogit.p11[i,j,k] <- B0.p11[k]+B7.p11.SVY*SVY[i,j,k]
          logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
          multilogit.p10[i,j,k] <- B0.p10[k]+B7.p10.SVY*SVY[i,j,k]
          logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
          multilogit.b[i,j,k] <- B0.b[k]+B7.b.SVY*SVY[i,j,k]
          logit(b[i,j,k]) <- multilogit.b[i,j,k]
    
          y[i,j,k] ~ dcat(p[,i,j,k,z[i,k]+1])

        }#k
      }#j
    }#i


#  Derived parameters

  #  Mean annual parameter for psi and number of occupied sites           
  
    for(k in 1:K){
      psik[k] <- mean(psi[,k])
      n.occ[k] <- sum(z[1:nSite, k])
    }#k

  #  Area occpupied indexed by year and region
    
    for (i in 1:nsite){
      for(k in 1:K){
        A[i,k] <- psi[i,k]*area[i,k]
      }
    }




#####################

# 2.2. Territory model 

# Input includes area occupied (A) indexed by year (k) and region (r) from occupancy model (2.1.)
# Area occupied is divided by territory size (T)
#   Simple model of the mean based on results from Rich et al. 2012

# Output is number of packs (P) indexed by year (k) and region (r)
    
####################



    P[k,r] <- A[k,r] / T

????????????????????????????????????????????????????????????????????????



#####################

# 2.3. Survival likelihood 
    
####################





#####################

# 2.4. Group level counts likelihood 

# Input data are group counts (y.group)
# Input estimates are survival (s) from survival model indexed by year (k) and region (r) and
#   recruitment (number of pups per pack, gamma) indexed by year, region, and group (i)

# Output is mean estimate of group size (G) which are indexed by year, region, and group
  
####################

# Ecological model/ system process

    for(k in 2:nyears){
      for(r in 1:nregions){
        for(i in 1:ngroups){
          g.mu[k,r,i] <- G[k-1,r,i] * (s[k-1,r] + delta + eps) + gamma[k-1,r,i]
          G[k,r,i] ~ dpois(g.mu[k,r,i])
        }
      }
    }

# Observation proccess

    for(k in 1:nyears){
      for(r in 1:nregions){
        for(i in 1:ngroups){
          y.group[k,r,i] ~ dnorm(G[k,r,i], tauy.group)
        }
      }
    }

# Derived parameters

    for(k in 1:nyears){
      for(r in 1:nregions){
        G.mean[k,r] <- mean(G[k,r,])
        gamma.mean[k,r] <- mean(gamma[k,r,])
      }
    }




#####################

# 2.5. Population level model 
    
####################

# Ecological model/ system process

    for(k in 1:nyears){
      for(r in 1:nregions){
        N.mu[k,r]<-P[k,r]*G.mean[k,r]
      }
    }

    for(k in 1:nyears){
      for(r in 1:nregions){
        N.mu[k,r] <- N.mu[k-1,r] * (s[k,r] + delta + eps) + gamma.mean[k-1,r] * P[k-1,r]
        N[k,r] <- dpois(N.mu[k,r])
      }
    }
    

# Derived parameters
    
    ????????????????????????????????????????????????????????????


#####################

# 2.6. Recruitment model 
    
####################


# Generalized linear model with log link function for recruitment

  for(k in 1:nyears){
    for(r in 1:nregions){
      for(i in 1:ngroups){
        log(mu.gamma[k,r,i]) <- b0.gam
        gamma[k,r,i] <- dpois(mu.gamma[k,r,i])
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

    