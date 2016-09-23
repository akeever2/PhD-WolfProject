##############################################################################

#                     Empirical recruitment model code            
# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2016

##############################################################################

# Pull in data and call for appropriate packages


library(R2jags)
library(mcmcplots)


################################################################################
#  Specify model in BUGS language

sink("EmpRecruitModel.txt")
cat("
    model {

############################################################

#             1. Priors

############################################################

#  Priors for OCCUPANCY parameters
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




#  Priors for ?????
  #  ?????







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
    







#####################

# 2.2. Pop level model 
    
####################






#####################

# 2.3. Group level counts likelihood 
  
####################



#####################

# 2.4. Survival likelihood 
    
####################




#####################

# 2.5. Recruitment model 
    
####################




############################################################

#             3. Derived parameters

############################################################






############################################################

#             4. Bugs requirements

############################################################


}", fill=TRUE)
sink()

    