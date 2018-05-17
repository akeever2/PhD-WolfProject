# Input data and then set up for analyses

setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results")
encounter <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/encounterhistory.csv")
sitecovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/sitecovariates.csv")
survcovs <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/surveycovariates.csv")

y.occ <- array(unlist(encounter[,2:41]), dim=c(695, 5, 8))+1
nsites <- nrow(y.occ)
noccs <- ncol(y.occ)
nyears <- 8

POM.Group <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/GroupSize.csv")

T.overlap <- c(1.12, 1.08, 1.13, 1.16, 1.26, 1.27, 1.33, 1.24, 1.26, 1.32)
mu.T <- 599.83
sd.T <- 368.21


sink("MTOccupancy.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Occupancy priors
    
    #  psi1 coefficient (occupancy in year 1)
    
    B0.psi1 ~ dnorm(0,0.001)	  
    
    
    #  Priors for transition probabilities (survival and colonization)
    
    for(k in 1:(nyears-1)){
      B0.phi[k] ~ dnorm(0,0.001)
      B0.gamma[k] ~ dnorm(0,0.001)
    }#k
    
    
    # Priors for detection probabilities (only varies by year)
    
    for(k in 1:nyears){
      B0.p11[k] ~ dnorm(0,0.001)
      B0.p10[k] ~ dnorm(0,0.001)
    }#k

    B0.b ~ dnorm(0,0.001)
    
    # Priors for covariates
    
    b.pc1.psi ~ dnorm(0,0.001) 
    b.pcks5yrs.psi ~ dnorm(0,0.001) 
    b.pc1.gam ~ dnorm(0,0.001) 
    b.pcks5yrs.gam ~ dnorm(0,0.001) 
    b.area ~ dnorm(0,0.001) 
    b.huntdays.p11 ~ dnorm(0,0.001) 
    b.acv ~ dnorm(0,0.001) 
    b.huntdays.p10 ~ dnorm(0,0.001) 
    b.nonfrrds.p10 ~ dnorm(0,0.001)
    b.nonfrrds.p11 ~ dnorm(0,0.001)
    b.frrds.p11 ~ dnorm(0,0.001) 
    b.frrds.p10 ~ dnorm(0,0.001) 

    
    
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
      logit.psi1[i] <- B0.psi1 + b.pc1.psi * PC1[i] + b.pcks5yrs.psi * pcks5yrs[i,1]       
      logit(psi1[i]) <- logit.psi1[i]                                     
      z[i,1] ~ dbern(psi1[i])
 
      for(k in 1:(nyears-1)){
        logit.phi[i,k] <- B0.phi[k] 
        logit.gamma[i,k] <- B0.gamma[k] + b.pc1.gam * PC1[i] + b.pcks5yrs.gam * pcks5yrs[i,k] 
        logit(phi[i,k]) <- logit.phi[i,k]
        logit(gamma[i,k]) <- logit.gamma[i,k]
      }#k
    
      for(k in 2:nyears){
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
          multilogit.p11[i,j,k] <- B0.p11[k] + b.area * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i]
          logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
          multilogit.p10[i,j,k] <- B0.p10[k] + b.acv * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
          logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
          multilogit.b[i,j,k] <- B0.b
          logit(b[i,j,k]) <- multilogit.b[i,j,k]
          
          y.occ[i,j,k] ~ dcat(p[,i,j,k,(z[i,k]+1)])
        }#k
      }#j
    }#i
    
    
    # Derived parameters
    
    for(i in 1:nsites){
      psi[i,1] <- psi1[i]
      growthr[i,1] <- 1

      for (k in 2:nyears){
        psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*gamma[i,k-1]
        growthr[i,k] <- psi[i,k]/psi[i,k-1]
        turnover[i,k-1] <- (1 - psi[i,k-1]) * gamma[i,k-1]/psi[i,k]
      }#k
    }#i
    
    #  Area occpupied indexed by year
    
    for(k in 1:nyears){
      A[k] <- sum(psi[,k] * area[])
    }
    

    # Pull in data for the mean for territory size, this will also incorporate error
    # into the estimates of pack size from the average territory size and area
    # area occupied by wolves
    
    T ~ dnorm(mu.T, 1/(sd.T * sd.T))
    
    # Estimate number of packs from area occupied (A), territory size (T), and 
    # territory overlap (T.overlap). Territory overlap is a fixed value for each
    # year that is supplied as data. 
    
    for(k in 1:nyears){
      P[k] <- (A[k] / T)*T.overlap[k]
    }

    
    # Pull in group count data to determine mean group size each year with error
    
    for(k in 1:nyears){
      G[k] ~ dnorm(mu.G[k], 1/(sd.G[k] * sd.G[k]))    
    }


    # Estimate abundance each year based on estimated # of packs (P) and mean group
    # size (G). Then, add on the lone wolves in the population

    for(k in 1:nyears){
      N.est[k] <- P[k] * G[k]
      N.total[k] <- N.est[k] * 1.125
    }




    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "area"=sitecovs$AREASAMP,
                 "noccs"=noccs, "y.occ"=y.occ, "PC1"=sitecovs$PC1, 
                 "pcks5yrs"=sitecovs[,24:31], 
                 "huntdays"=array(unlist(survcovs$HUNTDAYS), dim=c(695, 5, 8)),
                 "nonforrds"=sitecovs$LOWUSENONFORESTRDS, 
                 "forrds"=sitecovs$LOWUSEFORESTRDS, 
                 "acv"=array(unlist(survcovs$ACV), dim=c(695, 5, 8)),
                 "T.overlap" =T.overlap, "mu.T" = mu.T, "sd.T" = sd.T,
                 "mu.G" = POM.Group[,2], "sd.G" = POM.Group[,3])


#  Initial Values	
zst <- apply(y.occ,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 0
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(z=zst)}


# Parameters to keep track of and report
params <- c("P", "phi", "gamma", "psi", 
            "p11", "p10", "b", "b.pc1.psi", 
            "b.pcks5yrs.psi", "b.pc1.gam", 
            "b.pcks5yrs.gam", "b.area", 
            "b.huntdays.p11", "b.acv",  
            "b.huntdays.p10", "b.nonfrrds.p10",
            "b.nonfrrds.p11", "b.frrds.p11",  
            "b.frrds.p10", "B0.phi", "B0.psi1", 
            "B0.gamma", "B0.p11", "B0.p10", "B0.b", "N.est", "N.total") 


# MCMC Settings 
ni <- 10000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "MTOccupancy.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)

