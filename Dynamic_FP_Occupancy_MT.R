###############################################################################/


#                         POM ABUNDANCE ESTIMATION


# Allison C. Keever
# akeever1122@gmail.com
# github.com/akeever2
# Montana Cooperative Wildlife Research Unit
# May 2018

# Occupancy code adapted from Sarah B Bassing
# sarah.bassing@gmail.com
# github.com/SarahBassing
# Montana Cooperative Wildlife Research Unit
# August 2016

###############################################################################/


# This code uses the POM abundance estimation framework established by MFWP. It
# uses patch occupancy models (POM) to estimate area occupied by wolves and then
# divides by mean territory size to determine number of pack that fit into the
# area. We then multiply by average group size to get an abundance estiamte for
# each year. This code also accounts for territory overlap and lone wolves in
# the same way that MFWP has done in the past. 


# Inputs needed:
  # Encounter/detection histories - These should be 0 (no detection), 
  #                                 1 (uncertain detection), and 2 
  #                                 (certain detection). Rows are the sites and
  #                                 columns are the survey occasions. Occasions
  #                                 are 5 secondary for 10 years (50 total)
  # Site covariates - Site specific covariates. Rows are the sites and columns
  #                   are the different covariates
  # ACV - Survey covariate that changes for each site and occasion/year. The
  #       rows are the sites and the columns are the occasions and years (50)
  # HuntDays - Survey covariate that changes for each site and occasion/year. 
  #             The rows are the sites and the columns are the occasions and
  #             years (50 total)
  # MapPPN - Survey covariate that changes for each site and occasion/year.
  #           The rows are the sites and the columns are the occasions and 
  #           years (50 total)
  # POM.Group - This is average and SD of group sizes for each year. Group sizes
  #             were based off of good or moderate quality counts. The rows are
  #             years and the second column is average size and the third column
  #             is the SD for that year. 
  # Territory - Territory overlap, mean territory size, and SD territory size 
  #             are hard coded. Territory overlap numbers come from MFWP. Mean 
  #             and SD territory sizes were estimated from Lindsey Rich's 
  #             territory work (CITE). 




#### Input & Format Data and Load Packages ####


# Load the appropriate packages
library(R2jags)
library(mcmcplots)
library(ggplot2)
library(snowfall)


# Set working directory, which is where output files will be stored
setwd("C:/Users/allison/Documents/Project/WolfData/OccupancyData/RJagsResults2018_Keever")


# Set the memory to max so it can store the output file
memory.limit(size = 7500000) #just shy of 8 tb limit i think
memory.limit(size = NA)


# Pull in the encounter/detection histories, site covariates, group sizes, and 
# ACV and HuntDays covariates
encounter <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/DetectionHistoriesPC.csv", row.names=1)
sitecovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/SiteCovars.csv")
ACV <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/ACV.csv", row.names=1)
HuntDays <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/HuntDays.csv", row.names=1)
MapPPN <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/MappedPPN.csv", row.names=1)
POM.Group <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/GroupSize.csv")


# Reorganize the 2d matrix - site (row) by occasion*year (columns) - to a 3d 
# array - site (row; 695) by occasion (column; 5) by year (3d; 10). You add
# 1 to encounter histories because the data can't have 0s for JAGS. So the 
# detection data are now 1, 2 and 3 instead of 0, 1, and 2. 
y.occ <- array(unlist(encounter), dim=c(695, 5, 10))+1


# Set the number of sites, number of occasions, and number of years
nsites <- nrow(y.occ)
noccs <- ncol(y.occ)
nyears <- 10


# Hard code in territory overlap and the average and SD of territory sizes
T.overlap <- c(1.12, 1.08, 1.13, 1.16, 1.26, 1.27, 1.33, 1.24, 1.26, 1.32)
mu.T <- 599.83
sd.T <- 368.21




#### Model Code ####

sink("MTOccupancy.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    #  psi1 coefficient (occupancy in year 1)
    
    B0.psi1 ~ dnorm(0,0.001)	  
    
    
    #  Priors for transition probabilities (survival and colonization)
    
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
    
    
    
    ############################################################
    
    #             2. Likelihood
    
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
      logit.psi1[i] <- B0.psi1 + b.pc1.psi * PC1[i] + b.recPC.psi * recPC[i,1]       
      logit(psi1[i]) <- logit.psi1[i]                                     
      z[i,1] ~ dbern(psi1[i])
    
      for(k in 1:(nyears-1)){
        logit.phi[i,k] <- B0.phi + b.pc1.phi * PC1[i]  
        logit.gamma[i,k] <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1] 
        logit(phi[i,k]) <- logit.phi[i,k]
        logit(gamma[i,k]) <- logit.gamma[i,k]
      }#k
    
      for(k in 2:nyears){
        muZ[i,k] <- z[i,k-1]*(1-phi[i,k-1]) + (1-z[i,k-1])*gamma[i,k-1]
        z[i,k] ~ dbern(muZ[i,k])
      }#k
    }#i
    
    
    #  Observation process/submodel
    #  z is either 0 or 1 (unoccupied or occupied) but needs to be 1 or 2, so we add 1 in the y.occ
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
          multilogit.p11[i,j,k] <- B0.p11 + b.area.p11 * area[i] + b.huntdays.p11 * huntdays[i,j,k] + b.nonfrrds.p11 * nonforrds[i] + b.frrds.p11 * forrds[i] + b.acv.p11 * acv[i,j,k] + b.map.p11 * mapppn[i,j,k]
          logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
          multilogit.p10[i,j,k] <- B0.p10 + b.acv.p10 * acv[i,j,k] + b.huntdays.p10 * huntdays[i,j,k] + b.nonfrrds.p10 * nonforrds[i] + b.frrds.p10 * forrds[i]
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
        psi[i,k] <- psi[i,k-1]*(1-phi[i,k-1]) + (1-psi[i,k-1])*gamma[i,k-1]
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
    T2 ~ dgamma(3.157, 0.00526)
    T3 ~ dlnorm(6.22985815, 1/0.58728123)

    # Estimate number of packs from area occupied (A), territory size (T), and
    # territory overlap (T.overlap). Territory overlap is a fixed value for each
    # year that is supplied as data.

    for(k in 1:nyears){
      P[k] <- (A[k] / T)*T.overlap[k]
    }

    for(k in 1:nyears){
      P2[k] <- (A[k] / T2)*T.overlap[k]
    }

     for(k in 1:nyears){
      P3[k] <- (A[k] / T3)*T.overlap[k]
    }

    # Pull in group count data to determine mean group size each year with error

    for(k in 1:nyears){
      G[k] ~ dnorm(mu.G[k], 1/(sd.G[k] * sd.G[k]))
    }


    # Estimate abundance each year based on estimated # of packs (P) and mean group
    # size (G). Then, add on the lone wolves in the population

    for(k in 1:nyears){
      N.est[k] <- P3[k] * G[k]
      N.total[k] <- N.est[k] * 1.125
    }
    
    
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


#### Set-up and Run Model ####
# Bundle data for JAGs. Double check that recPC data and the mean group size
# and SD group size are referencing the right columns for the covariate. 
win.data <- list("nsites"=nsites, "nyears"=nyears, "area"=sitecovs$AREASAMP,
                 "noccs"=noccs, "y.occ"=y.occ, "PC1"=sitecovs$PC1, 
                 "recPC"=sitecovs[,27:36], 
                 "huntdays"=array(unlist(HuntDays), dim=c(695, noccs, nyears)),
                 "mapppn"=array(unlist(MapPPN), dim=c(695, noccs, nyears)),
                 "nonforrds"=sitecovs$LOWUSENONFORESTRDS, 
                 "forrds"=sitecovs$LOWUSEFORESTRDS, 
                 "acv"=array(unlist(ACV), dim=c(695, noccs, nyears)),
                 "T.overlap" =T.overlap, "mu.T" = mu.T, "sd.T" = sd.T,
                 "mu.G" = POM.Group[,2], "sd.G" = POM.Group[,3])


# Set initial values for true state, z. Use the encounter history data to set z for each
# site and year. 
zst <- apply(y.occ,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 0
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

# If you are having convergence issues can set up intial values for beta coefficients
# E.G., B0.psi1=runif(1,-1,1) for 1 initial value or B0.phi=runif(nyears,-1,1)
# for covariates that need more than 1 initial value
inits <- function(){list(z=zst, B0.colo=runif((nyears-1),-6,-3), b.pc1.colo=runif(1,-2,-1), b.recPC.colo=runif(1,1,2),
                         B0.psi1=runif(1,-5,-3), b.pc1.psi=runif(1,-1,1), b.recPC.psi=runif(1,1,2), B0.phi=runif(1,-1,1), b.pc1.phi=runif(1,1,2),
                         B0.p10=runif(1,-4,-3), b.huntdays.p10=runif(1,-1,1), b.nonfrrds.p10=runif(1,-1,1), b.frrds.p10=runif(1,-1,1),
                         B0.p11=runif(1,-8,-7))}

# inits <- function(){list(z=zst, B0.colo=runif((nyears-1),-6,-3), b.pc1.colo=runif(1,-2,-1), b.recPC.colo=runif(1,1,2),
#                          B0.psi1=runif(1,-5,-3), b.pc1.psi=runif(1,-1,1), b.recPC.psi=runif(1,1,2), B0.phi=runif(1,-1,1), b.pc1.phi=runif(1,1,2),
#                          B0.b=runif(1,5,7), B0.p10=runif(1,-4,-3), b.huntdays.p10=runif(1,-1,1), b.nonfrrds.p10=runif(1,-1,1), b.frrds.p10=runif(1,-1,1),
#                          b.acv.p10=runif(1,5,6), B0.p11=runif(1,-8,-7), b.area.p11=runif(1,-1,1), b.huntdays.p11=runif(1,-1,1), b.acv.p11=runif(1,1,2),
#                          b.map.p11=runif(1,3,4), b.nonfrrds.p11=runif(1,-1,1), b.frrds.p11=runif(1,-1,1))}

# Parameters to keep track of and report in the output file
params <- c("P", "P2", "phi", "gamma", "psi", 
            "p11", "p10", "b", "B0.phi",
            "B0.colo", "b.pc1.colo", 
            "b.recPC.colo", "B0.psi1", 
            "b.pc1.psi", "b.recPC.psi", 
            "b.pc1.phi", 
            "b.area.p11", 
            "b.huntdays.p11", 
            "b.acv.p11", 
            "b.map.p11", 
            "b.nonfrrds.p11",
            "b.frrds.p11",
            "b.huntdays.p10",
            "b.nonfrrds.p10",
            "b.frrds.p10",
            "b.acv.p10",
            "B0.p11", "B0.p10", "B0.b", "N.est", "N.total", "A", 
            "T", "T2", "P3") 


# MCMC Settings 
ni <- 5000
nt <- 2
nb <- 700
nc <- 3


# Call JAGS and run the model
out <- jags(win.data, inits, params, "MTOccupancy.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

# Show summary of the output
print(out, dig=2)

# Check convergence with MCMC plots
mcmcplot(out)


jag.sum <- out$BUGSoutput$summary
write.table(x=jag.sum, file="C:/Users/allison/Documents/Project/WolfData/OccupancyData/RJagsResults2018_Keever/MT_OccuResults2.txt", sep="\t")


#### Results Formatting and Plotting ####







