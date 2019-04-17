#####################################################################################I

#             Integrated Population Model to Estimate Recruitment of Wolves
#             POPULAITON LEVEL ONLY


# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2018

#             This model is adjusted to use only collar and occupancy data from MT
#             to estiamte recruitment from 2007-2016 using an IPM. The IPM model 
#             is outlined in word document ?????


#####################################################################################I

# Set the working directory I want to save files to

setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results")

#### Bring in data ####

# Set the memory to max so it can store the output file
library(snowfall)
memory.limit(size = 7500000) #just shy of 8 tb limit i think


# Pull in the encounter/detection histories, site covariates, and 
# ACV and HuntDays covariates
encounter <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/DetectionHistoriesPC.csv", row.names=1)
sitecovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/SiteCovars.csv")
ACV <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/ACV.csv", row.names=1)
HuntDays <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/HuntDays.csv", row.names=1)
MapPPN <-read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/MappedPPN.csv", row.names=1)


# collar data for survival model
y.surv <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/ysurv_subset2.csv")


# Mean and SD group sizes from good and moderate quality counts
POM.Group <- read.csv("C:/Users/allison/Documents/Project/WolfData/OccupancyData/GroupSize.csv")




#### Organize/format data for use in model ####

# Load packages for manipulating data
library(tidyr)
library(dplyr)

# Reorganize the 2d matrix - site (row) by occasion*year (columns) - to a 3d 
# array - site (row; 695) by occasion (column; 5) by year (3d; 10). You add
# 1 to encounter histories because the data can't have 0s for JAGS. So the 
# detection data are now 1, 2 and 3 instead of 0, 1, and 2. Then set the
# number of sites, occasions, and years
y.occ <- array(unlist(encounter), dim=c(695, 5, 10))+1
nsites <- nrow(y.occ)
noccs <- ncol(y.occ)
nyears <- 10


# You can double check that the encounter data was correctly set up by
# comparing the original data brought in and the y.occ data using the
# following code, remember it should be +1 for every value:
encounter[1:5, 1:5]
y.occ[1:5,,1]


# The survival data currently includes 2007-2017, however I only need until 2016
# (9 years {nyears - 1}) worth of data so I need to remove the data that I do not
# need
y.surv <- y.surv[y.surv$year < 2016,]


# The width interval is the number of months in each period. Because they are 
# not all the same, the hazard is adjusted by multiplying by the number of 
# months in each period to account for differences when calculating the 
# cumulative harzard (H) and survival
width.interval = c(2, 3, 3, 4)


# The number of periods (which is currently 4) and the number of observations
nperiods = length(unique(y.surv$period))
nobs =  nrow(y.surv)


# Set my constants: em.group is a 10% dispersal rate from a pack, 
# territory size is set to 600 sq km, and territory overlap is set to the 
# values determined by FWP
em.group <- 0.08
mu.T <- 599.83
sd.T <- 368.21
shape <- 3.15727
rate <- 0.0052635
T.overlap <- c(1.12, 1.08, 1.13, 1.16, 1.26, 1.27, 1.33, 1.24, 1.26, 1.32)





#### Model code ####


# Call for appropriate packages


library(R2jags)
library(mcmcplots)


################################################################################I
#  Specify model in BUGS language

sink("POMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
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
    T2 ~ dgamma(3.157, 0.00526)
    T3 ~ dlnorm(6.22985815, 1/0.58728123)
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    for(k in 1:nyears){
    P[k] <- (A[k] / (T3 + 0.000001)) * T.overlap[k]
    }
    
    # Pull in group count data to determine mean group size each year with error
    
    for(k in 1:nyears){
      G[k] ~ dnorm(mu.G[k], 1 / (sd.G[k] * sd.G[k] + 0.000001))T(0,)
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


# Data
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

# Parameters to keep track of and report
params <- c("P","N.est","phi", 
            "colo", "psi", "B0.phi", "B0.colo", "b.pc1.colo", "b.pcks5yrs.colo",
            "A", "p11", "p10", "b") 


# MCMC Settings 
ni <- 100000
nt <- 2
nb <- 10000
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "POMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)


#### format output for population level code ####

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$N.est
n.est[,2] <- out$BUGSoutput$sd$N.est
P2 <- array(NA, dim=c(nyears, 2))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$sd$P
phi2 <- array(NA, dim=c(nsites, nyears-1, 2))
phi2[,,1] <- out$BUGSoutput$mean$phi
phi2[,,2] <- out$BUGSoutput$sd$phi
colo2 <- array(NA, dim=c(nsites, nyears-1, 2))
colo2[,,1] <- out$BUGSoutput$mean$colo
colo2[,,2] <- out$BUGSoutput$sd$colo
b2 <- array(NA, dim=c(nsites, noccs, nyears, 2))
b2[,,,1] <- out$BUGSoutput$mean$b
b2[,,,2] <- out$BUGSoutput$sd$b
p102 <- array(NA, dim=c(nyears,1, 2))
p102[,,1] <- out$BUGSoutput$mean$p10
p102[,,2] <- out$BUGSoutput$sd$p10
p112 <- array(NA, dim=c(nyears,1, 2))
p112[,,1] <- out$BUGSoutput$mean$p11
p112[,,2] <- out$BUGSoutput$sd$p11
psi2 <- array(NA, dim=c(nsites, nyears, 2))
psi2[,,1] <- out$BUGSoutput$mean$psi
psi2[,,2] <- out$BUGSoutput$sd$psi
betas <- data.frame("B0.phi" = out$BUGSoutput$mean$B0.phi, 
                    "B0.phi.sd" = out$BUGSoutput$sd$B0.phi,
                    "B0.colo" = out$BUGSoutput$mean$B0.colo, 
                    "B0.colo.sd" = out$BUGSoutput$sd$B0.colo,
                    "b.pc1.colo" = out$BUGSoutput$mean$b.pc1.colo, 
                    "b.pc1.colo.sd" = out$BUGSoutput$sd$b.pc1.colo,
                    "b.pcks5yrs.colo" = out$BUGSoutput$mean$b.pcks5yrs.colo, 
                    "b.pcks5yrs.colo.sd" = out$BUGSoutput$sd$b.pcks5yrs.colo)

#### population level code ####

sink("PopLevel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## Population priors
    
    # Initial population size
    N.tot[1] ~ dnorm(600, 0.0001)I(0,)

    # Prior for mean recruitment
    B0.gam ~ dunif(-10,10)
    
    
    ## Bring in data P, colo, and phi
    
    for(k in 1:nyears){
      P[k] ~ dnorm(P2[k,1], 1 / (P2[k,2] * P2[k,2]+ 0.000001))
    }
    
    for(k in 1:(nyears-1)){
      B0.phi[k] ~ dnorm(betas[k,1], 1 / (betas[k,2] * betas[k,2]+ 0.000001))
      B0.colo[k] ~ dnorm(betas[k,3], 1 / (betas[k,4] * betas[k,4]+ 0.000001))
    }
    
    b.pc1.colo ~ dnorm(betas[1,5], 1 / (betas[1,6] * betas[1,6]+ 0.000001))
    b.pcks5yrs.colo ~ dnorm(betas[1,7], 1 / (betas[1,8] * betas[1,8]+ 0.000001))
    
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    #####################
    
    # 2.1. Survival likelihood 
    
    # Current model is 
    # Output is survival indexed by year (k)
    
    ####################
    
    # Estimate the harzard. 
    # This part transforms the linear predictor (mu.surv)
    # using the cloglog link and relates it to the data (event) for each 
    # observation
    
    for(i in 1:nobs){
    event[i] ~ dbern(mu.surv[i])
    cloglog(mu.surv[i]) <- b0.surv + b.period.surv[Period[i]] + eps.surv[Year[i]]
    }#i
    
    
    # Predicted values
    
    # Baseline hazard
    
    for(k in 1:(nyears-1)){
    for(p in 1:nperiods){
    cloglog(mu.pred[p,k]) <- b0.surv + b.period.surv[p] + eps.surv[k]
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
    
    
    #####################
    
    # 2.2. Population likelihood 
    
    ####################

    # Ecological model/ system process
    
    # First determine colonization and extinction 
    
    for(i in 1:nsites){
    for(k in 1:(nyears-1)){
    colo[i,k] <- B0.colo[k] + b.pc1.colo * PC1[i] + b.pcks5yrs.colo * pcks5yrs[i,k]
    }
    }
    
    for(k in 1:(nyears-1)){
    phi[k] <- B0.phi[k]      
    }
    
    for(k in 2:nyears){
    N.rec[k] ~ dpois(P[k-1] * gamma.mean[k-1])
    N.ps[k] ~ dpois((P[k-1] + sum(colo[,k-1]*((area[]*T.overlap[k])/T - P[k-1])) - P[k-1] * (1 - phi[k-1])) * G.mean[k-1])
    N.ad[k] ~ dbin(s[k-1], N.ps[k])
    N.tot[k] <- N.ad[k] + N.rec[k]
    }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    n.est[k,1] ~ dnorm(N.tot[k], (1 / (n.est[k,2]*n.est[k,2]+.001)))
    }
    
    # Recruitment
    for(k in 1:(nyears-1)){
      log(mu.gamma[k]) <- B0.gam
      gamma.mean[k] ~ dpois(mu.gamma[k])
    }
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data2 <- list("nyears"=nyears, "phi2"=phi2, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "s2"=s2, "gamma2"=gamma2, "P2"=P2, 
                  "colo2"=colo2, "T"=Tr, "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "pcks5yrs"=sitecovs[,24:31], "betas"=betas)


#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("P", "s", "N.tot", "gamma.mean") 


# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

print(out2, dig=2)
mcmcplot(out2)

