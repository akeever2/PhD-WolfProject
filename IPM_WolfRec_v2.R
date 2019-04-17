####################################################################I

#     Integrated Population Model to Estimate Recruitment of Wolves


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


# This code uses changes in abundance and group size and estimates of
# survival to estimate recruitment of wolves. Abundance estiatmes are 
# generated from the POM abundance estimation framework established by
# MFWP (CITE). Group sizes come from group count data collected by 
# MFWP. Survival is estimated from GPS and VHF collar data collected
# by MFWP and collaborating agencies. Data for this analysis is from
# 2007-2016. This model is described in further detail in (CITE)


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
# Territory - Territory overlap, mean territory size, and SD territory size 
#             are hard coded. Territory overlap numbers come from MFWP. Mean 
#             and SD territory sizes were estimated from Lindsey Rich's 
#             territory work (CITE).
# Dispersal - Dispersal at the pack level is set at a fixed rate of 10%,
#             which comes from (CITE). 
# Group count data - These data are organized in long format with pack
#                     counts and other information from packs for each
#                     year. In the code I organize to wide format.
# Pack covariates - These data are a organized by pack (row) and pack
#                   specific covariates, including FWP region, forest
#                   cover and road (4wd and 2wd) density. 
# Collar data - This is a subset of the collar data for only adults
#               and yearlings that weren't removed for control action.
#               These data are expanded, so each individual wolf has
#               an observation for each time period and year it was
#               alive. It has a start time, a stop time, and the event
#               (1 or 0) for whether it failed (died, 1) or not (0). 
#               If the wolf was censored during a time period, then 
#               the last observation of that wolf will be from the 
#               previous time period. It already has all covariate data
#               needed. 




#### Bring in data ####

# Set the working directory I want to save files to
setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results")


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


# Pull in group count data and pack level covariates
group <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/GroupCountData.csv")
packcovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/FinalPackCovs.csv")


# Pull in collar data for survival model
y.surv <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/ysurv_subset2.csv")




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


# For the group count data, pull out only "good" (G) and "moderate" (M)
# quality counts and call it group.g
group.g <- group[group$Count.Quality=="G" | group$Count.Quality=="M", ]


# Make a new dataframe, g2,  with only the year, pack, recreation area, 
# and count of groups from the good and moderate counts
#### NOTE: check to see if using EoY counts gives better/different results ####
g2 <- data.frame(year=group.g$YEAR, pack=group.g$Revised.Pack.Name, area=group.g$RECAREA,  
                 count=group.g$CountPlusRemovals)
#### end NOTE ####

# Change the data from long format to wide format
g2 <-g2 %>% spread(year, count)


# Some packs were only around prior to 2007, so get rid of packs/group 
# counts that do not have any counts in years 2007-2016. Double check 
# that the specified columns correspond to the correct years of data.
g3 <- g2[as.logical((rowSums(is.na(g2[,22:31]))-nyears)),]


# The final group count data only needs to be from 2007-2016, so 
# double check that the specified columns correspond to the correct
# years of data and use only those columns. Final group count data 
# frame where the number of packs is the number of rows and the 
# columns are the counts per year. Then set the number of groups.
# Then set the number of groups
y.group <- g3[,22:31]
ngroups <- nrow(y.group)


# Set up group covariate data to only include packs in the group counts
# and be in the correct order. 
groupcovs <- merge(g3[,1:2], packcovs, by.x="pack", by.y="Pack", sort=FALSE)


# Now get rid of all the other random dataframes that aren't needed
rm(g2, g3, group.g)


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


###################################################I
#  Specify model in BUGS language

sink("RecIPM_m1.txt")
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
    for(k in 1:(nyears-1)){
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
    # for(k in 1:nyears){
    #   eps.gam[k] ~ dnorm(0, tau.gam)
    # }
    # 
    # sigma.gam ~ dunif(0,100)
    # tau.gam <- pow(sigma.gam, -2)
    # var.gam <- pow(sigma.gam, 2)
    
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
    T2 ~ dgamma(3.157, 0.00526)
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
    
    for(k in 1:(nyears-1)){
      for(p in 1:nperiods){
        cloglog(mu.pred[p,k]) <- b.period.surv[p] + eps.surv[k]
        hazard[p,k] <- -log(1 - mu.pred[p,k])
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
      G.mean[k] <- mean(G[,k])
      gamma.mean[k] <- mean(gamma[,k])
      n.est[k] <- P[k] * G.mean[k]
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


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "ngroups"=ngroups, 
                 "area"=sitecovs$AREASAMP, "noccs"=noccs, "y.occ"=y.occ, 
                 "PC1"=sitecovs$PC1, "recPC"=sitecovs[,27:36], 
                 "huntdays"=array(unlist(HuntDays), dim=c(695, noccs, nyears)),
                 "nonforrds"=sitecovs$LOWUSENONFORESTRDS, 
                 "forrds"=sitecovs$LOWUSEFORESTRDS, 
                 "acv"=array(unlist(ACV), dim=c(695, noccs, nyears)),
                 "mapppn"=array(unlist(MapPPN), dim=c(695, noccs, nyears)),
                 "em.group"=em.group, "y.group"=y.group, 
                 "T.overlap"=T.overlap, "mu.T"=mu.T, "sd.T"=sd.T,
                 "event"=y.surv$Event, 
                 "nperiods"=length(unique(y.surv$period)),  
                 "Period" = y.surv$period, 
                 "Year" = as.factor(y.surv$year),
                 "width.interval"=width.interval, 
                 "nobs" =  nrow(y.surv))


#  Initial Values	
zst <- apply(y.occ,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 1
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(B0.gam=runif(1,-1,1), sd.proc=runif(1,0,10), 
                         sigma.group=runif(1,0,10), z=zst, B0.colo=runif((nyears-1),-6,-3), b.pc1.colo=runif(1,-2,-1), b.recPC.colo=runif(1,1,2),
                         B0.psi1=runif(1,-5,-3), b.pc1.psi=runif(1,-1,1), b.recPC.psi=runif(1,1,2), B0.phi=runif(1,-1,1), b.pc1.phi=runif(1,1,2),
                         B0.p10=runif(1,-4,-3), b.huntdays.p10=runif(1,-1,1), b.nonfrrds.p10=runif(1,-1,1), b.frrds.p10=runif(1,-1,1),
                         B0.p11=runif(1,-8,-7))}

#### !!!add initial values!!! ####

# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", 
            "B0.phi", "B0.colo", "b.pc1.colo", 
            "b.recPC.colo",  
            "b.pc1.colo", 
            "b.recPC.phi",
            "b.pc1.phi", 
            "B0.gam", "A", "gamma", "G", "B1.gam", "b.period.surv", 
            "eps.surv", 
            "var.surv", "annual.s") 


# MCMC Settings 
ni <- 10000
nt <- 2
nb <- 1000
nc <- 3


# Call JAGS 
out_m1 <- jags(win.data, inits, params, "RecIPM_m1.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
               n.burnin=nb, jags.module = c("glm", "dic"))

#print(out_m1, dig=2)

#mcmcplot(out)


#### format output for population level code ####

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out_m1$BUGSoutput$mean$n.est
n.est[,2] <- out_m1$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- out_m1$BUGSoutput$mean$G.mean
G.mean2[,2] <- out_m1$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- out_m1$BUGSoutput$mean$annual.s
s2[,2] <- out_m1$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 2))
P2[,1] <- out_m1$BUGSoutput$mean$P
P2[,2] <- out_m1$BUGSoutput$sd$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- out_m1$BUGSoutput$mean$gamma.mean
gamma2[,2] <- out_m1$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(out_m1$BUGSoutput$mean$B0.phi,9), 
                    "B0.phi.sd" = rep(out_m1$BUGSoutput$sd$B0.phi,9),
                    "B0.colo" = out_m1$BUGSoutput$mean$B0.colo, 
                    "B0.colo.sd" = out_m1$BUGSoutput$sd$B0.colo,
                    "b.pc1.colo" = rep(out_m1$BUGSoutput$mean$b.pc1.colo,9), 
                    "b.pc1.colo.sd" = rep(out_m1$BUGSoutput$sd$b.pc1.colo,9),
                    "b.recPC.colo" = rep(out_m1$BUGSoutput$mean$b.recPC.colo,9), 
                    "b.recPC.colo.sd" = rep(out_m1$BUGSoutput$sd$b.recPC.colo,9),
                    "b.pc1.phi" = rep(out_m1$BUGSoutput$mean$b.pc1.phi,9), 
                    "b.pc1.phi.sd" = rep(out_m1$BUGSoutput$sd$b.pc1.phi,9),
                    "b.period.surv" = c(out_m1$BUGSoutput$mean$b.period.surv, 1,2,3,4,5), 
                    "b.period.surv.sd" = c(out_m1$BUGSoutput$sd$b.period.surv, 1,2,3,4,5),
                    "eps.surv" = out_m1$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = out_m1$BUGSoutput$sd$eps.surv)


#### population level code ####

sink("Rec_PopLevel_m1.txt")
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
      G.mean[k] ~ dnorm(G.mean2[k,1], 1 / (G.mean2[k,2] * G.mean2[k,2]+ 0.0000001))
      gamma.mean[k] ~ dnorm(gamma2[k,1], 1 / (gamma2[k,2] * gamma2[k,2]+ 0.0000001))
    }
    
    for(k in 1:(nyears-1)){
      B0.colo[k] ~ dnorm(betas[k,3], 1 / (betas[k,4] * betas[k,4]+ 0.0000001))
      eps.surv[k] ~ dnorm(betas[k,13], 1 / (betas[k,14] * betas[k,14]+ 0.0000001))
    }
    
    B0.phi ~ dnorm(betas[1,1], 1 / (betas[1,2] * betas[1,2]+ 0.0000001))
    b.pc1.colo ~ dnorm(betas[1,5], 1 / (betas[1,6] * betas[1,6]+ 0.0000001))
    b.recPC.colo ~ dnorm(betas[1,7], 1 / (betas[1,8] * betas[1,8]+ 0.0000001))
    b.pc1.phi ~ dnorm(betas[1,9], 1 / (betas[1,10] * betas[1,10]+ 0.0000001))
    
    for(p in 1:nperiods){
      b.period.surv[p] ~ dnorm(betas[p,11], 1 / (betas[p,12] * betas[p,12]+ 0.0000001))
    }
    
    T ~ dlnorm(6.22985815, 1/0.58728123)

    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    
    # Ecological model/ system process
    
    # First determine colonization and extinction 
    
    for(i in 1:nsites){
      for(k in 1:(nyears-1)){
        colo[i,k] <- B0.colo[k] + b.pc1.colo * PC1[i] + b.recPC.colo * recPC[i,k+1]
        phi[i,k] <- B0.phi + b.pc1.phi * PC1[i]
      }
    }
    
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
      N.ps[k] ~ dpois((P[k-1] + sum(colo[,k-1]*((area[] * T.overlap[k])/T - P[k-1])) - P[k-1] * (1 - sum(phi[,k-1]))) * G.mean[k-1])
      N.ad[k] ~ dbin(annual.s[k-1], N.ps[k])
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


# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas)


#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("P", "annual.s", "N.tot", "gamma.mean", "G.mean", 
             "N.rec") 


# Call JAGS 
out2_m1 <- jags(win.data2, inits2, params2, "Rec_PopLevel_m1.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                n.burnin=nb, jags.module = c("glm", "dic"))

print(out2_m1, dig=2)
#mcmcplot(out2_m1)


output <- data.frame("gamma" = out2_m1$BUGSoutput$mean$gamma.mean,
                     "gamma.sd" = out2_m1$BUGSoutput$sd$gamma.mean,
                     "N.rec" = c(out2_m1$BUGSoutput$mean$N.rec, NA),
                     "N.rec.sd" = c(out2_m1$BUGSoutput$sd$N.rec, NA), 
                     "year" = c(2007:2016))

library(ggplot2)

ggplot(output, aes(x=year, y=gamma))+
  geom_point()+
  geom_errorbar(aes(ymin=gamma-gamma.sd, ymax=gamma+gamma.sd), colour="black", width=.1)+
  scale_y_continuous(limits = c(0,4))+
  scale_x_continuous(breaks = c(2007:2016))

ggplot(output, aes(x=year, y=N.rec))+
  geom_point()+
  geom_errorbar(aes(ymin=N.rec-N.rec.sd, ymax=N.rec+N.rec.sd), colour="black", width=.1)+
  scale_y_continuous(limits = c(0,450))+
  scale_x_continuous(breaks = c(2007:2016))
