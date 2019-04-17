####################################################################I

# Analyses of models for wolf recruitment in Montana

# For more details regarding the IPM see the rscript IPM_WolfRec_V2

# For models see the rscript Models_WolfRec_IPM


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

#### Bring in data ####

# Set the working directory I want to save files to
setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/ModelResultFiles")


# Set the memory to max so it can store the output file
library(snowfall)
memory.limit(size = 7500000) #just shy of 8 tb limit i think


# Pull in the encounter/detection histories, site covariates, and 
# ACV and HuntDays covariates
encounter <- read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/DetectionHistoriesPC.csv", row.names=1)
sitecovs <- read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/SiteCovars.csv")
ACV <-read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/ACV.csv", row.names=1)
HuntDays <-read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/HuntDays.csv", row.names=1)
MapPPN <-read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/MappedPPN.csv", row.names=1)


# Pull in group count data and pack level covariates
group <- read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/GroupCountData.csv")
packcovs <- read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/FinalPackCovs.csv")


# Pull in collar data for survival model
y.surv <- read.csv("C:/Users/allison/Documents/Project/Dissertation/Recruitment/ysurv_subset2.csv")




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

g2 <- data.frame(year=group.g$YEAR, pack=group.g$Revised.Pack.Name, area=group.g$RECAREA,  
                 count=group.g$CountPlusRemovals)


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
em.group <- 0.10
T.overlap <- c(1.12, 1.08, 1.13, 1.16, 1.26, 1.27, 1.33, 1.24, 1.26, 1.32)
Harv <- c(1,1,2,1,2,2,2,2,2,2)
LogN <- log(c(623,694,836,849,955,899,1065,878,961,851))
Winter <- c(22.55, 17.79, 25.16, 22.53, 19.86, 30.41, 20.54, 18.5, 26.6, 15.92)
meanG <- c(7.03,6.65,6.37,6.16,5.71,4.96,5.66,5.39,5.61,4.96)
sdG <- c(3.13,3.69,3.59,2.97,2.92,2.46,2.83,2.6,2.69,2.24)


#### CONSTANTS ####

library(R2jags)

# Data
Method <- c(1,1,2,1,2,3,3,3,3,3,3)
win.data <- list("nsites"=nsites, "nyears"=nyears, "ngroups"=ngroups, 
                 "area"=sitecovs$AREASAMP, "noccs"=noccs, "y.occ"=y.occ, 
                 "PC1"=sitecovs$PC1, "recPC"=sitecovs[,27:36], 
                 "huntdays"=array(unlist(HuntDays), dim=c(695, noccs, nyears)),
                 "nonforrds"=sitecovs$LOWUSENONFORESTRDS, 
                 "forrds"=sitecovs$LOWUSEFORESTRDS, 
                 "acv"=array(unlist(ACV), dim=c(695, noccs, nyears)),
                 "mapppn"=array(unlist(MapPPN), dim=c(695, noccs, nyears)),
                 "y.group"=y.group, 
                 "T.overlap"=T.overlap, 
                 "event"=y.surv$Event, 
                 "nperiods"=length(unique(y.surv$period)),  
                 "Period" = y.surv$period, 
                 "Year" = as.factor(y.surv$year),
                 "width.interval"=width.interval, 
                 "nobs" =  nrow(y.surv), 
                 "nregions"=length(unique(groupcovs$Region)), 
                 "GroupReg"=groupcovs$Region, 
                 "indicator"=ifelse(is.na(y.group) , 0 , 1), 
                 "Method" = Method, 
                 "Forest" = scale(groupcovs$ForProp)[,1],
                 "FourWD" = scale(groupcovs$X4wdDens)[,1],
                 "TwoWD" = scale(groupcovs$X2wdDens)[,1],
                 "Harv" = Harv, 
                 "LogN" = scale(LogN)[,1], 
                 "Winter" = scale(Winter)[,1], 
                 "mu.G" = meanG, 
                 "sd.G" =sdG)


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


# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", 
            "B0.gam","B1.gam", "b.period.surv", 
            "eps.surv", "var.surv", "annual.s",
            "eps.gam", "var.gam", "eps.reg", "var.reg", 
            "b.method", "b.dd", "b.harv", "fit.new", "fit",
            "b.4wd", "b.2wd", "b.forest", "pop.growth", 
            "b.winter", "n.est2", "G.mean.high") 


# MCMC Settings 
ni <- 70000
nt <- 3
nb <- 20000
nc <- 3


#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("P", "annual.s", "N.tot", "gamma.mean", 
             "N.rec") 

#### M1; R ~ pack size  ####
start1<-Sys.time()
# Call JAGS 
m1_group_out <- jags(win.data, inits, params, "M1_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))

m1_group_out <- autojags(m1_group_out, n.iter = 20000, n.update = 4, n.thin = 3)

# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m1_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m1_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m1_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m1_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m1_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m1_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m1_group_out$BUGSoutput$mean$P
P2[,2] <- m1_group_out$BUGSoutput$sd$P
P2[,3] <- m1_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m1_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m1_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10),
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m1_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m1_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m1_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m1_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m1_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m1_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m1_group_out)
save(datums, file ="m1_group_out.RData")
end1<-Sys.time()

start1-end1

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m1_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m1_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m1_pop_outSummary.txt", sep="\t")
save(m1_pop_out, file ="m1_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop
#### M2; R ~ pack size + ran year ####
start1<-Sys.time()
# Call JAGS 
m2_group_out <- jags(win.data, inits, params, "M2_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))

m2_group_out <- autojags(m2_group_out, n.iter = 20000, n.update = 4, n.thin = 3)

# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m2_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m2_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m2_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m2_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m2_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m2_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m2_group_out$BUGSoutput$mean$P
P2[,2] <- m2_group_out$BUGSoutput$sd$P
P2[,3] <- m2_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m2_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m2_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(1,9), 
                    "B0.phi.sd" = rep(1,9),
                    "B0.colo" = rep(1,9),
                    "B0.colo.sd" = rep(1,9),
                    "b.pc1.colo" = rep(1,9), 
                    "b.pc1.colo.sd" = rep(1,9),
                    "b.recPC.colo" = rep(1,9), 
                    "b.recPC.colo.sd" = rep(1,9),
                    "b.pc1.phi" = rep(1,9),
                    "b.pc1.phi.sd" = rep(1,9),
                    "b.period.surv" = c(m2_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5), 
                    "b.period.surv.sd" = c(m2_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5),
                    "eps.surv" = m2_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m2_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m2_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m2_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m2_group_out)
save(datums, file ="m2_group_out.RData")
end1<-Sys.time()

start1-end1

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m2_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m2_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m2_pop_outSummary.txt", sep="\t")
save(m2_pop_out, file ="m2_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M3; R ~ pack size + ran region ####
start1<-Sys.time()
# Call JAGS 
m3_group_out <- jags(win.data, inits, params, "M3_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))

m3_group_out <- autojags(m3_group_out, n.iter = 20000, n.update = 4, n.thin = 3)

# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m3_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m3_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m3_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m3_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m3_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m3_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m3_group_out$BUGSoutput$mean$P
P2[,2] <- m3_group_out$BUGSoutput$sd$P
P2[,3] <- m3_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m3_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m3_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10), 
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m3_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m3_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m3_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m3_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m3_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m3_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m3_group_out)
save(datums, file ="m3_group_out.RData")
end1<-Sys.time()

start1-end1

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m3_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m3_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m3_pop_outSummary.txt", sep="\t")
save(m3_pop_out, file ="m3_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop

#### M4; R ~ pack size + ran year + ran region ####
start1<-Sys.time()
# Call JAGS 
m4_group_out <- jags(win.data, inits, params, "M4_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))

m4_group_out <- autojags(m4_group_out, n.iter = 20000, n.update = 4, n.thin = 3)

# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m4_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m4_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m4_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m4_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m4_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m4_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 2))
P2[,1] <- m4_group_out$BUGSoutput$mean$P
P2[,2] <- m4_group_out$BUGSoutput$sd$P
P2[,3] <- m4_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m4_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m4_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10),
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m4_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m4_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m4_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m4_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m4_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m4_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m4_group_out)
save(datums, file ="m4_group_out.RData")
end1<-Sys.time()

start1-end1

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m4_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m4_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m4_pop_outSummary.txt", sep="\t")
save(m4_pop_out, file ="m4_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M5; R ~ pack size + ran year + ran region + method ####

start1<-Sys.time()
# Call JAGS 
m5_group_out <- jags(win.data, inits, params, "M5_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))
end1<-Sys.time()

start1-end1

start2<-Sys.time()
m5_group_out <- autojags(m5_group_out, n.iter = 20000, n.update = 4, n.thin = 3)



# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m5_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m5_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m5_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m5_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m5_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m5_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m5_group_out$BUGSoutput$mean$P
P2[,2] <- m5_group_out$BUGSoutput$sd$P
P2[,3] <- m5_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m5_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m5_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10), 
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m5_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m5_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m5_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m5_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m5_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m5_group_outSummary2.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m5_group_out)
save(datums, file ="m5_group_out2.RData")
end2<-Sys.time()

start2-end2

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m5_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m5_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m5_pop_outSummary2.txt", sep="\t")
save(m5_pop_out, file ="m5_pop_out2.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M6; R ~ pack size + ran year + ran region + 4wd + 2wd ####

start1<-Sys.time()
# Call JAGS 
m6_group_out <- jags(win.data, inits, params, "M6_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))
end1<-Sys.time()

start1-end1

start2<-Sys.time()
m6_group_out <- autojags(m6_group_out, n.iter = 20000, n.update = 4, n.thin = 3)



# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m6_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m6_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m6_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m6_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m6_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m6_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m6_group_out$BUGSoutput$mean$P
P2[,2] <- m6_group_out$BUGSoutput$sd$P
P2[,3] <- m6_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m6_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m6_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10), 
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m6_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m6_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m6_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m6_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m6_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m6_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m6_group_out)
save(datums, file ="m6_group_out.RData")
end2<-Sys.time()

start2-end2

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m6_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m6_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m6_pop_outSummary.txt", sep="\t")
save(m6_pop_out, file ="m6_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M7; R ~ pack size + ran year + ran region + forcov ####

start1<-Sys.time()
# Call JAGS 
m7_group_out <- jags(win.data, inits, params, "M7_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))
end1<-Sys.time()

start1-end1

start2<-Sys.time()
m7_group_out <- autojags(m7_group_out, n.iter = 20000, n.update = 4, n.thin = 3)



# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m7_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m7_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m7_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m7_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m7_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m7_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m7_group_out$BUGSoutput$mean$P
P2[,2] <- m7_group_out$BUGSoutput$sd$P
P2[,3] <- m7_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m7_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m7_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10), 
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m7_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m7_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m7_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m7_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m7_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m7_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m7_group_out)
save(datums, file ="m7_group_out.RData")
end2<-Sys.time()

start2-end2

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m7_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m7_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m7_pop_outSummary.txt", sep="\t")
save(m7_pop_out, file ="m7_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M8; R ~  ####

start1<-Sys.time()
# Call JAGS 
m8_group_out <- jags(win.data, inits, params, "M8_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))
end1<-Sys.time()

start1-end1

start2<-Sys.time()
m8_group_out <- autojags(m8_group_out, n.iter = 20000, n.update = 4, n.thin = 3)



# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m8_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m8_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m8_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m8_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m8_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m8_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m8_group_out$BUGSoutput$mean$P
P2[,2] <- m8_group_out$BUGSoutput$sd$P
P2[,3] <- m8_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m8_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m8_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,10), 
                    "B0.phi.sd" = rep(NA,10),
                    "B0.colo" = rep(NA,10), 
                    "B0.colo.sd" = rep(NA,10),
                    "b.pc1.colo" = rep(NA,10), 
                    "b.pc1.colo.sd" = rep(NA,10),
                    "b.recPC.colo" = rep(NA,10), 
                    "b.recPC.colo.sd" = rep(NA,10),
                    "b.pc1.phi" = rep(NA,10), 
                    "b.pc1.phi.sd" = rep(NA,10),
                    "b.period.surv" = c(m8_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m8_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m8_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m8_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m8_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m8_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m8_group_out)
save(datums, file ="m8_group_out.RData")
end2<-Sys.time()

start2-end2

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m8_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m8_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m8_pop_outSummary.txt", sep="\t")
save(m8_pop_out, file ="m8_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M9; R ~ pack size + ran year + ran region + dd ####

start1<-Sys.time()
# Call JAGS 
m9_group_out <- jags(win.data, inits, params, "M9_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))
end1<-Sys.time()

start1-end1

start2<-Sys.time()
m9_group_out <- autojags(m9_group_out, n.iter = 20000, n.update = 4, n.thin = 3)



# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m9_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m9_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m9_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m9_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears), 2))
s2[,1] <- m9_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m9_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m9_group_out$BUGSoutput$mean$P
P2[,2] <- m9_group_out$BUGSoutput$sd$P
P2[,3] <- m9_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m9_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m9_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(1,10), 
                    "B0.phi.sd" = rep(1,10),
                    "B0.colo" = rep(1,10), 
                    "B0.colo.sd" = rep(1,10),
                    "b.pc1.colo" = rep(1,10), 
                    "b.pc1.colo.sd" = rep(1,10),
                    "b.recPC.colo" = rep(1,10), 
                    "b.recPC.colo.sd" = rep(1,10),
                    "b.pc1.phi" = rep(1,10), 
                    "b.pc1.phi.sd" = rep(1,10),
                    "b.period.surv" = c(m9_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5,6), 
                    "b.period.surv.sd" = c(m9_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5,6),
                    "eps.surv" = m9_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m9_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m9_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m9_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m9_group_out)
save(datums, file ="m9_group_out.RData")
end2<-Sys.time()

start2-end2

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m9_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m9_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m9_pop_outSummary.txt", sep="\t")
save(m9_pop_out, file ="m9_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### M10; R ~ pack size + ran year + ran region + harvest ####

start1<-Sys.time()
# Call JAGS 
m10_group_out <- jags(win.data, inits, params, "M10_GroupRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                     n.burnin=nb, jags.module = c("glm", "dic"))
end1<-Sys.time()

start1-end1

start2<-Sys.time()
m10_group_out <- autojags(m10_group_out, n.iter = 20000, n.update = 4, n.thin = 3)



# Format output for population level code 

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- m10_group_out$BUGSoutput$mean$n.est
n.est[,2] <- m10_group_out$BUGSoutput$sd$n.est
G.mean2 <- array(NA, dim=c(nyears, 2))
G.mean2[,1] <- m10_group_out$BUGSoutput$mean$G.mean
G.mean2[,2] <- m10_group_out$BUGSoutput$sd$G.mean
s2 <- array(NA, dim=c((nyears-1), 2))
s2[,1] <- m10_group_out$BUGSoutput$mean$annual.s
s2[,2] <- m10_group_out$BUGSoutput$sd$annual.s
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- m10_group_out$BUGSoutput$mean$P
P2[,2] <- m10_group_out$BUGSoutput$sd$P
P2[,3] <- m10_group_out$BUGSoutput$median$P
gamma2 <- array(NA, dim=c(nyears, 2))
gamma2[,1] <- m10_group_out$BUGSoutput$mean$gamma.mean
gamma2[,2] <- m10_group_out$BUGSoutput$sd$gamma.mean
betas <- data.frame("B0.phi" = rep(NA,9), 
                    "B0.phi.sd" = rep(NA,9),
                    "B0.colo" = rep(NA,9), 
                    "B0.colo.sd" = rep(NA,9),
                    "b.pc1.colo" = rep(NA,9), 
                    "b.pc1.colo.sd" = rep(NA,9),
                    "b.recPC.colo" = rep(NA,9), 
                    "b.recPC.colo.sd" = rep(NA,9),
                    "b.pc1.phi" = rep(NA,9), 
                    "b.pc1.phi.sd" = rep(NA,9),
                    "b.period.surv" = c(m10_group_out$BUGSoutput$mean$b.period.surv, 1,2,3,4,5), 
                    "b.period.surv.sd" = c(m10_group_out$BUGSoutput$sd$b.period.surv, 1,2,3,4,5),
                    "eps.surv" = m10_group_out$BUGSoutput$mean$eps.surv, 
                    "eps.surv.sd" = m10_group_out$BUGSoutput$sd$eps.surv)

jag.sum <- m10_group_out$BUGSoutput$summary
write.table(x=jag.sum, file="m10_group_outSummary.txt", sep="\t")
datums <- list(betas, gamma2, G.mean2, s2, P2, n.est, m10_group_out)
save(datums, file ="m10_group_out.RData")
end2<-Sys.time()

start2-end2

# Data
win.data2 <- list("nyears"=nyears, "nsites"=nsites, "area"=sitecovs$AREASAMP,
                  "n.est"=n.est, "gamma2"=gamma2, "P2"=P2, 
                  "G.mean2"=G.mean2, "PC1"=sitecovs$PC1,
                  "recPC"=sitecovs[,27:36], "nperiods"=length(unique(y.surv$period)),
                  "width.interval"=width.interval,"betas"=betas, "T.overlap"=T.overlap)

start1pop<-Sys.time()
# Call JAGS 
m10_pop_out <- jags(win.data2, inits2, params2, "PopRecIPM.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
                   n.burnin=nb, jags.module = c("glm", "dic"))

jag.sum2 <- m10_pop_out$BUGSoutput$summary
write.table(x=jag.sum2, file="m10_pop_outSummary.txt", sep="\t")
save(m10_pop_out, file ="m10_pop_out.RData")

end1pop<-Sys.time()

start1pop-end1pop


#### Diagnostics ####
# Load packages
library(mcmcplots)
library(superdiag)
library(ggmcmc)


# First convert outputs to MCMC object
m10_group_mcmc<-as.mcmc(m10_group_out)
m10_pop_mcmc<-as.mcmc(m10_pop_out)


# Produce a PDF of diagnosic plots for the group level model
pdf("m10_group_diagnosticplots.pdf")

# Print trace and density plots
# plot(m3_group_mcmc, parms=c("B0.gam", "B1.gam", "b.period.surv", "var.surv", "var.reg","P", 
#                             "G.mean", "gamma.mean", "n.est", "annual.s", "var.group"))

# Autocorrelation plots
# autocorr.plot(m3_group_mcmc, parms=c("B0.gam", "B1.gam", "b.period.surv", "var.surv", "var.reg","P", 
#                                      "G.mean", "gamma.mean", "n.est", "annual.s", "var.group"))

# Density plots
denplot(m10_group_mcmc, parms=c("B0.gam", "B1.gam", "b.period.surv", "var.surv", "var.gam","P", 
                               "G.mean", "gamma.mean", "n.est", "annual.s", "var.group", "var.reg", 
                               "b.harv"))

# Traceplots
traplot(m10_group_mcmc, parms=c("B0.gam", "B1.gam", "b.period.surv", "var.surv", "var.gam","P", 
                               "G.mean", "gamma.mean", "n.est", "annual.s", "var.group", "var.reg", 
                               "b.harv"))

dev.off()


# Produce a PDF of diagnosic plots for the pop level model
pdf("m10_pop_diagnosticplots.pdf")

# Print trace and density plots
# plot(m3_pop_mcmc, parms=c("P", "G.mean", "gamma.mean", "N.tot", "annual.s", "N.rec"))

# Autocorrelation plots
# autocorr.plot(m3_pop_mcmc, parms=c("P", "G.mean", "gamma.mean", "N.tot", "annual.s", "N.rec"))

# Density plots
denplot(m10_pop_mcmc, parms=c("P", "G.mean", "gamma.mean", "N.tot", "annual.s", "N.rec"))

# Traceplots
traplot(m10_pop_mcmc, parms=c("P", "G.mean", "gamma.mean", "N.tot", "annual.s", "N.rec"))

dev.off()


# Get quick numerical representations of diagonstics with superdiag
m10_group_diags <- superdiag(m10_group_mcmc, burnin=100)
save(m10_group_diags, "m10_group_diags.RData")

m10_pop_diags <- superdiag(m10_pop_mcmc, burnin=100)
save(m10_pop_diags, file="m10_pop_diags.RData")


# Produce and save a coefficient dot plots for the group level
pdf("m10_group_caterplot.pdf")
caterplot(m10_group_mcmc, parms=c("B0.gam", "B1.gam", "b.period.surv", "var.surv", "var.gam", "var.reg", "b.harv"))
caterplot(m10_group_mcmc, parms=c("P"))
caterplot(m10_group_mcmc, parms=c("G.mean", "gamma.mean"))
caterplot(m10_group_mcmc, parms=c("n.est"))
caterplot(m10_group_mcmc, parms=c("annual.s"))
dev.off()


# Produce and save a coefficient dot plots for the pop level
pdf("m10_pop_caterplot.pdf")
caterplot(m10_pop_mcmc, parms=c("P"))
caterplot(m10_pop_mcmc, parms=c("G.mean", "gamma.mean"))
caterplot(m10_pop_mcmc, parms=c("N.tot", "N.rec"))
caterplot(m10_pop_mcmc, parms=c("annual.s"))
dev.off()


# Produce plots using ggmcmc for diagnostics
m10_group_gg <- ggs(m10_group_mcmc)
ggmcmc(m10_group_gg, file="m10_group_ggmcmc.pdf")

m10_pop_gg <- ggs(m10_pop_mcmc)
ggmcmc(m10_pop_gg, file="m10_pop_ggmcmc.pdf")
