######################################################################################I

#         Simulation for IPM recruitment model at the Population Level ONLY

######################################################################################I


##### Bring in the data from the full recruitment model simulation ####

library(snowfall)
memory.limit(size = 7500000) #just shy of 8 tb limit i think

setwd("C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults")

# Call for appropriate packages

library(R2jags)
library(mcmcplots)

mu.G <- apply(y.group.50, 2, mean)
sd.G <- apply(y.group.50, 2, sd)


#### POM MODEL ####
#  Specify model in BUGS language

sink("POMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    #  psi1 coefficient (occupancy in year 1)
    
    psi1 ~ dunif(0,1)
    
    
    #  Priors for transition probabilities (survival and colonization)
    
    for(k in 1:(nyears-1)){
      phi[k] ~ dunif(0,1)
      colo[k] ~ dunif(0,1)
    }#k
    
    
    # Priors for detection probabilities (only varies by year)
    
    for(k in 1:nyears){
      p11[k] ~ dunif(0,1)
      p10[k] ~ dunif(0,1)
      b[k] ~ dunif(0,1)
    }#k
    
    
    
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
      z[i,1] ~ dbern(psi1)
      
      for(k in 2:nyears){
        muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*colo[k-1]
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
          y.occ[i,j,k] ~ dcat(p[,i,j,k,z[i,k]+1])
        }#k
      }#j
    }#i
    
    
    # Derived parameters
    
    for(i in 1:nsites){
      psi[i,1] <- psi1
      growthr[i,1] <- 1  
    
      for (k in 2:nyears){                                          
        psi[i,k] <- psi[i,k-1]*phi[k-1] + (1-psi[i,k-1])*colo[k-1]
      }#k
    }#i
    
    #  Area occpupied indexed by year
    
    for(k in 1:nyears){
      A[k] <- sum(psi[,k] * area[,k])
    }

    #####################
    
    # 2.2. Territory model 
    
    ####################
    
    # Estimate number of packs from area occupied (A) and territory size (T)
    
    for(k in 1:nyears){
      P[k] <- A[k] / T
    }
    
    
    # Pull in group count data to determine mean group size each year with error
    
    for(k in 1:nyears){
      G[k] ~ dnorm(mu.G[k], 1 / (sd.G[k] * sd.G[k] + 0.000001))T(0,)
    }
    
    
    # Estimate abundance each year based on estimated # of packs (P) and mean group
    # size (G). Then, add on the lone wolves in the population
    
    for(k in 1:nyears){
      n.est[k] <- P[k] * G[k]
    }
    
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nyears"=nyears, 
                 "ngroups"=50, "nobs"=20, 
                 "em.group"=em.group, 
                 "y.group"=y.group.50, "event"=y.surv.20, 
                 "nsites"=nsites, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ,
                 "mu.G" = mu.G, "sd.G" = sd.G)


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
inits <- function(){list(z=zst, colo=runif(nyears-1, 0.05, 0.15), phi=runif(nyears-1, 0.7,0.8),
                         b=runif(nyears, 0.1, 0.3), p10=runif(nyears, 0.05, 0.15),
                         p11=runif(nyears, 0.4, 0.6))}

# Parameters to keep track of and report
params <- c("P", "n.est", "phi", "colo", "psi", "p11", "p10", "b") 


# MCMC Settings 
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "POMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

#print(out, dig=2)




#### format output for population level code ####

n.est <- array(NA, dim=c(nyears, 2))
n.est[,1] <- out$BUGSoutput$mean$n.est
n.est[,2] <- out$BUGSoutput$sd$n.est
P2 <- array(NA, dim=c(nyears, 3))
P2[,1] <- out$BUGSoutput$mean$P
P2[,2] <- out$BUGSoutput$median$P
P2[,3] <- out$BUGSoutput$sd$P

shapes.colo <- shapes.phi <-list()
colo.alpha <- colo.beta <- phi.alpha <- phi.beta <-NA
for(i in 1:(nyears-1)){
  shapes.colo[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$colo[i], out$BUGSoutput$sd$colo[i]^2, "Survival")
  colo.alpha[i] <- shapes.colo[[i]]$alpha
  colo.beta[i] <- shapes.colo[[i]]$beta
  shapes.phi[[i]] <- beta.MoM.fcn(out$BUGSoutput$mean$phi[i], out$BUGSoutput$sd$phi[i]^2, "Survival")
  phi.alpha[i] <- shapes.phi[[i]]$alpha
  phi.beta[i] <- shapes.phi[[i]]$beta
}

betas <- data.frame("colo.alpha" = colo.alpha,
                    "colo.beta" = colo.beta,
                    "phi.alpha" = phi.alpha,
                    "phi.beta" = phi.beta)

datums <-list(betas, n.est, P2, out)
save(datums, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out1_full.RData")


#### population level code ####

sink("PopLevelOnly_Sim.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## Population priors
    
    # Initial population size
    
    # N.ad[1] ~ dnorm(1000, 0.0001)I(0,)
    # N.rec[1] ~ dnorm(615, 0.0001)I(0,)
    N.tot[1] ~ dnorm(1000, 0.0001)I(0,)

    
    ## Bring in data s, G.mean, gamma.mean, P, colo, and phi
    
    for(k in 1:nyears){
      P[k] ~ dnorm(P2[k,1], 1 / (P2[k,3] * P2[k,3]+ 0.0000001))
    }
    


    # Priors for beta coefficients for recruitment
    
      for(k in 1:nyears){
        b0.gam[k] ~ dunif(-10,10)
      }

    # Survival priors
    
    # Intercept and covariates
    
    b0.surv ~ dnorm(0,0.001)

    # Random effect for year
    for(k in 1:nyears){
      eps.surv[k] ~ dnorm(0, tau.surv)
    }
    
    sigma.surv ~ dunif(0,100)
    tau.surv <- pow(sigma.surv, -2)
    var.surv <- pow(sigma.surv, 2)
    
    
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
    
    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    for(k in 1:nyears){
      for(i in 1:nobs){
        event[i,k] ~ dbern(mu.surv[i,k])
        cloglog(mu.surv[i,k]) <- b0.surv + eps.surv[k]
      }
    }
    
    
    # Predicted values
    
    for(k in 1:(nyears-1)){
      cloglog(mu.pred[k]) <- b0.surv + eps.surv[k]
      hazard[k] <- -log(1-mu.pred[k])
    }
    
    # Cumulative hazard and survival 
    
    for(k in 1:(nyears-1)){
      H[k] <- hazard[k]
      annual.s[k] <- exp(-H[k])
    }
    
    
    #####################
    
    # 2.2. Population likelihood 
    
    ####################
    
    # Ecological model/ system process
    
    # First determine colonization and extinction 
    
    for(k in 2:nyears){
      # N.rec[k] ~ dpois(P[k-1] * gamma[k-1])
      # N.ad[k] ~ dbin(annual.s[k-1], N.tot[k-1])
      mu.N[k] <- N.tot[k-1] * annual.s[k-1] + P[k-1] * gamma[k-1] * (1-em.group[k-1])
      N.tot[k] ~ dpois(mu.N[k])
    }

    # for(k in 1:nyears){
    #   N.tot[k] <- round(N.ad[k] + N.rec[k])
    # }
    
    # Linking pack size (P) and mean group size (G.mean) as data (n.est) to abundance (N.tot)
    
    for(k in 1:nyears){
    n.est[k,1] ~ dnorm(N.tot[k], (1 / (n.est[k,2]*n.est[k,2]+0.00001)))
    }
    
    # Recruitment
    for(k in 1:(nyears-1)){
      mu.gamma[k] <- exp(b0.gam[k])
      gamma[k] ~ dpois(mu.gamma[k])
    }#k
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data2 <- list("nyears"=nyears, "em.group"=em.group,
                  "n.est"=n.est, "P2"=P2, 
                  "nobs"=20, "event"=y.surv.20)


#  Initial Values	
inits2 <- function(){list(b0.gam=runif(nyears, 1.4, 1.7), 
                          b0.surv=runif(1, -1, 0),
                          N.tot=rpois(nyears, n.est[,1]),
                          N.ad=rpois(nyears, c(n.est[,1]-gamma.mean*P2[,1])), 
                          N.rec=rpois(nyears, c(gamma.mean*P2[,1])))}


# Parameters to keep track of and report
params2 <- c("P", "annual.s", "N.tot", "gamma", "var.surv",
             "N.rec", "b0.surv", "N.ad", "eps.surv") 

ni <- 500000
nb <- 250000
nt <- 5
nc <- 3

# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevelOnly_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

# out2.update <- autojags(out2, n.iter=100000, n.update=5,  Rhat=1.1)

# out2.update <- update(out2, n.iter = 100000)

#print(out2, digits=2)

full <- data.frame("n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[1:15,3],
                   "n.est.uci"=out2$BUGSoutput$summary[1:15,7],
                   "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[31:44,3],NA),
                   "s.uci"=c(out2$BUGSoutput$summary[31:44,7],NA),
                   "gamma"=c(out2$BUGSoutput$mean$gamma,NA), "gamma.lci"=c(out2$BUGSoutput$summary[62:75,3],NA),
                   "gamma.uci"=c(out2$BUGSoutput$summary[62:75,7],NA),"dataset"=rep("full", nyears))

write.csv(full, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_PopOnly_full.csv")
save(out2, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out2_full.RData")


#### S: 10 surv ####

# Data
win.data2 <- list("nyears"=nyears, "em.group"=em.group,
                  "n.est"=n.est, "P2"=P2, 
                  "nobs"=10, "event"=y.surv.10)

out2 <- jags(win.data2, inits2, params2, "PopLevelOnly_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

# out2.update <- autojags(out2, n.iter=100000, n.update=5,  Rhat=1.1)

# out2.update <- update(out2, n.iter = 100000)

#print(out2, digits=2)

S <- data.frame("n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[1:15,3],
                   "n.est.uci"=out2$BUGSoutput$summary[1:15,7],
                   "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[31:44,3],NA),
                   "s.uci"=c(out2$BUGSoutput$summary[31:44,7],NA),
                   "gamma"=c(out2$BUGSoutput$mean$gamma,NA), "gamma.lci"=c(out2$BUGSoutput$summary[62:75,3],NA),
                   "gamma.uci"=c(out2$BUGSoutput$summary[62:75,7],NA),"dataset"=rep("full", nyears))

write.csv(S, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_PopOnly_S.csv")
save(out2, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out2_S.RData")

rm(out2)


#### T: 20.2 surv ####

# Data
win.data2 <- list("nyears"=nyears, "em.group"=em.group,
                  "n.est"=n.est, "P2"=P2, 
                  "nobs"=20, "event"=y.surv.20.2)

out2 <- jags(win.data2, inits2, params2, "PopLevelOnly_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

# out2.update <- autojags(out2, n.iter=100000, n.update=5,  Rhat=1.1)

# out2.update <- update(out2, n.iter = 100000)

#print(out2, digits=2)

Tt <- data.frame("n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[1:15,3],
                "n.est.uci"=out2$BUGSoutput$summary[1:15,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[31:44,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[31:44,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma,NA), "gamma.lci"=c(out2$BUGSoutput$summary[62:75,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[62:75,7],NA),"dataset"=rep("full", nyears))

write.csv(Tt, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_PopOnly_T.csv")
save(out2, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out2_T.RData")

rm(out2)

#### U: 10.2 surv ####

# Data
win.data2 <- list("nyears"=nyears, "em.group"=em.group,
                  "n.est"=n.est, "P2"=P2, 
                  "nobs"=10, "event"=y.surv.10.2)

out2 <- jags(win.data2, inits2, params2, "PopLevelOnly_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

# out2.update <- autojags(out2, n.iter=100000, n.update=5,  Rhat=1.1)

# out2.update <- update(out2, n.iter = 100000)

#print(out2, digits=2)

U <- data.frame("n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[1:15,3],
                "n.est.uci"=out2$BUGSoutput$summary[1:15,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[31:44,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[31:44,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma,NA), "gamma.lci"=c(out2$BUGSoutput$summary[62:75,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[62:75,7],NA),"dataset"=rep("full", nyears))

write.csv(U, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_PopOnly_U.csv")
save(out2, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out2_U.RData")

rm(out2)

#### V: 20.5 surv ####

# Data
win.data2 <- list("nyears"=nyears, "em.group"=em.group,
                  "n.est"=n.est, "P2"=P2, 
                  "nobs"=20, "event"=y.surv.20.5)

out2 <- jags(win.data2, inits2, params2, "PopLevelOnly_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

# out2.update <- autojags(out2, n.iter=100000, n.update=5,  Rhat=1.1)

# out2.update <- update(out2, n.iter = 100000)

#print(out2, digits=2)

V <- data.frame("n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[1:15,3],
                "n.est.uci"=out2$BUGSoutput$summary[1:15,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[31:44,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[31:44,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma,NA), "gamma.lci"=c(out2$BUGSoutput$summary[62:75,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[62:75,7],NA),"dataset"=rep("full", nyears))

write.csv(V, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_PopOnly_V.csv")
save(out2, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out2_V.RData")

rm(out2)

#### W: 10.5 surv ####
# Data
win.data2 <- list("nyears"=nyears, "em.group"=em.group,
                  "n.est"=n.est, "P2"=P2, 
                  "nobs"=10, "event"=y.surv.10.5)

out2 <- jags(win.data2, inits2, params2, "PopLevelOnly_Sim.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
             n.burnin=nb, jags.module = c("glm", "dic"))

# out2.update <- autojags(out2, n.iter=100000, n.update=5,  Rhat=1.1)

# out2.update <- update(out2, n.iter = 100000)

#print(out2, digits=2)

W <- data.frame("n.est"=out2$BUGSoutput$mean$N.tot, "n.est.lci"=out2$BUGSoutput$summary[1:15,3],
                "n.est.uci"=out2$BUGSoutput$summary[1:15,7],
                "s"=c(out2$BUGSoutput$mean$annual.s,NA), "s.lci"=c(out2$BUGSoutput$summary[31:44,3],NA),
                "s.uci"=c(out2$BUGSoutput$summary[31:44,7],NA),
                "gamma"=c(out2$BUGSoutput$mean$gamma,NA), "gamma.lci"=c(out2$BUGSoutput$summary[62:75,3],NA),
                "gamma.uci"=c(out2$BUGSoutput$summary[62:75,7],NA),"dataset"=rep("full", nyears))

write.csv(W, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/final_PopOnly_W.csv")
save(out2, file ="C:/Users/allison/Documents/Project/Dissertation/Recruitment/Results/SimulationResults/PopOnly_out2_W.RData")

rm(out2)
