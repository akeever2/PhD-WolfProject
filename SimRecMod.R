#########################################################################################

###         MAKE FAKE DATA

#########################################################################################

# Load package for categorical distribution to create false-positive encounter history
library(LaplacesDemon)

data.fn <- function(nyears=4, ngroups=12, s.prob.range=c(0.62, 0.67), P.range=c(90,110), 
                    em.range=c(0.08, 0.18), b0.gam.range=c(1.4, 1.7), sd.proc=4, 
                    sigma.group=1, nobs=20){
  
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
  y.group <- array(dim=c(ngroups, nyears))
  G <- array(dim=c(ngroups, nyears))
  sigma.proc.group <- array(dim=c(nyears))
  y.temp <- array(dim=c(ngroups, nyears))
  G.mean <- array(dim=c(nyears))
  N.tot <- array(dim=c(nyears))
  N.rec <- array(dim=c(nyears))
  N.ad <- array(dim=c(nyears))
  N.alt <- array(dim=c(nyears))
  y.surv <- array(dim=c(nobs,nyears))

  # Recruitment and group size
  # Recruitment
  b0.gam <- runif(n=nyears, min=b0.gam.range[1], max=b0.gam.range[2])
  
  gamma <- exp(b0.gam)
  
  # Random noise for first year group size
  for(i in 1:ngroups){
    G[i,1] <- max(2, rpois(1, 4))
  }
  
  # Process noise for group counts
  # for(k in 1:nyears){
  #   for(r in 1:nregions){
  #     sigma.proc.group[k,r] <- rnorm(1, 0, sd.proc)
  #   }
  # }
  
  # Determine annual group immigration and emigration based on bounds
  
  s.prob <- runif(nyears, min=s.prob.range[1], max=s.prob.range[2])
  em.group <- runif(n=nyears-1, min=em.range[1], max=em.range[2])
  
  # Group size in successive years
  for(i in 1:ngroups){
    for(k in 2:nyears){
      G[i,k] <- ifelse(G[i,k-1]==0, 0, rpois(1, (G[i,k-1]*s.prob[k]*(1-em.group[k-1])+gamma[k]))) #+sigma.proc.group[k-1,r])))
    }
  }
  
  # Group counts
  for(i in 1:ngroups){
    for(k in 1:nyears){
      y.temp[i,k] <- round(rnorm(1, G[i,k], sigma.group))
      y.temp[i,k] <- ifelse(y.temp[i,k] < 0, 0, y.temp[i,k])
      y.group[i,k] <- ifelse(y.temp[i,k] > G[i,k], G[i,k], y.temp[i,k])
    }
  }
  
  # Observation error
  var.group <- sigma.group*sigma.group
  
  # Mean group size
  for(k in 1:nyears){
    G.mean[k] <- mean(G[,k])
  }

  # survival data
  for(k in 1:nyears){
    for(i in 1:nobs){
      y.surv[i,k] <- rbinom(1,1,(1-s.prob[k]))
    }
  }
  
  P <- round(runif(nyears, min=P.range[1], max=P.range[2]))
  
  # Population level
  # Total population in the first year
  N.tot <- P*G.mean


  
  
  return(list(nyears=nyears, ngroups=ngroups, nobs=nobs, P=P, 
              s.prob=s.prob, em.group=em.group, y.surv=y.surv,
              y.group=y.group, b0.gam=b0.gam, G.mean=G.mean,
              gamma.mean=gamma, N.tot=N.tot, var.group=var.group))
}



datum<-data.fn(nyears=10, ngroups=50, s.prob.range=c(.73, .77), P.range=c(90,110), 
               em.range=c(0.08, 0.18), b0.gam.range=c(1.4, 1.7), sd.proc=4, 
               sigma.group=1, nobs=20)

#detach()
attach(datum)
str(datum)





##############################################################################

#             Integrated population model to estimate recruitment

# Allison C. Keever
# Montana Cooperative Wildlife Research Unit
# 2017

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
    
      # Intercept and covariates
      
      for(k in 1:nyears){
        b0.surv[k] ~ dnorm(0,0.001)
      }
      
      
      # Prior for time periods
      
      # for(i in 1:nperiods){
      #   b.period[i] ~ dnorm(0, 0.001)
      # }
    


    ## 1.4 Group priors
    
      # Initial group sizes
    
      for(i in 1:ngroups){
        G[i,1] ~ dpois(5)I(1,)
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
    
      for(k in 1:nyears){
        b0.gam[k] ~ dunif(-10,10)
      }
    
    
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Occupancy likelihood 
    
    ####################


    
    #####################
    
    # 2.2. Territory model 
    
    ####################


    
    #####################
    
    # 2.3. Survival likelihood 
    
    ####################
    
    # Survivial probability based on informative prior, set up to be same for each year
    # and region
    
    # for(i in 1:nobs){
    #   event[i] ~ dbern(mu.surv[i])
    #   cloglog(mu.surv[i]) <- b0.surv[Year[i]] + b.period[Period[i]]
    # }#i

    for(k in 1:nyears){
      for(i in 1:nobs){
        event[i,k] ~ dbern(mu.surv[i,k])
        cloglog(mu.surv[i,k]) <- b0.surv[k]
      }
    }
    
    
    # Predicted values
    
    # Baseline hazard
    # for(k in 1:nyears){
    #   # for(p in 1:nperiods){
    #     cloglog(mu.pred[p,k]) <- b0.surv[k] + b.period[p]
    #     hazard[p,k] <- -log(1-mu.pred[p,k])
    #   # }#p
    # }#k

    for(k in 1:nyears){
      cloglog(mu.pred[k]) <- b0.surv[k]
      hazard[k] <- -log(1-mu.pred[k])
    }
    
    # Cumulative hazard and survival 
    
    # for(k in 1:nyears){
    #   base.H[1,k] <- hazard[1,k] * width.interval[1]
    #   for(p in 2:nperiods){
    #     base.H[p,k] <- base.H[p-1, k] + hazard[p,k] * width.interval[p]
    #   }#p
    # }#k
    # 
    # for(k in 1:nyears){
    #   for(p in 1:nperiods){
    #     base.s[p,k] <- exp(-base.H[p,k])
    #     annual.s[k] <- base.s[length(p), k]
    #   }#p
    # }#k

    for(k in 1:nyears){
      H[k] <- hazard[k]
      annual.s[k] <- exp(-H[k])
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
        g.mu[i,k] <- G[i,k-1] * annual.s[k-1] * (1 - em.group[k-1]) + gamma[i,k-1] # + sigma.proc.group[k-1])
        G[i,k] ~ dnorm(g.mu[i,k], 1/(g.mu[i,k]+.01))T(0,70)
      }#k
    }#i
    
    # Observation proccess
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        y.group[i,k] ~ dnorm(G[i,k], tauy.group)T(0,)
      }#k
    }#i
    
    # Derived parameters
    
    for(k in 1:nyears){
      G.mean[k] <- mean(G[,k])
      gamma.mean[k] <- mean(gamma[,k])
      n.est[k] <- P[k] * G.mean[k]
    }#k
    
    
    
    
    #####################
    
    # 2.5. Recruitment model 
    
    ####################
    
    
    # Generalized linear model with log link function for recruitment
    
    for(i in 1:ngroups){
      for(k in 1:nyears){
        log(mu.gamma[i,k]) <- b0.gam[k]
        gamma[i,k] ~ dpois(mu.gamma[i,k])
      }#k
    }#i
    
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nyears"=nyears, "P"=P, 
                 "ngroups"=15, "nobs"=5, 
                 "em.group"=em.group, "y.group"=y.group.15, "event"=y.surv.5)



inits <- function(){list(b0.gam=runif(nyears, 1.4, 1.7))} #, b0.surv=runif(nyears, 0.71, 0.74))}


# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", "annual.s") 


# MCMC Settings 
ni <- 200000
nt <- 2
nb <- 50000
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)

truth <- data.frame("G.mean"=G.mean,
                    "n.est"=N.tot, 
                    "s"=s.prob,
                    "gamma"=gamma.mean, 
                    "dataset"=rep("Truth", nyears), "year"=c(1:10))


full <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                  "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                  "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                  "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                  "dataset"=rep("full", nyears))

A <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                   "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                   "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                   "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                   "dataset"=rep("A", nyears), "year"=c(1:10))

B <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                 "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                 "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                 "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                 "dataset"=rep("B", nyears), "year"=c(1:10))

C <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                 "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                 "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                 "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                 "dataset"=rep("C", nyears), "year"=c(1:10))

D <- data.frame("G.mean"=out$BUGSoutput$mean$G.mean, "G.mean.sd"=out$BUGSoutput$sd$G.mean, 
                 "n.est"=out$BUGSoutput$mean$n.est, "n.est.sd"=out$BUGSoutput$sd$n.est, 
                 "s"=out$BUGSoutput$mean$annual.s, "s.sd"=out$BUGSoutput$sd$annual.s, 
                 "gamma"=out$BUGSoutput$mean$gamma.mean, "gamma.sd"=out$BUGSoutput$sd$gamma.mean,
                 "dataset"=rep("D", nyears), "year"=c(1:10))


output <- bind_rows(truth, full, A, B, C, D)
output$year <- rep(1:10, 6)
output[1:10, 7:10] <- 0
# output$G.mean.true <- rep(G.mean, 5)
# output$n.est.true <- rep(N.tot, 5)
# output$gamma.true <- rep(gamma.mean, 5)
# output$s.true <- rep(s.prob, 5)


y.group.25 <- sample_n(as.data.frame(y.group), 25)
y.group.15 <- sample_n(as.data.frame(y.group), 15)

y.surv.10 <- sample_n(as.data.frame(y.surv), 10)
y.surv.5 <- sample_n(as.data.frame(y.surv), 5)


ggplot(data=output2, aes(x=year, y=gamma, colour=dataset))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=(gamma-gamma.sd), ymax=(gamma+gamma.sd)), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)+
  scale_x_continuous(name="Year", breaks=c(1:10))+
  scale_y_continuous(name="Recruitment", breaks=c(1,2,3,4,5,6,7,8))+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14))


ggplot(data=output2, aes(x=year, y=s, colour=dataset))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=(s-s.sd), ymax=(s+s.sd)), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)+
  scale_x_continuous(name="Year", breaks=c(1:10))+
  scale_y_continuous(name="Survival probability", breaks=c(0.2, 0.4, 0.6, 0.8, 1))+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14))


ggplot(data=output2, aes(x=year, y=n.est, colour=dataset))+
  geom_point(position=position_dodge(width=0.4), size=3)+
  geom_errorbar(aes(ymin=(n.est-n.est.sd), ymax=(n.est+n.est.sd)), 
                width=0.2, position=position_dodge(width=0.4), 
                size=1)+
  scale_y_continuous(name="Abundance")+
  theme_bw()+
  theme(axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16))

write.csv(output, "simoutput.csv")
output2 <- read.csv("simoutput.csv")
