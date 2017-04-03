
#########################################################################################

###         MAKE FAKE DATA

#########################################################################################



data.fn <- function(nsites=20, nyears=4, nregions=3, ngroups=12, s.prob=0.8, Tr=600, 
                    imm.group=0.1, em.group=0.15, imm.p=0, em.p=0, b0.gam=1.5, sd.proc=4, 
                    sigma.group=1){
  
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
  for(k in 1:nyears){
    for(r in 1:nregions){
      sigma.proc.group[k,r] <- rnorm(1, 0, sd.proc)
    }
  }
  
  # Group size in successive years
  for(i in 1:ngroups){
    for(k in 2:nyears){
      for(r in 1:nregions){
        G[i,k,r] <- ifelse(G[i,k-1,r]==0, 0, rpois(1, (G[i,k-1,r]*s.prob*(1+imm.group-em.group)+gamma +sigma.proc.group[k-1,r])))
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
  # Psi
  for(i in 1:nsites){
    for(r in 1:nregions){
      psi[i,,r] <- runif(1, 0.35, 0.85)
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
      N.ad[k,r] <- rbinom(1, round(N.tot[k-1,r]*(1+imm.p-em.p)), s.prob)
      N.tot[k,r] <- N.rec[k,r]+N.ad[k,r]
    }
  }
  
  # Alt pop size
  for(k in 1:nyears){
    for(r in 1:nregions){
      N.alt[k,r] <- P[k,r]*G.mean[k,r]
    }
  }
  
  
  
  
  
  return(list(nsites=nsites, nyears=nyears, nregions=nregions, ngroups=ngroups, 
              s.prob=s.prob, Tr=Tr, imm.group=imm.group, em.group=em.group, imm.p=imm.p, 
              em.p=em.p, y.group=y.group, b0.gam=b0.gam, P=P, G.mean=G.mean, A=A,
              gamma.mean=gamma, N.tot=N.tot, psi=psi, area=area, var.group=var.group, 
              N.alt=N.alt))
}



datum<-data.fn(nsites=50, nyears=10, nregions=5, ngroups=50, s.prob=0.8, Tr=600, 
               imm.group=0.1, em.group=0.15, imm.p=0, em.p=0, b0.gam=1.5, sd.proc=1, 
               sigma.group=1)

detach()
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
    
    b0.gam ~ dunif(-10,10)
    
    ############################################################
    
    #             2. Likelihoods
    
    ############################################################
    
    
    #####################
    
    # 2.1. Occupancy likelihood 
    
    ####################
    
    # Simplified model to get area occupied based on fixed values of psi (occupancy probability)
    
    
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
    #   Simple model of the mean based on results from Rich et al. 2012
    
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
          g.mu[i,k,r] <- G[i,k-1,r] * s[k-1,r] * (1 + imm.group - em.group) + 
            gamma[i,k-1,r]# + sigma.proc.group[k-1,r])
          G[i,k,r] ~ dnorm(g.mu[i,k,r], 1/g.mu[i,k,r])
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
          log(mu.gamma[i,k,r]) <- b0.gam
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
                 "ngroups"=ngroups, "psi"=psi, "T"=Tr, "s.prob"=s.prob, "imm.group"=imm.group, 
                 "em.group"=em.group, "y.group"=y.group)


#  Initial Values	
inits <- function(){list(b0.gam=runif(1,-1,1), sd.proc=runif(1,0,10), sigma.group=runif(1,0,10))}


# Parameters to keep track of and report
params <- c("P","var.group", "G.mean", "gamma.mean", "n.est", "s") 


# MCMC Settings 
ni <- 10000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb)

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


    ## Bring in data s, G.mean, gamma.mean
    
    for(k in 1:nyears){
      for(r in 1:nregions){
        #s[k,r] ~ dnorm(s2[k,r,1], 1 / (s2[k,r,2] * s2[k,r,2]))
        s[k,r] <- 0.8
        #P[k,r] ~ dnorm(P2[k,r,1], 1 / (P2[k,r,2] * P2[k,r,2]))
        P[k,r] <- 30
        gamma.mean[k,r] ~ dnorm(gamma2[k,r,1], 1 / (gamma2[k,r,2] * gamma2[k,r,2]))
      }
    }
      
    ############################################################
    
    #             2. Likelihood
    
    ############################################################
      
    
    # Ecological model/ system process
    
    for(k in 2:nyears){
      for(r in 1:nregions){
        N.rec[k,r] ~ dpois(P[k-1,r] * gamma.mean[k-1,r])
        N.ad[k,r] ~ dbin(s[k-1,r], round(N.tot[k-1,r] * (1+imm.p-em.p)))
        N.tot[k,r] <- N.ad[k,r] + N.rec[k,r]
        #n.mu[k,r] <- N.tot[k-1,r]*(1+imm.p-em.p)*s[k-1,r]+P[k-1,r]*gamma.mean[k-1,r]
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
win.data2 <- list("nyears"=nyears, "nregions"=nregions, "em.p"=em.p,
                 "n.est"=n.est, "s2"=s2, "gamma2"=gamma2, "P2"=P2, "imm.p"=imm.p)


#  Initial Values	
inits2 <- function(){list()}


# Parameters to keep track of and report
params2 <- c("P", "s", "N.tot", "gamma.mean") 


# Call JAGS 
out2 <- jags(win.data2, inits2, params2, "PopLevel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb)

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






