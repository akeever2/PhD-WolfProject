
#########################################################################################

###         MAKE FAKE DATA

#########################################################################################

# Load package for categorical distribution to create false-positive encounter history
library(LaplacesDemon)

data.fn <- function(nsites=20, nyears=4, nregions=3, ngroups=12, s.prob=0.8, Tr=600, 
                    imm.range=c(0,0.1), em.range=c(0.1, 0.2), b0.gam=1.5, sd.proc=4, 
                    sigma.group=1, psi1=0.4, p11.range=c(0.4,0.6), p10.range=c(0,0.1), 
                    b.range=c(0.2,0.4), patch.surv.range=c(0.7,0.9), patch.col.range=c(0,0.1),
                    noccs=3){
  
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
  y.occ <- array(dim=c(nsites, noccs, nyears, nregions))
  z <- muZ <- array(dim=c(nsites, nyears, nregions))
  p.0 <- p.1 <- array(dim=c(nyears, 3))
  prob <- array(dim=c(nsites, 3, nyears, nregions))
  
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
  # for(k in 1:nyears){
  #   for(r in 1:nregions){
  #     sigma.proc.group[k,r] <- rnorm(1, 0, sd.proc)
  #   }
  # }
  
  # Determine annual group immigration and emigration based on bounds
  
  imm.group <- runif(n=nyears-1, min=imm.range[1], max=imm.range[2])
  em.group <- runif(n=nyears-1, min=em.range[1], max=em.range[2])
  
  # Group size in successive years
  for(i in 1:ngroups){
    for(k in 2:nyears){
      for(r in 1:nregions){
        G[i,k,r] <- ifelse(G[i,k-1,r]==0, 0, rpois(1, (G[i,k-1,r]*s.prob*(1+imm.group[k-1]-em.group[k-1])+gamma))) #+sigma.proc.group[k-1,r])))
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
  # Initial occupancy, psi
  for(i in 1:nsites){
    for(r in 1:nregions){
      psi[i,1,r] <- psi1
    }
  }
  
  # Determine detection and patch colinization and survivial
  p11 <- runif(n=nyears, min=p11.range[1], max=p11.range[2])
  p10 <- runif(n=nyears, min=p10.range[1], max=p10.range[2])
  b <- runif(n=nyears, min=b.range[1], max=b.range[2])
  patch.col <- runif(n=nyears-1, min=patch.col.range[1], max=patch.col.range[2])
  patch.surv <- runif(n=nyears-1, min=patch.surv.range[1], max=patch.surv.range[2])
  
  # Generate latent states of occurance
  # Year 1
  for(r in 1:nregions){
    z[,1,r] <- rbinom(nsites, 1, psi[1])
  }
  
  # Later years
  for(i in 1:nsites){
    for(r in 1:nregions){
      for(k in 2:nyears){
        muZ[k] <- z[i,k-1,r]*patch.surv[k-1]+(1-z[i,k-1,r])*patch.col[k-1]
        z[i,k,r] <- rbinom(1,1,muZ[k])
      }
    }
  }
  
  # Generate encounter histories using detection and true state of occurance
  # Set up detection vectors to be used in categorical distribution for when z=0 and z=1
  # p.0 is z=0, p.1 is z=1. The order is no detection (y=0), uncertain det (y=1) and 
  # certain det (y=2) for the probabilites
  for(k in 1:nyears){
    p.0[k,] <- c(1-p10[k], p10[k], 0)
    p.1[k,] <- c(1-p11[k], p11[k]*(1-b[k]), b[k]*p11[k])
  }
  
  # Create occupancy data
  for(i in 1:nsites){
    for(r in 1:nregions){
      for(k in 1:nyears){
        if(z[i,k,r]==0){
          prob[i,,k,r] <- p.0[k,]
        } else {
          prob[i,,k,r] <- p.1[k,]
        }
        for(j in 1:noccs){
          y.occ[i,j,k,r] <- rcat(1, prob[i,,k,r])
        }
      }
    }
  }
  
  # Compute annual occupancy so I can get area occupied
  for(i in 1:nsites){
    for(r in 1:nregions){
      for(k in 2:nyears){
        psi[i,k,r] <- psi[i,k-1,r]*patch.surv[k-1]+(1-psi[i,k-1,r])*patch.col[k-1] 
      }
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
      N.ad[k,r] <- rbinom(1, round(N.tot[k-1,r]*(1+patch.col[k-1]-(1-patch.surv[k-1]))), s.prob)
      N.tot[k,r] <- N.rec[k,r]+N.ad[k,r]
    }
  }
  
  # Alt pop size
  for(k in 1:nyears){
    for(r in 1:nregions){
      N.alt[k,r] <- P[k,r]*G.mean[k,r]
    }
  }
  
  
  
  
  
  return(list(nsites=nsites, nyears=nyears, nregions=nregions, ngroups=ngroups, noccs=noccs,
              s.prob=s.prob, Tr=Tr, imm.group=imm.group, em.group=em.group, patch.col=patch.col, 
              patch.surv=patch.surv, y.group=y.group, b0.gam=b0.gam, P=P, G.mean=G.mean, A=A,
              gamma.mean=gamma, N.tot=N.tot, psi=psi, area=area, var.group=var.group, 
              N.alt=N.alt, psi=psi, y.occ=y.occ, p10=p10, p11=p11, b=b, z=z))
}



datum<-data.fn(nsites=50, nyears=10, nregions=4, ngroups=30, s.prob=0.8, Tr=600, 
               imm.range=c(0,0.1), em.range=c(0.1, 0.2), b0.gam=1.5, sd.proc=4, 
               sigma.group=1, psi1=0.4, p11.range=c(0.4,0.6), p10.range=c(0,0.1), 
               b.range=c(0.1,0.2), patch.surv.range=c(0.7,0.9), patch.col.range=c(0.05,0.1),
               noccs=10)

#detach()
attach(datum)
str(datum)





sink("IPMmodel.txt")
cat("
    model {
    
    ############################################################
    
    #             1. Priors
    
    ############################################################
    
    ## 1.1 Occupancy priors
    
    #  psi1 coefficient (occupancy in year 1)
    
    #B0.psi1 ~ dnorm(0,0.001)	  
    psi1 ~ dunif(0,1)
    
    
    #  Priors for transition probabilities (survival and colonization)
    
    for(k in 1:(nyears-1)){
      #B0.phi[k] ~ dnorm(0,0.001)
      #B0.gamma[k] ~ dnorm(0,0.001)
      phi[k] ~ dunif(0,1)
      colo[k] ~ dunif(0,1)
    }#k
    
    
    # Priors for detection probabilities (only varies by year)
    
    for(k in 1:nyears){
      #B0.p11[k] ~ dnorm(0,0.001)
      #B0.p10[k] ~ dnorm(0,0.001)
      #B0.b[k] ~ dnorm(0,0.001)
      p11[k] ~ dunif(0,1)
      p10[k] ~ dunif(0,1)
      b[k] ~ dunif(0,1)
    }#k
    
    

    
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
      for(r in 1:nregions){
        # logit.psi1[i] <- B0.psi1       
        # logit(psi1[i]) <- logit.psi1[i]                                     
        # z[i,1] ~ dbern(psi1[i])
        z[i,1,r] ~ dbern(psi1)
        
        # for(k in 1:nyears-1){                                                    
        #   # logit.phi[i,k] <- B0.phi[k]                                       
        #   # logit.gamma[i,k] <- B0.gamma[k]
        #   # logit(phi[i,k]) <- logit.phi[i,k]                         
        #   # logit(gamma[i,k]) <- logit.gamma[i,k]
        # }#k
    
        for(k in 2:nyears){
          muZ[i,k,r] <- z[i,k-1,r]*phi[k-1] + (1-z[i,k-1,r])*colo[k-1]
          z[i,k,r] ~ dbern(muZ[i,k,r])
        }#k
      }#r
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
          for(r in 1:nregions){
            # multilogit.p11[i,j,k] <- B0.p11[k]
            # logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
            # multilogit.p10[i,j,k] <- B0.p10[k]
            # logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
            # multilogit.b[i,j,k] <- B0.b[k]
            # logit(b[i,j,k]) <- multilogit.b[i,j,k]
            
            y.occ[i,j,k,r] ~ dcat(p[,i,j,k,z[i,k,r]+1])
          }#r
        }#k
      }#j
    }#i
    
    
    # Derived parameters
    
    for(i in 1:nsites){
      for(r in 1:nregions){
        psi[i,1,r] <- psi1
        growthr[i,1,r] <- 1  
        
        for (k in 2:nyears){                                          
          psi[i,k,r] <- psi[i,k-1,r]*phi[k-1] + (1-psi[i,k-1,r])*colo[k-1]
          growthr[i,k,r] <- psi[i,k,r]/psi[i,k-1,r]
          turnover[i,k-1,r] <- (1 - psi[i,k-1,r]) * colo[k-1]/psi[i,k,r]
        }#k
      }#r
    }#i
    
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
    #   Territory size is fixed value
    
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
    
    
    ############################################################
    
    #             3. Bugs requirements
    
    ############################################################
    
    
    }", fill=TRUE)
sink()


# Data
win.data <- list("nsites"=nsites, "nyears"=nyears, "nregions"=nregions, "area"=area,
                 "noccs"=noccs, "T"=Tr,  "y.occ"=y.occ)


#  Initial Values	
zst <- apply(y.occ,c(1,3,4), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 1
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(z=zst, colo=runif(nyears-1, 0.05, 0.15), phi=runif(nyears-1, 0.7,0.8),
                         b=runif(nyears, 0.1, 0.3), p10=runif(nyears, 0.05, 0.15),
                         p11=runif(nyears, 0.4, 0.6))}


# Parameters to keep track of and report
params <- c("P", "phi", "colo", "psi", 
            "p11", "p10", "b") 


# MCMC Settings 
ni <- 8000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "IPMmodel.txt", n.chains=nc, n.thin=nt, n.iter=ni, 
            n.burnin=nb, jags.module = c("glm", "dic"))

print(out, dig=2)

mcmcplot(out)



b.m<-as.data.frame(out$BUGSoutput$mean$b)
rownames(b.m)<-c(1:nyears)
colnames(b.m)<-"b.model"
b.m<-b.m %>% mutate(Year=c(1:nyears))

p11.m<-as.data.frame(out$BUGSoutput$mean$p11)
rownames(p11.m)<-c(1:nyears)
colnames(p11.m)<-"p11.model"
p11.m<-p11.m %>% mutate(Year=c(1:nyears))

p10.m<-as.data.frame(out$BUGSoutput$mean$p10)
rownames(p10.m)<-c(1:nyears)
colnames(p10.m)<-"p10.model"
p10.m<-p10.m %>% mutate(Year=c(1:nyears))

phi.m<-as.data.frame(out$BUGSoutput$mean$phi)
rownames(phi.m)<-c(1:(nyears-1))
colnames(phi.m)<-"phi.model"
phi.m<-phi.m %>% mutate(Year=c(1:(nyears-1)))

colo.m<-as.data.frame(out$BUGSoutput$mean$colo)
rownames(colo.m)<-c(1:(nyears-1))
colnames(colo.m)<-"colo.model"
colo.m<-colo.m %>% mutate(Year=c(1:(nyears-1)))

output<- b.m
output<- left_join(output, p11.m, by=c("Year"))
output<- left_join(output, p10.m, by=c("Year"))
output<- left_join(output, phi.m, by=c("Year"))
output<- left_join(output, colo.m, by=c("Year"))

phi.t<-as.data.frame(patch.surv)
rownames(phi.t)<-c(1:(nyears-1))
colnames(phi.t)<-"phi.true"
phi.t<-phi.t %>% mutate(Year=c(1:(nyears-1)))

colo.t<-as.data.frame(patch.col)
rownames(colo.t)<-c(1:(nyears-1))
colnames(colo.t)<-"colo.true"
colo.t<-colo.t %>% mutate(Year=c(1:(nyears-1)))

output <- output %>% mutate(b.true=b, p11.true=p11, p10.true=p10)
output <- left_join(output, colo.t, by=c("Year"))
output <- left_join(output, phi.t, by=c("Year"))




ggplot(data=output, aes(x=Year, y=b.model, color="blue"))+
  geom_point(data=output, aes(x=Year, y=b.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=b.model, color="blue"))

ggplot(data=output, aes(x=Year, y=p11.model, color="blue"))+
  geom_point(data=output, aes(x=Year, y=p11.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=p11.model, color="blue"))


ggplot(data=output, aes(x=Year, y=p10.model, color="blue"))+
  geom_point(data=output, aes(x=Year, y=p10.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=p10.model, color="blue"))

ggplot(data=output, aes(x=Year, y=phi.model, color="blue"))+
  geom_point(data=output, aes(x=Year, y=phi.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=phi.model, color="blue"))

ggplot(data=output, aes(x=Year, y=colo.model, color="blue"))+
  geom_point(data=output, aes(x=Year, y=colo.true, color="red"))+
  geom_point(data=output, aes(x=Year, y=colo.model, color="blue"))





