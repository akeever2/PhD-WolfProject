N.surv = 100
round=2
x1 = rbinom(N.surv,size=1,prob=0.5)
x = t(rbind(x1))
censortime = runif(N.surv,0,1)
survtime= rexp(N.surv,rate=exp(-.5*x1))
survtime = round(survtime,digits=round)
event = as.numeric(censortime>survtime)
y = survtime; 
y[event==0] = censortime[event==0]
t=sort(unique(y[event==1]))
t=c(t,max(censortime))
bigt=length(t)-1






##########################################################

# Cox proportional hazards model code with a counting process. 
# Adapted from Luek Winbugs example

##########################################################





library(R2jags)
library(mcmcplots)
################################################################################
#  Specify model in BUGS language

sink("survival.txt")
cat("
    
    data {
    
    # Set up the data
    
    for(i in 1:N.surv) {
    for(j in 1:T) {
    # Create risk set; Yij = 1 if obs.t >= t and start.t <= t
    Y[i,j] <- step(obs.t[i] - t[j] + eps) #* step(t[j] - start.t[i] + eps) 
    
    # Counting process jump; dN = 1 if obs.t is in interval [t[j], t[j+1) 
    # and fail = 1
    dN[i,j] <- Y[i,j] * step(t[j+1] - obs.t[i] - eps) * fail[i] 
    }
    }
    
    
    }    
    
    
    model {
    
    # Now that the data are set up, run the model
    
    
    # Specify priors
    
    # for(j in 1:T){
    # Prior for baseline hazard if using the Poisson trick
    # b0.surv[j] ~ dnorm(0, 0.001)
    # dL0[j] <- exp(b0.surv[j])
    # }
    
    # Prior for beta coefficient
    beta ~ dnorm(0, 0.0001)
    
    
    # Likelihood and Intensity process...What is the intensity process?
    
    for(j in 1:T) {
    for(i in 1:N.surv) {
    
    dN[i,j] ~ dpois(Idt[i,j])
    Idt[i,j] <- Y[i,j] * exp(beta * x[i]) * dL0[j]
    
    # Intensity process if using the Poisson trick - independent log-normal
    # increments which means drop dL0, c, r, and mu from model
    #Idt[i,j] <- Y[i,j] * exp(b0.surv[j] + beta * x[i])
    }
    
    # Prior mean hazard
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c
    
    # Survival function for both groups of fake data; multiply beta coefficient
    # by covariate value (1 or 0 in this example)
    S.treat[j] <- pow(exp(-sum(dL0[1:j])), exp(beta * 1))
    S.no[j] <- pow(exp(-sum(dL0[1:j])), exp(beta * 0))
    
    }
    
    # Confidence in guess for dL0
    c <- 0.001
    
    # Prior guess at failure rate
    r <- 0.1
    
    for(j in 1:T){
    dL0.star[j] <- r * (t[j+1] - t[j])
    }
    
    # Annual survival in this example because observations end after 1 year
    S.annual.treat <- S.treat[T]
    S.annual.no <- S.no[T]
    
    
    }
    ", fill=TRUE)
sink()


win.data <- list(x=x1, obs.t=y, t=t, T=bigt, N.surv=N.surv,
                 fail=event, eps=1E-10)



#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("S.annual.treat", "S.annual.no","S.no", "S.treat", "dL0", "beta") 


# MCMC Settings 
ni <- 5000
nt <- 2
nb <- 2500
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)



