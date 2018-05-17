    #  Multinomial example
    #  Josh Nowak
    #  10/2017
################################################################################
    require(R2jags)
    require(dplyr)
################################################################################
    #  Simulate data
    
    n_observers <- 6
    
    #  Make some observer data, each observer views 20 papers
    dat <- tibble::tibble(
      observer = rep(1:n_observers, each = 100),
      y = as.integer(NA)
    )
    
    #  The mean for each of 3 response categories are defined as...must sum to 1
    p <- numeric(length = 3)
    p[1] <- 0.5
    p[2] <- 0.2
    p[3] <- 1 - (p[1] + p[2])
    
    #  Convert from link to real scale
    real_p <- numeric(length = 3)
    real_p[1] <- log(p[1]/p[3])
    real_p[2] <- log(p[2]/p[3])
    real_p[3] <- log(p[3]/p[3])    
    
    stopifnot(real_p[3] == 0)

    #  Add random effect of observer using a multivariate normal b/c as one
    #   value goes up something else must go down
    ind_p <- MASS::mvrnorm(n_observers, real_p, diag(1, 3))
    
    #  Check that all rows sum to one after the link function is applied
    #  Function for link
    mlogit_link <- function(x){
      #  Takes numeric vector
      exp_x <- exp(x)
      sapply(exp_x, function(z) z/sum(exp_x))
      
    }
    
    #  Apply link function for each observer (rows)
    p_obs <- apply(ind_p, 1, mlogit_link)
    
    #  Check that all columns sum to 1, rounding to 9 decimal places
    stopifnot(all(round(colSums(p_obs), 9) == 1))

    #  Create observations
    for(i in 1:nrow(dat)){
      dat$y[i] <- sample(1:length(p), 1, prob = p_obs[,dat$observer[i]])
    }
    
    #  Relative to our probabilities how did we do?
    dat %>%
      count(y) %>%
      mutate(
        prob = round(n/sum(n), 2)
      )
      
    #  Pretty good
################################################################################
    sink("multinom.txt")
    cat("
     model{
      #  Priors
      
      #  Intercept
      rp[1] ~ dnorm(0, 0.001)
      rp[2] ~ dnorm(0, 0.001)
      rp[3] <- 0
      
      #  Random effects for individual observers
      mu_eps[1] <- 0
      mu_eps[2] <- 0
      
      #  Diagonal matrix for Sigma prior
      V[1, 1] <- 1
      V[2, 2] <- 1
      V[1, 2] <- 0
      V[2, 1] <- 0      
      
      #  Prior on precision matrix, with 3 degrees of freedom
      omega[1:2, 1:2] ~ dwish(V[,], 3)
      
      #  Sigma matrix
      sigma2[1:2, 1:2] <- inverse(omega[,])
      
      #  Random effects
      for(i in 1:n_observers){
        eps[i, 1:2] ~ dmnorm(mu_eps[], omega[1:2, 1:2])
      }

      #  Linear predictor
      for(i in 1:n_observers){
        lp[i, 1] <- exp(rp[1] + eps[i, 1])
        lp[i, 2] <- exp(rp[2] + eps[i,2])
        lp[i, 3] <- exp(0)
        sum_lp[i] <- lp[i,1] + lp[i,2] + lp[i,3]
        p[i, 1] <- lp[i,1]/sum_lp[i]
        p[i, 2] <- lp[i,2]/sum_lp[i]
        p[i, 3] <- lp[i,3]/sum_lp[i]
      }
      
      #  Likelihood
      for(i in 1:nobs){
        y[i] ~ dcat(p[observer[i],])
      }
      
      #  Derived
      for(i in 1:3){
        global_lp[i] <- exp(rp[i])
        global_p[i] <- global_lp[i]/global_sum
      }
      global_sum <- sum(global_lp[])
      
      
    }    
    ", fill = T)
    sink()
################################################################################
    #  Jags data
    jagsd <- list(
      y = dat$y,
      nobs = nrow(dat),
      n_observers = n_distinct(dat$observer),
      observer = dat$observer
    )
    
    #  Initial values
    inits <- function(){
      list(
        rp = c(rnorm(length(p)-1), NA)
      )
    }
    
    #  Params
    parms <- c("rp", "p", "global_p")
    
    #  Call
    out <- jags(
      jagsd,
      inits,
      parms,
      
      "multinom.txt",
      n.iter = 3000,
      n.burn = 500
    )
    
    #  Compare for observer #1
    tibble::tibble(
      Truth = p_obs[,1],
      Est = out$BUGS$mean$p[1,],
      LCL = apply(out$BUGS$sims.list$p[,1,], 2, quantile, 0.025),
      UCL = apply(out$BUGS$sims.list$p[,1,], 2, quantile, 0.975)
    ) %>%
    mutate_all(round, digits = 2)
    
    #  Compare for global
    tibble::tibble(
      Truth = p,
      Est = out$BUGS$mean$global_p,
      LCL = apply(out$BUGS$sims.list$global_p, 2, quantile, 0.025),
      UCL = apply(out$BUGS$sims.list$global_p, 2, quantile, 0.975)
    ) %>%
    mutate_all(round, digits = 2)    
################################################################################
    #  End