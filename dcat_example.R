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
      observer = rep(1:n_observers, each = 20)
    )
    
    #  The mean for each of 3 response categories are defined as...must sum to 1
    p <- numeric(length = 3)
    p[1] <- 0.5
    p[2] <- 0.2
    p[3] <- 1 - (p[1] + p[2])
    
    #  Convert from link to real scale to for fun
    real_p <- numeric(length = 3)
    real_p[1] <- log(p[1]/p[3])
    real_p[2] <- log(p[2]/p[3])
    real_p[3] <- log(p[3]/p[3])    
    
    stopifnot(real_p[3] == 0)

    #  For fun apply multinomial logit to get back probabilities
    sum_rp <- sum(exp(real_p))
    exp(real_p[1])/sum_rp
    exp(real_p[2])/sum_rp
    #  Magic happens here!!!
    exp(real_p[3])/sum_rp
    
    #  Hopefully we get the link function now and understand why setting one
    #   category to 0 doesn't matter.
    
    #  Create observations
    dat$y <- sample(1:length(p), nrow(dat), replace = T, prob = p)
    
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
      rp[1] ~ dnorm(0, 0.001)
      rp[2] ~ dnorm(0, 0.001)
      
      #  Linear predictor - verbose
      lp[1] <- exp(rp[1])
      lp[2] <- exp(rp[2])
      lp[3] <- exp(0)
      
      #  Sum of all values
      sum_lp <- sum(lp[])
      
      #  Probability...'prop of total made up of each'
      p[1] <- lp[1]/sum_lp
      p[2] <- lp[2]/sum_lp      
      p[3] <- lp[3]/sum_lp
      
      #  Likelihood
      for(i in 1:nobs){
        y[i] ~ dcat(p[])
      }
    }    
    ", fill = T)
    sink()
################################################################################
    #  Jags data
    jagsd <- list(
      y = dat$y,
      ncat = length(p),
      nobs = nrow(dat)
    )
    
    #  Initial values
    inits <- function(){
      list(
        rp = rnorm(length(p)-1)
      )
    }
    
    #  Params
    parms <- c("rp", "p")
    
    #  Call
    out <- jags(
      jagsd,
      inits,
      parms,
      
      "multinom.txt",
      n.iter = 3000,
      n.burn = 500
    )
    
    #  Compare
    tibble::tibble(
      Truth = p,
      Est = out$BUGS$mean$p,
      LCL = apply(out$BUGS$sims.list$p, 2, quantile, 0.025),
      UCL = apply(out$BUGS$sims.list$p, 2, quantile, 0.975)
    ) %>%
    mutate_all(round, digits = 2)
    
    #  Prove that at each iteration the sum of p is 1
    #  Use unique b/c machine rounding and all...
    unique(apply(out$BUGS$sims.list$p, 1, sum))
    
    
    
    
    
    
    
    
    