# Example from http://data.princeton.edu/pop509/recid3.html

# Let's have another look at the recidivism data. We will split duration into single years with an open-ended category at 5+ and fit a piecewise exponential model with the same covariates as Wooldridge.
# 
# We will then treat the data as discrete, assuming that all we know is that recidivism occured somewhere in the year. We will fit a binary data model with a logit link, which corresponds to the discrete time model, and using a complementary-log-log link, which corresponds to a grouped continuous time model.

library(foreign)
library(dplyr)

recid <- read.dta("http://www.stata.com/data/jwooldridge/eacsap/recid.dta")

library(survival)

recid$fail <- 1 - recid$cens

recidx <- survSplit(recid, cut = seq(12, 60, 12), start = "t0", end = "durat", event = "fail", episode = "interval")

labels <- paste("(",seq(0,60,12),",",c(seq(12,60,12),81), "]",sep="")

recidx <- mutate(recidx, exposure = durat - t0, interval = factor(interval + 1, labels = labels))

mf <-  fail ~  interval + workprg + priors + tserved + felon + alcohol + drugs + black + married + educ + age

pwe <- glm(mf, offset = log(exposure), data = recidx, family = poisson)

coef(summary(pwe))


# Finally we use a complementary log-log link

cloglog <- glm(mf, data = recidx, family = binomial(link = cloglog))

coef(summary(cloglog))







# Example from http://data.princeton.edu/wws509/R/recidivism.html, German Rodriguez 
# at Princeton University GLM class. 


library(foreign)
recid <- read.dta("http://www.stata.com/data/jwooldridge/eacsap/recid.dta")
nrow(recid)



recid <- mutate(recid, fail = 1 - cens, id = row_number())
filter(recid, id == 9) %>% select(id, durat, fail)


# To create pseudo-observations for survival analysis we will use the survSplit() 
# function in the `survival` package. We will split the data into single-year 
# intervals of duration from 0-12 to 48-60 with an open-ended category 60+.

library(survival)
breaks <- seq(12, 60, by=12)
recidx <- survSplit(Surv(durat, fail) ~ ., data = recid, cut = breaks, 
                    episode = "interval", start = "start")


# The function codes the interval variable using integer codes, and we turn that 
# into a factor for convenience We calculate exposure time for each episode as the 
# difference between duration at the start and end.

recidx <- mutate(recidx, exposure = durat - start, 
                 interval = factor(interval,  
                                   labels = paste("(", c(0,breaks), ",", 
                                                  c(breaks,100), "]", sep=""))) %>%
  rename(events = fail)
nrow(recidx)


# Finally we show how the observation for subject 9 above becomes five 
# pseudo-observations, with 12 months of exposure in years one to four with no 
# events, and 6 months of exposure in year five with one event.

filter(recidx, id == 9) %>% select(id, start, durat, interval, events, exposure)



# We are now ready to fit a proportional hazards model with a piecewise exponential 
# baseline where the hazard changes from year to year. We use the same model as 
# Wooldridge(2002), involving ten predictors, all fixed covariates.


fit=glm(events~interval+workprg+priors+tserved+felon+alcohol+drugs+black+married+
          educ+age+offset(log(exposure)), data=recidx, family=poisson)
summary(fit)




# The results are exactly the same as in the sister page using Stata. We see that the 
# risk of recidivism is about the same in the first two years, but then decreases 
# substantially with duration since release. At any given duration felons have 25% 
# lower risk of recidivism than non-felons with the same observed characteristics. 
# Subjects imprisoned for alcohol or drug related offenses have much higher risk of 
# recidivism, everything else being equal





# We now illustrate the calculation of survival probabilities. We start with the 
# baseline hazard, which we obtain by adding the constant to the interval 
# coefficients using zero for (0-12] to obtain log-hazards, and exponentiating to 
# get hazards. To obtain the cumulative or integrated hazard we multiply each hazard 
# by the width of the interval, which happens to be 12 for all intervals, and sum. 
# We then obtain the survival as the exponential of the negative cumulative hazard. 
# Note that we only need the first five years.

b = coef(fit)
h = exp( b[1] + c(0,b[2:6]) )
H = cumsum( 12*h)
S = exp(-H)
S

# These calculations apply to the reference cell and are not very meaningful because 
# they set age to zero (and age, by the way, is measured in months).
# 
# We will now estimate the probability of staying out of prison for five years given 
# average values of the predictors. In calculating the mean of each predictor we have 
# to be careful to include only one observation per person, so we restrict the 
# calculation to the first interval, which is "(0,12]". To do this I extract the 
# first row for each person and the columns corresponding to fixed predictors into a 
# matrix X. The relevant coefficients are in slots 7 to 16.

xvars = c("workprg","priors","tserved","felon","alcohol","drugs","black","married",
          "educ","age")
X = recidx[recidx$interval=="(0,12]", xvars]
xbar = colMeans(X)
bx = b[7:16]
xb = sum(xbar * bx)
exp(-H[5] * exp(xb))



# Thus, the probability of staying out of prison for the average person is 65.7%. 
# We can easily calculate this probability for felons and non-felons keeping all 
# other variables at their means. Note that felon is the 4-th predictor in X

xb0 = sum(xbar[-4] * bx[-4])  
exp(-H[5] * exp(xb0))

xb1 = xb0 + bx[4]
exp(-H[5] * exp(xb1))


# The predicted probability is 70.8% for felons and 63.2% for non-felons when all 
# other characteristics are set to the mean, a difference of 7.6 percentage points.
# 
# An alternative calculation sets every person to be a felon or non-felon leaving 
# all other characteristics as they are, and then averages the predicted probability 
# of surviving five years without returning to prison.

felon = X[,"felon"]
X[,"felon"]= 0
xb = as.matrix(X) %*% bx
mean(exp(-(H[5] * exp(xb))))

X[,"felon"]= 1
xb = as.matrix(X) %*% bx
mean(exp(-(H[5] * exp(xb))))

X[,"felon"] = felon

# The average probability of staying out of prison for five years is 68.6% if a 
# felon and 61.2% if not, a difference of 7.4 percentage points.




## Another example: http://data.princeton.edu/wws509/R/c7s1.html

# The datasets page has the original tabulation of children by sex, cohort, age and 
# survival status (dead or still alive at interview), as analyzed by Somoza (1980).
# 
# As is often the case with survival data, a good part of the effort is to convert 
# the raw data into the counts of events and exposure needed for analysis.
# 
# Data Preparation
# 
# We will start by reading the data and collapsing over sex, and will then compute 
# events and exposure to reproduce Table 7.1 in the lecture notes.

somoza <- read.dta("http://data.princeton.edu/wws509/datasets/somoza.dta")

s <- aggregate(dead ~ cohort + age, data=somoza, FUN=sum)

s$alive <- aggregate(alive ~ cohort + age, data=somoza, FUN=sum)[,"alive"]

s <- s[order(s$cohort, s$age),]

# The next step is to compute exposure time. We proceed by cohort, sum all dead and 
# alive to get the total number in the cohort, and then proceed by age subtract 
# those dying or still alive in each each group to obtain the number entering and 
# exiting each age group. Exposure is then the number in the mid-point times the 
# width of the age group, and I chose to express it in years.

w <- c(1,2,3,6,12,36,60,0)/12

s$exposure <- 0

for(cohort in levels(s$cohort)) {
     i <- which(s$cohort == cohort)
     data <- s[i,]
     n <- sum(data$alive + data$dead)
     exit <- n - cumsum(data$alive + data$dead)
     enter <- c(n, exit[-length(exit)])
     s[i,"exposure"] <- w*(enter+exit)/2
   }

co <- subset(s, age != "10+ years")

co$age <- factor(co$age) # dropping level 10+

names(co)[3] <- "deaths"


# After calculating exposure I dropped kids older than ten, as we are only interested
# in survival to age ten. I also renamed "dead" to "deaths", which makes more sense. 
# The data are now ready for analysis.

# In preparation for model fitting I calculate the offset or log of exposure and add 
# it to the data frame.

co$os <- log(co$exposure)

# Now on to the additive model with main effects of age and cohort, which is 
# equivalent to a proportional hazards model:

ph <- glm(deaths ~ age + cohort + offset(os), family=poisson, data=co)

summary(ph)$coefficients


# Let us calculate the fitted life table shown in Table 7.4 of the lecture notes. 
# The predict command following a Poisson regression has and argument type="response" 
# to calculate the expected number of events, which we then divide by exposure to 
# obtain fitted rates. (An alternative is to leave out the typeargument to obtain 
# the linear predictor, substract the offset, and exponentiate to obtain the rate. 
# I'll let you verify that you get the same prediction.)

co$hazard <- predict(ph, type="response")/co$exposure

# At this point recall that the age intervals have different widths. We just don't 
# need the width for 10+. We then loop by cohort, calculate the cumulative hazard 
# by summing the product of the hazard in each age by the width of the interval, and 
# then exponentiate the negative cumulative hazard to get the survival function:

w <- w[-8]

co$survival <- 0

for(cohort in levels(co$cohort)) {
     i <- which(co$cohort == cohort)
     Hazard <- cumsum(co$hazard[i] * w)
     co$survival[i] <- exp(-Hazard)
   }

# Rather than list the data frame I will tabulate the survival by age and cohort:
  
S <- tapply(co$survival,list(co$age,co$cohort),function(x) round(x,4))
S












##########################################################

#         Make fake data

##########################################################



# N = sample size    
# lambda = scale parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C

simulCoxExp <- function(N, h0, beta, rateC)
{
  # covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))
  
  # Simulate survival times using exonential latent event times
  u <- runif(n=N, 0, 1)
  Tlat <- -(log(u) / (h0 * exp(x * beta)))
  
  # censoring times
  C <- rexp(n=N, rate=rateC)
  
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  range2 <- function (x){ (x-min(x))/(max(x)-min(x))*(365-0)+0}
  time2 <- range2(x=time)
  
  # data set
  data.frame(id=1:N,
             time=time,
             time2=time2,
             status=status,
             x=x)
}


datum <- simulCoxExp(100, 0.001, -.6, .001)






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

obs.t <- datum$time
fail <- datum$status
x <- datum$x
t <- c(unique(obs.t),max(obs.t))
T<- as.numeric(length(unique(obs.t)))
N.surv<-as.numeric(length(unique(datum$id)))
eps <-0.000001
c <- 0.001  #confidence in guess for dL0
r <- 0.1    #prior guess at failure rate

win.data <- list("obs.t"=obs.t, "T"=T, "t"=t, "eps"=eps, "x"=x,
                 "fail"=fail, "N.surv"=N.surv)


#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("S.annual.treat", "S.annual.no","S", "dL0", "beta") 


# MCMC Settings 
ni <- 500
nt <- 2
nb <- 250
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)

























































fit <- coxph(Surv(time, status)~x, data=datum)
survs <- survfit(fit, newdata=data.frame(x=1))
plot(survs)

fit <- coxph(Surv(y, event)~x)
survs <- survfit(fit, newdata=data.frame(x=1))
survs2 <- survfit(fit, newdata=data.frame(x=0))
plot(survs)
lines(survs2)


n = 100
round=2
x1 = rbinom(n,size=1,prob=0.5)
x2 = rbinom(n,size=1,prob=0.5)
x3 = rbinom(n,size=1,prob=0.5)
x = t(rbind(x1,x2,x3))
censortime = runif(n,0,1)
survtime= rexp(n,rate=exp(x1+2*x2+3*x3))
survtime = round(survtime,digits=round)
event = as.numeric(censortime>survtime)
y = survtime; 
y[event==0] = censortime[event==0]
t=sort(unique(y[event==1]))
t=c(t,max(censortime))
bigt=length(t)-1




#######

# Cox proportional hazards model code with a counting process. 
# Adapted from Luek Winbugs example

#######


#data for example

start.t = c(276,285,130,1,136,1,143,1,1,147,1,149,201,1,164,29,1,1,1,1,85,1,326,1,1,111,1,1,1,1,248,1,1,1,1,312,1,305,1,1,352,1,1,1,1,1,27,48,1,1,73,1,1,40,1,1,1,313,1,95,268,1,28,114,86,326,1,1,1,1,156,1,1,1,1,1,41,1,133,130,1,1,135,1,136,139,1,140,281,1,134,134,142,1,142,1,144,224,1,132,1,138,1,133,1,1,209,1,135,1,130,132,133,146,133,1,1,1,1,1,145,1,140,1,1,1,1,204,1,1,206,1,1,145,126,1,286,1,1,1,1,136,1,208,1,203,1,1,201,1,1,1,1,1,1,221,1,142,1,1,1,1,218,1,1,141,134,194,1,195,1,200,204,1,1,1,1,1,134,1,1,233,1,1,234,1,1,158,1,129,1,1,137,1,1,158,1,201,126,1,1,1,1,179,187,195,196,132,124,1,1,240,1,150,125,174,1,1,1,160,211,1,1,1,132,1,1,1,141,1,210,1,1,1,1,167,1,1,1,131,1,1,129,209,1,144,162,1,296,1,1,1,1,211,1,219,1,1,220,1,1,1,241,1,152,1,1,1,139,1,184,1,127,1,1,253,1,1,1,1,145,1,1,1,1,1,176,1,190,1,1,1,1,1,1,1,1,248,1,1,1,1,133,1,1,133,133,141,1,1,248,248,124,1,1,1,1,141,1,1,157,172,1,1,324,1,1,1,1,1,1,179,1,195,204,1,1,1,1,211,1,1,1,1,101,1,1,1,1,1,180,1,1,1,1,1,298,253,1,1,143,1,162,163,1,1,236,1,246,1,1,304,1,157,1,165,1,1,197,1,1,240,197,1,1,207,1,147,220,1,1,1,1,241,1,1,1,1,1,241,1,1,1,1,1,1,1,1,1,260,1,1,241,1,1,132,1,1,143,1,161,1,250,1,133,1,133,1,1,102,130,1,1,1,235,1,236,1,1,1,1,1,248,1,1,1,1,301,142,156,1,1,1,133,178,1,1,1,1,141,1,1,1,143,1,1,143,1,1,145,1,1,1,148,1,1,1,153,1,1,1,1,159,1,1,1,181,138,1,262,266,1,319,1,1,1,1,259,264,1,1,1,1,1,1,118,211,247,1,293,1,1,1,132,1,131,1,1,1,144,144,154,1,1,1,199,1,209,1,1,1,1,1,217,1,1,1,151,1,1,1,228,163,1,1,143,169,1,1,1,256,1,179,1,1,1,1,146,1,193,176,1,1,181,1,1,163,1,1,1,206,1,1,1,212,1,1,1,1,217,1,254,1,1,208,1,209,1,215,1,222,1,229,1,204,217,1,1,1,1,1,252,223,225,183,1,1,187,154,1,1,138,152,1,206,1,149,1,207,1,194,141,1,1,142,1,1,1,1,159,1,184,187,1,204,1,176,1,1,1,1,164,1,1,1,1,1,1,172,1,1,1,140,1,199,1,1,201,1,203,1,1,1,200,1,1,1,1,1,199,193,1,1,1,253,1,247,1,299,1,290,1,152,1,1,136,347,1,1,183,1,1,218,1,247,311,1,180,1,1,1,1,1,138,1,1,1,176,1,1,1,1,1,1,172,176,1,1,1,1,247,1,235,1,241,1,188,1,1,161,167,237,1,1,250,1,1,1,249,1,1,1,1,235,176,316,1,179,1,1,1,1,253,1,1,59,160,1,1,1,1,1,167,1,1,134,1,1,1,147,1,1,1,1,1,1,1,287,1,1,1,238,233,1,233,1,176,1,151,1,1,194,1,1,1,1,1,1,1,1,1,300,1,77,1,1,179,1,1,1,1,251,1,1,1,217,1,136,1,1,1,151,154,164,1,1,1,1,305,1,1,1,147,1,144,1,1,144,1,172,222,1,1,138,271,1,138,164,1,1,1,1,146,1,1,176,1,1,1,224,1,1,1,177,158,1,1,162,239,1,139,1,1,1,1,1,125,1,131,171,168,166,1,1,133,1,1,1,270,1,144,1,1,1,177,1,1,1,1,1,162,1,142,1,364,1,156,1,131,1,1,1,1,1,1,138,144,1,1,134,122,1,175,1,154,1,194,1,270,155,1,161,1,1,1,1,1,135,1,1,1,1,1,302,1,1,129,1,147,139,164,149,1,1,1,1,152,178,1,1,178,1,1,1,1,1,179,1,1,1,169,1,1,1,1,168,1,1,136,296,1,168,1,1,324,1,1,1,1,182,1,1,1,293,1,168,301,1,
            162,130,145,153,1,182,1,1,365,1,1,157,1,1,1,1,172,1,1,1,228,1,1,1,1,250,159,1,1,1,51,168,1,314,1,1,1,156,1,202,315,247,1,174,1,307,1,1,1,313,1,1,1,1,150,1,1,156,1,1,228,1,1,1,179,117,1,1,1,171,1,1,1,173,1,1,182,1,1,1,1,301,1,250,1,169,158,168,1,180,36,138,1,265,150,1,1,181,175,1,1,1,180,1,1,1,1,301,1,1,322,1,1,304,1,1,1,171,1,1,192,1,1,180,1,180,1,1,313,313,300,1,236,1,148,1,1,1,1,1,1,64,176,1,1,291,1,1,179,1,142,1,1,1,306,1,302,1,1,1,126,1,1,1,127,1,1,1,307,1,1,1,158,154,1,171,1,319,1,320,1,131,1,1,132,265,255,1,139,1,1,1,176,1,158,181,1,1,1,243,210,1,212,1,1,1,272,195,161,1,306,181,1,132,1,42,291,355,1,1,294,1,1,207,298,1,1,1,316,1,179,1,178,1,2,137,134,1,153,176,1,165,1,42,259,1,1,194,194,1,300,312,1,107,1,311,1,170,1,167,1,1,27,159,1,130,107,1,312,1,1,135)

obs.t = c(333,327,365,48,365,204,365,365,341,365,239,266,365,209,230,365,365,365,365,194,365,5,365,365,149,365,365,365,365,4,365,365,365,365,161,365,17,365,365,326,365,365,365,365,365,178,35,365,365,142,365,365,39,365,365,365,71,365,190,334,365,79,122,125,99,365,365,365,365,301,365,365,365,365,365,273,365,325,344,365,365,327,365,26,208,365,363,251,365,341,299,318,365,87,365,39,347,365,174,365,153,365,59,365,365,274,365,11,365,35,350,214,250,229,365,365,365,365,365,135,365,14,365,365,365,365,135,365,365,55,365,365,29,334,365,11,365,365,365,365,45,365,170,365,69,365,365,341,365,365,365,365,365,365,340,365,29,365,365,365,365,239,365,365,102,313,313,365,38,365,16,282,365,365,365,365,365,107,365,365,26,365,365,27,365,365,120,365,236,365,365,211,365,365,353,365,183,361,365,365,365,365,239,267,344,220,256,289,365,365,271,365,76,246,313,365,365,365,79,331,365,365,365,206,365,365,365,3,365,182,365,365,365,365,44,365,365,365,142,365,365,240,306,365,236,166,365,250,365,365,365,365,72,365,13,365,365,346,365,365,365,118,365,346,365,365,365,310,365,91,365,76,365,365,149,365,365,365,365,66,365,365,365,365,365,41,365,263,365,365,365,365,365,365,365,365,45,365,365,365,365,217,365,365,52,235,336,365,365,10,353,353,365,365,365,365,244,365,365,121,354,365,365,202,365,365,365,365,365,365,218,365,344,351,365,365,365,365,163,365,365,365,365,324,365,365,365,365,365,205,365,365,365,365,365,206,344,365,365,30,365,24,354,365,365,267,365,60,365,365,118,365,82,365,14,365,365,70,365,365,198,337,365,365,98,365,21,308,365,365,365,365,304,365,365,365,365,365,92,365,365,365,365,365,365,365,365,365,4,365,365,323,365,365,266,365,365,165,365,17,365,92,365,155,365,144,365,365,361,247,365,365,365,45,365,130,365,365,365,365,365,76,365,365,365,365,38,313,281,365,365,365,305,308,365,365,365,365,224,365,365,365,149,365,365,45,365,365,284,365,365,365,175,365,365,365,140,365,365,365,365,6,365,365,365,224,231,365,331,319,365,88,365,365,365,365,195,319,365,365,365,365,365,365,23,136,317,365,110,365,365,365,8,365,145,365,365,365,329,213,363,365,365,365,323,365,10,365,365,365,365,365,303,365,365,365,9,365,365,365,9,323,365,365,64,323,365,365,365,353,365,73,365,365,365,365,145,365,2,290,365,365,78,365,365,119,365,365,365,174,365,365,365,313,365,365,365,365,213,365,16,365,365,13,365,15,365,15,365,15,365,15,365,15,329,365,365,365,365,365,95,310,280,284,365,365,245,343,365,365,117,184,365,228,365,295,365,69,365,44,273,365,365,34,365,365,365,365,135,365,13,336,365,64,365,7,365,365,365,365,179,365,365,365,365,365,365,10,365,365,365,354,365,356,365,365,192,365,76,365,365,365,61,365,365,365,365,365,183,273,365,365,365,27,365,345,365,353,365,97,365,98,365,365,166,351,365,365,308,365,365,11,365,118,357,365,234,365,365,365,365,365,204,365,365,365,193,365,365,365,365,365,365,287,237,365,365,365,365,28,365,18,365,286,365,56,365,365,216,318,343,365,365,335,365,365,365,240,365,365,365,365,244,329,279,365,44,365,365,365,365,169,365,365,52,76,365,365,365,365,365,12,365,365,333,365,365,365,228,365,365,365,365,365,365,365,100,365,365,365,324,357,365,74,365,74,365,4,365,365,123,365,365,365,365,365,365,365,365,365,84,365,73,365,365,108,365,365,365,365,20,365,365,365,73,365,353,365,365,365,301,363,321,365,365,365,365,201,365,365,365,129,365,324,365,365,131,365,24,321,365,365,38,305,365,17,329,365,365,365,365,99,365,365,253,365,365,365,248,365,365,365,364,245,365,365,283,343,365,40,365,365,365,365,365,191,365,190,324,324,324,365,365,42,365,365,365,242,365,18,365,365,365,47,365,365,365,365,365,179,365,178,365,50,365,361,365,303,365,365,365,365,365,365,102,324,365,365,9,138,365,114,365,113,365,275,365,85,353,365,325,365,365,365,365,365,324,365,365,365,365,365,215,365,365,273,365,231,276,262,248,365,365,365,365,271,344,365,365,5,365,365,365,365,365,30,365,365,365,229,365,365,365,365,99,365,365,81,325,365,79,365,365,123,365,365,365,365,347,365,365,365,59,365,263,351,365,350,
          283,189,325,365,293,365,365,4,365,365,71,365,365,365,365,4,365,365,365,257,365,365,365,365,314,325,365,365,365,236,57,365,13,365,365,365,69,365,305,298,330,365,55,365,356,365,365,365,89,365,365,365,365,219,365,365,100,365,365,99,365,365,365,115,282,365,365,365,331,365,365,365,234,365,365,321,365,365,365,365,255,365,20,365,293,341,317,365,35,314,278,365,214,322,365,365,110,334,365,365,365,51,365,365,365,365,59,365,365,67,365,365,236,365,365,365,100,365,365,87,365,365,325,365,5,365,365,164,334,321,365,144,365,53,365,365,365,365,365,365,21,165,365,365,142,365,365,10,365,108,365,365,365,100,365,10,365,365,365,317,365,365,365,99,365,365,365,100,365,365,365,35,229,365,332,365,33,365,184,365,262,365,365,99,297,327,365,32,365,365,365,100,365,40,164,365,365,365,100,301,365,137,365,365,365,84,333,327,365,59,333,365,38,365,93,51,325,365,365,46,365,365,93,283,365,365,365,100,365,129,365,108,365,79,93,248,365,32,324,365,7,365,45,222,365,365,93,305,365,253,332,365,24,365,100,365,79,365,100,365,365,84,262,365,100,336,365,99,365,365,39,354)

fail = c(0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,
         0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0)


t <- c(unique(obs.t),max(obs.t))
T<- 238 #length(unique(obs.t))
N.surv<-1270
eps <-0.000001
c <- 0.001  #confidence in guess for dL0
r <- 0.1    #prior guess at failure rate




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
    
        Y[i,j] <- step(obs.t[i] - t[j] + eps)*step(t[j] - start.t[i] + eps) # risk set = 1 if obs.t >= t and start.t <= t
        dN[i,j] <- Y[i,j] * step(t[j+1] - obs.t[i] - eps) * fail[i] 
        # counting process jump if obs.t is in interval and fail = 1
      }
    }


}    


model {

# Now that the data are set up, run the model

    
# Specify priors
  
  # Prior for independent log-normal hazard increments

    #for(j in 1:T){
      #b0.surv[j] ~ dnorm(0, 0.001)
      #dL0[j] <- exp(b0.surv[j])
    #}

  # If you don't want to use the above for hazard, use this adds dL0, c, r, mu to model

    # for(j in 1:T){
    #   dL0[j] ~ dgamma(mu[j], c)
    #   mu[j] <- dL0.star[j]*c
    # }
    


# Likelihood and Intensity process...What is the intensity process?
    
      for(j in 1:T) {
        for(i in 1:N.surv) {

          dN[i,j] ~ dpois(Idt[i,j]) 
          #Idt[i,j] <- Y[i,j]*exp(b0.surv[j])  #Intensity process, this is where you could add in covariates
          Idt[i,j] <- Y[i,j] * dL0[j]        #Intensity process for the other modeling framework
        }

        dL0[j] ~ dgamma(mu[j], c)
        mu[j] <- dL0.star[j]*c


        S[j] <- exp(-sum(dL0[1:j]))
        #S[j] <- pow(exp(-sum(b0.surv[1:j])), exp(b1*X)) #Survival if you have coefficients
    
      }

      #c <- 0.001  #confidence in guess for dL0
      #r <- 0.1    #prior guess at failure rate
    
      for(j in 1:T){
        dL0.star[j] <- r*(t[j+1]-t[j])
      }

      S.annual <- S[T]
    

    }
", fill=TRUE)
sink()




win.data <- list("start.t"=start.t, "obs.t"=obs.t, "T"=T, "t"=t, "eps"=eps, 
                 "fail"=fail, "N.surv"=N.surv, "c"=c, "r"=r)


#  Initial Values	
inits <- function(){list()}


# Parameters to keep track of and report
params <- c("S.annual", "S", "dL0") 


# MCMC Settings 
ni <- 500
nt <- 2
nb <- 250
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)


















n = 100
round=2
x1 = rbinom(n,size=1,prob=0.5)
x = t(rbind(x1))
censortime = runif(n,0,1)
survtime= rexp(n,rate=exp(-.5*x1))
survtime = round(survtime,digits=round)
event = as.numeric(censortime>survtime)
y = survtime; 
y[event==0] = censortime[event==0]
t=sort(unique(y[event==1]))
t=c(t,max(censortime))
bigt=length(t)-1







sink("survival.txt")
cat("
    
    # Set up data

    data {

      for(i in 1:N) {
        for(j in 1:T) {
          Y[i,j] <- step(obs.t[i] - t[j] + eps)
          dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i]
        }
      }

    }


    # Specify model

    model{

      for(i in 1:N){
        betax[i,1] <- 0
      
        for(k in 2:(p+1)){
          betax[i,k] <- betax[i,k-1] + beta[k-1]*x[i,k-1]
        }
      }

      for(j in 1:T) {
        for(i in 1:N) {
          dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
          Idt[i, j] <- Y[i, j] * exp(betax[i,p+1]) * dL0[j] # Intensity
        }
      
        dL0[j] ~ dgamma(mu[j], c)
        mu[j] <- dL0.star[j] * c # prior mean hazard
      }
      
      c <- 0.001
      r <- 0.1
      
      for (j in 1 : T) {
        dL0.star[j] <- r * (t[j + 1] - t[j])
      } 
      
      for(k in 1:p){
        beta[k] ~ dnorm(0.0,0.000001)
      }

    }


", fill=TRUE)
sink()





model=c(1)
x <- x[,model==1]
p <- sum(model) # models have betas with 1
params <- c("beta","dL0")
win.data <- list(x=x,obs.t=y,t=t,T=bigt,N=n,fail=event,eps=1E-10,p=p)
inits <-  function(){list( beta = rep(0,p), dL0 = rep(0.0001,bigt))}


# MCMC Settings 
ni <- 500
nt <- 2
nb <- 250
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "survival.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)

print(out, dig=2)

out2 <- out
mcmcplot(out2)





# Adapted from Princeton prof. http://data.princeton.edu/wws509/R/c7s1.html


require(foreign)

somoza <- read.dta("http://data.princeton.edu/wws509/datasets/somoza.dta")

s <- aggregate(dead ~ cohort + age, data=somoza, FUN=sum)

s$alive <- aggregate(alive ~ cohort + age, data=somoza, FUN=sum)[,"alive"]

s <- s[order(s$cohort, s$age),]


w <- c(1,2,3,6,12,36,60,0)/12

s$exposure <- 0

for(cohort in levels(s$cohort)) {
     i <- which(s$cohort == cohort)
     data <- s[i,]
     n <- sum(data$alive + data$dead)
     exit <- n - cumsum(data$alive + data$dead)
     enter <- c(n, exit[-length(exit)])
     s[i,"exposure"] <- w*(enter+exit)/2
   }

co <- subset(s, age != "10+ years")

co$age <- factor(co$age) # dropping level 10+

names(co)[3] <- "deaths"



































#Survival from IPM example

data {
  for(i in 1:N.surv) {
    for(j in 1:T) {
      Y[i,j] <- step(obs.t[i] - t[j] + eps)*step(t[j] - start.t[i] + eps) # risk set = 1 if obs.t >= t and start.t <= t
      dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i] 
      # counting process jump if obs.t is in interval and fail = 1
    }
  }
}
# Integrated population model 
model {
  # survival
  for(j in 1:T) {
    for(i in 1:N.surv) {
      dN[i, j]   ~ dpois(Idt[i, j]) 
      Idt[i, j] <- Y[i, j] * dL0[j] # Intensity process
    }    
    dL0[j] ~ dgamma(mu2[j], c)
    mu2[j] <- dL0.star[j] * c    
    S[j] <- exp(-sum(dL0[1 : j]))	
  }
  c <- 0.001 # confidence in guess for dL0
  r <- 0.1 # prior guess at failure rate
  for (j in 1 : T) {  dL0.star[j] <- r * (t[j + 1] - t[j])  } 
  # prior guess at the hazard function
  S.annual <- S[T]
  




#data
  
  start.t = c(276,285,130,1,136,1,143,1,1,147,1,149,201,1,164,29,1,1,1,1,85,1,326,1,1,111,1,1,1,1,248,1,1,1,1,312,1,305,1,1,352,1,1,1,1,1,27,48,1,1,73,1,1,40,1,1,1,313,1,95,268,1,28,114,86,326,1,1,1,1,156,1,1,1,1,1,41,1,133,130,1,1,135,1,136,139,1,140,281,1,134,134,142,1,142,1,144,224,1,132,1,138,1,133,1,1,209,1,135,1,130,132,133,146,133,1,1,1,1,1,145,1,140,1,1,1,1,204,1,1,206,1,1,145,126,1,286,1,1,1,1,136,1,208,1,203,1,1,201,1,1,1,1,1,1,221,1,142,1,1,1,1,218,1,1,141,134,194,1,195,1,200,204,1,1,1,1,1,134,1,1,233,1,1,234,1,1,158,1,129,1,1,137,1,1,158,1,201,126,1,1,1,1,179,187,195,196,132,124,1,1,240,1,150,125,174,1,1,1,160,211,1,1,1,132,1,1,1,141,1,210,1,1,1,1,167,1,1,1,131,1,1,129,209,1,144,162,1,296,1,1,1,1,211,1,219,1,1,220,1,1,1,241,1,152,1,1,1,139,1,184,1,127,1,1,253,1,1,1,1,145,1,1,1,1,1,176,1,190,1,1,1,1,1,1,1,1,248,1,1,1,1,133,1,1,133,133,141,1,1,248,248,124,1,1,1,1,141,1,1,157,172,1,1,324,1,1,1,1,1,1,179,1,195,204,1,1,1,1,211,1,1,1,1,101,1,1,1,1,1,180,1,1,1,1,1,298,253,1,1,143,1,162,163,1,1,236,1,246,1,1,304,1,157,1,165,1,1,197,1,1,240,197,1,1,207,1,147,220,1,1,1,1,241,1,1,1,1,1,241,1,1,1,1,1,1,1,1,1,260,1,1,241,1,1,132,1,1,143,1,161,1,250,1,133,1,133,1,1,102,130,1,1,1,235,1,236,1,1,1,1,1,248,1,1,1,1,301,142,156,1,1,1,133,178,1,1,1,1,141,1,1,1,143,1,1,143,1,1,145,1,1,1,148,1,1,1,153,1,1,1,1,159,1,1,1,181,138,1,262,266,1,319,1,1,1,1,259,264,1,1,1,1,1,1,118,211,247,1,293,1,1,1,132,1,131,1,1,1,144,144,154,1,1,1,199,1,209,1,1,1,1,1,217,1,1,1,151,1,1,1,228,163,1,1,143,169,1,1,1,256,1,179,1,1,1,1,146,1,193,176,1,1,181,1,1,163,1,1,1,206,1,1,1,212,1,1,1,1,217,1,254,1,1,208,1,209,1,215,1,222,1,229,1,204,217,1,1,1,1,1,252,223,225,183,1,1,187,154,1,1,138,152,1,206,1,149,1,207,1,194,141,1,1,142,1,1,1,1,159,1,184,187,1,204,1,176,1,1,1,1,164,1,1,1,1,1,1,172,1,1,1,140,1,199,1,1,201,1,203,1,1,1,200,1,1,1,1,1,199,193,1,1,1,253,1,247,1,299,1,290,1,152,1,1,136,347,1,1,183,1,1,218,1,247,311,1,180,1,1,1,1,1,138,1,1,1,176,1,1,1,1,1,1,172,176,1,1,1,1,247,1,235,1,241,1,188,1,1,161,167,237,1,1,250,1,1,1,249,1,1,1,1,235,176,316,1,179,1,1,1,1,253,1,1,59,160,1,1,1,1,1,167,1,1,134,1,1,1,147,1,1,1,1,1,1,1,287,1,1,1,238,233,1,233,1,176,1,151,1,1,194,1,1,1,1,1,1,1,1,1,300,1,77,1,1,179,1,1,1,1,251,1,1,1,217,1,136,1,1,1,151,154,164,1,1,1,1,305,1,1,1,147,1,144,1,1,144,1,172,222,1,1,138,271,1,138,164,1,1,1,1,146,1,1,176,1,1,1,224,1,1,1,177,158,1,1,162,239,1,139,1,1,1,1,1,125,1,131,171,168,166,1,1,133,1,1,1,270,1,144,1,1,1,177,1,1,1,1,1,162,1,142,1,364,1,156,1,131,1,1,1,1,1,1,138,144,1,1,134,122,1,175,1,154,1,194,1,270,155,1,161,1,1,1,1,1,135,1,1,1,1,1,302,1,1,129,1,147,139,164,149,1,1,1,1,152,178,1,1,178,1,1,1,1,1,179,1,1,1,169,1,1,1,1,168,1,1,136,296,1,168,1,1,324,1,1,1,1,182,1,1,1,293,1,168,301,1,
              162,130,145,153,1,182,1,1,365,1,1,157,1,1,1,1,172,1,1,1,228,1,1,1,1,250,159,1,1,1,51,168,1,314,1,1,1,156,1,202,315,247,1,174,1,307,1,1,1,313,1,1,1,1,150,1,1,156,1,1,228,1,1,1,179,117,1,1,1,171,1,1,1,173,1,1,182,1,1,1,1,301,1,250,1,169,158,168,1,180,36,138,1,265,150,1,1,181,175,1,1,1,180,1,1,1,1,301,1,1,322,1,1,304,1,1,1,171,1,1,192,1,1,180,1,180,1,1,313,313,300,1,236,1,148,1,1,1,1,1,1,64,176,1,1,291,1,1,179,1,142,1,1,1,306,1,302,1,1,1,126,1,1,1,127,1,1,1,307,1,1,1,158,154,1,171,1,319,1,320,1,131,1,1,132,265,255,1,139,1,1,1,176,1,158,181,1,1,1,243,210,1,212,1,1,1,272,195,161,1,306,181,1,132,1,42,291,355,1,1,294,1,1,207,298,1,1,1,316,1,179,1,178,1,2,137,134,1,153,176,1,165,1,42,259,1,1,194,194,1,300,312,1,107,1,311,1,170,1,167,1,1,27,159,1,130,107,1,312,1,1,135)
              
  obs.t = c(333,327,365,48,365,204,365,365,341,365,239,266,365,209,230,365,365,365,365,194,365,5,365,365,149,365,365,365,365,4,365,365,365,365,161,365,17,365,365,326,365,365,365,365,365,178,35,365,365,142,365,365,39,365,365,365,71,365,190,334,365,79,122,125,99,365,365,365,365,301,365,365,365,365,365,273,365,325,344,365,365,327,365,26,208,365,363,251,365,341,299,318,365,87,365,39,347,365,174,365,153,365,59,365,365,274,365,11,365,35,350,214,250,229,365,365,365,365,365,135,365,14,365,365,365,365,135,365,365,55,365,365,29,334,365,11,365,365,365,365,45,365,170,365,69,365,365,341,365,365,365,365,365,365,340,365,29,365,365,365,365,239,365,365,102,313,313,365,38,365,16,282,365,365,365,365,365,107,365,365,26,365,365,27,365,365,120,365,236,365,365,211,365,365,353,365,183,361,365,365,365,365,239,267,344,220,256,289,365,365,271,365,76,246,313,365,365,365,79,331,365,365,365,206,365,365,365,3,365,182,365,365,365,365,44,365,365,365,142,365,365,240,306,365,236,166,365,250,365,365,365,365,72,365,13,365,365,346,365,365,365,118,365,346,365,365,365,310,365,91,365,76,365,365,149,365,365,365,365,66,365,365,365,365,365,41,365,263,365,365,365,365,365,365,365,365,45,365,365,365,365,217,365,365,52,235,336,365,365,10,353,353,365,365,365,365,244,365,365,121,354,365,365,202,365,365,365,365,365,365,218,365,344,351,365,365,365,365,163,365,365,365,365,324,365,365,365,365,365,205,365,365,365,365,365,206,344,365,365,30,365,24,354,365,365,267,365,60,365,365,118,365,82,365,14,365,365,70,365,365,198,337,365,365,98,365,21,308,365,365,365,365,304,365,365,365,365,365,92,365,365,365,365,365,365,365,365,365,4,365,365,323,365,365,266,365,365,165,365,17,365,92,365,155,365,144,365,365,361,247,365,365,365,45,365,130,365,365,365,365,365,76,365,365,365,365,38,313,281,365,365,365,305,308,365,365,365,365,224,365,365,365,149,365,365,45,365,365,284,365,365,365,175,365,365,365,140,365,365,365,365,6,365,365,365,224,231,365,331,319,365,88,365,365,365,365,195,319,365,365,365,365,365,365,23,136,317,365,110,365,365,365,8,365,145,365,365,365,329,213,363,365,365,365,323,365,10,365,365,365,365,365,303,365,365,365,9,365,365,365,9,323,365,365,64,323,365,365,365,353,365,73,365,365,365,365,145,365,2,290,365,365,78,365,365,119,365,365,365,174,365,365,365,313,365,365,365,365,213,365,16,365,365,13,365,15,365,15,365,15,365,15,365,15,329,365,365,365,365,365,95,310,280,284,365,365,245,343,365,365,117,184,365,228,365,295,365,69,365,44,273,365,365,34,365,365,365,365,135,365,13,336,365,64,365,7,365,365,365,365,179,365,365,365,365,365,365,10,365,365,365,354,365,356,365,365,192,365,76,365,365,365,61,365,365,365,365,365,183,273,365,365,365,27,365,345,365,353,365,97,365,98,365,365,166,351,365,365,308,365,365,11,365,118,357,365,234,365,365,365,365,365,204,365,365,365,193,365,365,365,365,365,365,287,237,365,365,365,365,28,365,18,365,286,365,56,365,365,216,318,343,365,365,335,365,365,365,240,365,365,365,365,244,329,279,365,44,365,365,365,365,169,365,365,52,76,365,365,365,365,365,12,365,365,333,365,365,365,228,365,365,365,365,365,365,365,100,365,365,365,324,357,365,74,365,74,365,4,365,365,123,365,365,365,365,365,365,365,365,365,84,365,73,365,365,108,365,365,365,365,20,365,365,365,73,365,353,365,365,365,301,363,321,365,365,365,365,201,365,365,365,129,365,324,365,365,131,365,24,321,365,365,38,305,365,17,329,365,365,365,365,99,365,365,253,365,365,365,248,365,365,365,364,245,365,365,283,343,365,40,365,365,365,365,365,191,365,190,324,324,324,365,365,42,365,365,365,242,365,18,365,365,365,47,365,365,365,365,365,179,365,178,365,50,365,361,365,303,365,365,365,365,365,365,102,324,365,365,9,138,365,114,365,113,365,275,365,85,353,365,325,365,365,365,365,365,324,365,365,365,365,365,215,365,365,273,365,231,276,262,248,365,365,365,365,271,344,365,365,5,365,365,365,365,365,30,365,365,365,229,365,365,365,365,99,365,365,81,325,365,79,365,365,123,365,365,365,365,347,365,365,365,59,365,263,351,365,350,
            283,189,325,365,293,365,365,4,365,365,71,365,365,365,365,4,365,365,365,257,365,365,365,365,314,325,365,365,365,236,57,365,13,365,365,365,69,365,305,298,330,365,55,365,356,365,365,365,89,365,365,365,365,219,365,365,100,365,365,99,365,365,365,115,282,365,365,365,331,365,365,365,234,365,365,321,365,365,365,365,255,365,20,365,293,341,317,365,35,314,278,365,214,322,365,365,110,334,365,365,365,51,365,365,365,365,59,365,365,67,365,365,236,365,365,365,100,365,365,87,365,365,325,365,5,365,365,164,334,321,365,144,365,53,365,365,365,365,365,365,21,165,365,365,142,365,365,10,365,108,365,365,365,100,365,10,365,365,365,317,365,365,365,99,365,365,365,100,365,365,365,35,229,365,332,365,33,365,184,365,262,365,365,99,297,327,365,32,365,365,365,100,365,40,164,365,365,365,100,301,365,137,365,365,365,84,333,327,365,59,333,365,38,365,93,51,325,365,365,46,365,365,93,283,365,365,365,100,365,129,365,108,365,79,93,248,365,32,324,365,7,365,45,222,365,365,93,305,365,253,332,365,24,365,100,365,79,365,100,365,365,84,262,365,100,336,365,99,365,365,39,354)
                        
  fail = c(0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,
           0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0)
                                 

t <- c(unique(obs.t),max(obs.t))
T<- 238 #length(unique(obs.t))
N.surv<-1270
eps <-0.000001









#Survival code for Cox proportional hazards, using the Leuk example

  
  
#data
  
  

model leuk;

const
N = 42,    # number of patients
T = 17,    # number of unique failure times
eps = 0.000001;  # used to guard against numerical 
# imprecision in step function  

var
obs.t[N],  # observed failure or censoring time for each patient
t[T+1],  # unique failure times + maximum censoring time
dN[N,T],  # counting process increment
Y[N,T],  # 1=subject observed; 0=not observed 
Idt[N,T],  # intensity process
Z[N],  # covariate     
beta,  # regression coefficient
dL0[T],  # increment in unknown hazard function
beta0[T],  # log(increment in unknown hazard function)
dL0.star[T],  # prior guess at hazard function
c,  # degree of confidence in prior guess for dL0
mu[T],  # location parameter for Gamma (= c * dL0.star)
r,  # prior guess at failure rate
fail[N],  # failure = 1; censored = 0
S.treat[T],  # survivor function for treatment group
S.placebo[T];  # survivor function for placebo group

data obs.t, fail, Z in "leuk.dat", t in "failtime.dat";
inits in "leuk.in";



{
  # Set up data
  
  for(i in 1:N) {
    for(j in 1:T) {
      
      # risk set = 1 if obs.t >= t
      Y[i,j] <- step(obs.t[i] - t[j] + eps);  
      
      # counting process jump = 1 if obs.t in [ t[j], t[j+1] )
      #                      i.e. if t[j] <= obs.t < t[j+1]
      dN[i,j] <- Y[i,j]*step(t[j+1] - obs.t[i] - eps)*fail[i]; 
      
    }
  }
  # Model 
  
  for(j in 1:T) {
    
    #   beta0[j] ~ dnorm(0,0.001); # include this when using Poisson trick
    
    for(i in 1:N) {
      
      dN[i,j]   ~ dpois(Idt[i,j]);              # Likelihood
      Idt[i,j] <- Y[i,j]*exp(beta*Z[i])*dL0[j]; # Intensity 
      
      #  Try Poisson trick - independent log-normal hazard increments
      #                    - enables dL0, c, r, mu to be dropped from model
      #  Idt[i,j] <- Y[i,j]*exp(beta0[j]+beta*Z[i]); # Intensity    
      
    }     
    
    dL0[j] ~ dgamma(mu[j], c);
    mu[j] <- dL0.star[j] * c;    # prior mean hazard
    
    # Survivor function = exp(-Integral{l0(u)du})^exp(beta*z)    
    S.treat[j] <- pow(exp(-sum(dL0[1:j])), exp(beta * -0.5));
    S.placebo[j] <- pow(exp(-sum(dL0[1:j])), exp(beta * 0.5));	
    
  }
  
  c <- 0.001; r <- 0.1; 
  for (j in 1:T) {  
    dL0.star[j] <- r * (t[j+1]-t[j])  
  } 
  
  beta ~ dnorm(0.0,0.000001);                 
}