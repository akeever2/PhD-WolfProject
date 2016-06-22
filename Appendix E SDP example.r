## Example from Conroy and Peterson 2013 SDP example in Appendix E
## Stochastic dynamic programming requires this library
library(MDPtoolbox)
# With a decision specific transition matrices matrix
P <- array(0, c(3,3,3))
## decision 1 harvest = 0.1
P[,,1] <- matrix(c(0.2, 0.5, 0.3, 0.2,.3,0.5,0.1,0.3,0.6), 3, 3, byrow=TRUE)
## decision 1 harvest = 0.2
P[,,2] <- matrix(c(0.5,0.3,0.2,0.2,0.5,0.3,0.2,0.4,0.4), 3, 3, byrow=TRUE)
## decision 1 harvest = 0.3
P[,,3] <- matrix(c(0.7,0.3,0.0,0.6,0.3,0.1,0.3,0.5,0.2), 3, 3, byrow=TRUE)
##Reward matrix
R <- matrix(c(0.5,1.0,1.5,1,2,3,1.5,3.0,4.5), 3, 3, byrow=TRUE)

## step through each iteration and examine the policies
ini<-c(0,0,0)
val<- mdp_bellman_operator(P, R, 1, ini)
val
ini<-val$V
val<- mdp_bellman_operator(P, R, 1, ini)
val
ini<-val$V
val<- mdp_bellman_operator(P, R, 1, ini)
val
ini<-val$V
val<- mdp_bellman_operator(P, R, 1, ini)
ini<-val$V
val<- mdp_bellman_operator(P, R, 1, ini)
val
ini<-val$V
val<- mdp_bellman_operator(P, R, 1, ini)
val

## heres the better way if automatically 
# iertates and stops when the policy is stable
best<-mdp_policy_iteration(P, R, discount=.99999)
best
# extract best policy
the.best.dec<-best$policy

##########################################################################
##########################################################################
## Now, forward simulation with optimal state dependent policy for 50 yrs.
## We don't need to change much. We can use the TMP we created above
## BUT BUT BUT they are row stochastic (rows sum to one) We need to transpose
## "t()" to make column stochastic.

## Lets see, this is the transposed TMP for decision 1
t(P[,,1])

## This is the transposed TMP for decision 2 etc
t(P[,,2])

## Harvest rate matrix for each decision
h.rate<- c(0.1,0.2,0.3)

#' Let's say that the estimated current population is 8 million with 
#' a sd of 3 the cutoff points for the 3 population states are 
#' 0 to 7.5 mil, 7.5 to 12.5 mil, and 12.5 plus mil (midpoint popn sizes of 5,10,15)
#' Let's use a normal CFD to get the state specfic probabilities for 
#' each state:

  lo <- round(pnorm(7.5, mean = 8, sd = 3),3)
  med <-round(pnorm(12.5, mean = 8, sd = 3) - lo,3)
  hi <- 1-(lo+med)

## create population state vector
popn<-c(lo,med,hi)

# E.G., what is population state?
max(popn) == popn

# place to put total 50 year harvest
T.harvest <- 0

# lets keep track of the population each year
all.popn <- NULL

### enough set up lets run and start year loop 
for(year in 1:50){
  
# find best decision using decision vector above
# this means the decision will change depending on last years
# population size
dec<- sum((max(popn) == popn)* the.best.dec)

# use best decision number to specify best decision TPM
# and calculate best population
pop.new<-t(P[,,dec]) %*% popn

# calculate annual harvest note we estimate population by
# probability weighting the
harv<- sum(popn*c(5,10,15))*h.rate[dec]

## Calculate total harvest
T.harvest<- T.harvest + harv

## add population to time series
all.popn<- c(all.popn,sum(popn*c(5,10,15)))

## current population becomes future population next time step
popn<- pop.new
}

## total harvest over 50 years
T.harvest

## plot population size each year
plot(all.popn~c(1:50))


