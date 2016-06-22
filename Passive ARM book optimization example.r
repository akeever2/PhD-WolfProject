## Example from Conroy and Peterson 2013 Passive ARM SDP example in Appendix E
## Stochastic dynamic programming requires this library
library(MDPtoolbox)
# With a decision specific transition matrices matrix
P <- array(0, c(3,3,3))

#### weight of model 1
mod1.wt<- 0.5

## decision 1 harvest = 0.1
Model1.1 <- matrix(c(0.2, 0.5, 0.3, 0.2,.3,0.5,0.1,0.3,0.6), 3, 3, byrow=TRUE)
Model2.1 <- matrix(c(0.3, 0.5, 0.2, 0.3,0.3,0.4,0.2,0.3,0.5), 3, 3, byrow=TRUE)
## model averaged transition matrix
P[,,1] <- Model1.1*mod1.wt + Model2.1*(1-mod1.wt)

## decision 1 harvest = 0.2
Model1.2<- matrix(c(0.5,0.3,0.2,0.2,0.5,0.3,0.2,0.4,0.4), 3, 3, byrow=TRUE)
Model2.2<- matrix(c(0.7,0.3,0.0,0.4,0.5,0.1,0.5,0.4,0.1), 3, 3, byrow=TRUE)
P[,,2] <-  Model1.2*mod1.wt + Model2.2*(1-mod1.wt)

## decision 1 harvest = 0.3
Model1.3<- matrix(c(0.7,0.3,0.0,0.6,0.3,0.1,0.3,0.5,0.2), 3, 3, byrow=TRUE)
Model2.3<- matrix(c(0.9,0.1,0.0,0.7,0.3,0.0,0.4,0.6,0.0), 3, 3, byrow=TRUE)
P[,,3] <-  Model1.3*mod1.wt + Model2.3*(1-mod1.wt)

##Reward matrix
R <- matrix(c(0.5,1.0,1.5,1,2,3,1.5,3.0,4.5), 3, 3, byrow=TRUE)


## heres the better way if automatically 
# iertates and stops when the policy is stable
mdp_policy_iteration(P, R, discount=.9999999)


