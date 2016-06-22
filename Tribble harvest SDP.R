
#' back to our tribble harvest function
#' here we modified it to handle vectors as inputs
#' that is multiple N.t and H.rate

popn.harv<-function(N.t,H.rate){
  ## maximum population growth rate note use of length function
  rmax<-rgamma(length(N.t),36,1/0.0083) ##mean 0.3 sd 0.05
  ## carrying capacity
  K<-round(rgamma(length(N.t),25,1/5)) ##mean 125 sd 25
  ## pop dynamics
  N.t2<-round(N.t+rmax*N.t*(1-N.t/K))
  ## harvest some tribbles
  harvest<- rbinom(length(N.t),N.t2,H.rate)  
  ## surviving tribbles
  N.t <- N.t2 - harvest
  # return a matrix "cbind" with harvest and next population size
  return(cbind(harvest,N.t))
}

#############################################################
#' now it is time to build the tramsition probability matrices for
#' to solve the tribble harvest MDP problem. we will be constructing a
#' TMP and return (utility) vector for each decision alternative
#' 

  ## matrix each row is alpha and beta of beta distribution
  ## for each decision for 4 harvest levels and cv 30%
  harv=matrix(c(9.9, 89.1, 8.69,34.76, 7.48, 
                17.45,6.27,9.4),ncol=2,byrow = T)
R<- NULL
#for(decis in 1:4){
 decis = 1

#' it is generally easier to simulate the population dymanics and use those to
#' create the TMP and return vector lets do that for a decision

  #'Initial population size, we'd like to catch the full range and replicate
  #'that multi times 
    N.tm1<- rep(seq(3,200,by=5),200)

    ## generate random harvest for a decision each year
    H.rate<-rbeta(length(N.tm1),harv[decis,1],harv[decis,2])

    ## run the dynamics function
    ret.n.fut<-popn.harv(N.tm1,H.rate)

    #combine into single matrix
    sims<-cbind(ret.n.fut,N.tm1)

    # lets break into 5 groups. We would like data to be relatively
    # evenly distributed among the groups so we use quantile function 
    # to figure out the breaks between groups
    # note there is a ggplot function that does this but be careful
    breaks<-quantile(c(sims[,2:3]), probs=seq(0,1, by= 0.2))

   ## now break up into groups using cut function in Hmisc package
require(Hmisc)
  # population size in t
  Ntgrp<- cut(sims[,2],breaks=as.numeric(breaks),labels = F)
  # population size in t
  Ntm1grp<- cut(sims[,3],breaks=as.numeric(breaks),labels = F)

## create a table of transition frequencies that will be turned into 
## state transition probabilities one for each decision alternative
## calculate frequencies of initial population and future population
TM<- table(Ntm1grp,Ntgrp)

### These create the transition matrices using the above table you will want one for each 
## harvest decision alternative note procuces ROW STOCHASTIC  
TM <- prop.table(TM,1)
rowSums(TM)

## calculate the average (expected) return for each population
## state syntax is calculate average of things in sims[,1] by 
## Ntm1 group (initial popn group), calculate mean
Return<-tapply(sims[,1],Ntm1grp, mean)

# Set up arrays for solving 
## dimensions are #initial state, # final states, no decisions
if(decis == 1) P <- array(0, c(5,5,4))
P[,,decis] <- TM

# each column is a decision and row is a system state
R <- cbind(R,Return)

#}

# require(MDPtoolbox)
# 
# ### now find optimal state dependent harvest 
# mdp_policy_iteration(P, R, discount=.99999)

 

