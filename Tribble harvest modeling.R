
##' Now lets think of a dynamic decision here we want to determine the 
##' optimal harvest rate of tribbles, tribble life history is such that
##' we can use logistic population model to cover the dynamics of the
##' whole population below we specify the population parameters and evaluate
##' the **cumulative harvest** over a ten year period. Note that the harvest 
##' rate remains constant across the years and the model is deterministic. 

## maximum population growth rate
rmax<-0.3
## carrying capacity
K<-125
#population dynamics function
tada<-c()
for(H.rate in 0:90*0.01){
  #initial population
  N.t<- 100  
  ## initial cumulative harvest
  T.harvest<-0
  ## year loop
  for(i in 1:10){
    ## pop dynamics
    N.t2<-N.t+rmax*N.t*(1-N.t/K)
    ## harvest some tribbles
    harvest<- N.t2*H.rate  
    ## surviving tribbles
    N.t <- N.t2 - harvest
    ## accumulate the harvest
    T.harvest<- T.harvest + harvest
  }
  tada<- rbind(tada,c(H.rate,T.harvest))
}

plot(tada[,2]~tada[,1],xlab = "Harvest rate", ylab = "Total harvest")

############################################################
## lets add a little annual stochasisity to the model 


tada<-c()
for(H.rate in 0:90*0.01){
  #initial population
  N.t<- 100
  ## initial cumulative harvest
  T.harvest<-0
  ## year loop
  for(i in 1:10){
    ## maximum population growth rate
    rmax<-rgamma(1,36,1/0.0083) ##mean 0.3 sd 0.05
    ## carrying capacity
    K<-round(rgamma(1,25,1/5)) ##mean 125 sd 25
    ## pop dynamics
    N.t2<-round(N.t+rmax*N.t*(1-N.t/K))
    ## harvest some tribbles
    harvest<- rbinom(1,N.t2,H.rate)  
    ## surviving tribbles
    N.t <- N.t2 - harvest
    ## accumulate the harvest
    T.harvest<- T.harvest + harvest
  }
  tada<- rbind(tada,c(H.rate,T.harvest))
}

plot(tada[,2]~tada[,1],xlab = "Harvest rate", ylab = "Total harvest")

#' we could solve this using simulation but first lets create a function
#' that does the dynamics in 1 year
#' here the inputs are the harvest rate and the current tribble 
#' population size. Remember this is our system state and decision 
#' (for later think state dependent decision making)

popn.harv<-function(N.t,H.rate){
  ## maximum population growth rate
  rmax<-rgamma(1,36,1/0.0083) ##mean 0.3 sd 0.05
  ## carrying capacity
  K<-round(rgamma(1,25,1/5)) ##mean 125 sd 25
  ## pop dynamics
  N.t2<-round(N.t+rmax*N.t*(1-N.t/K))
  ## harvest some tribbles
  harvest<- rbinom(1,N.t2,H.rate)  
  ## surviving tribbles
  N.t <- N.t2 - harvest
  # return a vector with harvest and next population size
  return(c(harvest,N.t))
}


################ using the function we now have


tada<-c()
for(H.rate in 0:90*0.01){
  #initial population
  N.t<- 100
  ## initial cumulative harvest
  T.harvest<-0
  ## year loop
  for(i in 1:10){
    ret.n.fut<-popn.harv(N.t,H.rate)
    ## add yearly harvest to total
    T.harvest<- T.harvest + ret.n.fut[1]
    ## adult tribbles surviving to reproduce
    N.t<-ret.n.fut[2]
  }
  tada<- rbind(tada,c(H.rate,T.harvest))
}

plot(tada[,2]~tada[,1],,xlab = "Harvest rate", ylab = "Total harvest")


#' lets add a replicate loop so we can run multiple
#' sumulations and calculate the expected cumulative harvest
#' (i.e., the average cumulative harvest)

tada<-c()
for(H.rate in 0:90*0.01){
  all.harv<- NULL
  for(rep in 1:500){
    ## initial population
    N.t<- 100
    ## initial cumulative harvest
    T.harvest<-0
    ## year loop
    for(i in 1:10){
      ret.n.fut<-popn.harv(N.t,H.rate)
      ## add yearly harvest to total
      T.harvest<- T.harvest + ret.n.fut[1]
      ## adult tribbles surviving to reproduce
      N.t<-ret.n.fut[2]
    } ## end of annual loop
    
    all.harv<-c(all.harv,T.harvest)
  } #end of replicate loop
  
  tada<- rbind(tada,c(H.rate,mean(all.harv)))
} # end of decision loop

plot(tada[,2]~tada[,1],xlab = "Harvest rate", ylab = "Total harvest")

## where is maximum harvest rate?
tada[tada[,2]== max(tada[,2]),1]




