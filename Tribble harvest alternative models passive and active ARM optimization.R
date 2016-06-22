
require(rgenoud)

#' we are now ready to try to find the optimal solution using 
#' something other than brute force
#' here the harvest function as before
#' 
popn.harv<-function(N.t,H.rate,mod1.wt){
  ## maximum population growth rate
  rmax<-rgamma(1,36,1/0.0083) ##mean 0.3 sd 0.05
  ## carrying capacity
  K<-round(rgamma(1,25,1/5)) ##mean 125 sd 25
  ## pop dynamics
  ### carrying capacity
  mod1.N.t2<-round(N.t+rmax*N.t*(1-N.t/K))
  ### no carrying capacity
  mod2.N.t2<-round(N.t+rmax*N.t)
  ## model averaged estimate
  N.t2<-round(mod1.N.t2*mod1.wt+mod2.N.t2*(1-mod1.wt))
  ## harvest some tribbles
  harvest<- rbinom(1,N.t2,H.rate)  
  ## surviving tribbles
  N.t <- N.t2 - harvest
  # return a vector with harvest and next population size
  return(c(harvest,N.t))
}


## decis in now a vector
find.best<-function(decis){
  ## matrix each row is alpha and beta of beta distribution
  ## for each decision for 4 harvest levels and cv 30%
  harv=matrix(c(9.9, 89.1, 8.69,34.76, 7.48, 
                17.45,6.27,9.4),ncol=2,byrow = T)
  N.t<- initial.N
  mod1.wt<-model.weight
  T.harvest<-0
  for(i in 1:length(decis)){
    ## generate random harvest for a decision each year
    H.rate<-rbeta(1,harv[decis[i],1],harv[decis[i],2])
    ret.n.fut<-popn.harv(N.t,H.rate,mod1.wt)
    ## add yearly harvest to total
    T.harvest<- T.harvest + ret.n.fut[1]
    ## adult tribbles surviving to reproduce
    N.t<-ret.n.fut[2]
  }
  ## notice NO MORE negative sign
  return(T.harvest)
}

## does it work?
model.weight<- 0.5
initial.N <- 10
find.best(rep(1,5))

## lets get the optimal harvest under passive adaptive management over a range of current popn states
# and model weights (aka information states)
state.dep<-NULL

yrs<-10
# we need 10 upper and lower bounds
min.max<-cbind(rep(1,yrs),rep(4,yrs))
### this will take some time because of low tolerance and large pop.size for the genetic algorithm
for(initial.N in c(50,150)){
  for(model.weight in c(0.33,0.66)){
  cols<-NULL
  for(i in 1:5){
  xx<- genoud(
    fn = find.best,             
    nvars = yrs,                 
    max = TRUE, 
    pop.size = 2000,
    data.type.int = TRUE,         
    starting.values=rep(1,yrs),   
    Domains = min.max,
    boundary.enforcement = 2,
    solution.tolerance = 0.0000000001) 
 cols<-rbind(cols,xx$par)
}
  state.dep<- rbind(state.dep,c(initial.N,model.weight,round(colMeans(cols))))
  }
}

## lets take a look
state.dep

#' now lets allow the model weights to change through time, this will take a few more
#' code changes 
#' 

popn.harv<-function(N.t,decis,mod.wt){
  
 ## beta parameters for harvest rates for each decision
  harv=matrix(c(9.9, 89.1, 8.69,34.76, 7.48, 
                17.45,6.27,9.4),ncol=2,byrow = T)
  ## set total harvest to zero
  T.harvest<- 0
  for(i in 1:length(decis)){
  H.rate<-rbeta(1,harv[decis[i],1],harv[decis[i],2])
  ## maximum population growth rate
  rmax<-rgamma(1,36,1/0.0083) ##mean 0.3 sd 0.05
  ## carrying capacity
  K<-round(rgamma(1,25,1/5)) ##mean 125 sd 25
  ## pop dynamics
  ### carrying capacity
  mod1.N.t2<-round(N.t+rmax*N.t*(1-N.t/K))
  ### no carrying capacity
  mod2.N.t2<-round(N.t+rmax*N.t)
  ## model averaged estimate
  N.t2<-round(mod1.N.t2*mod.wt+mod2.N.t2*(1-mod.wt))
  ## harvest some tribbles
  harvest<- rbinom(1,N.t2,H.rate) 
  ## surviving tribbles
  N.t <- N.t2 - harvest
  # Bayes Rule to calculate posterior. Lets say I ran enough sinulations to know the the standard deviation of the estimate
  # is normally dist with a CV = 0.32
  mod.wt<- (dnorm(x=N.t2, mean =mod1.N.t2, sd = mod1.N.t2*0.32)*mod.wt)/(dnorm(x=N.t2, mean =mod1.N.t2, sd = mod1.N.t2*0.32)*mod.wt + dnorm(x=N.t2, mean =mod2.N.t2, sd = mod2.N.t2*0.32)*(1-mod.wt))
  
  T.harvest<- T.harvest + harvest
  N.t<- N.t2
  }
  return(T.harvest)
}

## decis in now a vector
find.best<-function(decis){
  N.t<- initial.N
  mod.wt<-model.weight
  T.harv<-popn.harv(N.t,decis,mod.wt)
  return(T.harv)
}
model.weight<- 0.5
initial.N <- 10
find.best(rep(1,5))

## lets get the optimal harvest over a range of current popn
state.dep<-NULL

yrs<-10
# we need 10 upper and lower bounds
min.max<-cbind(rep(1,yrs),rep(4,yrs))

for(initial.N in c(50,150)){
  for(model.weight in c(0.33,0.66)){
    cols<-NULL
    for(i in 1:5){
      xx<- genoud(
        fn = find.best,             
        nvars = yrs,                 
        max = TRUE, 
        pop.size = 1000,
        data.type.int = TRUE,         
        starting.values=rep(1,yrs),   
        Domains = min.max,
        boundary.enforcement = 2,
        solution.tolerance = 0.00001) 
      cols<-rbind(cols,xx$par)
    }
    state.dep<- rbind(state.dep,c(initial.N,model.weight,round(colMeans(cols))))
  }
}

## lets take a look
state.dep
