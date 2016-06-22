
#' we are now ready to try to find the optimal solution using 
#' something other than brute force
#' here the harvest function as before

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
### now we need to build a function to use with an optimizer


 
find.best<-function(H.rate){
  N.t<- 100
  T.harvest<-0
  for(i in 1:10){
    ret.n.fut<-popn.harv(N.t,H.rate)
    ## add yearly harvest to total
    T.harvest<- T.harvest + ret.n.fut[1]
    ## adult tribbles surviving to reproduce
    N.t<-ret.n.fut[2]
  }
  ## notice the negative sign
  return(-T.harvest)
}

find.best(.1)

# first we'll use the optim function 
# with the bounded quasi newton
fit1<- optim(
  par=c(0.5),        
  fn=find.best,        
  method="L-BFGS-B", 
  lower=c(0),  
  upper=c(1))   
fit1$par


## next we'll try the genetic algorithm
library(rgenoud)
min.max<-t(c(0,1)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = find.best,             
  nvars = 1,                 
  max = FALSE,               
  pop.size = 1000,           
  starting.values=c(0.5),   
  Domains = min.max,        
  boundary.enforcement = 2) 

xx$par

#######################################################
#######################################################
##' often managers implement a decision but the decision 
##' doesn't quite hit the mark. for example, managers may intend 
##' to implement a 30% tribble harvest rate, but due to
##' unforseen cirmcumstances the true harvest rate is 10%
##' this is known as ..... come on..... i'm waiting...

#' to incorporate this source of uncertainty, lets consider 
#' 4 levels of harvest 10%, 20%, 30%, 40% each with CV 30% 

find.best<-function(decis){
  ## matrix each row is alpha and beta of beta distribution
  ## for each decision for 4 harvest levels and cv 30%
  harv=matrix(c(9.9, 89.1, 8.69,34.76, 7.48, 
                17.45,6.27,9.4),ncol=2,byrow = T)
  N.t<- 100
  T.harvest<-0
  for(i in 1:10){
    ## generate random harvest for a decision each year
    H.rate<-rbeta(1,harv[decis,1],harv[decis,2])
    ret.n.fut<-popn.harv(N.t,H.rate)
    ## add yearly harvest to total
    T.harvest<- T.harvest + ret.n.fut[1]
    ## adult tribbles surviving to reproduce
    N.t<-ret.n.fut[2]
  }
  ## notice NO MORE negative sign
  return(T.harvest)
}

find.best(3)

## because the decisions are now discrete we can no longer use the optimx
## function we can use the genetic algorithm, with a little modification

min.max<-t(c(1,4)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = find.best,             
  nvars = 1,                 
  max = TRUE, 
  pop.size = 1000,
  data.type.int = TRUE,         
  starting.values=c(1),   
  Domains = min.max,        
  boundary.enforcement = 2) 

xx$par

################################################################
################################################################
#' so far the current state (size) of the population has remained
#' the same. MDP, however, assumes that the system state is known 
#' (maybe imperfectly) we'd also like to evaluate if the optimal
#' decision changes defending on the system state
#' lets change the function so that we can change the initial population


find.best<-function(decis){
  ## matrix each row is alpha and beta of beta distribution
  ## for each decision for 4 harvest levels and cv 30%
  harv=matrix(c(9.9, 89.1, 8.69,34.76, 7.48, 
                17.45,6.27,9.4),ncol=2,byrow = T)
  N.t<- initial.N
  T.harvest<-0
  for(i in 1:10){
    ## generate random harvest for a decision each year
    H.rate<-rbeta(1,harv[decis,1],harv[decis,2])
    ret.n.fut<-popn.harv(N.t,H.rate)
    ## add yearly harvest to total
    T.harvest<- T.harvest + ret.n.fut[1]
    ## adult tribbles surviving to reproduce
    N.t<-ret.n.fut[2]
  }
  ## notice NO MORE negative sign
  return(T.harvest)
}

initial.N <- 10
find.best(3)

## lets get the optimal harvest over a range of current popn
state.dep<-NULL
for(initial.N in seq(25,150,by=25)){

xx<- genoud(
  fn = find.best,             
  nvars = 1,                 
  max = TRUE, 
  pop.size = 1500,
  data.type.int = TRUE,         
  starting.values=c(1),   
  Domains = min.max,        
  boundary.enforcement = 2) 

state.dep<- rbind(state.dep,c(initial.N,xx$par))
}

## lets take a look
state.dep

#'generally takes lots of populations, optimal 
#'is actually below
#     [,1] [,2]
#[1,]   25    2
#[2,]   50    2
#[3,]   75    3
#[4,]  100    3
#[5,]  125    3
#[6,]  150    3

#'the results suggest that the current system state
#'affects tribble harvest decision making even 10 years out
#'maybe we should consider the population size each year after 
#'but we also want to make sure that we don't make a decision that
#'may be good for the current time period but really restrict 
#'future options, in this sense we would like to find the optimal
#'*series* of decisions now we will modify the function
#'to read in a vector of decisions

## decis in now a vector
find.best<-function(decis){
  ## matrix each row is alpha and beta of beta distribution
  ## for each decision for 4 harvest levels and cv 30%
  harv=matrix(c(9.9, 89.1, 8.69,34.76, 7.48, 
                17.45,6.27,9.4),ncol=2,byrow = T)
  N.t<- initial.N
  T.harvest<-0
  for(i in 1:length(decis)){
    ## generate random harvest for a decision each year
    H.rate<-rbeta(1,harv[decis[i],1],harv[decis[i],2])
    ret.n.fut<-popn.harv(N.t,H.rate)
    ## add yearly harvest to total
    T.harvest<- T.harvest + ret.n.fut[1]
    ## adult tribbles surviving to reproduce
    N.t<-ret.n.fut[2]
  }
  ## notice NO MORE negative sign
  return(T.harvest)
}

initial.N <- 10
find.best(rep(1,5))

## lets get the optimal harvest over a range of current popn
state.dep<-NULL

yrs<-10
# we need 10 upper and lower bounds
min.max<-cbind(rep(1,yrs),rep(4,yrs))

for(initial.N in seq(25,150,by=25)){
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
  state.dep<- rbind(state.dep,c(initial.N,round(colMeans(cols))))
}

## lets take a look
state.dep
###############################################################
"
[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
[1,]   25    1    1    1    2    1    2    2    2     3     4
[2,]   50    1    2    2    2    1    2    2    3     2     4
[3,]   75    2    2    2    3    1    2    3    3     3     4
[4,]  100    3    2    2    2    2    1    3    2     4     4
[5,]  125    3    2    3    2    3    2    2    3     3     4
[6,]  150    4    2    2    2    2    2    2    2     3     4"
