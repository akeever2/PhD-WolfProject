#Fixed parameter values for population dynamics
litter.size<-4
pup.surv<-0.6
nb.surv<-0.7
b.surv<-0.7
prop.stay<-0.7

#Empty vector to store infomation
N.t2=c()

#Function to grow the population. It requires N, a vector of 3 numbers with the first
#being pups, the second being non-breeders, and the third being breeders, harvest rate,
#litter size, survival rates for pups, non-breeders, and breeders, and the proportion
#of non-breeders that remain as non-breeders. The function returns a vector of pop
#size in order and the number harvested
pop.growth<-function(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv){
  N.t2[1]<-round(N[3]*litter.size*0.5*(1-harv))
  N.t2[2]<-round(N[1]*pup.surv*(1-harv)+N[2]*nb.surv*prop.stay*(1-harv))
  N.t2[3]<-round(N[2]*nb.surv*(1-prop.stay)*(1-harv)+N[3]*b.surv*(1-harv))
  harvest<-sum(N)*harv
    
  return(c(N.t2, harvest))
}

pop.growth(N=c(400, 250, 150), harv=0.2, litter.size, pup.surv, nb.surv, prop.stay, b.surv)


#The function to feed to the optimization methods. Only takes harvest rate and uses
#the pop.growth function. Returns the negative total cumulative harvest

########   WHY DOES IT RETURN THE NEGATIVE TOTAL HARVEST????????    ########

best.policy<-function(harv){
  N<-c(400, 250, 150)
  total.harvest<-0
  for(i in 1:10){
    fut.pop<-pop.growth(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv)
    total.harvest<-total.harvest+fut.pop[4]
    N<-fut.pop[1:3]
  }
  return(-total.harvest)
}

best.policy(0.2)


#The R optim function
fit1<- optim(
  par=c(0.5),        
  fn=best.policy,        
  method="L-BFGS-B", 
  lower=c(0),  
  upper=c(1))   
fit1$par



#The genetic algorithm optimization method
library(rgenoud)
min.max<-t(c(0,1)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = best.policy,             
  nvars = 1,                 
  max = FALSE,               
  pop.size = 1000,           
  starting.values=c(0.5),   
  Domains = min.max,        
  boundary.enforcement = 2) 

xx$par





######################################################################################

##       Now I am going to add in livestock loss as a simple linear model     ##

#####################################################################################


#Fixed parameter values for population dynamics for wolves and livestock
litter.size<-4
pup.surv<-0.6
nb.surv<-0.7
b.surv<-0.7
prop.stay<-0.7

stock.r<-0.3
stock.k<-50

#Empty vector to store infomation for wolves
N.t2=c()
N.stock.t2=c()

#Function to grow the population of wolves and stock. It requires N, a vector of 3 numbers with the first
#being pups, the second being non-breeders, and the third being breeders, harvest rate,
#litter size, survival rates for pups, non-breeders, and breeders, and the proportion
#of non-breeders that remain as non-breeders. The function returns a vector of pop
#size in order and the number harvested
pop.growth<-function(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock, stock.r, stock.k){
  loss<-(sum(N)/40000)*N.stock
  N.stock.t2<-N.stock+N.stock*stock.r*(1-N.stock/stock.k)-loss
  N.t2[1]<-round(N[3]*litter.size*0.5*(1-harv))
  N.t2[2]<-round(N[1]*pup.surv*(1-harv)+N[2]*nb.surv*prop.stay*(1-harv))
  N.t2[3]<-round(N[2]*nb.surv*(1-prop.stay)*(1-harv)+N[3]*b.surv*(1-harv))
  harvest<-sum(N)*harv
  
  return(c(N.t2, harvest, N.stock.t2, loss))
}

pop.growth(N=c(400, 250, 150), harv=0.2, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock=30, stock.r, stock.k)


#The function to feed to the optimization methods. Only takes harvest rate and uses
#the pop.growth function. Returns the negative total cumulative harvest

########   WHY DOES IT RETURN THE NEGATIVE TOTAL HARVEST????????    ########

best.policy<-function(harv){
  N<-c(400, 250, 150)
  N.stock<-30
  total.harvest<-0
  avg.loss<-0
  for(i in 1:10){
    fut.pop<-pop.growth(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock, stock.r, stock.k)
    total.harvest<-total.harvest+fut.pop[4]
    N<-fut.pop[1:3]
    N.stock<-fut.pop[5]
    avg.loss<-(avg.loss+fut.pop[6])/i
  }
  stock.util<-ifelse(avg.loss>0.5, 0, -2*avg.loss+1)
  harv.util<-0.001*total.harvest
  util<-stock.util*0.6+harv.util*0.4
  
  return(-util)
}

best.policy(0.2)


#The R optim function
fit1<- optim(
  par=c(0.5),        
  fn=best.policy,        
  method="L-BFGS-B", 
  lower=c(0),  
  upper=c(1))   
fit1$par



#The genetic algorithm optimization method
library(rgenoud)
min.max<-t(c(0,1)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = best.policy,             
  nvars = 1,                 
  max = FALSE,               
  pop.size = 1000,           
  starting.values=c(0.5),   
  Domains = min.max,        
  boundary.enforcement = 2) 

xx$par




#hydroPSO function name and package name
library(hydroPSO)

hydroPSO(fn=best.policy, 
         lower=0, 
         upper=1 
         )






######################################################################################

##       Now I am going to add stochasticity in the model     ##

#####################################################################################


#Fixed parameter values for population dynamics for wolves and livestock
litter.size<-4
prop.stay<-0.7

#for stochasticy, the alpha and beta shape parameters for a beta distribution
pup.surv<-c(28.2, 18.8)
nb.surv<-c(146.3, 62.7)
b.surv<-c(146.3, 62.7)

stock.r<-0.3
stock.k<-50

#Empty vector to store infomation for wolves
N.t2=c()
N.stock.t2=c()

#Function to grow the population of wolves and stock. It requires N, a vector of 3 numbers with the first
#being pups, the second being non-breeders, and the third being breeders, harvest rate,
#litter size, survival rates for pups, non-breeders, and breeders, and the proportion
#of non-breeders that remain as non-breeders. The function returns a vector of pop
#size in order and the number harvested
pop.growth<-function(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock, stock.r, stock.k){
  loss<-(sum(N)/40000)*N.stock
  N.stock.t2<-N.stock+N.stock*stock.r*(1-N.stock/stock.k)-loss
  N.t2[1]<-round(sum(rpois(N[3], litter.size))*0.5*(1-harv))
  N.t2[2]<-rbinom(1, N[1], rbeta(1, pup.surv[1], pup.surv[2])*(1-harv))+rbinom(1, N[2], rbeta(1, nb.surv[1], nb.surv[2])*prop.stay*(1-harv))
  N.t2[3]<-rbinom(1, N[2], rbeta(1, nb.surv[1], nb.surv[2])*(1-prop.stay)*(1-harv))+rbinom(1, N[3], rbeta(1, b.surv[1], b.surv[2])*(1-harv))
  harvest<-sum(N)*harv
  
  return(c(N.t2, harvest, N.stock.t2, loss))
}

pop.growth(N=c(400, 250, 150), harv=0.2, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock=30, stock.r, stock.k)


#The function to feed to the optimization methods. Only takes harvest rate and uses
#the pop.growth function. Returns the negative total cumulative harvest

best.policy<-function(harv){
  N<-c(400, 250, 150)
  N.stock<-30
  total.harvest<-0
  avg.loss<-0
  for(i in 1:10){
    fut.pop<-pop.growth(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock, stock.r, stock.k)
    total.harvest<-total.harvest+fut.pop[4]
    N<-fut.pop[1:3]
    N.stock<-fut.pop[5]
    avg.loss<-(avg.loss+fut.pop[6])/i
  }
  stock.util<-ifelse(avg.loss>0.5, 0, -2*avg.loss+1)
  harv.util<-0.001*total.harvest
  util<-stock.util*0.6+harv.util*0.4
  
  return(-util)
}

best.policy(0.2)


#The genetic algorithm optimization method
library(rgenoud)
min.max<-t(c(0,1)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = best.policy,             
  nvars = 1,                 
  max = FALSE,               
  pop.size = 1000,           
  starting.values=c(0.5),   
  Domains = min.max,        
  boundary.enforcement = 2) 

xx$par






######################################################################################

##       Now I am going to make it state dependent and dynamic     ##

#####################################################################################


#Fixed parameter values for population dynamics for wolves and livestock
litter.size<-4
prop.stay<-0.7

#for stochasticy, the alpha and beta shape parameters for a beta distribution
pup.surv<-c(28.2, 18.8)
nb.surv<-c(146.3, 62.7)
b.surv<-c(146.3, 62.7)

stock.r<-0.3
stock.k<-50

#Empty vector to store infomation for wolves
N.t2=c()
N.stock.t2=c()

#Function to grow the population of wolves and stock. It requires N, a vector of 3 numbers with the first
#being pups, the second being non-breeders, and the third being breeders, harvest rate,
#litter size, survival rates for pups, non-breeders, and breeders, and the proportion
#of non-breeders that remain as non-breeders. The function returns a vector of pop
#size in order and the number harvested
pop.growth<-function(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock, stock.r, stock.k){
  loss<-(sum(N)/40000)*N.stock
  N.stock.t2<-N.stock+N.stock*stock.r*(1-N.stock/stock.k)-loss
  N.t2[1]<-round(sum(rpois(N[3], litter.size))*0.5*(1-harv))
  N.t2[2]<-rbinom(1, N[1], rbeta(1, pup.surv[1], pup.surv[2])*(1-harv))+rbinom(1, N[2], rbeta(1, nb.surv[1], nb.surv[2])*prop.stay*(1-harv))
  N.t2[3]<-rbinom(1, N[2], rbeta(1, nb.surv[1], nb.surv[2])*(1-prop.stay)*(1-harv))+rbinom(1, N[3], rbeta(1, b.surv[1], b.surv[2])*(1-harv))
  harvest<-sum(N)*harv
  
  return(c(N.t2, harvest, N.stock.t2, loss))
}

pop.growth(N=c(400, 250, 150), harv=0.2, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock=30, stock.r, stock.k)


#The function to feed to the optimization methods. Only takes harvest rate and uses
#the pop.growth function. Returns the negative total cumulative harvest

best.policy<-function(harv){
  N<-initial.N
  N.stock<-30
  total.harvest<-0
  avg.loss<-0
  for(i in 1:10){
    fut.pop<-pop.growth(N, harv, litter.size, pup.surv, nb.surv, prop.stay, b.surv, N.stock, stock.r, stock.k)
    total.harvest<-total.harvest+fut.pop[4]
    N<-fut.pop[1:3]
    N.stock<-fut.pop[5]
    avg.loss<-(avg.loss+fut.pop[6])/i
  }
  stock.util<-ifelse(avg.loss>0.5, 0, -2*avg.loss+1)
  harv.util<-0.001*total.harvest
  util<-stock.util*0.6+harv.util*0.4
  
  return(-util)
}

best.policy(0.2)


#The genetic algorithm optimization method
library(rgenoud)
min.max<-t(c(0,1)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = best.policy,             
  nvars = 1,                 
  max = FALSE,               
  pop.size = 1000,           
  starting.values=c(0.5),   
  Domains = min.max,        
  boundary.enforcement = 2) 

xx$par



