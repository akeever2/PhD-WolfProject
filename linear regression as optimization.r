## lets try to estimate the parameters of a simple linear regression using 
## two optimization routines
## First create a linear combination of variables.
x<-rnorm(100,10,4)
y<-rnorm(100,-3 + 0.5*x,3)

## if we had to estimate a slope and intercept from simple linear
## regression by hand (or using a spreadsheet), what would be do?
## use least squares regression maybe

##################################################################
## create a function that calculates sum of squares
## given some vector of x and y variables and 2 parameters

ss<-function(x,y,parm){
  ## sum of squares
  a<-sum((y - (parm[1] + parm[2]*x))^2)
  
  return(a)
}

## try different combinations of intercept and slope
## see how the sum of squares changes, the intercet and slope values
## that results in the smallest sum of squares is our MLE
ss(x,y,c(-0.5,1.5))

###### How de we find the minimum usning the ss function?
####### We could be dim and try multiple combinations till we found the minimum
## lets create a matrix that has all combos from -5 to 5 by .5
trys<-merge(seq(-5,5,by=0.5),seq(-5,5,by=0.5))

ss.val<-NULL
for(i in 1:nrow(trys)){ss.val<-c(ss.val,ss(x,y,trys[i,]))}

## which row contains the minimum?
which(ss.val == min(ss.val))
## what are the slope and intercept for the row 
trys[253,]
## or we could do it this way
trys[ss.val == min(ss.val),]


## Optimization routines require functions that
## input only those parameters that it needs to find
## here we create a (w)rapper function- Yo! That reduces the
## input to a single object (here a vector with two objects, the
## intercept and slope)
## wrapper for function
regress<-function(parm){
  z<-ss(x,y,parm)
  return(z)
}

# first we'll use the optim function 
# with the bounded quasi newton
fit1<- optim(
    par=c(1,1),        ## initial values of intercept,slope
    fn=regress,        ## name of function
    method="L-BFGS-B", ## specifies bounded quasi newton
    lower=c(-15,-15),  ## lower bound of intercept,slope
    upper=c(15,15))    ## uppes bound of intercept,slope
fit1

# MLE are in fit1$par, lets plug into regress function
regress(fit$par)

# then we'll use the optim function 
# with simulated annealing
fit2<- optim(
  par=c(-2,-5),
  fn=regress,
  method="SANN")
fit2


## next we'll try the genetic algorithm
library(rgenoud)
min.max<-cbind(c(-15,-15),c(15,15)) ## lower and upper bound of intercept,slope

xx<- genoud(
  fn = regress,              ### name of function
  nvars = 2,                 ### number of parameters to estimate
  max = FALSE,               ### find the minimum
  pop.size = 1000,           ### population size (number of searchers)
  starting.values=c(5,5),    ### initial values of intercept,slope
  Domains = min.max,         ### matrix of boundary values created above
  boundary.enforcement = 2)  ### specifies stay within upper and lower bound
xx


### now lets try GLM
glm(y~x)

### compare all 
fit1$par; fit2$par;  xx$par


