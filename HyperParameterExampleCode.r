it <- 100 #Number of iterations or simulations
yrs <- 10 #Number of years you run each simulation for

#First I made some empty matrices to store everything
female.harvest.beta <- matrix(,it,yrs) #creates an empty matrix for the female harvest rate using the beta distribution
female.harvest.unif <- matrix(,it,yrs) #creates an empty matrix for the female harvest rate using the uniform distribution

#Then I defined the hyper-parameter distributions. This will be the distribution to draw from and determines the value for each iteration/simulation

##Uniform distribution example
#Create an empty matrix to hold the values for each iteration
harvest.iter <-matrix(,it,1)
for (i in 1:it){
#Determine the mean value for each iteration
harvest.iter[i] <- runif(1, 0.10, 0.20)
	for (t in 1:yrs){
	#Determine the value to be used for each year withing an iteration by adding or subtracting 0.02 from the value for each iteration. You can add variation in here however you see fit, this just seemed the easiest way. 
	female.harvest.unif[i,t] <- runif(1, harvest.iter[i]-0.02, harvest.iter[i]+0.02)
	}
}

#Plot it
LC.unif <- apply(female.harvest.unif, 2, quantile, probs=0.025, na.rm=TRUE)
UC.unif <- apply(female.harvest.unif, 2, quantile, probs=0.975, na.rm=TRUE)
plot(1:yrs, colMeans(female.harvest.unif, na.rm=TRUE), xlab="Year", ylab="Mean Female Harvest Rate/ Iteration", main="Female Harvest from Uniform Distribution", ylim=c(0,1))
arrows(1:yrs, colMeans(female.harvest.unif, na.rm=TRUE), 1:yrs, LC.unif, length=0.05, angle=90)
arrows(1:yrs, colMeans(female.harvest.unif, na.rm=TRUE), 1:yrs, UC.unif, length=0.05, angle=90)








##Beta distribution example
#But first, a function to calculate the alpha a beta shape parameters from the mean and variance and display a plot for the beta distribution using the Method of Moments calculation
beta.MoM.fcn <- function (mean, var){
  alpha <- mean*((mean*(1-mean)/var)-1)
  beta <- (1-mean)*((mean*(1-mean)/var)-1)
  y <- rbeta(1:1000, alpha, beta)
  beta.plot <- plot(y, dbeta(y, alpha, beta), xlab="Female Harvest", ylab="Frequency")
  return(c(list(alpha=alpha, beta=beta), beta.plot))
  }
shape.params <- beta.MoM.fcn(0.15, 0.001) #For the function you put the mean value, here I used 0.15 for female harvest, and the variance which I played around with until I got a good looking distribution. The function stores the alpha and beta shape parameters so you can access them later. 

#I used the beta distribution to determine the mean for each iteration/simulation. 
fhb.mean <- rbeta (1:it, shape.params$alpha, shape.params$beta)
fhb.var <- 0.001  #This is the variance the iteration loop will use to create the next distribution which will determine values for each year within an iteration/simulation. You can change this to which ever value you want, the smaller the number the tighter the distribution

#Now that we have a mean and variance for each iteration we need to create empty matrices for the alpha and beta shape parameters for each iteration which creates a new distribution for each iteration that will be used to determine the value each year.
alpha.iter <- matrix(,it,1)
beta.iter <- matrix(,it,1) 

for (i in 1:it){
#Calculate the alpha and beta shape parameters for each iteration from the mean and variance using the Method of Moments calculation. 
alpha.iter[i] <- fhb.mean[i]*((fhb.mean[i]*(1-fhb.mean[i])/fhb.var)-1)
beta.iter[i] <- (1-fhb.mean[i])*((fhb.mean[i]*(1-fhb.mean[i])/fhb.var)-1)

	for (t in 1:yrs){
	#Calculate female harvest rate each year within an iteration
	female.harvest.beta[i,t] <- rbeta (1, alpha.iter, beta.iter)
	}

}

#Plot female harvest rate
LC.beta <- apply(female.harvest.beta, 2, quantile, probs=0.025, na.rm=TRUE)
UC.beta <- apply(female.harvest.beta, 2, quantile, probs=0.975, na.rm=TRUE)
plot(1:yrs, colMeans(female.harvest.beta, na.rm=TRUE), xlab="Year", ylab="Mean Female Harvest Rate/ Iteration", main="Female Harvest from Beta Distribution", ylim=c(0,1))
arrows(1:yrs, colMeans(female.harvest.beta, na.rm=TRUE), 1:yrs, LC.beta, length=0.05, angle=90)
arrows(1:yrs, colMeans(female.harvest.beta, na.rm=TRUE), 1:yrs, UC.beta, length=0.05, angle=90)



