beta.MoM.fcn <- function (mean, var, param){
  alpha <-mean*((mean*(1-mean)/var)-1)
  beta <-(1-mean)*((mean*(1-mean)/var)-1)
  y <-rbeta(1:1000, alpha, beta)
  beta.plot <-plot(y, dbeta(y, alpha, beta), ylab="Frequency", xlab=param)
  return(c(list(alpha=alpha, beta=beta), beta.plot))
  }
beta.MoM.fcn(0.15, 0.001, "Survival")


lognorm.MoM.fcn <- function (mean, var, param){
  alpha <-log((mean^2)/sqrt(var+(mean^2)))
  beta <-sqrt(log(1+(var/(mean^2))))
  y <-rlnorm(1:1000, alpha, beta)
  beta.plot <-plot(y, dlnorm(y, alpha, beta), ylab="Frequency", xlab=param)
  return(c(list(alpha=alpha, beta=beta), beta.plot))
}