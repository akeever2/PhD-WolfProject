
## maximum population growth rate
rmax<-0.3

## carrying capacity
K<-125
#population dynamics function
N.t<- 100
tada<-c()
for(H.rate in 0:90*0.01){
T.harvest<-0

for(i in 1:10){
N.t2<-N.t+rmax*N.t*(1-N.t/K)

harvest<- N.t2*H.rate  

N.t <- N.t2 - harvest
 
T.harvest<- T.harvest + harvest
}
tada<- rbind(tada,c(H.rate,T.harvest))
}

plot(tada[,2]~tada[,1])