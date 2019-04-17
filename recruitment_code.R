################-The effects of HC on recruitment-##################

###b

HC<-seq(0,1, by=0.02) #harvest rate
prop.pup<-0.45         #proportion of population that is pups
prop.breed<-0.4      #Proportion of adult population that are breeders
c<-5                 #Some constant that determines steepness of slope
breed<-1              #Background proportion of packs that breed

#b for feb-mar which is determined by harvest rate (HC) * the proportion 
#of populaiton that are adults (1-prop.pup) * the proportion of adults
#that are breeders (prop.breed)* the slope constant (c) + the intercept
#proportion of packs that breed (breed)
b.feb<--HC*(1-prop.pup)*prop.breed*c+breed

#b for apr-jan which is just the background proportion of packs that breed
b.ap<-rep(breed, length(HC))

plot(HC, b.feb, type="l", ylim=c(0,1), ylab="Probability a pack successfully breeds and reporudces (b)", xlab="Human-caused mortality")
lines(HC, b.ap, lty=2)
legend("bottomleft", c("Feb-Mar", "Apr-Jan"), lty=c(1,2))


###litter size

litter.size<-5
l=rep(litter.size, length(HC))
plot(HC, l, type="l", ylim=c(0,7), ylab="Litter size (l)", xlab="Human-caused mortality")

###pup survival

s.p<-0.8      #background pup survival

#direct effects of harvest on pup survival. It is the proportion of the 
#population that are pups (prop.pup)*harvest rate*slope constant (h)+ the
#intercept backgroup pup survival (s.p)
s.p.direct<-s.p*(1-(HC*prop.pup))
plot(HC, s.p.direct, type="l", ylim=c(0,1))

a<-1.25       #Some constant that determines slope

#Indirect effects of HC on pup survival. It is some constant that
#determines slope (a) *harvest rate^2 + intercept background survivial
s.p.indirect<--a*(1-prop.pup)*(HC)^2+s.p
plot(HC, s.p.indirect, type="l", ylim=c(0,1))


plot(HC, s.p.indirect, type="l", ylim=c(0,1), xlab="Human-caused mortality", ylab="S0")
lines(HC, s.p.direct, lty=2)
legend("topright", c("Indirect", "Direct"), lty=c(1,2))

#geometric mean to determine total survival
#pup.surv<-sqrt(s.p.direct*s.p.indirect)
#pup.surv<-sqrt(s.p*(s.p.direct+s.p.indirect-1))
#pup.surv<-sqrt(s.p*s.p.direct)+sqrt(s.p*s.p.indirect)-s.p
pup.surv<-s.p.direct+s.p.indirect-s.p
plot(HC, pup.surv, type="l", ylim=c(0,1), lwd=2)
lines(HC, s.p.direct, lty=1)
lines(HC, s.p.indirect, lty=2)
legend("bottomleft", c("Total", "Direct", "Indirect"), lwd=c(2,1,1), lty=c(1,1,2))


############
#recruitment
rec.feb<-pup.surv*litter.size*b.feb
rec.ap<-pup.surv*litter.size*b.ap
wgtrc<-rec.feb*0.25+rec.ap*0.75
plot(HC, rec.feb, type="l", ylim=c(0,5),xlab="Human-caused mortality", ylab="Recruitment")
lines(HC, rec.ap, lty=2)
lines(HC, wgtrc, lwd=2)
legend("topright", c("Feb-Mar", "Apr-Jan"), lty=c(1,2))

par(mfrow=c(2,2), mar=c(2,4,0,0), oma=c(1,0,1,1))
plot(HC, b.feb, type="l", ylim=c(0,1), xlab="",ylab="b", xaxt="n")
lines(HC, b.ap, lty=2)
legend("bottomleft", c("Feb-Mar", "Apr-Jan"), lty=c(1,2))
plot(HC, l, type="l", ylim=c(0,7), xlab="", ylab="l", xaxt="n")
plot(HC, pup.surv, type="l", ylim=c(0,1), lwd=2, ylab="So", xlab="Human-caused mortality")
lines(HC, s.p.direct, lty=2)
lines(HC, s.p.indirect, lty=3)
legend("bottomleft", c("Total", "Direct", "Indirect"), lwd=c(2,1,1), lty=c(1,2,3))
plot(HC, rec.feb, type="l", ylim=c(0,5), xlab="Human-caused mortality", ylab="Recruitment")
lines(HC, rec.ap, lty=2)
lines(HC, wgtrc, lwd=2)
legend("topright", c("Feb-Mar", "Apr-Jan", "Weighted"), lty=c(1,2, 1), lwd=c(1,1,2))



################-The effects of density on recruitment-##################


###b
dens<-seq(0,40, by=0.5)
b<-rep(1, length(dens))
plot(dens, b, type="l")

###Litter size
m<--0.03
litter.size<-5

l<-m*dens+litter.size
plot(dens, l, type="l", ylim=c(0,10))

###Pup survival
k<-2
s.p<-0.8
j<--.0003

pup.surv<-j*(dens-k)^2+s.p
plot(dens, pup.surv, type="l", ylim=c(0,1))

#############
#recruitment
rec<-b*pup.surv*l
plot(dens, rec, type="l", ylim=c(0,5))


################-The effects of group size on recruitment-##################


gs<-seq(2,20, by=1)
###b
b<-rep(1, length(gs))

###Litter size
l<-rep(5, length(gs))


###Pup survival
h<-15
s.p<-0.8
a<--.0025

pup.surv<-a*(gs-h)^2+s.p
plot(gs, pup.surv, type="l", ylim=c(0,1))


#############
#recruitment
rec<-b*pup.surv*litter.size
plot(gs, rec, type="l", ylim=c(0,5))



################-The effects of prey on recruitment-##################

prey<-seq(0,50, by=1)
###b
b<-rep(1, length(prey))

###Litter size
l<-rep(5, length(prey))


###Pup survival
h<-15
s.p<-0.8
a<--.0025

pup.surv<-a*(gs-h)^2+s.p
plot(gs, pup.surv, type="l", ylim=c(0,1))


#############
#recruitment
rec<-b*pup.surv*litter.size
plot(gs, rec, type="l", ylim=c(0,5))




