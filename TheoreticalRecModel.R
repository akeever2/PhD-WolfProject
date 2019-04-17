###################################################:

##        Theoretical recruitment model code

###################################################:


#### Base/null model of averages ####

library(dplyr)
library(tidyr)
library(ggplot2)


# Function to probalistically estimate recruitment
prob.rec<-function (npacks, prob.breed, litter.size, pup.survival){
  
  # Number of packs that breed
  packs.breed <- rbinom(1, npacks, prob.breed)
  
  # Number of pups produced in breeding packs
  pups.pack <- rpois(packs.breed, litter.size)
  
  # Number of pups that survive
  surv.pups <- rbinom(packs.breed, pups.pack, pup.survival)
  
  # Mean pups per pack for all packs
  mean.rec <- sum(surv.pups)/npacks
  
  # Mean pups per pack for breeding packs
  mean.rec.breed <- sum(surv.pups)/packs.breed

return(list(mean.rec=mean.rec, mean.rec.breed=mean.rec.breed))

}


# Run function once 
prob.rec(npacks=100, prob.breed=0.85, litter.size=5, pup.survival=0.75)


# Run function for 10000 iterations to get variation
simoutput <- lapply(rep(100, 1000), FUN = prob.rec, prob.breed=0.85, 
                    litter.size=5, pup.survival=0.75)
rec.sim.datum <- as.data.frame(bind_rows(simoutput))

mean(rec.sim.datum[,1])
var(rec.sim.datum[,1])

mean(rec.sim.datum[,2])
var(rec.sim.datum[,2])


# Plot the mean and variation with group size varying
rec.sim.datum$group.size <- rep(seq(1,20, by=1), nrow(rec.sim.datum)/20)
rec.sim.datum$group.rec.mean <- rep(mean(rec.sim.datum$mean.rec), nrow(rec.sim.datum))

mean.dat <- rec.sim.datum %>% group_by(group.size) %>% summarise(group.mean=mean(mean.rec))
rec.sim.datum <- rec.sim.datum %>% nest(-group.size) %>% 
  mutate(group.mean=mean.dat$group.mean) %>% unnest

ggplot(data=rec.sim.datum, aes(x=group.size, y=mean.rec))+
  geom_point(alpha=0.3, color="darkgray")+
  geom_point(y=rec.sim.datum$group.mean, size=4)+
  scale_x_continuous(name="Group Size", breaks=c(2,4,6,8,10,12,14,16,18,20))+
  scale_y_continuous(name="Pups per Pack", limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text=element_text(size=12), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.title=element_blank(),  
        legend.key=element_blank())

rec.non.prob <- 0.85 * 5 * 0.75



#### H: Pup survival increases linearly with pack size ####

# Explore relationship
# Group size
gs<-seq(2,20, by=1)

# proportion breeding
b<-rep(1, length(gs))

# Litter size
l<-rep(5, length(gs))


###Pup survival
h<-12
s.p<-0.8
a<--.0045

pup.surv<-0.45 + .02*gs
plot(gs, pup.surv, type="l", ylim=c(0,1))

#recruitment
rec<-b*pup.surv*l
plot(gs, rec, type="l", ylim=c(0,5))

# Function to probalistically estimate recruitment
prob.rec.1<-function (npacks, prob.breed, litter.size,  
                      pack.size, min.surv){
  
  # Survival of pups
  pup.survival <- min.surv + 0.02 * pack.size
  
  # Number of packs that breed
  packs.breed <- rbinom(1, npacks, prob.breed)
  
  # Number of pups produced in breeding packs
  pups.pack <- rpois(packs.breed, litter.size)
  
  # Number of pups that survive
  surv.pups <- rbinom(packs.breed, pups.pack, pup.survival)
  
  # Mean pups per pack for all packs
  mean.rec <- sum(surv.pups)/npacks
  
  # Mean pups per pack for breeding packs
  mean.rec.breed <- sum(surv.pups)/packs.breed
  
  return(list(mean.rec=mean.rec, mean.rec.breed=mean.rec.breed, pack.size=pack.size))
  
}


# Run function once 
prob.rec.1(npacks=100, prob.breed=0.85, litter.size=5, pack.size=10, min.surv=0.45)


# Run function for 10000 iterations to get variation
simoutput <- lapply(rep(seq(2, 20, by=1), 500), FUN = prob.rec.1, npacks=100, prob.breed=0.85, 
                    litter.size=5, min.surv=0.45)

rec.sim.datum.1 <- as.data.frame(bind_rows(simoutput))



# Plot the mean and variation with group size varying
mean.dat.1 <- rec.sim.datum.1 %>% group_by(pack.size) %>% summarise(group.mean=mean(mean.rec))
rec.sim.datum.1 <- rec.sim.datum.1 %>% nest(-pack.size) %>% 
  mutate(group.mean=mean.dat.1$group.mean) %>% unnest

ggplot(data=rec.sim.datum.1, aes(x=pack.size, y=mean.rec))+
  geom_point(alpha=0.3, color="darkgray")+
  geom_point(y=rec.sim.datum.1$group.mean, size=4)+
  scale_x_continuous(name="Group Size", breaks=c(2,4,6,8,10,12,14,16,18,20))+
  scale_y_continuous(name="Pups per Pack", limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text=element_text(size=12), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.title=element_blank(),  
        legend.key=element_blank())


#### H: Direct effects of harvest only ####

prob.rec.2<-function (npacks, prob.breed, litter.size, max.survival, HC){
  
  # Survival with harvest
  pup.survival <- max.surv * (1-HC)
  
  # Number of packs that breed
  packs.breed <- rbinom(1, npacks, prob.breed)
  
  # Number of pups produced in breeding packs
  pups.pack <- rpois(packs.breed, litter.size)
  
  # Number of pups that survive
  surv.pups <- rbinom(packs.breed, pups.pack, pup.survival)
  
  # Mean pups per pack for all packs
  mean.rec <- sum(surv.pups)/npacks
  
  # Mean pups per pack for breeding packs
  mean.rec.breed <- sum(surv.pups)/packs.breed
  
  return(list(mean.rec=mean.rec, mean.rec.breed=mean.rec.breed, HC=HC))
  
}


# Run function once 
prob.rec.2(npacks=100, prob.breed=0.85, litter.size=5, max.survival=0.8, HC=0.3)


# Run function for 10000 iterations to get variation
simoutput <- lapply(rep(seq(0,1,by=.05), 500), FUN = prob.rec.2, npacks=100, prob.breed=0.85, 
                    litter.size=5, max.survival=0.8)
rec.sim.datum.2 <- as.data.frame(bind_rows(simoutput))


# Plot the mean and variation with group size varying
mean.dat.2 <- rec.sim.datum.2 %>% group_by(HC) %>% summarise(HC.mean=mean(mean.rec))
rec.sim.datum.2 <- rec.sim.datum.2 %>% nest(-HC) %>% 
  mutate(HC.mean=mean.dat.2$HC.mean) %>% unnest

ggplot(data=rec.sim.datum.2, aes(x=HC, y=mean.rec))+
  geom_point(alpha=0.3, color="darkgray")+
  geom_point(y=rec.sim.datum.2$HC.mean, size=4)+
  scale_x_continuous(name="Human-caused mortality")+
  scale_y_continuous(name="Pups per Pack", limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text=element_text(size=12), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.title=element_blank(),  
        legend.key=element_blank())



#### H: Direct and indirect effects ####

prob.rec.3<-function (npacks, prob.breed, litter.size, max.survival, HC){
  
  # Survival with harvest
  direct <- max.survival * (1-HC)
  indirect <- -.95*(HC^2) + max.survival
  pup.survival <- max(0,direct + indirect - max.survival)
  
  # Number of packs that breed
  packs.breed <- rbinom(1, npacks, prob.breed)
  
  # Number of pups produced in breeding packs
  pups.pack <- rpois(packs.breed, litter.size)
  
  # Number of pups that survive
  surv.pups <- rbinom(packs.breed, pups.pack, pup.survival)
  
  # Mean pups per pack for all packs
  mean.rec <- sum(surv.pups)/npacks
  
  # Mean pups per pack for breeding packs
  mean.rec.breed <- sum(surv.pups)/packs.breed
  
  return(list(mean.rec=mean.rec, mean.rec.breed=mean.rec.breed, HC=HC))
  
}


# Run function once 
prob.rec.3(npacks=100, prob.breed=0.85, litter.size=5, max.surv=0.8, HC=0.3)


# Run function for 10000 iterations to get variation
simoutput <- lapply(rep(seq(0,1,by=.05), 500), FUN = prob.rec.3, npacks=100, prob.breed=0.85, 
                    litter.size=5, max.surv=0.8)
rec.sim.datum.3 <- as.data.frame(bind_rows(simoutput))


# Plot the mean and variation with group size varying
mean.dat.3 <- rec.sim.datum.3 %>% group_by(HC) %>% summarise(HC.mean=mean(mean.rec))
rec.sim.datum.3 <- rec.sim.datum.3 %>% nest(-HC) %>% 
  mutate(HC.mean=mean.dat.3$HC.mean) %>% unnest

ggplot(data=rec.sim.datum.3, aes(x=HC, y=mean.rec))+
  geom_point(alpha=0.3, color="darkgray")+
  geom_point(y=rec.sim.datum.3$HC.mean, size=4)+
  scale_x_continuous(name="Human-caused mortality")+
  scale_y_continuous(name="Pups per Pack", limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text=element_text(size=12), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.title=element_blank(),  
        legend.key=element_blank())



#### Combined graphing ####
  
ggplot(data=rec.sim.datum.1, aes(x=pack.size, y=mean.rec))+
  geom_point(alpha=0.3, color="darkgray")+
  geom_point(y=rec.sim.datum.1$group.mean, size=4)+
  geom_point(data=rec.sim.datum, aes(x=group.size, y=mean.rec), alpha=0.3, 
             color="lightblue")+
  geom_point(data=rec.sim.datum, aes(x=group.size, y=group.mean), size=4,  
             color="darkblue")+
  scale_x_continuous(name="Group Size", breaks=c(2,4,6,8,10,12,14,16,18,20), limits=c(2,20))+
  scale_y_continuous(name="Pups per Pack", limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.justification=c(1,0), legend.position=c(1,0),
        legend.title=element_blank(),
        legend.key=element_blank())




ggplot(data=rec.sim.datum.3, aes(x=HC, y=mean.rec))+
  geom_point(alpha=0.3, color="darkgray")+
  geom_point(y=rec.sim.datum.3$HC.mean, size=4)+
  geom_point(data=rec.sim.datum.2, aes(x=HC, y=mean.rec), alpha=0.3, color="lightblue")+
  geom_point(data=rec.sim.datum.2, aes(x=HC, y=HC.mean), size=4, color="darkblue")+
  scale_x_continuous(name="Human-caused mortality")+
  scale_y_continuous(name="Pups per Pack", limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text=element_text(size=12), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"),
        legend.title=element_blank(),  
        legend.key=element_blank())
