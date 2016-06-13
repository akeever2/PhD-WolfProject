#Code to simulate an individual-based wolf population. Since the number of packs 
#is estimated by POM, I will start with the number of packs on the landscape. 
#Using the mean pack size I will determine the number of indivudals that make up 
#each pack. This code was adopted from Josh Nowak's code. 
#############################################################################################


#First I need to create a starting IB wolf population where indiviudals are 
#characterized by individual ID, pack ID, sex, and stage (Breeder or Non-breeder). 
#######################################
require(dplyr)
require(tidyr)

#Number of packs
n.packs<-20

#Determine pack size from mean pack size. Could also use negative binomial or 
#lognormal
pack.size<-rpois(n.packs, 5.5)

#Fix the cases where there are less than 2 individuals per pack
pack.size[pack.size<2]<-2

#Determine sex of the individuals in each pack and assign pack id number
wolf.pop<-bind_rows(lapply(1:n.packs, function(i){
  data.frame(sex=c(rep("F", floor(pack.size[i]/2)), 
                   rep("M", ceiling(pack.size[i]/2))), pack=i, stringsAsFactors=F)
}))

#Add individual id and choose breeders by sex and pack
wolf.pop<- wolf.pop %>% mutate(id=1:n()) %>% group_by(pack, sex) %>% 
  mutate(stage = rmultinom(1, size=1, prob=rep(1/n(), n()))[,])
wolf.pop$stage[wolf.pop$stage==1]<-"B"
wolf.pop$stage[wolf.pop$stage==0]<-"NB"





#Now that I have the starting wolf population set up the functions I will use to 
#simulate wolf population dynamics
#####################################


#Reproduction function
##########################

#This function produces pups for each specfic pack based on a mean litter size. 
#Eventually I could also make the number of pups born per pack be a function of 
#pack attributes. The function creates a new data frame for the pack pups and 
#then adds it to the wolf population.

reproduction<-function(wolf.pop, litter.size){
  n.pups<-matrix(NA, nrow=length(unique(wolf.pop$pack)), ncol=2)
  colnames(n.pups)<-c("pups", "pack")
  #prob.breed<-rbinom(length(unique(wolf.pop$pack)), 1, prob=0.95)
  for(i in 1:length(unique(wolf.pop$pack))){
    #Determine number of pups using mean litter size and then assign the pack 
    #number
    n.pups[i,1]<-rpois(1, litter.size)
    n.pups[,2]<-unique(wolf.pop$pack)
  }
  
  #Make data frame called pack.pups that repeats the pack id for each pup in the 
  #pack, then assign stage, sex, and id number to each pup. 
  n.pups<-as.data.frame(n.pups)
  pack.pups<-n.pups[rep(row.names(n.pups), n.pups$pups), ] %>% mutate(stage="P") %>% 
    group_by(pack) %>% 
    mutate(sex=ifelse(rbinom(n(), size=1, prob=0.5)==1, "M", "F"))
  pack.pups<-pack.pups %>% ungroup() %>%
    mutate(id=seq(max(wolf.pop$id)+1,length.out=nrow(pack.pups)))
  
  #Get rid of the pups column and add the pups to the wolf population. 
  pack.pups$pups<-NULL
  wolf.pop<-bind_rows(wolf.pop, pack.pups) %>% arrange(pack, stage)
}


#Harvest function
######################

#Harvest function takes the wolf population which will now include adults and pups
#and applys randomly based on proportion of each class in the population plus
#an added vulnerability for pups to be harvested. It also takes the argument
#harvest rate, which for now is just the percentage of the population that is 
#removed. 

harvest<-function(wolf.pop, harv.rate){
  #The probability of being a breeder given you are an adult. This will be used 
  #later to determine the probability of a breeder being harvested. This assumes 
  #that breeders and non breeders have an equal probability of being harvested 
  #conditional on the proportion of breeders in the adult population.
  prob.breed<-length(wolf.pop$stage[wolf.pop$stage=="B"])/
    length(wolf.pop$stage[wolf.pop$stage!="P"])
  
  #The probability of harvesting a pup is determined by the propotion of pups in 
  #the population and a stupidity factor. Right now the stupidity factor is a 
  #constant 0.10. 
  prob.pup<-length(wolf.pop$stage[wolf.pop$stage=="P"])/
    length(wolf.pop$stage[wolf.pop$stage!="P"]) + 0.10
  
  #Harvest rate for any one stage depends on thier frequency in the population 
  #and total harvest rate. As stated above, pups have an additional stupidity
  #factor. 
  prob.harv.pup<-prob.pup*harv.rate
  prob.harv.nb<-(1-prob.pup)*(1-prob.breed)*harv.rate
  prob.harv.b<-(1-prob.pup)*prob.breed*harv.rate
  
  
  #This is making harvest probabilistic, for right now I want total control, but 
  #to randomly select which individuals get harvested. 
  
  wolf.pop<-wolf.pop %>% group_by(stage) %>% 
    mutate(harvest=rbinom(n(), size=1, 
                          prob=ifelse(wolf.pop$stage=="P", 
                                      prob.harv.pup, 
                                      ifelse(wolf.pop$stage=="NB",
                                             prob.harv.nb, prob.harv.b))))
  
}



#Transition and recruitment
#####################################


recruit<-function(wolf.pop){
  
  wolf.packs<-gather(wolf.pop, pack, id)
  
    
    
    array(NA, dim=c(1, ncol(wolf.pop), length(unique(wolf.pop$pack))))
  for(i in 1:nrow(wolf.pop)){
    for(j in 1:ncol(wolf.pop)){
      for(k in 1:unique(length(wolf.pop$pack))){
        packs[,,k]<- wolf.pop %>% filter(pack==wolf.pop$pack[,,i])
      }
    }
  }
  
}

x<-seq(0,1, by=0.02)
m<-5
b<-10
rec<--m*x+b
plot(x, rec)

a<-10
adrec<-(a+x)/(x+1)
plot(x, adrec)

for(t in 2:Time){
  N.ad[t]<-round(N.ad[t-1]*(S.ad-(1-prob.harv.pups[t-1])*HC))+
    round(N.pups[t-1]*(S.pups-prob.harv.pups[t-1]*HC))
  N.pups[t]<-round(packs*l)
  N.breed[t]<-packs*2
  N.help[t]<-N.ad[t]-N.breed[t]
  prob.harv.pups[t]<-N.pups[t-1]/(N.ad[t-1]+N.pups[t-1])+0.25
  Pb[t]<-N.breed[t]/N.ad[t]


























###########################################################################################

###Failed harvest code###
#wolf.pop<-wolf.pop %>% group_by(stage) %>% 
# mutate(harvest=0)


#wolf.pop$harvest[wolf.pop$stage=="P"][sample.int
#   (length(wolf.pop$stage[wolf.pop$stage=="P"]), 
#   length(wolf.pop$stage[wolf.pop$stage=="P"])*prob.harv.pup, 
#  prob=rep(prob.harv.pup, length(wolf.pop$stage[wolf.pop$stage=="P"])))]<-1
#wolf.pop$harvest[wolf.pop$stage=="NB"][sample.int
#   (length(wolf.pop$stage[wolf.pop$stage=="NB"]), 
#   length(wolf.pop$stage[wolf.pop$stage=="NB"])*prob.harv.nb, 
#  prob=rep(prob.harv.nb, length(wolf.pop$stage[wolf.pop$stage=="NB"])))]<-1
#wolf.pop$harvest[wolf.pop$stage=="B"][sample.int
#   (length(wolf.pop$stage[wolf.pop$stage=="B"]), 
#   length(wolf.pop$stage[wolf.pop$stage=="B"])*prob.harv.b, 
#  prob=rep(prob.harv.b, length(wolf.pop$stage[wolf.pop$stage=="B"])))]<-1









############################################################################################
#Simulate an individual-based wolf population: Josh's code
##########################################################
require(dplyr)
#  Number of individuals and packs known, pack size unknown

#  Number of individuals
nind <- 100

#  Number of packs
npacks <- 10

#  Expected mean pack size 
nind/npacks

#  Create a data.frame to hold all of the individuals and their attributes
#  Then assign each animal to a pack
#  Then assign sex with equal probability of M and F
#  Then choose breeder status if both sexes occur in the same pack
df <- data.frame(id = 1:nind) %>%
  mutate(pack = sample(1:npacks, nind, replace = T),
         sex = sample(c("M", "F"), n(), replace = T)) %>%
  group_by(pack, sex) %>%
  mutate(breeder = rmultinom(1, size = 1, prob = rep(1/n(), n()))) %>%
  #  Fix the cases where only one sex occurs in a pack
  group_by(pack) %>%
  mutate(breeder = replace(breeder, length(unique(sex)) < 2, 0),
         breeder = as.logical(breeder))

#  Pack sizes - increase npacks to reduce?  
#  Or build a recursive function that describes the probability of pack 
#   membership as a function of existing packs on the landscape and current pack membership
table(df$pack)		

#  View one pack, the third just cuz
p3 <- filter(df, pack == 3)
print(p3)

#  Summarize pack composition
table(p3$sex, p3$breeder)

#  Table breeder status by sex
table(df$sex, df$breeder)

#####################################################




#One biological argument would suggest that we should begin with 100 individuals and create pairs of breeders.  OTOH, if the population has been going for a while then we can put a number of packs on the ground, populate the packs with a mean number of individuals and then assign two of the animals breeder status.  Something like....

######################################################
npacks <- 10

#  Grow packs by mean pack size, negative binomial or lognormal would also work,
#   I am being a little lazy
pack_size <- rpois(npacks, 6)

#  Fake truncation so pack has at least two members
pack_size[pack_size < 2] <- 2

#  Populate each pack
df <- bind_rows(lapply(1:npacks, function(i){ 
  #  I assigned sex a different way here, it will make sex ratio within a pack
  #  nearly equal, use the other method for random assignment, but beware the 
  #  pack that is all one sex...I would start by assign the first 2 members 
  #  deterministically and then assign the remainder stochastically
  
  data.frame(sex = c(rep("F", floor(pack_size[i]/2)), 
                     rep("M", ceiling(pack_size[i]/2))),
             pack = i, stringsAsFactors = F) 
  
}))

#  Add individual id and then choose breeder by sex and pack 
df <- df %>% 
  mutate(id = 1:n()) %>%
  group_by(pack, sex) %>%
  mutate(breeder = rmultinom(1, size = 1, prob = rep(1/n(), n())))

#  Emphasizing the equal sex ratio comment
table(df$pack, df$sex)

#  Pack sizes
table(df$pack)

#  2 breeders in each pack
tapply(df$breeder, df$pack, function(x) sum(x) == 2)	
###########################################################################################

















x<-seq(0,1, by=0.02)
m<-5
b<-10
rec<--m*x+b
plot(x, rec)

a<-10
adrec<-(a+x)/(x+1)
plot(x, adrec)


#Time is the number of years
Time<-10

#N.breed is the number of breeders, N.help is the number of helpers, 
  #and N.pups is the number of pups. N.ad is the total number of adults. 
N.ad<-rep(NA, Time, by=1)
N.pups<-rep(NA, Time, by=1)
packs<-10
l<-5
S.ad<-0.85
S.pups<-0.6


N.ad[1]<-50
N.pups[1]<-50
N.breed[1]<-packs*2
N.help<-N.ad-N.breed
Pb<-N.breed/N.ad
prob.harv.pups<-N.pups/(N.ad+N.pups) + 0.25
HC<-0.3
  
  
for(t in 2:Time){
  N.ad[t]<-round(N.ad[t-1]*(S.ad-(1-prob.harv.pups[t-1])*HC))+
    round(N.pups[t-1]*(S.pups-prob.harv.pups[t-1]*HC))
  N.pups[t]<-round(packs*l)
  N.breed[t]<-packs*2
  N.help[t]<-N.ad[t]-N.breed[t]
  prob.harv.pups[t]<-N.pups[t-1]/(N.ad[t-1]+N.pups[t-1])+0.25
  Pb[t]<-N.breed[t]/N.ad[t]
  
}

