###           Functions to simulate population dynamics          ###

#Reproduction function
##########################
#This function produces pups for each specfic pack based on a mean litter size. 
#Eventually I could also make the number of pups born per pack be a function of 
#pack attributes. The function creates a new data frame for the pack pups and 
#then adds it to the wolf population.

reproduction<-function(wolf.pop, litter.size){
  n.pups<-matrix(NA, nrow=length(unique(wolf.pop$pack)), ncol=2)
  colnames(n.pups)<-c("pups", "pack")
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
#and applys. It also takes the argument harvest rate, which for now is just the 
#percentage of the population that is removed. 

harvest<-function(wolf.pop, harv.rate){
  wolf.pop<-wolf.pop %>% ungroup() %>% 
    mutate(harvest=rbinom(nrow(wolf.pop), 1, harv.rate))
}



#Transition and recruitment
#####################################
recruit<-function(wolf.pop){
  wolf.pop<-wolf.pop[wolf.pop$harvest==0,]
  wolf.pop$stage<-ifelse
}



#wolf.pop<-wolf.pop %>% group_by(pack)
#if(wolf.pop$harvest[wolf.pop$harvest==1] && grep("B", wolf.pop$stage))
  #{
  #wolf.pop$harvest<-1
#}











require(dplyr)
require(tidyr)

#Number of packs
n.packs<-20

#Determine pack size from mean pack size. Could also use negative binomial or 
#lognormal
pack.size<-rpois(n.packs, 5)

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







pop.grow<-function(n.start, harv){
  
}