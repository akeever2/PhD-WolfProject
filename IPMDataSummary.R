library(dplyr)
library(tidyverse)


#### Bring in data and set it up ####
# Pull in group count data and pack level covariates
group <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/GroupCountData.csv")
packcovs <- read.csv("C:/Users/allison/Documents/Project/WolfData/GroupCount/FinalPackCovs.csv")


# Pull in collar data for survival model
y.surv <- read.csv("C:/Users/allison/Documents/Project/WolfData/Keever_SurvivalData/ysurv_subset2.csv")

# For the group count data, pull out only "good" (G) and "moderate" (M)
# quality counts and call it group.g
group.g <- group[group$Count.Quality=="G" | group$Count.Quality=="M", ]


# Make a new dataframe, g2,  with only the year, pack, recreation area, 
# and count of groups from the good and moderate counts
g2 <- data.frame(year=group.g$YEAR, pack=group.g$Revised.Pack.Name, area=group.g$RECAREA,  
                 count=group.g$CountPlusRemovals)


# Change the data from long format to wide format
g2 <-g2 %>% spread(year, count)


# Some packs were only around prior to 2007, so get rid of packs/group 
# counts that do not have any counts in years 2007-2016. Double check 
# that the specified columns correspond to the correct years of data.
g3 <- g2[as.logical((rowSums(is.na(g2[,22:31]))-nyears)),]


# The final group count data only needs to be from 2007-2016, so 
# double check that the specified columns correspond to the correct
# years of data and use only those columns. Final group count data 
# frame where the number of packs is the number of rows and the 
# columns are the counts per year. Then set the number of groups.
# Then set the number of groups
y.group <- g3[,22:31]
ngroups <- nrow(y.group)


# Set up group covariate data to only include packs in the group counts
# and be in the correct order. 
groupcovs <- merge(g3[,1:2], packcovs, by.x="pack", by.y="Pack", sort=FALSE)


# Now get rid of all the other random dataframes that aren't needed
rm(g2, g3, group.g)


# The survival data currently includes 2007-2017, however I only need until 2016
# (9 years {nyears - 1}) worth of data so I need to remove the data that I do not
# need
y.surv <- y.surv[y.surv$year < 2017,]


#### Now summarize survival data for results ####

# Number of collared wolves
length(unique(y.surv$WOLFNO))

# Years of observations
unique(y.surv$year)

# Get rid of duplicated wolves for simplicity, but keep track of different years
y.surv2 <- y.surv %>% distinct(WOLFNO, year, .keep_all=TRUE)

# Number of collared wolves by year
group_by(y.surv2, year) %>% summarise(avg=mean(n()), sd=sd(n()))

# Make dataframe for just unique wolves with no replications for year
y.surv3 <- y.surv %>% distinct(WOLFNO, .keep_all=TRUE)

# Get the number of unique causes of death
COD<-group_by(y.surv3, COD) %>% summarise(count=n())
COD<-COD[order(COD$count),]

# Get the number of days wolves survived
y.surv3$days<-as.Date(y.surv3$stop, format = "%m/%d/%Y")-as.Date(y.surv3$start, format = "%m/%d/%Y")
mean(y.surv3$days)
sd(y.surv3$days)


# % in yearling and adult stage
group_by(y.surv3, CAPTURE_AGE_CLASS) %>% summarise(n=n())

# Number of packs represented and number from same pack
group_by(y.surv3, PACK) %>% summarise(n=n()) %>% summarise(avg=mean(n), sd=mean(n), max=max(n), min=min(n))

# Regions represented
group_by(y.surv3, Region) %>% summarise(n=n())

# Sex
group_by(y.surv3, SEX) %>% summarise(n=n())



#### Now group data summaries ####

# Make a group count dataframe with the group name included
y.group2 <- g3[,c(1,2,22:31)]

# Alter other group dataframes to only be the years needed
group2<-group[group$YEAR > 2006,]
group.g2<-group.g[group.g$YEAR > 2006,]

# Double check number of observations
length(which(!is.na(y.group)))

# Observations per year
group_by(group.g2, YEAR) %>% summarise(N=n()) %>% summarise(avg=mean(N), sd=sd(N))

# Obs from each pack
group_by(group.g2, Revised.Pack.Name) %>% summarise(N=n()) %>% summarise(avg=mean(N), sd=sd(N), max=max(N), min=min(N))

# Average pack size total
group_by(group.g2, YEAR) %>% summarise(avg=mean(EoY.Count), sd=sd(EoY.Count))
mean(group.g2$EoY.Count)
sd(group.g2$EoY.Count)

# Average pack size durning non harvest
g.noharv <- group.g2[group.g2$YEAR == 2007 | 
                       group.g2$YEAR == 2008 |
                       group.g2$YEAR == 2010,]
mean(g.noharv$EoY.Count)
sd(g.noharv$EoY.Count)

# Average pack size during harvest
g.harv <- group.g2[group.g2$YEAR == 2009 | 
                       group.g2$YEAR == 2011 |
                       group.g2$YEAR == 2012 | 
                     group.g2$YEAR == 2013 |
                     group.g2$YEAR == 2014 | 
                     group.g2$YEAR == 2015 |
                     group.g2$YEAR == 2016,]
mean(g.harv$EoY.Count)
sd(g.harv$EoY.Count)


