###########################################################################################

#                   Analyses for breeding pair estimator

###########################################################################################

# Set the working directory so all output files go to correct 
# folder
setwd("C:/Users/allison/Documents/SideProjects/BreedingPairs")


# Pull in the data (CSV file) and name it bp based on filepath
bp <- read.csv("C:/Users/allison/Documents/SideProjects/BreedingPairs/Data for UM 9-7-17.csv")


# Pull in the list of pack names that includes the FWP region
packs <- read.csv("C:/Users/allison/Documents/SideProjects/BreedingPairs/PackList.csv")


# Merge the bp dataframe and packs dataframe so that FWP region 
# is in the new bp2 dataframe
bp2 <- merge(bp, packs, sort=FALSE)


# Check out the structure of the dataframe to get a feel for it
str(bp2)


# The structure of the data shows us a few issues. First, 
# the breeding pair column is considered a factor with 3 levels, 
# "N", "Y", and "U". We want to get rid of the U's and make it 
# 1 or 0 for Y or N. We will call the variable BP and retain the
# original Breeding.Pair variable and get rid of the observations
# with U. 
bp2 <- bp2[bp2$Breeding.Pair == "N" | bp2$Breeding.Pair == "Y",]
bp2$Breeding.Pair <- factor(bp2$Breeding.Pair)
bp2$BP <- ifelse(bp2$Breeding.Pair == "N", 0, 1)


# Make sure the numbers are still correct
table(bp2$Breeding.Pair)
table(bp2$BP)


# Another issues when looking at the structure is that count 
# quality is a factor with 6 levels. It should only be 3, "G", 
# "M", and "P". We need to see what the factors are and then 
# fix the issue and get rid of unclassified ones. 
table(bp2$Count.Quality)


# The table shows us there are some classified as "CENSOR" and 
# some classified as "U", the 6th factor was an empty space. 
# We will get rid of the unnecessary factor levels. 
bp2 <- bp2[bp2$Count.Quality == "G" | bp2$Count.Quality == "M" | 
             bp2$Count.Quality == "P",]
bp2$Count.Quality <- factor(bp2$Count.Quality)
table(bp2$Count.Quality)


# Now look at distribution of the data
hist(bp2$EoY.Count)
hist(bp2$YEAR)
hist(bp2$MORT_CONTROL)
hist(bp2$FWP.Region)
hist(bp2$MORT_HARVEST)
hist(bp2$TOTAL_REMOVED)




# Add a harvest covariate dependent on year. Harvest is 1 if the year is greater than 
# 2010 or if the year is 2009, else it is 0
bp2$harvest <- factor(ifelse(bp2$YEAR > 2010 | bp2$YEAR== 2009, 1, 0))

# Make sure the harvest covariate is right
table(bp2$YEAR, bp2$harvest)

# Check the structure of the data
str(bp2)


# Subset the data for good counts only and prepare to run analyses

# Create a subset of the data
bp.good <- bp2[bp2$Count.Quality == "G",]
bp.good$Count.Quality <- factor(bp.good$Count.Quality)

# Now only packs of 4+ wolves
bp.good2 <- bp.good[bp.good$EoY.Count >=4,]



# Bring in the covariate data and merge with the BP data
covs <- read.csv("./CovData.csv")
colnames(covs) <- c("YEAR", "PopHarvMort", "PopContMort", "PopOtherMort", "PopPerHarv", 
                    "PopPerCont", "PopPerOther", "Abundance", "Lambda", "P.Abundance", 
                    "P.Lambda")
BP<-merge(bp.good2, covs)

# Create density covariate based off of estimated abundance by dividing by the max area occupied 
# by wolves from POM and then multiplying by 1000 to put into wolves per 1000 km^2
BP$Density <- BP$P.Abundance/76215*1000

# Create the pack level percent mortality based off of count + removals and the number removed
# by harvest, control, or other human mortality
BP$PackPerHarv <- BP$MORT_HARVEST/BP$CountPlusRemovals
BP$PackPerCont <- BP$MORT_CONTROL/BP$CountPlusRemovals
BP$PackPerOther <- BP$MORT_HUMAN/BP$CountPlusRemovals



######################################################################################

# Final analyses

######################################################################################


library(lme4)

m1 <- glmer(BP~EoY.Count+P.Lambda+(1|Revised.Pack.Name), data=BP, family="binomial")
m2 <- glmer(BP~EoY.Count+P.Lambda+PopPerHarv+(1|Revised.Pack.Name), data=BP, family="binomial")
m3 <- glmer(BP~EoY.Count+P.Lambda+PopPerCont+(1|Revised.Pack.Name), data=BP, family="binomial")
m4 <- glmer(BP~EoY.Count+P.Lambda+PackPerHarv+(1|Revised.Pack.Name), data=BP, family="binomial")
m5 <- glmer(BP~EoY.Count+P.Lambda+PackPerCont+(1|Revised.Pack.Name), data=BP, family="binomial")
m6 <- glmer(BP~EoY.Count+P.Lambda+I(PopPerHarv+PopPerCont+PopPerOther)+(1|Revised.Pack.Name), data=BP, family="binomial")
m7 <- glmer(BP~EoY.Count+P.Lambda+I(PackPerHarv+PackPerCont+PackPerOther)+(1|Revised.Pack.Name), data=BP, family="binomial")
m8 <- glmer(BP~EoY.Count*I(PopPerHarv+PopPerCont+PopPerOther)+P.Lambda+(1|Revised.Pack.Name), data=BP, family="binomial")
m9 <- glmer(BP~EoY.Count*I(PackPerHarv+PackPerCont+PackPerOther)+P.Lambda+(1|Revised.Pack.Name), data=BP, family="binomial")
m10 <- glmer(BP~EoY.Count+(1|Revised.Pack.Name), data=BP, family="binomial")
m11 <- glmer(BP~EoY.Count+Density+(1|Revised.Pack.Name), data=BP, family="binomial")

m <- list()
head(m)

m[[1]]=m1
m[[2]]=m2
m[[3]]=m3
m[[4]]=m4
m[[5]]=m5
m[[6]]=m6
m[[7]]=m7
m[[8]]=m8
m[[9]]=m9
m[[10]]=m10
m[[11]]=m11


head(m)
## then name our models .

model.names <-c("Count+Lambda","Count+Lambda+PopHarvest", "Count+Lambda+PopControl",
                "Count+Lambda+PackHarvest", "Count+Lambda+PackControl", 
                "Count+Lambda+PopMort", "Count+Lambda+PackMort",
                "Count+Lambda+PopMort+Count*PopMort", 
                "Count+Lambda+PackMort+Count*PackMort","Count", "Count+Density")

library(AICcmodavg)
aictab(cand.set = m, modnames = model.names)


# Look at top model
summary(m6)

# The random effect doesn't seem to be important. It is 0. After looking at all top models, it
# seems we can lose the random effect. So now I will run all the models without the random effect
# and use that for our inference
m1g <- glm(BP~EoY.Count+P.Lambda, data=BP, family="binomial")
m2g <- glm(BP~EoY.Count+P.Lambda+PopPerHarv, data=BP, family="binomial")
m3g <- glm(BP~EoY.Count+P.Lambda+PopPerCont, data=BP, family="binomial")
m4g <- glm(BP~EoY.Count+P.Lambda+PackPerHarv, data=BP, family="binomial")
m5g <- glm(BP~EoY.Count+P.Lambda+PackPerCont, data=BP, family="binomial")
m6g <- glm(BP~EoY.Count+P.Lambda+I(PopPerHarv+PopPerCont+PopPerOther), data=BP, family="binomial")
m7g <- glm(BP~EoY.Count+P.Lambda+I(PackPerHarv+PackPerCont+PackPerOther), data=BP, family="binomial")
m8g <- glm(BP~EoY.Count*I(PopPerHarv+PopPerCont+PopPerOther)+P.Lambda, data=BP, family="binomial")
m9g <- glm(BP~EoY.Count*I(PackPerHarv+PackPerCont+PackPerOther)+P.Lambda, data=BP, family="binomial")
m10g <- glm(BP~EoY.Count, data=BP, family="binomial")
m11g <- glm(BP~EoY.Count+Density, data=BP, family="binomial")
m12g <- glm(BP~EoY.Count+harvest, data=BP, family="binomial")

mg <- list()
head(mg)

mg[[1]]=m1g
mg[[2]]=m2g
mg[[3]]=m3g
mg[[4]]=m4g
mg[[5]]=m5g
mg[[6]]=m6g
mg[[7]]=m7g
mg[[8]]=m8g
mg[[9]]=m9g
mg[[10]]=m10g
mg[[11]]=m11g
mg[[12]]=m12g


head(mg)
## then name our models .

model.names <-c("Count+Lambda","Count+Lambda+PopHarvest", "Count+Lambda+PopControl",
                "Count+Lambda+PackHarvest", "Count+Lambda+PackControl", 
                "Count+Lambda+PopMort", "Count+Lambda+PackMort",
                "Count+Lambda+PopMort+Count*PopMort", 
                "Count+Lambda+PackMort+Count*PackMort","Count", "Count+Density", "Count+harv")

library(AICcmodavg)
aictab(cand.set = mg, modnames = model.names)



# Check for overdispersion. If overdispered then the Laplance Approximiation in not appropriate: 
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}

overdisp_fun(m6g)


# Plot the fit to see how well the model does at describing data and to see if we are missing
# something that is important in the variation

plot(fitted(m6), residuals(m6), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(m6), residuals(m6)))


pchisq(sum(residuals(m6g, type="deviance")^2), df=m6g$df.residual, lower.tail=FALSE)

sum(residuals(m6g, type="pearson")^2)/m6g$df.residual




modavg(cand.set=mg, parm="EoY.Count", modnames=model.names, 
       exclude=list("EoY.Count:I(PopPerHarv + PopPerCont + PopPerOther)", 
                    "EoY.Count:I(PackPerHarv + PackPerCont + PackPerOther)"))
modavg(cand.set=mg, parm="P.Lambda", modnames=model.names)
modavg(cand.set=mg, parm="(Intercept)", modnames=model.names)




# Determine AUC and graph ROC curve. 
# Get fitted values from the model
mod1fitted <- fitted(m11g)

# Make ROCR prediction and performance classes, but first open up the package
library(ROCR)
predClass <- prediction(mod1fitted, BP$BP)
str(predClass)

# Make ROCR performance class from the prediction class
prefClass <- performance(predClass, "tpr","fpr") # change 2nd and/or 3rd arguments for other metrics

# Calculate the Area Under Curve
mod1auc <- performance(predClass, measure="auc") 
str(mod1auc)
auc <- as.numeric(mod1auc@y.values)
auc

mod1sens <- performance(predClass, measure="sens") 
str(mod1sens)
sens <- as.numeric(mod1sens@y.values)
sens

# Max the sum of sensitivity and specificity 
fpr <- prefClass@x.values[[1]]
tpr <- prefClass@y.values[[1]]
sum <- tpr + (1-fpr)
index <- which.max(sum)
cutoff <- prefClass@alpha.values[[1]][[index]]
cutoff

# Plot ROC Curve with cut point and AUC
plot(prefClass, colorize = T, lwd = 5, print.cutoffs.at=seq(0,1,by=0.1),
     text.adj=c(1.2,1.2),
     main = "ROC Curve")
text(0.5, 0.5, "AUC = 0.889")
abline(v=cutoff, col = "red", lwd = 3)



newdat <- data.frame("EoY.Count"=rep(seq(4,15,by=1), 3), 
                     "P.Lambda"=rep(1.05, 12*3), 
                     "PopPerHarv"=c(rep(.084, 12), rep(.26, 12),rep(.40, 12)),
                     "PopPerCont"=rep(0,12*3),
                     "PopPerOther"=rep(0,12*3))
pred.m6fin <- predict(m6g, newdat, se.fit=TRUE)

# Transform to probabilities for plotting
newdat$prob <- with(pred.m6fin, exp(fit)/(1+exp(fit)))
newdat$UL <- with(pred.m6fin, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
newdat$LL <- with(pred.m6fin, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))

ggplot(newdat, aes(x=EoY.Count, y=prob, linetype=as.factor(PopPerHarv), colours=as.factor(PopPerHarv))) +
  geom_line(size=1.25) +
  geom_ribbon(aes(ymin=LL,ymax=UL, fill=as.factor(PopPerHarv)),alpha=0.1)+
  theme_classic()+
  scale_color_grey()+
  scale_x_continuous(name="Pack size", breaks=seq(4,15, by=2))+
  scale_y_continuous(name="Probability pack contains a successful breeding pair")+
  labs(fill='% Mortality', linetype="% Mortality")





# Plot data stuff

ggplot(covs, aes(x=YEAR, y=P.Abundance))+
  geom_line()+
  theme_classic()+
  scale_color_grey()+
  scale_x_continuous(name="Year", breaks=seq(1986,2016, by=4))+
  scale_y_continuous(name="Abundance")
























# First test to see if any of the mortalities at the pack or population level can be combineded, 
# meaning no difference, using model  selection

# Run models. 
m1 <- glm(BP~EoY.Count+PopPerHarv+PopPerCont+PopPerOther, data=BP, family=(binomial))
m2 <- glm(BP~EoY.Count+I(PopPerHarv+PopPerCont)+PopPerOther, data=BP, family=(binomial))
m3 <- glm(BP~EoY.Count+PopPerHarv+I(PopPerCont+PopPerOther), data=BP, family=(binomial))
m4 <- glm(BP~EoY.Count+I(PopPerHarv+PopPerOther)+PopPerCont, data=BP, family=(binomial))
m5 <- glm(BP~EoY.Count+I(PopPerHarv+PopPerCont+PopPerOther), data=BP, family=(binomial))
m6 <- glm(BP~EoY.Count+PopPerHarv, data=BP, family=(binomial))
m7 <- glm(BP~EoY.Count+PopPerCont, data=BP, family=(binomial))
m8 <- glm(BP~EoY.Count+PopPerOther, data=BP, family=(binomial))

# Model selection
m <- list()
head(m)

m[[1]]=m1
m[[2]]=m2
m[[3]]=m3
m[[4]]=m4
m[[5]]=m5


head(m)
## then name our models .

model.names <-c("Separate","Harvest and Control", 
                "Control and Other", "Harvest and Other", 
                "All")

library(AICcmodavg)
aictab(cand.set = m, modnames = model.names)

# Now at pack level
m1p <- glm(BP~EoY.Count+PackPerHarv+PackPerCont+PackPerOther, data=BP, family=(binomial))
m2p <- glm(BP~EoY.Count+I(PackPerHarv+PackPerCont)+PackPerOther, data=BP, family=(binomial))
m3p <- glm(BP~EoY.Count+PackPerHarv+I(PackPerCont+PackPerOther), data=BP, family=(binomial))
m4p <- glm(BP~EoY.Count+I(PackPerHarv+PackPerOther)+PackPerCont, data=BP, family=(binomial))
m5p <- glm(BP~EoY.Count+I(PackPerHarv+PackPerCont+PackPerOther), data=BP, family=(binomial))
m6p <- glm(BP~EoY.Count+PackPerHarv, data=BP, family=(binomial))
m7p <- glm(BP~EoY.Count+PackPerCont, data=BP, family=(binomial))
m8p <- glm(BP~EoY.Count+PackPerOther, data=BP, family=(binomial))

# Model selection
mp <- list()
head(mp)

mp[[1]]=m1p
mp[[2]]=m2p
mp[[3]]=m3p
mp[[4]]=m4p
mp[[5]]=m5p


head(mp)
## then name our models .

model.names <-c("Separate","Harvest and Control", 
                "Control and Other", "Harvest and Other", 
                "All")

library(AICcmodavg)
aictab(cand.set = mp, modnames = model.names)

# Not much support either way... seems like better to just keep separate for both. 

# Now I am going to run combinations of density and lambda to see which base model is better, but first
# I need to test for collinearity

cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X, use="complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R
}

cor.prob(as.matrix(BP[,c("Density","P.Lambda","PopPerHarv","PopPerCont","PopPerOther",
                         "PackPerHarv","PackPerCont","PackPerOther", "PopMort", "PackMort")]))


# From this we see some collinearity issues.  We can fix this by combinging the human caused 
# mortalities into 1 covariate, which is fine according to models and model selectin up above. 
# Also, we should not run pop level morts and pack level morts in the same model because they are
# collinear. Also, density and lambda have high collinearity and should not be included in the same
# model. 

m1b <- glm(BP~EoY.Count, data=BP, family=(binomial))
m2b <- glm(BP~EoY.Count+Lambda, data=BP, family=(binomial))
m3b <- glm(BP~EoY.Count+Density, data=BP, family=(binomial))


mb <- list()
head(mb)

mb[[1]]=m1b
mb[[2]]=m2b
mb[[3]]=m3b


head(mb)
## then name our models .

model.namesb <-c("Count","Count+Lambda", "Count+Density")

aictab(cand.set = mb, modnames = model.namesb)


# The base model with most support is Count + Lambda, so now I will test hypotheses with effects
# of harvest
m1fin <- glm(BP~EoY.Count+Lambda, data=BP, family=(binomial))
m2fin <- glm(BP~EoY.Count+Lambda+PopPerHarv+PopPerCont+PopPerOther, data=BP, family=(binomial))
m3fin <- glm(BP~EoY.Count+Lambda+I(PopPerHarv+PopPerCont+PopPerOther), data=BP, family=(binomial))
m4fin <- glm(BP~EoY.Count*I(PopPerHarv+PopPerCont+PopPerOther)+Lambda, data=BP, family=(binomial))
m5fin <- glm(BP~EoY.Count+Lambda+PackPerHarv+PackPerCont+PackPerOther, data=BP, family=(binomial))
m6fin <- glm(BP~EoY.Count+Lambda+I(PackPerHarv+PackPerCont+PackPerOther), data=BP, family=(binomial))
m7fin <- glm(BP~EoY.Count*I(PackPerHarv+PackPerCont+PackPerOther)+Lambda, data=BP, family=(binomial))

m8fin <- glm(BP~EoY.Count+I(EoY.Count^2)+Lambda, data=BP, family=(binomial))
mfin <- list()
head(mfin)

mfin[[1]]=m1fin
mfin[[2]]=m2fin
mfin[[3]]=m3fin
mfin[[4]]=m4fin
mfin[[5]]=m5fin
mfin[[6]]=m6fin
mfin[[7]]=m7fin

head(mfin)
## then name our models .

model.namesfin <-c("Count+Lambda","Count+Lambda+PopHarvest+PopControl+PopOther","Count+Lambda+PopMort", 
                   "Count*PopMort+Lambda","Count+Lambda+PackHarvest+PackControl+PackOther",
                   "Count+Lambda+PackMort", "Count*PackMort+Lambda")

aictab(cand.set = mfin, modnames = model.namesfin)



# Now I am going to make things for plotting
newdat <- data.frame("EoY.Count"=rep(seq(4, 23, 1), 2), 
                     "Lambda"=c(rep(max(BP$Lambda, na.rm=TRUE),20), rep(1,20)), 
                     "PopPerHarv"=c(rep(max(BP$PopPerHarv), 20),rep(0,20)), 
                     "PopPerCont"=c(rep(max(BP$PopPerCont), 20),rep(0,20)), 
                     "PopPerOther"=c(rep(max(BP$PopPerOther), 20),rep(0,20)))
pred.m2fin <- predict(m2fin, newdat, se.fit=TRUE)

# Transform to probabilities for plotting
newdat$prob <- with(pred.m2fin, exp(fit)/(1+exp(fit)))
newdat$UL <- with(pred.m2fin, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
newdat$LL <- with(pred.m2fin, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))

ggplot(newdat, aes(x=EoY.Count, y=prob, group=PopPerHarv, colour=PopPerHarv)) +
  geom_line(size=1.25) +
  geom_ribbon(aes(ymin=LL,ymax=UL, fill=PopPerHarv),alpha=0.1)+
  labs(x="Count", y="Probability a pack has a breeding pair")










# Now test which forms of harvest are necessary
m1p <- glm(BP~EoY.Count+PackPerHarv+PackPerCont+PackPerOther, data=BP, family=(binomial))
m1p <- glm(BP~EoY.Count+PackPerHarv+PackPerCont, data=BP, family=(binomial))
m1p <- glm(BP~EoY.Count+PackPerHarv+PackPerOther, data=BP, family=(binomial))
m1p <- glm(BP~EoY.Count+PackPerCont+PackPerOther, data=BP, family=(binomial))
m1p <- glm(BP~EoY.Count+PackPerHarv, data=BP, family=(binomial))
m1p <- glm(BP~EoY.Count+PackPerHarv+PackPerCont+PackPerOther, data=BP, family=(binomial))



















m1 <- glm(BP~EoY.Count+PopHarvMort+PopContMort+PopOtherMort, data=BP, family=(binomial))
m2 <- glm(BP~EoY.Count+I(PopHarvMort+PopContMort+PopOtherMort), data=BP, family=(binomial))
m3 <- glm(BP~EoY.Count+PopHarvMort+PopContMort, data=BP, family=(binomial))
m4 <- glm(BP~EoY.Count*PopHarvMort+PopContMort, data=BP, family=(binomial))
m5 <- glm(BP~EoY.Count*PopContMort+PopHarvMort, data=BP, family=(binomial))
m6 <- glm(BP~EoY.Count, data=BP, family=(binomial))
m7 <- glm(BP~EoY.Count+MORT_HUMAN+MORT_HARVEST+MORT_CONTROL, data=BP, family=(binomial))
m8 <- glm(BP~EoY.Count+I(MORT_HUMAN+MORT_HARVEST+MORT_CONTROL), data=BP, family=(binomial))
m9 <- glm(BP~EoY.Count+Lambda+Est.pop.size, data=BP, family=(binomial))
m10 <- glm(BP~EoY.Count+PopPerHarv+PopPerCont, data=BP, family=(binomial))
m11 <- glm(BP~EoY.Count+I(MORT_HUMAN/CountPlusRemovals)+I(MORT_HARVEST/CountPlusRemovals)+I(MORT_CONTROL/CountPlusRemovals), data=BP, family=(binomial))










# Model 1
m1 <- glm(BP~EoY.Count, data=bp.good2, family=(binomial))
m2 <- glm(BP~EoY.Count+factor(harvest), data=bp.good2, family=(binomial))
m3 <- glm(BP~EoY.Count*harvest, data=bp.good2, family=(binomial))
m4 <- glm(BP~EoY.Count+RECAREA, data=bp.good2, family=(binomial))
m5 <- glm(BP~EoY.Count+RECAREA+harvest, data=bp.good2, family=(binomial))
m6 <- glm(BP~EoY.Count+RECAREA*harvest, data=bp.good2, family=(binomial))
m7 <- glm(BP~EoY.Count+YEAR, data=bp.good2, family=(binomial))
m8 <- glm(BP~EoY.Count*RECAREA, data=bp.good2, family=(binomial))

# load the package for random effects and run random effects model
library(lme4)

mr1 <- glmer(BP~EoY.Count+(1|YEAR), data=bp.good, family="binomial")
mr2 <- glmer(BP~EoY.Count+harvest+(1|YEAR), data=bp.good, family="binomial")



m <- list()
head(m)

m[[1]]=m1
m[[2]]=m2
m[[3]]=m3
m[[4]]=m4
m[[5]]=m5
m[[6]]=m6
#m[[7]]=m7
m[[7]]=m8


head(m)
## then name our models .

model.names <-c("Count","Count+harvest", 
                "Count+harvest+Count*harvest", "Count+area", 
                "Count+area+harvest", 
                "Count+area+harvest+area*harvest", 
                #"Count+year", 
                "Count+area+count*area")

library(AICcmodavg)
aictab(cand.set = m, modnames = model.names)


cor.prob(as.matrix(bp.good[,c("YEAR","EoY.Count","harvest")]))

summary(m3)

n.sim<-1000
simu<-sim(mr1, n.sims=n.sim)
head(simu@fixef)
apply(simu@fixef,2,quantile,prob=c(0.025,0.5,0.975))

nsim <- 1000
bsim <- sim(mr1, n.sim=nsim)
newdat <- data.frame("EoY.Count"=seq(4, 15, 0.2))
Xmat <- model.matrix(~EoY.Count, data=newdat)
predmat <- matrix(ncol=nsim, nrow=nrow(newdat))
predmat<-apply(bsim@fixef,1,function(x) exp(Xmat%*%x)/(1+exp(Xmat%*%x)))
newdat$lower <- apply(predmat, 1, quantile, prob=0.025)
newdat$upper <- apply(predmat, 1, quantile, prob=0.975)
newdat$med<-apply(predmat, 1, quantile, prob=0.5)

ggplot(newdat, aes(x=EoY.Count, y=med)) +
  geom_line(size=1.25) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.1)+
  labs(x="Count", y="Probability a pack has a breeding pair")



# # First grab all of the estimates and standard errors
# models = rbind(summary(m1)$coefficients[,1:2], 
#                summary(m2)$coefficients[,1:2], 
#                summary(m3)$coefficients[,1:2], 
#                summary(m4)$coefficients[,1:2], 
#                summary(m5)$coefficients[,1:2], 
#                summary(m6)$coefficients[,1:2], 
#                summary(m7)$coefficients[,1:2])
# # Name your models
# modelnames = c("Count","Count+harvest", 
#                "Count+harvest+Count*harvest", "Count+area", 
#                "Count+area+harvest", 
#                "Count+area+harvest+area*harvest", 
#                "Random effect")
# # Now put all of your estimates in a pretty table with names that you'll remember!
# estimates.all = matrix(models, nrow=2*length(modelnames), ncol=2, dimnames = list(paste(rep(modelnames, each=2),c("intercept", "coefficient")), c("B", "SE")))
# estimates.all














# Make BP a number and not a factor
bp.good$BP2 <- ifelse(bp.good$BP == "1", 1, 0)

# First model: harvest
harv.good <- glm(BP~harvest, data=bp.good, family=(binomial))
summary(harv.good)

# Second model: study area
area.good <- glm(BP~AREA, data=bp.good, family=(binomial))
summary(area.good)

# Third model: End of year count
eoy.good <- glm(BP~EoY.Count, data=bp.good, family=(binomial))
summary(eoy.good)



# Now run the single model equivalent to analyses Mike did earlier. 
mod1 <- glm(BP~EoY.Count+factor(harvest)+AREA+(factor(harvest):AREA)+(factor(harvest):EoY.Count), data=bp.good, family=binomial())
summary(mod1)

# Determine AUC and graph ROC curve. 
  # Get fitted values from the model
  mod1fitted <- fitted(mr1)
  
  # Make ROCR prediction and performance classes, but first open up the package
  library(ROCR)
  predClass <- prediction(mod1fitted, bp.good$BP)
  str(predClass)
  
  # Make ROCR performance class from the prediction class
  prefClass <- performance(predClass, "tpr","fpr") # change 2nd and/or 3rd arguments for other metrics
  
  # Calculate the Area Under Curve
  mod1auc <- performance(predClass, measure="auc") 
  str(mod1auc)
  auc <- as.numeric(mod1auc@y.values)
  auc
  
  mod1sens <- performance(predClass, measure="sens") 
  str(mod1sens)
  sens <- as.numeric(mod1sens@y.values)
  sens
  
  # Max the sum of sensitivity and specificity 
  fpr <- prefClass@x.values[[1]]
  tpr <- prefClass@y.values[[1]]
  sum <- tpr + (1-fpr)
  index <- which.max(sum)
  cutoff <- prefClass@alpha.values[[1]][[index]]
  cutoff
  
  # Plot ROC Curve with cut point and AUC
  plot(prefClass, colorize = T, lwd = 5, print.cutoffs.at=seq(0,1,by=0.1),
       text.adj=c(1.2,1.2),
       main = "ROC Curve")
  text(0.5, 0.5, "AUC = 0.889")
  abline(v=cutoff, col = "red", lwd = 3)


# Now plot results from model
  
  datum <- data.frame("EoY.Count"=rep(rep(seq(4, 15, 0.2),3),2), "harvest" = factor(rep(c(rep(0, 56), rep(1,56)),3)),
                      "RECAREA"=rep(c(rep("NWMT",56), rep("CID",56), rep("GYA",56)),2))
  
  datum <- data.frame("EoY.Count"=rep(seq(4, 15, 0.2),2), "harvest" = factor(c(rep(0, 56), rep(1,56))))
  
  # First create new data frames to use for predictions
  dat0CID <-data.frame("EoY.Count"=seq(4,15,0.2), "harvest"=factor(rep(0, 56)), 
                    "AREA"=rep("CID", 56))
  dat1CID <-data.frame("EoY.Count"=seq(4,15,0.2), "harvest"=factor(rep(1, 56)), 
                       "AREA"=rep("CID", 56))
  dat0GYA <-data.frame("EoY.Count"=seq(4,15,0.2), "harvest"=factor(rep(0, 56)), 
                       "AREA"=rep("GYA", 56))
  dat1GYA <-data.frame("EoY.Count"=seq(4,15,0.2), "harvest"=factor(rep(1, 56)), 
                       "AREA"=rep("GYA", 56))
  dat0NWMT <-data.frame("EoY.Count"=seq(4,15,0.2), "harvest"=factor(rep(0, 56)), 
                       "AREA"=rep("NWMT", 56))
  dat1NWMT <-data.frame("EoY.Count"=seq(4,15,0.2), "harvest"=factor(rep(1, 56)), 
                       "AREA"=rep("NWMT", 56))
  
  # Now make predictions based on new data
  pred0CID <- predict(m6, datum, se.fit=TRUE)
  pred1CID <- predict(mod1, dat1CID, se.fit=TRUE)
  pred0GYA <- predict(mod1, dat0GYA, se.fit=TRUE)
  pred1GYA <- predict(mod1, dat1GYA, se.fit=TRUE)
  pred0NWMT <- predict(mod1, dat0NWMT, se.fit=TRUE)
  pred1NWMT <- predict(mod1, dat1NWMT, se.fit=TRUE)
  
  # Transform to probabilities for plotting
  datum$prob <- with(pred0CID, exp(fit)/(1+exp(fit)))
  dat1CID$prob <- with(pred1CID, exp(fit)/(1+exp(fit)))
  dat0GYA$prob <- with(pred0GYA, exp(fit)/(1+exp(fit)))
  dat1GYA$prob <- with(pred1GYA, exp(fit)/(1+exp(fit)))
  dat0NWMT$prob <- with(pred0NWMT, exp(fit)/(1+exp(fit)))
  dat1NWMT$prob <- with(pred1NWMT, exp(fit)/(1+exp(fit)))
  
  # Transform the standared errors to CIs
  datum$UL <- with(pred0CID, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  datum$LL <- with(pred0CID, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  dat1CID$UL <- with(pred1CID, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  dat1CID$LL <- with(pred1CID, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  dat0GYA$UL <- with(pred0GYA, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  dat0GYA$LL <- with(pred0GYA, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  dat1GYA$UL <- with(pred1GYA, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  dat1GYA$LL <- with(pred1GYA, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  dat0NWMT$UL <- with(pred0NWMT, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  dat0NWMT$LL <- with(pred0NWMT, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  dat1NWMT$UL <- with(pred1NWMT, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  dat1NWMT$LL <- with(pred1NWMT, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  
  # Make a dataframe that contains all the data for easy plotting
  datum <- rbind.data.frame(dat0CID, dat0GYA, dat0NWMT, dat1CID, dat1GYA, dat1NWMT)
  
  # Now finally plot it
  ggplot(datum, aes(x=EoY.Count, y=prob, group=harvest, colour=harvest)) +
    geom_line(size=1.25) +
    geom_ribbon(aes(ymin=LL,ymax=UL, fill=harvest),alpha=0.1)+
    labs(x="Group Count", y="Probability of breeding pair")+
    facet_grid(~RECAREA)


# Run the single models to compare with Mike's results
  
  cid.good <- bp.good[bp.good$AREA == "CID",]
  gya.good <- bp.good[bp.good$AREA == "GYA",]
  nwmt.good <- bp.good[bp.good$AREA == "NWMT",]
  
  cid.mod <- glm(BP~EoY.Count+factor(harvest), data=cid.good, family=(binomial))
  gya.mod <- glm(BP~EoY.Count+factor(harvest), data=gya.good, family=(binomial))
  nwmt.mod <- glm(BP~EoY.Count+factor(harvest), data=nwmt.good, family=(binomial))
  
  summary(cid.mod)
  summary(gya.mod)
  summary(nwmt.mod)
  
  

# model 2
  
  mod2 <- glm(BP~EoY.Count*factor(harvest), data=bp.good, family=binomial())
  summary(mod2)
  
  
  datum2 <- data.frame("EoY.Count"=rep(seq(4,15,.2),2), "harvest"=c(rep(0,56),rep(1,56)))
  pred <- predict(mod2, datum2, se.fit=TRUE)
  datum2$prob <- with(pred, exp(fit)/(1+exp(fit)))
  datum2$UL <- with(pred, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  datum2$LL <- with(pred, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  
  
  ggplot(datum2, aes(x=EoY.Count, y=prob, group=harvest, colour=harvest)) +
    geom_line(size=1.25) +
    geom_ribbon(aes(ymin=LL,ymax=UL, fill=harvest),alpha=0.1)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# Other method for ROC

predicted <- predict(mod1)
library(ROCR)
prob <- prediction(predicted, bp.good$BP)
tprfpr <- performance(prob, "tpr", "fpr")
tpr <- unlist(slot(tprfpr, "y.values"))
fpr <- unlist(slot(tprfpr, "x.values"))
roc <- data.frame(tpr, fpr)
ggplot(roc) + geom_line(aes(x = fpr, y = tpr)) + 
  geom_abline(intercept = 0, slope = 1, colour = "gray") + 
  ylab("Sensitivity") + 
  xlab("1 - Specificity")




