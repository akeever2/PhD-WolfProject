###########################################################################################

#                   Analyses for breeding pair estimator

###########################################################################################


# Pull in the data (CSV file) and name it bp. Function opens up a browser to find the file. 
bp <- read.csv(file.choose())

# Look at distribution of the data to get a feel for it
hist(bp$EoY.Count)
hist(bp$YEAR)
hist(bp$MORT_CONTROL)
hist(bp$TOTAL_REMOVED)

# Add a harvest covariate dependent on year. Harvest is 1 if the year is greater than 
# 2010 or if the year is 2009, else it is 0
bp$harvest <- ifelse(bp$YEAR > 2010 | bp$YEAR== 2009, 1, 0)

# Make sure the harvest covariate is right
table(bp$YEAR, bp$harvest)

# Check the structure of the data
str(bp)

# It is treating BP as a factor of 3 becacuse of ".", so I am removing the . rows and 
# also the censored rows from the count quality 

bp <- bp[!((bp$BP == ".") | (bp$Count.Quality == "CENSOR") | bp$Count.Quality=="U"),]
bp$BP <- factor(bp$BP)
bp$Count.Quality <- factor(bp$Count.Quality)

str(bp)


# Now run some basic analyses

# First model: harvest
harv <- glm(BP~harvest, data=bp, family=(binomial))
summary(harv)

# Second model: study area
area <- glm(BP~AREA, data=bp, family=(binomial))
summary(area)

# Third model: End of year count
eoy <- glm(BP~EoY.Count, data=bp, family=(binomial))
summary(eoy)

# Fourth model: Count quality
qual <- glm(BP~Count.Quality, data=bp, family=(binomial))
summary(qual)


# Having a poor or medium quality count means the pack is less likely to have a breeding
# pair. Could this potentially bias results...

# Subset the data for good counts only, and rerun the above simple analyses

# Create a subset of the data
bp.good <- bp[bp$Count.Quality == "G",]
bp.good$Count.Quality <- factor(bp.good$Count.Quality)

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
  mod1fitted <- fitted(mod1)
  
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
  text(0.5, 0.5, "AUC = 0.872")
  abline(v=cutoff, col = "red", lwd = 3)


# Now plot results from model
  
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
  pred0CID <- predict(mod1, dat0CID, se.fit=TRUE)
  pred1CID <- predict(mod1, dat1CID, se.fit=TRUE)
  pred0GYA <- predict(mod1, dat0GYA, se.fit=TRUE)
  pred1GYA <- predict(mod1, dat1GYA, se.fit=TRUE)
  pred0NWMT <- predict(mod1, dat0NWMT, se.fit=TRUE)
  pred1NWMT <- predict(mod1, dat1NWMT, se.fit=TRUE)
  
  # Transform to probabilities for plotting
  dat0CID$prob <- with(pred0CID, exp(fit)/(1+exp(fit)))
  dat1CID$prob <- with(pred1CID, exp(fit)/(1+exp(fit)))
  dat0GYA$prob <- with(pred0GYA, exp(fit)/(1+exp(fit)))
  dat1GYA$prob <- with(pred1GYA, exp(fit)/(1+exp(fit)))
  dat0NWMT$prob <- with(pred0NWMT, exp(fit)/(1+exp(fit)))
  dat1NWMT$prob <- with(pred1NWMT, exp(fit)/(1+exp(fit)))
  
  # Transform the standared errors to CIs
  dat0CID$UL <- with(pred0CID, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  dat0CID$LL <- with(pred0CID, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
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
    geom_ribbon(aes(ymin=LL,ymax=UL, fill=harvest),alpha=0.1) +
    facet_grid(AREA~.)


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




