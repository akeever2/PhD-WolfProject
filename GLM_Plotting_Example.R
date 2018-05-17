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

