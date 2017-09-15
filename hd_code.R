####################################################################################

#             Code for TWS HD to get proportions (means) and plot results

####################################################################################

# Pull in the clean data from out dropbox folder, I did not specify location. 
# This will open a browser where you can select the file. 

hd <- read.csv(file.choose())


# Manipulate data and create new variables

# hd$HD<-hd$NoHD
# hd2$HD[hd2$HD==0]<-2


# Open packages needed to analyze and plot

library(nlme)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)


# Run linear models (lm) to get the means, i.e. proportion in each group for:
# H-D multiple hypotheses (hdm or Hdmult), H-D single hypothesis (hds or Hdsing)
# H-D vague or benefit of doubt (hdv or VagueHD), no H-D (nhd or NoHD), and 
# whether an expirment was conducted (ex or Experiment)

hdm <- lm(Hdmult~1, data=hd)
hds <- lm(Hdsing~1, data=hd)
hdv <- lm(VagueHD~1, data=hd)
nhd <- lm(NoHD~1, data=hd)
ex <- lm(Experiment~1, data=hd)


# Store the means in a data frame for plotting

results <- data.frame(Type=c("Hdmult", "Hdsing", "VagueHD", "NoHD", "Experiment"),
                      HD=c(1,1,1,0,NA),
                      Means=c(coef(hdm), coef(hds), coef(hdv), coef(nhd), coef(ex)),
                      SEs=c(coef(summary(hdm))[,2], coef(summary(hds))[,2], 
                            coef(summary(hdv))[,2], coef(summary(nhd))[,2], 
                            coef(summary(ex))[,2]))


# Plot the data in a pie graph... mmm pie...

  # Create a blank theme to simplify the pie graph

  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x=element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )


  # Make a bar plot of the data, without experiment in there
  
  bp <- ggplot(data=results[1:4,], aes(x="", y=Means[1:4], fill=Type))+
    geom_bar(stat="identity", width=1)
  
  bp
  
  
  #Make a pie graph of the data
  
  pie <- bp + coord_polar("y") + blank_theme + #scale_fill_grey() +
    geom_text(aes(y=Means[c(3:4,1:2)]/2 + c(0,cumsum(Means[c(3:4,1:2)])[-length(Means[c(3:4,1:2)])]), 
                  label=percent(Means[c(3:4,1:2)])), size=4.5) +
    scale_fill_grey(name="Category", 
                    breaks=c("Hdmult","Hdsing","VagueHD","NoHD"),
                    labels=c("Multiple Hypotheses", "Single Hypothesis", "Vague Hypotheses", 
                                 "Unreliable"),
                    start=0.4, end=0.9)

  pie



  
  # Export image
  
  
  png(file="piegraph.png", width=5.5, height=5, units="in", res=300, type="cairo")
  par(mar=c(5,3,2,2))
  pie
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
# Run random effects models for funsies

hdm2 <- lme(Hdmult~1, data=hd, random=(~1|Observer))
hds2 <- lme(Hdsing~1, data=hd, random=(~1|Observer))
hdv2 <- lme(VagueHD~1, data=hd, random=(~1|Observer))
nhd2 <- lme(NoHD~1, data=hd, random=(~1|Observer))
ex2 <- lme(Experiment~1, data=hd, random=(~1|Observer))

results <- results %>% 
  mutate(MeansR=c(fixef(hdm2), fixef(hds2), fixef(hdv2), fixef(nhd2), fixef(ex2)))#,
         # SEs=c(fixef(summary(hdm2))[,2], fixef(summary(hds2))[,2], 
         #       fixef(summary(hdv2))[,2], fixef(summary(nhd2))[,2], 
         #       fixef(summary(ex2))[,2]))



