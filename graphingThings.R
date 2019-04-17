#######################

# Graphing and things

#######################


datum <- read.csv(file.choose())

library(ggplot2)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





ggplot(data=datum, aes(x=Year, y=Rec))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Rec_L, ymax=Rec_U), 
                width=0.2, size=1)+
  scale_x_continuous(name="Year", breaks=c(2007:2016))+
  scale_y_continuous(name="Recruitment", breaks=c(1,2,3,4,5), limits=c(0,5))+
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text=element_text(size=12), 
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black"), 
        strip.background = element_rect(fill="#999999", color="black"), 
        strip.text = element_text(face="bold"),
        panel.spacing = unit(0, "lines"),
        legend.title=element_blank(),  
        legend.key=element_blank())

  
ggsave(filename="RecRate.png")




 Pop <- ggplot(data=datum, aes(x=Year))+
  geom_line(aes(y=Ntot), color="black")+
  geom_ribbon(aes(ymin=Ntot_L, ymax=Ntot_U), 
                color="grey", alpha=0.3)+
   geom_line(aes(y=FWP), color="blue")+
   geom_ribbon(aes(ymin=FWP_L, ymax=FWP_U), color="blue", 
               fill="blue", alpha=0.3)+
   scale_x_continuous(name="Year", breaks=c(2007:2016))+
  scale_y_continuous(name="Abundance", limits=c(250,1350), 
                     breaks=c(250, 500, 750, 1000, 1250))+
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
 
 Pack <- ggplot(data=datum, aes(x=Year, y=P))+
   geom_point(size=3)+
   geom_errorbar(aes(ymin=P_L, ymax=P_U), 
                 width=0.2, size=1)+
   scale_x_continuous(name="Year", breaks=c(2007:2016))+
   scale_y_continuous(name="# Packs")+
   theme(axis.text.x=element_text(size=12), 
         axis.text.y=element_text(size=12), 
         axis.title.x=element_text(size=14), 
         axis.title.y=element_text(size=14), 
         legend.text=element_text(size=12), 
         panel.background=element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_rect(fill=NA, color="black"), 
         strip.background = element_rect(fill="#999999", color="black"), 
         strip.text = element_text(face="bold"),
         panel.spacing = unit(0, "lines"),
         legend.title=element_blank(),  
         legend.key=element_blank())
 
 Surv <- ggplot(data=datum, aes(x=Year, y=S))+
   geom_point(size=3)+
   geom_errorbar(aes(ymin=S_L, ymax=S_U), 
                 width=0.2, size=1)+
   scale_x_continuous(name="Year", breaks=c(2007:2016))+
   scale_y_continuous(name="Survival",  limits=c(0.25,1))+
   theme(axis.text.x=element_text(size=12), 
         axis.text.y=element_text(size=12), 
         axis.title.x=element_text(size=14), 
         axis.title.y=element_text(size=14), 
         legend.text=element_text(size=12), 
         panel.background=element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_rect(fill=NA, color="black"), 
         strip.background = element_rect(fill="#999999", color="black"), 
         strip.text = element_text(face="bold"),
         panel.spacing = unit(0, "lines"),
         legend.title=element_blank(),  
         legend.key=element_blank())
 
 Group <- ggplot(data=datum, aes(x=Year, y=G))+
   geom_point(size=3)+
   geom_errorbar(aes(ymin=G_L, ymax=G_U), 
                 width=0.2, size=1)+
   scale_x_continuous(name="Year", breaks=c(2007:2016))+
   scale_y_continuous(name="Group size", limits=c(4,8))+
   theme(axis.text.x=element_text(size=12), 
         axis.text.y=element_text(size=12), 
         axis.title.x=element_text(size=14), 
         axis.title.y=element_text(size=14), 
         legend.text=element_text(size=12), 
         panel.background=element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_rect(fill=NA, color="black"), 
         strip.background = element_rect(fill="#999999", color="black"), 
         strip.text = element_text(face="bold"),
         panel.spacing = unit(0, "lines"),
         legend.title=element_blank(),  
         legend.key=element_blank())
 
 
ggarrange(Pop, Group, Surv, labels=c("A","B","C"),  ncol=1, nrow=3, align="v")

ggsave(filename="N_G_S.png")


