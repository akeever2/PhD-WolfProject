# % difference line fitting to estimate abundance

x<-1986:2016
y<-c(rep(0, 10), rep(NA,11), .317817,.279539, .37201,.328622,.319372,.304783,.413146,.373576,.453694,.441833)
m1<-lm(y~x)
m2<-lm(y~poly(x,2))
m3<-lm(y~poly(x,3))
m4<-lm(y~poly(x,4))


library("AICcmodavg", lib.loc="~/R/win-library/3.4")
 
m<-list()
m[[1]]<-m1
m[[2]]<-m2
m[[3]]<-m3
m[[4]]<-m4
model.names<-c("linear", "2nd", "3rd", "4th")
aictab(cand.set=m, modnames=model.names)

anova(m1,m2)
anova(m2,m3)
anova(m3,m4)

library("ggplot2", lib.loc="~/R/win-library/3.4")

ggplot(datum, aes(x=x, y=y))+
  geom_point(size=2)+
  geom_smooth(method=lm, se=TRUE, 
              formula=y~poly(x,3))+
  theme_classic()+
  scale_color_grey()+
  scale_x_continuous(name="Year", breaks=seq(1986,2016, by=4))+
  scale_y_continuous(name="% Difference")





