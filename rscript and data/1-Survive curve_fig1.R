#survive curve

library(tidyverse)
library(MASS)


aphid<-read.csv("aphid_surv.csv", sep =',',header = T) #read in survival data

###calculate lt50 
species = unique(as.character(aphid$spec))
temps=unique(aphid$temp)

lt <- data.frame(species=character(),temp=integer(),lt50=numeric(),lt75=numeric(),se=numeric())

for(i in 1:length(species)){
  for(k in 1:length(temps)){
    tmp = aphid%>%filter(spec == species[i]&temp == temps[k])
    if(dim(tmp)[1]!=0){
      LT50<-dose.p(glm(cbind(surv,total-surv) ~ dur,data=tmp, family = binomial(link=logit)),p=0.5)
      LT75<-dose.p(glm(cbind(surv,total-surv) ~ dur,data=tmp, family = binomial(link=logit)),p=0.75)
      df=data.frame(species[i],temps[k],lt50=unname(LT50[1:1]),lt75=unname(LT75[1:1]),SE=unname(attributes(LT50)$SE))
      names(df)=c("species","temp","lt50","lt75","se")
      lt<-rbind(lt,df)}
  }
}

###plot survival curve
pdf("figure1_surv_34.pdf", width = 3.5, height = 3.5)
ggplot(aphid,aes(dur, survival_rate,colour=spec))+
  geom_point(size=2.0,alpha = 0.4)+
  stat_smooth(aes(group=spec),size=1.3,se=F,method = glm,fullrange=T,method.args = list(family = binomial))+
  scale_x_log10(breaks=c(10, 60 ,180,600,1440),labels=c("10min","1h","3h","10h","1d"))+
  scale_color_discrete(labels = c("M. dirhodum", "R. padi", "S. avenae"))+
  labs(x="Exposure time",y="Survival proportion",color = "Species")+
  geom_hline(yintercept = 0.5, linetype='dotted',size= 0.5)+
  theme( panel.background=element_rect(fill="white"),
         legend.text = element_text(face = "italic",size = 8),
         legend.position = c(0.84,0.85),
         legend.title = element_text(size = 8),
         axis.line=element_line(colour="black",size=1.0),
         axis.text.y=element_text(size=8,colour = "black"),
         axis.title =element_text(size=8),
         axis.text.x = element_text(size= 8,colour = "black",angle = 40, hjust=1))

dev.off()
