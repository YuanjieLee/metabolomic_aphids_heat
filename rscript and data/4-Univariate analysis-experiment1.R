library(scales)
library(tidyverse)

#Univariate analysis-experiment1

metabo_foldchange<-read.csv("metabo_foldchange.csv")

spec = unique(as.character(metabo_foldchange$species))
treatment = unique(metabo_foldchange$TRT)
temps = unique(metabo_foldchange$temp)
metabolite = unique(metabo_foldchange$metabolites)

p.ttest <- data.frame(species=character(),TRT = character(),temp=integer(),metabolites = character(),pvaule=numeric())

for(i in 1:length(spec)){
  for(j in 1:length(treatment)){
    for(k in 1:length(temps)){
      for(h in 1:length(metabolite)){
        data = metabo_foldchange%>%filter(species == spec[i]&TRT == treatment[j]&temp == temps[k]&metabolites == metabolite[h])
        if(dim(data)[1]!=0){
          test<-t.test(data$log_foldchange,mu = 0) ## we test log fold change
          p = data.frame(spec[i],treatment[j],temps[k],metabolite[h],round(test$p.value,4))
          names(p)=c("species","treatments","temp","metabolites","pvalue")
          p.ttest = rbind(p,p.ttest)}
      }
    }
  }
}
#write.csv(p.ttest,"p.t-test_metabolite.csv")

#### compare pairwise the interspecific difference for each metabolite
metabo_foldchange<-metabo_foldchange%>%filter(injurylevel!="C")
treatment = unique(metabo_foldchange$TRT)
temps = unique(metabo_foldchange$temp)
metabolite = unique(metabo_foldchange$metabolites)

pair.test <- data.frame(treatments = character(),temp=integer(),metabolites = character(),pvalue = numeric())

for(j in 1:length(treatment)){
  for(k in 1:length(temps)){
    for(h in 1:length(metabolite)){
      data = metabo_foldchange%>%filter(TRT == treatment[j]&temp == temps[k]&metabolites == metabolite[h])
      if(dim(data)[1]!=0){
        ano<-aov(log_foldchange ~species, data)
        a<-as.data.frame(TukeyHSD(ano)$species)
        p = data.frame(treatment[j],temps[k],metabolite[h],a[4])
        names(p)=c("treatments","temp","metabolites","pvalue")
        pair.test = rbind(p,pair.test)}
    }
  }
}


####figure 2. univariate analysis of experiment 1

###filter the metabolites that different from control using p-value <0.05 for treatment 2/2 D
p.test_metabo<-p.ttest%>%filter(pvalue<0.05&treatments=="2/2 D"&temp=="34")

###filter the metabolites that showed significant interspecific differences  
p.metabo<-pair.test%>%filter(temp=="34"&pvalue < 0.05)%>%
  filter(treatments=="2/2 D")%>%
  subset(.,metabolites%in%p.test_metabo$metabolites)

### get the data for these metabolites 
metabo_foldchange<-read.csv("metabo_foldchange.csv", header = T)

metabo_sj<-metabo_foldchange%>%filter(temp=='34')%>%
  filter(TRT=="2/2 D")%>%
  subset(.,metabolites %in% p.metabo$metabolites)
#metabo_sj<-metabo_sj[order(metabo_sj$fold_change, decreasing = T),]
#metabo_sj$metabolites<-factor(metabo_sj$metabolites,levels = unique(metabo_sj$metabolites))
metabo_sj$species<-factor(metabo_sj$species,levels=c("MD","SA","RP"))
metabo_sj$group<-factor(metabo_sj$group, levels=c( "Polyols","Carbohydrates","Amino acids","Amines", "TCA intermediates", "Organic acids",  "Phosphorylated compounds"))
metabo_sj$metabolites_ordered <- factor(metabo_sj$metabolites, 
                                        levels = unique(metabo_sj$metabolites[order(metabo_sj$group)]))

pdf("figure2_uni_sj.pdf", width = 7, height = 6)
ggplot(metabo_sj, aes(x = species, y = fold_change,fill = species)) + 
  geom_boxplot(outlier.size = 0.3,lwd =0.3,fatten =0.8)+
  facet_wrap(~metabolites_ordered, scales = "free",ncol = 6)+
  scale_fill_manual(name= "Species", values = c("#F8766D" ,"#619CFF","#00BA38"), 
                    labels = c("M. dirhodum", "S. avenae", "R. padi"))+
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Fold change in concentration of metabolites ")+
  theme(panel.background=element_rect(fill="white"),
        panel.spacing.y = unit(1.5, "lines"),
        legend.position = c(0.83,0.1),
        legend.text = element_text(face = "italic", size = 10),
        legend.direction = "vertical",
        legend.spacing.x = unit(0.3, 'cm'),
        strip.background = element_rect(fill="white"),
        strip.text.x = element_text(margin = margin(0.3,0,0.5,0, "cm")),
        strip.text = element_text(size = 9, family ="Arial"),
        axis.line=element_line(size=0.4),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
dev.off()

##3.2 figure s1.s2 individual metabolites in experiment 1

metabo_average<-read.csv("metabo_average.csv", header = T)
heat_sj1<-metabo_average%>%
  filter(treatment=='1/2 D'|treatment=='1/2 R'|treatment == 'C D'|treatment == 'C R')
heat_sj1$treatment[heat_sj1$treatment=='C D']<-'C'
heat_sj1$treatment[heat_sj1$treatment=='C R']<-'C'
heat_sj1$treatment[heat_sj1$treatment=='1/2 D']<-'D'
heat_sj1$treatment[heat_sj1$treatment=='1/2 R']<-'R'
heat_sj1$treatment<-factor(heat_sj1$treatment,levels = c("C","D","R"))
#heat_sj1<-heat_sj1[order(heat_sj1$fold_change.mean, decreasing = T),]
#heat_sj1$metabolites<-factor(heat_sj1$metabolites,levels = unique(heat_sj1$metabolites))

pdf("figureS3_heat_sj1.pdf", width = 10, height = 6)
ggplot(data = heat_sj1, aes(x = treatment, y = fold_change.mean, group = species, colour = species)) + geom_line()+
  geom_point() + geom_errorbar(aes(x= treatment, ymin = fold_change.mean-fold_change.se,ymax = fold_change.mean + fold_change.se), stat="identity", width = 0.1)+
  facet_wrap(~metabolites_ordered,scales = "free",nrow = 4, ncol = 8)+
  scale_colour_discrete(name = "Species", labels = c("M. dirhodum","R. padi","S. avenae"))+
  ylab("Fold change in metabolites concentration")+
  xlab("Treatment")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.text = element_text(face = "italic"),
        axis.text.x  = element_text(size = 8),
        panel.background=element_rect(fill="white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 7),
        axis.line=element_line(colour="grey",size=0.8))

dev.off()

### 2/2 Direct & Recovery
heat_sj2<-metabo_average%>%
  filter(treatment=='2/2 D'|treatment=='2/2 R'|treatment == 'C D'|treatment == 'C R')
heat_sj2$treatment[heat_sj2$treatment=='C D']<-'C'
heat_sj2$treatment[heat_sj2$treatment=='C R']<-'C'
heat_sj2$treatment[heat_sj2$treatment=='2/2 D']<-'D'
heat_sj2$treatment[heat_sj2$treatment=='2/2 R']<-'R'
heat_sj2$treatment<-factor(heat_sj2$treatment,levels = c("C","D","R"))
#heat_sj2<-heat_sj2[order(heat_sj2$fold_change.mean, decreasing = T),]
#heat_sj2$metabolites<-factor(heat_sj2$metabolites,levels = unique(heat_sj2$metabolites))

pdf("figureS4_heat_sj2.pdf", width = 10, height = 6)
ggplot(data = heat_sj2, aes(x = treatment, y = fold_change.mean, group = species, colour = species)) + geom_line(linetype = 2)+
  geom_point() + geom_errorbar(aes(x= treatment, ymin = fold_change.mean-fold_change.se,ymax = fold_change.mean + fold_change.se), stat="identity", width = 0.1)+
  facet_wrap(~metabolites_ordered,scales = "free",nrow = 4, ncol = 8)+
  scale_colour_discrete(name = "Species", labels = c("M. dirhodum","R. padi","S. avenae"))+
  ylab("Fold change in metabolites concentration")+
  xlab("Treatment")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.text = element_text(face = "italic"),
        axis.text.x  = element_text(size = 8),
        panel.background=element_rect(fill="white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 7),
        axis.line=element_line(colour="grey",size=0.8))

dev.off()
