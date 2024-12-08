
metabo_foldchange<-read.csv("metabo_foldchange.csv")

### exp2: same exposure time
heat_se<-metabo_average%>%
  filter(treatment=='F D'|treatment=='F R'|treatment == 'C D'|treatment == 'C R')

heat_se$treatment[heat_se$treatment=='C D']<-'C'
heat_se$treatment[heat_se$treatment=='C R']<-'C'
heat_se$treatment[heat_se$treatment=='F D']<-'D'
heat_se$treatment[heat_se$treatment=='F R']<-'R'
heat_se$treatment<-factor(heat_se$treatment,levels = c("C","D","R"))
#heat_se<-heat_se[order(heat_se$fold_change.mean, decreasing = T),]
#heat_se$metabolites<-factor(heat_se$metabolites,levels = unique(heat_se$metabolites))

pdf("figure6_heat_se.1.pdf", width = 8.2, height = 5)
ggplot(data = heat_se, aes(x = treatment, y = fold_change.mean, group = species, colour = species)) + geom_line()+
  geom_point() + geom_errorbar(aes(x= treatment, ymin = fold_change.mean-fold_change.se,ymax = fold_change.mean + fold_change.se), stat="identity", width = 0.1)+
  facet_wrap(~metabolites_ordered,scales = "free",nrow = 4, ncol = 8)+
  scale_colour_discrete(name = "Species", labels = c("M. dirhodum","R. padi","S. avenae"))+
  ylab("Fold change in metabolites concentration")+
  xlab("Treatment")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.text = element_text(face = "italic"),
        axis.text  = element_text(size = 8),
        panel.background=element_rect(fill="white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 8),
        axis.line=element_line(colour="grey",size=0.5))
dev.off()