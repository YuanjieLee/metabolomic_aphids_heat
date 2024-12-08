library(dplyr)
library(ade4)


metabo_exp1<-read.delim("metabo_exp1.txt",header = T,stringsAsFactors = T)
head(metabo_exp1)
###split data according to species##
data_MD_34<-metabo_exp1%>%filter(species == "MD")
data_RP_34<-metabo_exp1%>%filter(species == "RP")
data_SA_34<-metabo_exp1%>%filter(species == "SA")
data_MD_34<-data_MD_34[,c(7:39)] 
data_RP_34<-data_RP_34[,c(7:39)]
data_SA_34<-data_SA_34[,c(7:39)]
data_MD_34$TRT<-as.factor(data_MD_34$TRT)
data_RP_34$TRT<-as.factor(data_RP_34$TRT)
data_SA_34$TRT<-as.factor(data_SA_34$TRT)

#### Colour setting for treatments of experiment 1
mycols<-c("#ffa836","#98bF64","#990F02","#6485d1","#000000")
##colour setting for metabolites that correlated to PC
redorange<-colorRampPalette(c("orange","goldenrod","darkred"))
barplot(1:32,pch =16, col = redorange(32))

##experiment1: MD at 34C
mol <- data_MD_34[,c(2:33)]
summary(data_MD_34$TRT)
levels(data_MD_34$TRT)
pca = dudi.pca(log(mol), scannf=F, nf=31)

### Between PCA ----> reveal  the  differences  between groups.

be=bca(pca, fac = data_MD_34$TRT, scannf=F, nf=31)
valeur_propre =be$eig*100/sum(be$eig)
barplot (valeur_propre, main="% inertia per component")

## monte carlo test which indicates if the discrimination between group is significant or not
rand1 <- randtest(bca(pca, data_MD_34$TRT, scan = FALSE), 1000)
rand1
plot(rand1, main = "Monte-Carlo test", 
     sub="Simulated p-value = XXX (1000 permutations)", cex.sub=0.6)


#PC1 vs PC2
pdf(file="md_pca.pdf",width = 4.5, height = 4)
par(tck = -0.015, mgp = c(3, 0.5, 0))
plot(be$ls[,1:2], cex=0, main= NULL, xlab="", ylab = "",cex.axis=0.8)
abline (v=0, col='gray')
abline (h=0, col='gray')
s.class(be$ls[,1:2], xax=1, yax=2, fac=data_MD_34$TRT, cellipse=1, axesell=F, cstar=1, clabel=0.8, pch=10, cpoint=0.5, add.plot=T,col=mycols)
text(x = 3.3, y = 4.3, labels = substitute(italic('M. dirhodum')), cex = 1.1) 
# Add custom x-axis title with specific distance
mtext(side = 1, line = 1.5, text = paste("PC1 (", paste(round((be$eig*100/sum(be$eig))[1], 2), "%"),")"),cex = 1)
# Add custom y-axis title
mtext(side = 2, line = 1.5, text = paste("PC2 (", round((be$eig*100/sum(be$eig))[2], 2), "%", ")"), cex = 1)

dev.off()

###plot the correlations of each metabolite with PC1 and PC2
s.arrow(be$co, clab = 0.6, xax=1, yax=2, sub="correlations between-PCA [PC1vsPC2]", csub=0.7)

## extract correlation values to PC1 and PC2 of each metabolite 
md_co<-be$co[,1:2]
colnames(md_co)<-c("md_pc1","md_pc2")

###extract pc1 vaules of each sample 
md_pc1<-cbind(rep("M. dirhodum",38),as.character(data_MD_34$TRT), be$ls[,1])


##experiment1: RP at 34C

mol <- data_RP_34[,c(2:33)]
summary(data_RP_34$TRT)

pca=dudi.pca(log(mol), scannf=F, nf=31)

### Between PCA ----> reveal  the  differences  between groups.
be=bca(pca, fac = data_RP_34$TRT, scannf=F, nf=31)
## graphs: explaination of each PC
valeur_propre =be$eig*100/sum(be$eig)
barplot (valeur_propre, main="% inertia per component")

## monte carlo test which indicates if the discrimination between group is significant or not
rand1 <- randtest(bca(pca, data_RP_34$TRT, scan = FALSE), 1000)
rand1
plot(rand1, main = "Monte-Carlo test", 
     sub="Simulated p-value = 0.0009 (1000 permutations)", cex.sub=0.6)

###plot between pca
pdf(file="rp_pca.pdf",width = 4.5, height = 4)
par(tck = -0.015, mgp = c(3, 0.5, 0))
plot(be$ls[,1:2], cex=0, main= NULL, xlab="", ylab = "",cex.axis=0.8)
abline (v=0, col='gray')
abline (h=0, col='gray')
s.class(be$ls[,1:2], xax=1, yax=2, fac=data_RP_34$TRT, cellipse=1, axesell=F, cstar=1, clabel=0.8, pch=10, cpoint=0.5, add.plot=T,col=mycols)
text(x = 4, y = 3.3, labels = substitute(italic('R. padi')), cex = 1.1) 
# Add custom x-axis title with specific distance
mtext(side = 1, line = 1.5, text = paste("PC1 (", paste(round((be$eig*100/sum(be$eig))[1], 2), "%"),")"),cex = 1)
# Add custom y-axis title
mtext(side = 2, line = 1.5, text = paste("PC2 (", round((be$eig*100/sum(be$eig))[2], 2), "%", ")"), cex = 1)
dev.off()

s.arrow(be$co, clab = 0.8, xax=1, yax=2, sub="correlations between-PCA [PC1vsPC2]", csub=0.7)

## extract correlation values to PC1 and PC2 of each metabolite
rp_co<-be$co[,1:2]

colnames(rp_co)<-c("rp_pc1","rp_pc2")

### extract PC1 values for each sample
rp_pc1<-cbind(rep("R. padi",40),as.character(data_RP_34$TRT),be$ls[,1])

# experiment1: SA at 34C

mol <- data_SA_34[,c(2:33)]
summary(data_SA_34$TRT)

pca=dudi.pca(log(mol), scannf=F, nf=31)

### Between PCA ----> reveal  the  differences  between groups.

be=bca(pca, fac = data_SA_34$TRT, scannf=F, nf=31)

valeur_propre =be$eig*100/sum(be$eig)
barplot (valeur_propre, main="% inertia per component")

## monte carlo test which indicates if the discrimination between group is significant or not

rand1 <- randtest(bca(pca, data_SA_34$TRT, scan = FALSE), 1000)
rand1
plot(rand1, main = "Monte-Carlo test", 
     sub="Simulated p-value = .00099 (1000 permutations)", cex.sub=0.6)

###plot between pca
pdf(file="sa_pca.pdf",width = 4.5, height = 4)
par(tck = -0.015, mgp = c(3, 0.5, 0))
plot(be$ls[,1:2], cex=0, main= NULL, xlab="", ylab = "",cex.axis = 0.8)
abline (v=0, col='gray')
abline (h=0, col='gray')
s.class(be$ls[,1:2], xax=1, yax=2, fac=data_SA_34$TRT, cellipse=1, axesell=F, cstar=1, clabel=0.8, pch=10, cpoint=0.5, add.plot=T,col=mycols)
text(x = 4, y = 4.5, labels = substitute(italic('S. avenae')), cex = 1.1) 
# Add custom x-axis title with specific distance
mtext(side = 1, line = 1.5, text = paste("PC1 (", paste(round((be$eig*100/sum(be$eig))[1], 2), "%"),")"),cex = 1)
# Add custom y-axis title
mtext(side = 2, line = 1.5, text = paste("PC2 (", round((be$eig*100/sum(be$eig))[2], 2), "%", ")"), cex = 1)
dev.off()


s.arrow(be$co, clab = 0.8, xax=1, yax=2, sub="correlations between-PCA [PC1vsPC2]", csub=0.7)

## extract correlation values to PC1 and PC2 of each metabolite
sa_co<-be$co[,1:2]
colnames(sa_co)<-c("sa_pc1","sa_pc2")

## extract pc1 values of each sample
sa_pc1<-cbind(rep("S.avenae",39),as.character(data_SA_34$TRT),be$ls[,1])


#### figure1b. PC1 values- experiment 1: summary of PC1 values of each sample
#pc1<-rbind(md_pc1,sa_pc1,rp_pc1)
#colnames(pc1)<-c("species","treatment","pc1")
#pc1<-as.data.frame(pc1)
#write.csv(pc1,"pc1.1.csv")
##in this csv, i add a new column "injury level" manually
pc1<-read.csv("pc1.csv",header=T)
pc1<-pc1%>%group_by(species,injurylevel,treatment)%>%summarise(mean=mean(pc1),se= std.error(pc1))

pc1$species<-factor(pc1$species,levels=c("M. dirhodum","S. avenae","R. padi"))
df<-data.frame(species = c("M. dirhodum","S. avenae","R. padi"))
df$species<-factor(df$species,levels=c("M. dirhodum","S. avenae","R. padi"))

pdf(file="plot_pc1.pdf",width = 2.7, height = 5.5)
ggplot(as.data.frame(pc1)) + 
  geom_line(aes(x = treatment, y = mean,group= injurylevel,linetype=injurylevel), size = 0.4)+
  geom_point(aes(x = treatment, y = mean,group= injurylevel),size = 0.7) + 
  geom_errorbar(aes(x= treatment, ymin = mean-se,ymax = mean + se), stat="identity", width = 0.05)+
  facet_grid2(species~.,scales = "free_y",
              axes = "x", space = "free_y", switch = "y")+
  coord_cartesian(clip="off", ylim=c(-4.7, 5))+
  scale_linetype_manual(name = "Injury level", values = c(1,2),labels = c("1/2","2/2"))+
  ylab("PC1")+ xlab("Treatment")+
  guides(fill = F)+
  theme(legend.title = element_text(size = 8),
        strip.background = element_rect(fill=NA),
        legend.text = element_text(size = 8),
        panel.background=element_rect(fill="white"),
        legend.position=c(0.5,0.99),
        legend.direction = "horizontal",
        axis.line=element_line(colour="black",linewidth = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8,color="black"),
        strip.text = element_blank())

dev.off()




### figure 1b.statistics of PC1: experiment1: two-way anova}
#compare the control, direct
pc1<-read.csv("pc1.csv",header=T)
data.md<-pc1%>%filter(species == "M. dirhodum") 
data.sa<-pc1%>%filter(species == "S. avenae")
data.rp<-pc1%>%filter(species == "R. padi")

md.aov<-aov(pc1 ~ treatment + injurylevel + treatment*injurylevel,data = data.md)
sa.aov<-aov(pc1 ~ treatment + injurylevel + treatment*injurylevel,data = data.sa)
rp.aov<-aov(pc1 ~ treatment + injurylevel + treatment*injurylevel,data = data.rp)
summary(md.aov)
summary(sa.aov)
summary(rp.aov)##for rp, the interaction is not significant

TukeyHSD(md.aov)
TukeyHSD(sa.aov)
TukeyHSD(rp.aov)


### figure 1c.correlation heatmap

co_pc1<-cbind(sa_co,rp_co,md_co,metabolites=rownames(sa_co))%>%
  dplyr::select(sa_pc1,md_pc1,rp_pc1,metabolites)%>%mutate(mean_pc1= (sa_pc1+md_pc1+rp_pc1)/3)%>%gather(.,species,correlation_value,sa_pc1:rp_pc1)

co_pc.1<-co_pc1[order(co_pc1$mean_pc,decreasing = F),] ##get the order of metabolites by the value
co_pc1$metabolites<-factor(co_pc1$metabolites,levels = unique(co_pc.1$metabolites))
co_pc1$species<-factor(co_pc1$species,levels = c("md_pc1","sa_pc1","rp_pc1"))

pdf(file="pc1_correlation_heatmap.pdf",width = 2.8, height = 5.6)
ggplot(data = co_pc1, aes(x=species, y=metabolites, fill=correlation_value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation to PC1")+
  scale_x_discrete(labels=c("M. dirhodum","S. avenae","R. padi"))+
  geom_text(aes(species,metabolites, label = round(correlation_value,2)),
            color = "black", size =2)+
  theme(legend.position = "top",
        axis.text.x = element_text(face = "italic",size = 8, colour = "black"),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title= element_blank())

dev.off()








