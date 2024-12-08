library(dplyr)
library(ade4)
library(plotrix)
library(ggrepel)

## experiment 2: same exposure time with Control
metabo_exp2<-read.delim("metabo_exp2.txt",header = T,stringsAsFactors = T)

## color setting for each treatment of experiment 2
cols<-c ("#F8766D","#F8766D","#F8766D", "#00BA38","#00BA38","#00BA38","#619CFF","#619CFF","#619CFF")

##colour setting for metabolites that correlated to PC
redorange<-colorRampPalette(c("orange","goldenrod","darkred"))
barplot(1:32,pch =16, col = redorange(32))

### 2.2.1 PCA (figure 3) experiment2:same exposure time
### set colors

F_34_C<- metabo_exp2[,7:39]
F_34_C$TRT<-as.factor(F_34_C$TRT)
summary(F_34_C$TRT)

mol <- F_34_C[,c(2:33)]

###between PCA analysis
pca=dudi.pca(log(mol), scannf=F, nf=31)
be=bca(pca, fac = F_34_C$TRT, scannf=F, nf=31)

#PC1 vs PC2
pcdata<-cbind(be$ls[,1:2],str_sub(F_34_C$TRT,1,2),str_sub(F_34_C$TRT,-3,5))
colnames(pcdata)<-c("PC1","PC2","species","treatment")
pcdata$treatment[pcdata$treatment=="20D"]<-"C"
pcdata$treatment[pcdata$treatment=="20R"]<-"C"
pcdata$treatment[pcdata$treatment=="34R"]<-"R"
pcdata$treatment[pcdata$treatment=="34D"]<-"D"

###plot figure 3A:  between PCA 
pdf(file="fig3_se_pca.pdf",width = 6, height = 5)
par(mar=c(4,4,4,4),tck = -0.015, mgp = c(3, 0.5, 0))
plot(PC2~PC1,pch = c(16,17,18)[as.numeric(as.factor(treatment))], col=c ("#F8766D" ,"#00BA38","#619CFF")[as.numeric(as.factor(species))] , data= pcdata, cex=1, main= NULL, xlab ="", ylab = "")

legend("topright",legend = c("C", "D", "R"), title ="Treatment",pch = c(16,17,18))
abline (v=0, col='gray')
abline (h=0, col='gray')

s.class(be$ls[,1:2], xax=1, yax=2, fac=F_34_C$TRT, cellipse=1, axesell=F, cstar=1, clabel=0.5, pch="10", cpoint=0.5, add.plot=T,col= cols)
# Add custom x-axis title with specific distance
mtext(side = 1, line = 1.5, text = paste("PC1 (", paste(round((be$eig*100/sum(be$eig))[1], 2), "%"),")"),cex = 1)
# Add custom y-axis title
mtext(side = 2, line = 1.5, text = paste("PC2 (", round((be$eig*100/sum(be$eig))[2], 2), "%", ")"), cex = 1)
dev.off()

### plot figure 3B: PC2 values for treatments of each species
pc2data<-pcdata%>%group_by(species,treatment)%>%summarize(mean.pc2=mean(PC2),pc2.se=std.error(PC2))
pc2data$species[pc2data$species=="MD"]<-"M. dirhodum"
pc2data$species[pc2data$species=="SA"]<-"S. avenae"
pc2data$species[pc2data$species=="RP"]<-"R. padi"
pc2data$species<-factor(pc2data$species,levels = c("M. dirhodum","S. avenae","R. padi"))

df<-data.frame(species = c("M. dirhodum","S. avenae","R. padi"))
df$species<-factor(df$species,levels=c("M. dirhodum","S. avenae","R. padi"))

pdf("fig3b_pc2_se.pdf", width = 5, height = 4)
ggplot(pc2data) + 
  geom_line(aes(x = treatment, y = mean.pc2, group = species, color = species), size = 1.2)+
  geom_point(aes(x = treatment, y = mean.pc2, color = species),size=2) + 
  geom_errorbar(aes(x= treatment, ymin = mean.pc2-pc2.se,ymax = mean.pc2 + pc2.se, color = species), stat="identity", width = 0.15, size =0.7)+
  facet_wrap(~species, ncol = 3)+
  scale_color_manual(values = c("#F8766D" ,"#619CFF","#00BA38"))+
  ylab("PC2")+ xlab("Treatment")+
  ylim(-3,4.8)+
  theme(legend.position = "none",
        strip.background = element_rect(fill=NA),
        #legend.text = element_text(size = 14),
        panel.background=element_rect(fill="white"),
        axis.line=element_line(colour="black",linewidth = 1),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12,color="black"),
        strip.text = element_text(size = 14,face="italic")
        #axis.text.x = element_text(face = "italic",colour = Colors),
        #plot.title=element_text(hjust = 0.5, size = 18)
  )
dev.off()

### 2.2.1 correlations with pc1 and pc2(figure s2)
# Extract scores and loadings
loadings <- as.data.frame(be$co)

pdf("figs2a_correlation_se_1.pdf", width = 4, height =3,units = "in",res = 500)

ggplot() +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = Comp1, yend = Comp2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = loadings, aes(x = Comp1, y = Comp2, label = rownames(loadings)), 
                  size = 3, color = "red") +
  xlim(-1,1) + ylim(-1,1)+
  geom_segment(aes(x = 0, y = -1, xend = 0, yend = 1), color = "darkgrey",linetype="dashed" )+
  geom_segment(aes(x = -1, y = 0, xend = 1, yend = 0), color = "darkgrey",linetype="dashed" )+
  geom_segment(aes(x = -1, y = -1, xend = 1, yend = -1), color = "black") +  # Bottom line
  geom_segment(aes(x = -1, y = 1, xend = 1, yend = 1), color = "black") +    # Top line
  geom_segment(aes(x = -1, y = -1, xend = -1, yend = 1), color = "black") +  # Left line
  geom_segment(aes(x = 1, y = -1, xend = 1, yend = 1), color = "black") +  
  scale_x_continuous(sec.axis = dup_axis(),expand = c(0, 0)) +
  scale_y_continuous(sec.axis = dup_axis(),expand = c(0, 0)) +
  theme( panel.grid = element_blank(),
         panel.background = element_blank(),
         legend.position = "top", legend.title = element_text(size = 12), 
         legend.text = element_text(size = 12), 
         axis.text = element_text(size = 12),
         axis.title = element_blank())

dev.off()

#barplot: correlation value of PC1 and PC2 
pdf(file="figs2b_correlation PC1_SE.pdf",width = 2, height = 3,units = "in",res = 500)
par(mar=c(3,8,2,1))
barplot(sort(be$co[,1]),names.arg = row.names(be$co[order(be$co[,1]),]),horiz = TRUE, col=redorange(32), xlim=c(-1,1), ylab="",las=1.5, cex.names=0.8)
abline(v = 0.5 ,col = "grey", lty=2)
abline(v = -0.5 ,col = "grey", lty=2)
dev.off()

pdf(file="figs2c_correlation PC2_SE.pdf",width = 2, height = 3,units = "in",res = 500)
par(mar=c(3,8,2,1))
barplot(sort(be$co[,2]),names.arg = row.names(be$co[order(be$co[,2]),]),horiz = TRUE, col=redorange(32), xlim=c(-1,1), ylab="",las=1.5, cex.names=0.8)
abline(v = 0.5 ,col = "grey", lty=2)
abline(v = -0.5 ,col = "grey", lty=2)
dev.off()
