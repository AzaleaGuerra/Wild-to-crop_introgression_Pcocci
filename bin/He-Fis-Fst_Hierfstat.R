library(adegenet)
library(hierfstat)
library(ggplot2)

#Load metadata
meta.cocci <- read.delim("../meta/cocci_meta.txt")
meta.cocci <- meta.cocci[meta.cocci$Subset=="Included",]

#All snps
cocci <-read.delim("../data/cocci_subset.012", sep = "\t", header = FALSE, na.strings = -1)
cocci <- cocci[,-1]
cocci[cocci==0] <- 55
cocci[cocci==1] <- 56
cocci[cocci==2] <- 66

cocci <- cbind(meta.cocci$Locacion, meta.cocci$Population, cocci)
class(cocci) #Debe ser una dataframe

#Pair Fst W&C 1984 according to Cluster-Pop
fst.cocci <- pairwise.WCfst(cocci[,-1], diploid = TRUE)
fst.c <- as.matrix(fst.cocci)
write.table(fst.cocci, file = "../out/Fst_pop_Hierfstat_subset.txt")

############################  HETEROZYGOSYS AND Fis ######################
stats <- basic.stats(cocci[,-2])

#Ho
Ho.locus <-as.data.frame(stats$Ho)
Ho.locus[Ho.locus=="NaN"] <- NA
Ho <-numeric(0)
for (i in 1:ncol(Ho.locus)) {
  H <- cbind(colnames(Ho.locus)[i], mean(Ho.locus[,i], na.rm = TRUE), 
             sd(Ho.locus[,i], na.rm = TRUE), var(Ho.locus[,i], na.rm = TRUE))
  Ho <-rbind(H, Ho)
}
Ho <- Ho[order(Ho[,1]),]

#Hs
Hs.locus <- as.data.frame(stats$Hs)
Hs.locus[Hs.locus=="NaN"] <- NA
Hs <-numeric(0)
for (i in 1:ncol(Hs.locus)) {
  h <- cbind(colnames(Hs.locus)[i], mean(Hs.locus[,i], na.rm = TRUE),
             sd(Hs.locus[,i], na.rm = TRUE), var(Hs.locus[,i], na.rm = TRUE))
  Hs <-rbind(h, Hs)
}
Hs <- Hs[order(Hs[,1]),]

#Inbreeding coefficient
Fis <- boot.ppfis(cocci[,-1], nboot = 100, quant=c(0.025,0.975))

#Bind data
resum <- cbind(Ho, Hs[,-1], Fis$fis.ci)
for (i in 2:ncol(resum)) {
  resum[,i] <- as.numeric(as.character(resum[,i]))
}
colnames(resum) <- c("Pop", "Ho","sd.Ho","var.Ho", "He","sd.He", "var.He", "Fis.low", "Fis.high")
write.table(resum, file = "../out/Het&Fis_subset_location.txt", row.names = FALSE) 

####Perform a Mann-Whitney test for Het
H <- numeric(0)
for (i in 1:ncol(Hs.locus)) {
  H <- rbind(H, cbind(as.numeric(Hs.locus[,i]), (colnames(Hs.locus)[i])))
}
H <- as.data.frame(H)
H$V1 <- as.numeric(as.character(H$V1))

kruskal.hs <- kruskal.test(x=H[,1], g= as.factor(H[,2]))
posthoc.hs <- capture.output(pairwise.wilcox.test(x=H[,1], g=as.factor(H[,2]), 
              paired = TRUE, p.adjust.method="bonferroni"), file = "../data/H_subset.posthoc")
boxplot(V1 ~ V2, data = H)


############################### PLOT #########################
###################### PLOT ##################
#Bind het and pi data
resum <-read.delim("../out/Het&Fis_Hierfstat.txt", sep = "\t")

#Plot Fis of all snps
ggplot(data=resum) +
  geom_point(aes(x=Pop, y= He), shape = 16, size=6, colour= "#323dce") +
  geom_errorbar(aes(x=Pop, ymin = He-var.He, ymax = He+var.He, width = 0.6), colour= "#323dce", size=1.9) +
  geom_point(aes(x=Pop, y= Ho), shape = 16, size=6, colour= "#999ee7") +
  geom_errorbar(aes(x=Pop, ymin = Ho-var.Ho, ymax = Ho+var.Ho, width = 0.6), colour= "#999ee7", size=1.9) +
  geom_errorbar(aes(x=Pop, ymin = Fis.low, ymax = Fis.high, width = 0.6), colour= "#ca5670", size=1.9) +
  theme_bw(base_size = 33) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank()) 

#Heatmap
library(gplots)
rownames(resum) <- resum$Pop
resum <- as.matrix(resum[,c(5,2,10,10)])

my_palette <- colorRampPalette(c("#006ae8", "#fffecd", "#cb0000"))(n = 299)
heatmap.2(resum[,-c(3,4)],
          dendrogram = "none",
          #na.rm = TRUE,
          #key = TRUE,
          keysize = 1.7,
          #key.xlab = expression(italic('F'[ST])),
          #key.par=list(mgp=c(1.5, 0.5, 0),
          #             mar=c(7, 9, 4, 0)),
          #scale = "none",
          Rowv = FALSE,
          Colv = FALSE,
          revC = FALSE,
          cexRow = 1.9,
          cexCol = 1.5,
          cellnote = resum[,-c(3,4)],
          notecol = "black",
          notecex = 1.5,
          #adjRow="right",
          density.info= "none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,20),     # widens margins around plot
          col=my_palette)       # use on color palette defined earlier

my_palette <- colorRampPalette(c("#fff875", "#02a8ad", "#540070"))(n = 299)
heatmap.2(resum[,c(3,4)],
          dendrogram = "none",
          keysize = 1.7,
          Rowv = FALSE,
          Colv = FALSE,
          revC = FALSE,
          cexRow = 1.7,
          cexCol = 1.5,
          cellnote = resum[,c(3,4)],
          notecol = "black",
          notecex = 1.5,
          #adjRow="right",
          density.info= "none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,20),     # widens margins around plot
          col=my_palette)       # use on color palette defined earlier

