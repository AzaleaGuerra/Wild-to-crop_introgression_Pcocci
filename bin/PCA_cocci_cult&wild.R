############### PCA de Phaseolus coccineus################
library(SNPRelate)

###Tranform from vcf to gds
snpgdsVCF2GDS("../data/cult.recode.vcf",
              out.fn="../data/cult.gds")
snpgdsSummary("../data/cult.gds")

#Load gds file
cult <- snpgdsOpen("../data/cult.gds", allow.duplicate = TRUE)

#Get samples IDs
sample.id.cult <- read.gdsn(index.gdsn(cult, "sample.id"))

###Load metadata
meta <- read.delim("../meta/cocci_meta.txt")
meta.cult <- meta[meta$Estatus2 == "Cultivar" ,]

###PCA
pca.cult<- snpgdsPCA(cult)
pc.perc.cu<- pca.cult$varprop*100
pc.perc.cu<- (round(pc.perc.cu, 2))

## Order the results
pca.c <- data.frame(sample.id = pca.cult$sample.id,
                         EV1 = pca.cult$eigenvect[,1],    
                         EV2 = pca.cult$eigenvect[,2],    
                         EV3 = pca.cult$eigenvect[,3],
                         Pob = meta.cult$Population,   
                         stringsAsFactors = FALSE)
#Plot
cPalette <- c("#ffe082", "#fc8014", "#ffb39c", "#e50000","#c79a00", "#b76900", "#757575")
##According to genetic group
ggplot(pca.c, aes(x=EV1, y=EV2)) +
  geom_point(aes(col= pca.c$Pob), size=9) + 
  xlab(paste0("Eigenvector 1 explains ", pc.perc.cu[1], "%")) +
  ylab(paste0("Eigenvector 2 explains ", pc.perc.cu[2], "%")) +
  scale_colour_manual(values = alpha(cPalette, 0.7)) + 
  theme_bw(base_size = 29) + #quita el fondo gris
  theme(axis.text=element_text(size=29), axis.title=element_text(size=29), 
          legend.text=element_text(size=29), legend.title = element_blank())+
  guides(size=FALSE)

################### PCA WILD ######################################
snpgdsVCF2GDS("../data/wild.recode.vcf",
              out.fn="../data/wild.gds")
snpgdsSummary("../data/wild.gds")

#Load gds
wild <- snpgdsOpen("../data/wild.gds", allow.duplicate = TRUE)

# Samples IDs
sample.id.wild <- read.gdsn(index.gdsn(wild, "sample.id"))

###Metadata
meta<- read.delim("../meta/cocci_meta.txt")
meta.wild <- meta[meta$Estatus2 == "Wild" ,]

###PCA
pca.wild<- snpgdsPCA(wild)
pc.perc.w<- pca.wild$varprop*100
pc.perc.w<- (round(pc.perc.w, 2))

## Order the results
pca.w <- data.frame(sample.id = pca.wild$sample.id,
             EV1 = pca.wild$eigenvect[,1],    
             EV2 = pca.wild$eigenvect[,2],    
             EV3 = pca.wild$eigenvect[,3],
             Pob = meta.wild$Population,   
             stringsAsFactors = FALSE)

#Plot
wPalette <- c("#80deea", "#2196f3", "#ba68c8", "#7349bd", "#9da1db", "#a2cf6e", "#009700", "#30709e")
ggplot(pca.w, aes(x=EV1, y=EV2)) +
  geom_point(aes(col= pca.w$Pob), size=9) + 
  xlab(paste0("Eigenvector 1 explains ", pc.perc.w[1], "%")) +
  ylab(paste0("Eigenvector 2 explains ", pc.perc.w[2], "%")) +
  scale_colour_manual(values= alpha(wPalette, 0.7)) + 
  theme_bw(base_size = 29) + #quita el fondo gris
  theme(axis.text=element_text(size=29), axis.title=element_text(size=29), 
        legend.text=element_text(size=29), legend.title = element_blank(), legend.key = element_blank())+
  guides(size=FALSE)

