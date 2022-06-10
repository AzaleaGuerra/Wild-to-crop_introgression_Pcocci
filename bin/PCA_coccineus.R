############### PCA - Phaseolus coccineus################
library(SNPRelate)
library(ggplot2)

###Tranform from vcf to gds
snpgdsVCF2GDS("../data/coccineus.recode.vcf", 
              out.fn="../data/pca_cocci.gds")
snpgdsSummary("../data/pca_cocci.gds")

#Loading gds file
cocci <- snpgdsOpen("../data/pca_cocci.gds", allow.duplicate = TRUE)

# Get samples IDs
sample.id.cocci <- read.gdsn(index.gdsn(cocci, "sample.id"))

###Load metadata
meta.coccineus<- read.delim("../meta/cocci_meta.txt", sep = ",")

###PCA
pca.cocci<- snpgdsPCA(cocci)

# % of explaned variation
pc.perc.coc<- pca.cocci$varprop*100
pc.perc.coc<- (round(pc.perc.coc, 2))

## Order the results
pca.coccin <- data.frame(sample.id = pca.cocci$sample.id,
                  EV1 = pca.cocci$eigenvect[,1],    
                  EV2 = pca.cocci$eigenvect[,2],    
                  EV3 = pca.cocci$eigenvect[,3],
                  EV4 = pca.cocci$eigenvect[,4],
                  EV5 = pca.cocci$eigenvect[,5],
                  GenGroup = meta.coccineus$Cluster.Admix,
                  ClusPop = meta.coccineus$Population,
                  Est = meta.coccineus$Estatus2,
                  stringsAsFactors = FALSE)

#Plots
cbPalette.9 <- (c("#bc0021", "#e05384", "#00a4f9", "#edad2f", "#8d360e", "#2567d3", "#774ddd", "#e56329", "#2cb25d"))
ggplot(pca.coccin, aes(x=EV1, y=EV2)) +
  geom_point(aes(col= pca.coccin$GenGroup, size=5)) + 
  xlab(paste0("Eigenvector 1 explains ", pc.perc.coc[1], "%")) +
  ylab(paste0("Eigenvector 2 explains ", pc.perc.coc[2], "%")) +
  scale_colour_manual(values=cbPalette.9) + 
  theme_bw() + 
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20), 
        legend.text=element_text(size=18), legend.title = element_blank(), legend.key = element_blank())+
  guides(size=FALSE)

### Colored by genetic group and population
cbPalette.Pop <- (c("#ffe082", "#fc8014", "#ffb39c", "#e50000","#c79a00", "#b76900", "#757575",
                    "#80deea", "#2196f3", "#ba68c8", "#7349bd", "#9da1db", "#a2cf6e", "#009700", "#30709e"))
ggplot(pca.coccin, aes(x=EV1, y=EV2)) +
  geom_point(aes(col= pca.coccin$ClusPop), size=9) + 
  xlab(paste0("Eigenvector 1 explains ", pc.perc.coc[1], "%")) +
  ylab(paste0("Eigenvector 2 explains ", pc.perc.coc[2], "%")) +
  scale_colour_manual(values= alpha(cbPalette.Pop, 0.7)) + 
  theme_bw(base_size = 27) + #quita el fondo gris
  theme(axis.text=element_text(size=27), axis.title=element_text(size=27), 
        legend.text=element_text(size=27), legend.title = element_blank(), legend.key = element_blank())+
  guides(size=FALSE)


