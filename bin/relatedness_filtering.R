#Reltedness was estimated with vcftools:
#vcftools --vcf ../data/coccineus.vcf --relatedness2 --out ../data/cocci

#Load relatedness data estimated with vcftools
relat <- read.delim("../data/cocci.relatedness2")

#Load metadata
meta <- read.delim("../meta/cocci_meta.txt", sep = "\t")
meta <- meta[meta$Subset=="Included",]

#Merge data according to sample ID
data <- merge(relat, meta[,c(2,13)], by.x = "INDV1", by.y = "ID")
data <- merge(data, meta[,c(2,13)], by.x = "INDV2", by.y = "ID")

#Keep only rows where both ind are from the same pop
data <- data[data$Population.x == data$Population.y,]

#Remove self comparations
data <- data[data$INDV2 != data$INDV1,]

#Keep the ones with high relatedness
data.high <- data[data$RELATEDNESS_PHI > 0.05,]



