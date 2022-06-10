##Private alleles Cultivars, wild and feral populations
library(hierfstat)

#Load data. Transform to FSTAT format (.dat) using PGSspider. Population file= ../meta/pop_cocci_bayescan.txt
cocci <- read.fstat("../data/cocci_status.dat", na.s = c(0,00,000,0000,00000,000000,NA))

#Change populations names
cocci[cocci$Pop==1,1]<- "Feral"
cocci[cocci$Pop==2,1]<- "Wild"
cocci[cocci$Pop==3,1]<- "Cultivar"

#Estimate allele frequencies
freq <- pop.freq(cocci)
#Find frequencies=0
zeros <- lapply(freq, function(x) which(x==0))
#Find loci with two zeros
two.zeros <- lapply(zeros, function(x) which(length(x) == 2))
length(which(two.zeros==1)) #Amount of private alleles
private_SNPs <- freq[which(two.zeros==1)]

#Find private SNPs of Ferals
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][3] != 0 & private_SNPs[[i]][3] != 1
  aa <- rbind(aa, a)
}
priv.fer <- private_SNPs[which(aa==1)]

#Find private SNPs in wild populations
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][5] != 0 & private_SNPs[[i]][5] != 1
  aa <- rbind(aa, a)
}
priv.wild <- private_SNPs[which(aa==1)]

#Find private SNPs in cultivars
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][1] != 0 & private_SNPs[[i]][1] != 1
  aa <- rbind(aa, a)
}
priv.cult <- private_SNPs[which(aa==1)]

##Extract IDs
write.table(names(priv.fer), "../meta/IDs_priv_fer.txt", row.names = F, col.names = F, quote = F)
write.table(names(priv.wild), "../meta/IDs_priv_wild.txt", row.names = F, col.names = F, quote = F)
write.table(names(priv.cult), "../meta/IDs_priv_cult.txt", row.names = F, col.names = F, quote = F)

