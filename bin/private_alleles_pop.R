##Private alleles Cultivars, wild and feral populations
library(hierfstat)

#Load data. Transform to FSTAT format (.dat) using PGSspider. Population file= ../meta/pop_cocci_pgdspider.txt
cocci <- read.fstat("../data/cocci_pop.dat", na.s = c(0,00,000,0000,00000,000000,NA))

#Change populations names
cocci[cocci$Pop==1,1]<- "Feral"
cocci[cocci$Pop==2,1]<- "W-SUR-O"
cocci[cocci$Pop==3,1]<- "C-TMVB"
cocci[cocci$Pop==4,1]<- "C-Spa"
cocci[cocci$Pop==5,1]<- "W-CDMX"
cocci[cocci$Pop==6,1]<- "W-SaJo"
cocci[cocci$Pop==7,1]<- "W-Tepo"
cocci[cocci$Pop==8,1]<- "W-striatus"
cocci[cocci$Pop==9,1]<- "C-BlaTla"
cocci[cocci$Pop==10,1]<- "C-OV"
cocci[cocci$Pop==11,1]<- "C-SMOCC"
cocci[cocci$Pop==12,1]<- "W-EsDi"
cocci[cocci$Pop==13,1]<- "C-Sur"
cocci[cocci$Pop==14,1]<- "W-Rego"
cocci[cocci$Pop==15,1]<- "W-Sur-Ch"

#Estimate allele frequencies
freq <- pop.freq(cocci)

#Find frequencies=0
zeros <- lapply(freq, function(x) which(x==0))
#Find loci with 14 zeros, which means that are polymprphic in only one population
zeros <- lapply(zeros, function(x) which(length(x) == 14))
length(which(zeros==1)) #Amount of private alleles
private_SNPs <- freq[which(zeros==1)]

#Find private SNPs of Ferals
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][14] != 0 & private_SNPs[[i]][14] != 1
  aa <- rbind(aa, a)
}
fer <- private_SNPs[which(aa==1)]

#C-BlaTla
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][1] != 0 & private_SNPs[[i]][1] != 1
  aa <- rbind(aa, a)
}
C.BlaTla <- private_SNPs[which(aa==1)]

#C-OV
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][3] != 0 & private_SNPs[[i]][3] != 1
  aa <- rbind(aa, a)
}
C.OV <- private_SNPs[which(aa==1)]

#C-SMOCC
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][5] != 0 & private_SNPs[[i]][5] != 1
  aa <- rbind(aa, a)
}
C.smocc <- private_SNPs[which(aa==1)]

#C-Spa
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][7] != 0 & private_SNPs[[i]][7] != 1
  aa <- rbind(aa, a)
}
C.Spa <- private_SNPs[which(aa==1)]

#C-Sur
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][9] != 0 & private_SNPs[[i]][9] != 1
  aa <- rbind(aa, a)
}
C.Sur <- private_SNPs[which(aa==1)]

#C-TMVB
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][11] != 0 & private_SNPs[[i]][11] != 1
  aa <- rbind(aa, a)
}
C.tmvb <- private_SNPs[which(aa==1)]

#W-cdmx
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][15] != 0 & private_SNPs[[i]][15] != 1
  aa <- rbind(aa, a)
}
w.cdmx <- private_SNPs[which(aa==1)]

#W-EspDia
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][17] != 0 & private_SNPs[[i]][17] != 1
  aa <- rbind(aa, a)
}
w.edia <- private_SNPs[which(aa==1)]

#W-Rego
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][19] != 0 & private_SNPs[[i]][19] != 1
  aa <- rbind(aa, a)
}
w.rego <- private_SNPs[which(aa==1)]

#W-SanJoa
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][21] != 0 & private_SNPs[[i]][21] != 1
  aa <- rbind(aa, a)
}
w.sanjoa <- private_SNPs[which(aa==1)]

#W-striatus
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][23] != 0 & private_SNPs[[i]][23] != 1
  aa <- rbind(aa, a)
}
w.str <- private_SNPs[which(aa==1)]

#W-Sur-Ch
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][25] != 0 & private_SNPs[[i]][25] != 1
  aa <- rbind(aa, a)
}
w.sur.ch <- private_SNPs[which(aa==1)]

#W-Sur-O
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][27] != 0 & private_SNPs[[i]][27] != 1
  aa <- rbind(aa, a)
}
w.sur.o <- private_SNPs[which(aa==1)]

#W-Tepo
aa <- numeric(0)
for (i in 1:length(private_SNPs)) {
  a <- private_SNPs[[i]][29] != 0 & private_SNPs[[i]][29] != 1
  aa <- rbind(aa, a)
}
w.tepo <- private_SNPs[which(aa==1)]

