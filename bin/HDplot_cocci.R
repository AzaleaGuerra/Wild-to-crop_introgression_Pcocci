library(vcfR)
library(ggplot2)
library(dplyr)
library(stringr)

#HDplot function
HDPlot_fast<-function(vcfData){
  #set up results table
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=9))
  colnames(HDplotTable)<-c("Locus_ID","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z")
  #format allele reads from vcf data into matrix of comma separated values
  genos<-apply(vcfData@gt[,2:dim(vcfData@gt)[2]], 2, function(x) str_split_fixed(x,":",3)[,1])
  rownames(genos)<-vcfData@fix[,3]
  reads<-apply(vcfData@gt[,2:dim(vcfData@gt)[2]], 2, function(x) str_split_fixed(x,":",4)[,2])
  rownames(reads)<-vcfData@fix[,3] 
  #replace . with 0
  reads<-gsub("\\.","0",reads)
  alleleReads_1<-apply(reads,2,function(x) str_split_fixed(x,",",2)[,1])
  alleleReads_2<-apply(reads,2,function(x) str_split_fixed(x,",",2)[,2])
  #convert to numeric format
  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))
  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))
  rownames(alleleReads_1)<-vcfData@fix[,3]
  rownames(alleleReads_2)<-vcfData@fix[,3]
  #subset to heterozygous genotypes
  #make genotype matrix where heterozygotes are 1 and other genotypes are 0
  hetMatrix<-genos
  hetMatrix<-apply(hetMatrix,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))
  #multiply read count matrices by heterozygote matrix to get read counts for heterozygotes
  alleleReads_1_het<-alleleReads_1*hetMatrix
  alleleReads_2_het<-alleleReads_2*hetMatrix
  #rows are loci and columns are samples
  #sum reads per allele per locus for heterozygous samples
  A_reads<-apply(alleleReads_1_het,1,sum)
  B_reads<-apply(alleleReads_2_het,1,sum)
  totalReads<-A_reads+B_reads
  ratio<-A_reads/totalReads
  std<-sqrt(totalReads*0.5*0.5)
  z<- -(totalReads/2-A_reads)/std
  #get percent heterozygosity for each locus
  numHets<-apply(hetMatrix,1,sum)
  hetPerc<-numHets/dim(hetMatrix)[2]

  #assign results to HDplotTable
  HDplotTable$Locus_ID<-vcfData@fix[,3]
  HDplotTable$depth_a<-A_reads
  HDplotTable$depth_b<-B_reads
  HDplotTable$ratio<-ratio
  HDplotTable$num_hets<-numHets
  HDplotTable$num_samples<-dim(hetMatrix)[2]
  HDplotTable$het_perc<-hetPerc
  HDplotTable$std<-std
  HDplotTable$z<-z
  return(HDplotTable)
}

##################################### ALL SAMPLES (wild, cultivated and feral samples) ###########################
#load VCF file
vcf.cocci<-read.vcfR("../data/Prefiltering/cocci_site_filt.vcf")
#Run HDplot and store results
HDPlot.cocci<-HDPlot_fast(vcf.cocci)

#plot results of HDplot without any thresholds 
#traditional HDplot
ggplot()+geom_point(data=HDPlot.cocci,aes(x=het_perc,y=z),alpha=0.2)+xlab("H")+ylab("D")+theme_bw()
#HDplot with ratio instead of deviation
ggplot()+geom_point(data=HDPlot.cocci,aes(x=het_perc,y=ratio),alpha=0.2)+xlab("H")+ylab("Ratio")+theme_bw()

#Thresholds at this point are determined by visual examination of the plot.
thresh_H<-0.60
thresh_H_divDup<-0.87
thresh_D<-4
thresh_Dneg<--4

#Add paralog status to HDplot results table based on thresholds above
paralogStatus<-function(data,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup){
  #if(data$het_perc>thresh_H|data$z>thresh_D|data$z<thresh_Dneg){
  H<-as.numeric(data[7])
  D<-as.numeric(data[9])
  if(is.na(D)){
    paralog<-NA
  }else if(H>=thresh_H_divDup){
    paralog<-2
  }else if(H>=thresh_H|D>thresh_D|D<thresh_Dneg){
    paralog<-1
  }else{
    paralog<-0
  }
  return(paralog)
}
HDPlot.cocci$paralog<-apply(HDPlot.cocci,1,paralogStatus,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup)
#replot results color coded by paralog status
#traditional HDplot
ggplot()+geom_point(data=HDPlot.cocci,aes(x=het_perc,y=z,color=paralog),alpha=0.2)+xlab("H")+ylab("D")+theme_bw()
#HDplot with ratio instead of deviation
ggplot()+geom_point(data=HDPlot.cocci,aes(x=het_perc,y=ratio,color=paralog),alpha=0.2)+xlab("H")+ylab("Ratio")+theme_bw()

###Extract info and export it
HDPlot.cocci <- na.omit(HDPlot.cocci)
paral.cocci <- na.omit(HDPlot.cocci[HDPlot.cocci$paralog == 1,])
paral.div.cocci <- na.omit(HDPlot.cocci[HDPlot.cocci$paralog == 2,])

write.table(paral.cocci, file = "../data/paral_cocci.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(paral.div.cocci, file = "../data/paral_div_cocci.txt", sep = "\t", quote = F, row.names = F, col.names = T)


##################################### WILD SAMPLES ###########################
#load VCF file
vcf.wild<-read.vcfR("../data/Prefiltering/wild_site_filt.recode.vcf")
#Run HDplot and store results
HDPlot.wild<-HDPlot_fast(vcf.wild)

#plot results of HDplot without any thresholds 
#traditional HDplot
ggplot()+geom_point(data=HDPlot.wild,aes(x=het_perc,y=z),alpha=0.2)+xlab("H")+ylab("D")+theme_bw()
#HDplot with ratio instead of deviation
ggplot()+geom_point(data=HDPlot.wild,aes(x=het_perc,y=ratio),alpha=0.2)+xlab("H")+ylab("Ratio")+theme_bw()

#Thresholds at this point are determined by visual examination of the plot.
thresh_H<-0.62
thresh_H_divDup<-0.95
thresh_D<-7
thresh_Dneg<--7

#Add paralog status to HDplot results table based on thresholds above
paralogStatus<-function(data,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup){
  #if(data$het_perc>thresh_H|data$z>thresh_D|data$z<thresh_Dneg){
  H<-as.numeric(data[7])
  D<-as.numeric(data[9])
  if(is.na(D)){
    paralog<-NA
  }else if(H>=thresh_H_divDup){
    paralog<-2
  }else if(H>=thresh_H|D>thresh_D|D<thresh_Dneg){
    paralog<-1
  }else{
    paralog<-0
  }
  return(paralog)
}
HDPlot.wild$paralog<-apply(HDPlot.wild,1,paralogStatus,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup)
#replot results color coded by paralog status
#traditional HDplot
ggplot()+geom_point(data=HDPlot.wild,aes(x=het_perc,y=z,color=paralog),alpha=0.2)+xlab("H")+ylab("D")+theme_bw()
#HDplot with ratio instead of deviation
ggplot()+geom_point(data=HDPlot.wild,aes(x=het_perc,y=ratio,color=paralog),alpha=0.2)+xlab("H")+ylab("Ratio")+theme_bw()

###Extract info and export it
HDPlot.wild <- na.omit(HDPlot.wild)
paral.wild <- na.omit(HDPlot.wild[HDPlot.wild$paralog == 1,])
paral.div.wild <- na.omit(HDPlot.wild[HDPlot.wild$paralog == 2,])

write.table(paral.wild, file = "../data/paral_wild.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(paral.div.wild, file = "../data/paral_div_wild.txt", sep = "\t", quote = F, row.names = F, col.names = T)


##################################### CULTIVATED SAMPLES ###########################
#load VCF file
vcf.cult<-read.vcfR("../data/Prefiltering/cult_site_filt.recode.vcf")
#Run HDplot and store results
HDPlot.cult<-HDPlot_fast(vcf.cult)

#plot results of HDplot without any thresholds 
#traditional HDplot
ggplot()+geom_point(data=HDPlot.cult,aes(x=het_perc,y=z),alpha=0.2)+xlab("H")+ylab("D")+theme_bw()
#HDplot with ratio instead of deviation
ggplot()+geom_point(data=HDPlot.cult,aes(x=het_perc,y=ratio),alpha=0.2)+xlab("H")+ylab("Ratio")+theme_bw()

#Thresholds at this point are determined by visual examination of the plot.
thresh_H<-0.62
thresh_H_divDup<-0.93
thresh_D<-7
thresh_Dneg<--7

#Add paralog status to HDplot results table based on thresholds above
paralogStatus<-function(data,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup){
  #if(data$het_perc>thresh_H|data$z>thresh_D|data$z<thresh_Dneg){
  H<-as.numeric(data[7])
  D<-as.numeric(data[9])
  if(is.na(D)){
    paralog<-NA
  }else if(H>=thresh_H_divDup){
    paralog<-2
  }else if(H>=thresh_H|D>thresh_D|D<thresh_Dneg){
    paralog<-1
  }else{
    paralog<-0
  }
  return(paralog)
}
HDPlot.cult$paralog<-apply(HDPlot.cult,1,paralogStatus,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup)
#replot results color coded by paralog status
#traditional HDplot
ggplot()+geom_point(data=HDPlot.cult,aes(x=het_perc,y=z,color=paralog),alpha=0.2)+xlab("H")+ylab("D")+theme_bw()
#HDplot with ratio instead of deviation
ggplot()+geom_point(data=HDPlot.cult,aes(x=het_perc,y=ratio,color=paralog),alpha=0.2)+xlab("H")+ylab("Ratio")+theme_bw()

###Extract info and export it
HDPlot.cult <- na.omit(HDPlot.cult)
paral.cult <- na.omit(HDPlot.cult[HDPlot.cult$paralog == 1,])
paral.div.cult <- na.omit(HDPlot.cult[HDPlot.cult$paralog == 2,])

write.table(paral.cult, file = "../data/paral_cult.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(paral.div.cult, file = "../data/paral_div_cult.txt", sep = "\t", quote = F, row.names = F, col.names = T)

