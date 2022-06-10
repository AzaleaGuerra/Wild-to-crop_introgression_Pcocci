#########  Plot Admixture  ########
library(ggplot2)

###Load K error data
k.error<- read.delim("../data/admixture/coccineus_Kerror.txt", header = F, sep = ":")
rownames(k.error)<- c("k=1", "k=2", "k=3", "k=4", "k=5", "k=6", "k=7", "k=8", "k=9", "k=10", "k=11", "k=12", "k=13", "k=14", "k=15", "k=16", "k=17", "k=18", "k=19", "k=20", "k=21", "k=22")

#plot error
e.plot<- ggplot(data=k.error[c(1:15),], aes(x=1:15, y=V2)) + geom_point() + geom_line()
e.plot + xlab("k") + ylab("Error") + theme_bw(base_size = 23)

###Load meta data
meta<- read.delim("../meta/cocci_meta.txt")
pop<- meta$Poblacion

###load admixture data (K value with the lowest error)
admix.8<- read.delim("../admixture_linux-1.3.0/coccineus.8.Q", sep = " ", header = F)
admix.8 <- cbind(pop, admix.8) #Bind pop info
#write.table(admix.8, file = "../meta/cocci_admix_8.txt", row.names = F, quote = FALSE, sep = "\t")
#Order the samples according to genetic groups
admix.8.or <- read.delim("../meta/cocci_admix_8.txt")


##### PLOTS #####
mis.col <- palette(c("#ffe082", "#ba68c8", "#7d80af", "#7ec2f8", "#fdaa61", "#f28b00", "#00c99e", "#ff7070"))
##Plot admixture K=8
#par(mar=c(5.5, 4.1, 4.1, 11.5), xpd=TRUE)
barplot(t(as.matrix(admix.8.or[,4:11])), col= mis.col, las=2, cex.names = .5, names.arg = admix.8.or$Cluster.Admix,
      ylab="", border=NA, cex.axis = 2.7, axisnames = T)
