library(ggplot2)
#Load meta data
meta <- read.delim("../meta/cocci_meta.txt")

#Load data. Fragments of at least 300Kb
runs.300 <- read.delim("../data/LRHomoz/cocci_300kb.hom.indiv")
#Add meta info
runs.300$pop <- meta$Population
runs.300$Status <- meta$Estatus2
runs.300 <- droplevels(runs.300[runs.300$pop!="Excluded",])

#Violin plot of the number of segments
ggplot(runs.300, aes(x=pop, y=NSEG, fill=pop)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE) +
  scale_fill_manual(values=c("#ffe082", "#fc8014", "#c9bc1f", "#ffb39c", "#e50000","#c79a00", "#b76900",
                             "#80deea", "#2196f3", "#ba68c8", "#7349bd", "#9da1db", "#a2cf6e", "#009700", "#30709e")) +
  ylab("Number of fragments") + ggtitle("300 Kb")+ 
  theme_bw(base_size = 17)+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Violin plot of the total lenght of ROH
ggplot(runs.300, aes(x=pop, y=KB, fill=pop)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE) +
  scale_fill_manual(values=c("#ffe082", "#fc8014", "#ffb39c", "#e50000","#c79a00", "#b76900",
                             "#757575",
                             "#80deea", "#2196f3", "#ba68c8", "#7349bd", "#9da1db", "#a2cf6e", "#009700", "#30709e")) +
  ylab("Total lenght of ROH (Kb)") + ggtitle("300 Kb")+ 
  theme_bw(base_size = 17)+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 150000))

##############################################################
#Load data. Fragments of at least 500Kb
runs.500 <- read.delim("../data/LRHomoz/cocci_500kb.hom.indiv")
#Add meta info
runs.500$pop <- meta$Population
runs.500$Status <- meta$Estatus2
runs.500 <- droplevels(runs.500[runs.500$pop!="Excluded",])

#Violin plot of the number of segments
ggplot(runs.500, aes(x=pop, y=NSEG, fill=Status)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE) +
  ylab("Number of fragments") + ggtitle("500 Kb")+ 
  theme_bw(base_size = 17)+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Violin plot of the total lenght of ROH
ggplot(runs.500, aes(x=pop, y=KB, fill=pop)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE)+
  scale_fill_manual(values=c("#ffe082", "#fc8014", "#ffb39c", "#e50000", "#b76900","#c79a00",
                             "#5d4037",
                             "#80deea", "#2196f3", "#ba68c8", "#7349bd", "#9da1db", "#a2cf6e", "#009700", "#30709e")) +
  ylab("Total lenght of ROH (Kb)") + ggtitle("500 Kb")+ 
  theme_bw(base_size = 25)+ theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + ylim(c(0, 150000))

##############################################################
#Load data. Fragments of at least 1000Kb
runs.1000 <- read.delim("../data/LRHomoz/cocci_1000kb.hom.indiv")
#Add meta info
runs.1000$pop <- meta$Population
runs.1000$Status <- meta$Estatus2
runs.1000 <- droplevels(runs.1000[runs.1000$pop!="Excluded",])

#Violin plot of the number of segments
ggplot(runs.1000, aes(x=pop, y=NSEG, fill=Status)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE) +
  ylab("Number of fragments") + ggtitle("1000 Kb")+ 
  theme_bw(base_size = 17)+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Violin plot of the total lenght of ROH
ggplot(runs.1000, aes(x=pop, y=KB, fill=pop)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE)+
  scale_fill_manual(values=c("#ffe082", "#fc8014", "#ffb39c", "#e50000","#c79a00", "#b76900",
                             "#757575",
                             "#80deea", "#2196f3", "#ba68c8", "#7349bd", "#9da1db", "#a2cf6e", "#009700", "#30709e")) +
  ylab("Total lenght of ROH (Kb)") + ggtitle("1000 Kb")+ 
  theme_bw(base_size = 17)+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 150000))

