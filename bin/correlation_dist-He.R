library(ggplot2)
library(geosphere)
library(ggpubr)
library(stringr)

#Load data
meta <- read.delim("../meta/data-map-het.csv", sep = ",")

#Estimate centroid of Cult-TMVB locations
c.tmvb <- na.omit(meta[meta$Population == "Cult-TMVB", c(1,4,3,7:10)])
cult <- na.omit(meta[meta$Estatus2 == "Cultivar", c(1,4,3,7:12)])
cent.c.tmvb <- centroid(c.tmvb[,c(2,3)])
#Estimates ditances to centroid
dista <- distm(cent.c.tmvb, cult[,c(2,3)])
cult$dist.cent <- as.vector(dista/1000)

###Correlation between distance to centroid and Pi
ggscatter(cult, y = "He.Loc_subset",  x = "dist.cent", 
          #color = "Pop",
          #fill = "Status",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE,
          cor.method = "spearman",
          ylab = "He", xlab = "Distance",
          cor.coef.size = 6,
          font.tickslab = c(16,"plain","black"),
          font.x = c(20, "plain", "black"),
          font.y = c(20, "plain", "black"))
cor.test(cult$He.Loc_subset, cult$dist.cent,  method = c("spearman", "kendall"))
cor.test(cult$He.Loc, cult$dist.cent, method = "spearman")

#Colored by population
scatter_plot <- ggplot(cult, aes(y=He.Loc_subset, x=dist.cent))
scatter_plot + geom_point(aes(color=Population), size = 12, alpha=0.8) + labs(y = "He", x = "Distance (Km)") + 
  geom_smooth(method="lm", alpha =0.15) + theme_bw(base_size = 33) +
  annotate("text", x = 300, y = 0.036, label = "italic(r) == -0.68", parse = TRUE, size = 13) +
  annotate("text", x = 550, y = 0.0357, label = "italic(p) == 0.05", parse = TRUE, size = 13) +
  scale_color_manual(values=c("#ffe082", "#fc8014", "#ffb39c", "#e50000","#c79a00", "#b76900", "#757575"))

scatter_plot <- ggplot(cult, aes(y=He.Loc, x=dist.cent))
scatter_plot + geom_point(aes(color=Population), size = 12, alpha=0.8) + labs(y = "He", x = "Distance (Km)") + 
  geom_smooth(method="lm", alpha =0.15) + theme_bw(base_size = 33) +
  annotate("text", x = 300, y = 0.036, label = "italic(r) == -0.67", parse = TRUE, size = 13) +
  annotate("text", x = 550, y = 0.0357, label = "italic(p) == 0.05", parse = TRUE, size = 13) +
  scale_color_manual(values=c("#ffe082", "#fc8014", "#ffb39c", "#e50000","#c79a00", "#b76900", "#757575"))

#Colored by location
scatter_plot <- ggplot(cult, aes(y=He.Loc_subset, x=dist.cent))
scatter_plot + geom_point(aes(color=Locacion), size = 12, alpha=0.75) + labs(y = "He", x = "Distance (Km)") + 
  geom_smooth(method="lm", alpha =0.15) + theme_bw(base_size = 33) +
  annotate("text", x = 300, y = 0.036, label = "italic(r) == -0.68", parse = TRUE, size = 13) +
  annotate("text", x = 550, y = 0.0357, label = "italic(p) == 0.05", parse = TRUE, size = 13) 
  

