rm(list=ls(all=T))
library(ggplot2)
library(gridExtra)
library(FSA)
#*********************************************************************
#Upload data

data <- read.table("./data/climbing_assay.txt", sep= "\t", header=T)
data$line <- as.factor(data$line)
class(data$line)

data$haplotype <- factor(data$haplotype, levels=c('VNR/VNR', 'ANR/ANR','VSK/VSK', 'ASK/ASK','ANR/VNR',
                                                  'VNK/VNR','VSR/VNR','ASR/VNR','VSK/VNR','ANK/VNR','ASK/VNR'), order=T)


data <- subset(data, data$age %in% c("7-9","7","8","9","7-10","10","11"))
data <- subset(data, data$number >4)

data$normalized <- data$mean_5s/mean(data$mean_5s[data$haplotype=='ASK/ASK'])

n_fun <- function(data){
  return(data.frame(y = 0, label = paste0("n = ",length(data))))
}
#Plot locomotion normalized by ASK/ASK
p1 <- ggplot(data, aes(x=haplotype, y=normalized))
p1 <- p1 + geom_point(position=position_jitter(h=0, w=0.03), shape=1)
p1 <- p1 + stat_summary(fun = "mean", col="red", alpha=0.7, position= position_nudge(x=0.15))
p1 <- p1 + stat_summary(fun.data = mean_se, geom = "errorbar", col="red", position= position_nudge(x=0.15), size=1, width=0)
p1 <- p1 + ggtitle(label="climbing assay 10s") + xlab("percent reaching the line") + ylab(NULL)
p1 <- p1 + scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1,1.2), limits=c(0,1.2)) + theme_classic()
p1 <- p1 + stat_summary(fun.data = n_fun, geom = "text", size=3.5) 
p1 <- p1 +
  theme(
    panel.background = element_blank(),        # Remove gray rectangle
    panel.grid.major.x = element_blank(),      # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),  # Keep horizontal grid lines
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black")  # Keep x and y axes
  )

p1


hetero <- p1
homo <- p1
hetero

grid.arrange(homo, hetero, nrow=1)

#Figure article
data$haplotype <- as.factor(data$haplotype)
##perform Kruskal-Wallis Test 
kruskal.test(mean_5s ~ haplotype, data = data) 
# Perform Dunn's test
library(dunn.test)
dunn.test(x = data$mean_5s, g = data$haplotype, method = "bh")
