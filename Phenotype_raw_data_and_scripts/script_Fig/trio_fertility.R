rm(list=ls(all=T))

library(ggplot2)
library(gridExtra)
#*********************************************************************
#Upload data

data <- read.table("trio_fertility_summer2025.txt", sep= "\t", header=T)



n_fun <- function(data){
  return(data.frame(y = 25, label = paste0("n = ",length(data))))
}
#25 and 29 homo
p1 <- ggplot(data, aes(x=haplotype, y=number))
p1 <- p1 + geom_point(position=position_jitter(h=0, w=0.02), size=2, shape=1)
p1 <- p1 + stat_summary(fun = "mean", alpha=0.7, position= position_nudge(x=0.1), color="red")
p1 <- p1 + stat_summary(fun.data = mean_se, geom = "errorbar", position= position_nudge(x=0.1), size=1, width=0, color="red")
p1 <- p1 + ggtitle(label="fertility") + ylab("total number of progeny") + xlab(NULL)
p1 <- p1 + theme_classic()
p1 <- p1 + stat_summary(fun.data = n_fun, geom = "text", size=3.5) 
p1 <- p1 + theme(panel.background = element_rect(fill = "white", colour = "grey50"))
#p1 <- p1 + facet_wrap(~count_progeny, ncol=2)
p1



wilcox.test(data$number[data$haplotype=='ASK'], data$number[data$haplotype=='VNR'])


