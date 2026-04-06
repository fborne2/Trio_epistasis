rm(list=ls(all=T))
setwd("C:/Users/borne/Documents/PostDoc/PostDoc_Project/trio_epistasis/trio_gene/")
library(ggplot2)
library(gridExtra)
#*********************************************************************
#Upload data

data <- read.table("./data/trio_fertility_summer2025.txt", sep= "\t", header=T)
data <- subset(data, data$haplotype %in% c("ASK","VNR"))



n_fun <- function(data){
  return(data.frame(y = 25, label = paste0("n = ",length(data))))
}
#25 and 29 homo
p1 <- ggplot(data, aes(x=haplotype, y=number))
p1 <- p1 + geom_point(position=position_jitter(h=0, w=0.02), size=2, shape=1)
p1 <- p1 + stat_summary(fun = "mean", alpha=0.7, position= position_nudge(x=0.1), color="red")
p1 <- p1 + stat_summary(fun.data = mean_se, geom = "errorbar", position= position_nudge(x=0.1), size=1, width=0, color="red")
p1 <- p1 + ggtitle(label="Number of progeny per fertile cross") + ylab("Number of progeny") + xlab(NULL)
p1 <- p1 + theme_classic()
p1 <- p1 + stat_summary(fun.data = n_fun, geom = "text", size=3.5) 
p1 <- p1 + theme(panel.background = element_rect(fill = "white", colour = "grey50"))
#p1 <- p1 + facet_wrap(~count_progeny, ncol=2)
p1


t.test(data$number[data$haplotype=='ASK'], data$number[data$haplotype=='VNR'])


##############################################""""
##Check Power 

#install.packages("pwr")
library(pwr)
#Summary statistics
library(Rmisc)
datastat <- summarySE(data, measurevar = "number", groupvars = c("haplotype"))

#
n1 <- datastat$N[datastat$haplotype=="ASK"]
mean1 <- datastat$number[datastat$haplotype=="ASK"]
sd1 <- datastat$sd[datastat$haplotype=="ASK"]
  
n2 <- datastat$N[datastat$haplotype=="VNR"]
mean2 <- datastat$number[datastat$haplotype=="VNR"]
sd2 <- datastat$sd[datastat$haplotype=="VNR"]

#Observed difference between the two haplotypes
delta <- mean2 - mean1
delta
#In terms of percent differences
delta_per = (mean2-mean1)/mean1*100
delta_per
#VNR is 4.8% increase

#Compute pooled standard deviation
sp <- sqrt(
  ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) /
    (n1 + n2 - 2)
)

sp

#Compute Cohen's d for the data
d <- delta / sp
d
#d=0.09 corresponds to about 5% difference increase in mean


############################################################
# Required flies vs % difference (80% power)
############################################################

# significance level and power
alpha <- 0.05
power_target <- 0.8

# Range of percent differences to evaluate
percent_seq <- seq(1, 50, by = 1)

# Function: compute required n for each % difference
n_required <- sapply(percent_seq, function(p) {
  
  delta <- mean1 * (p / 100)   # convert % to absolute difference
  d <- delta / sp              # Cohen's d
  
  pwr.t.test(
    d = d,
    power = power_target,
    sig.level = alpha,
    type = "two.sample"
  )$n
})

############################################################
# Plot
############################################################

df <- data.frame(
  percent = percent_seq,
  n_required = n_required
)

p2<- ggplot(df, aes(x = percent, y = n_required)) +
  geom_line(linewidth = 1.2) +
  
  # axis limits
  coord_cartesian(ylim = c(0, 2000)) +
  
  # vertical lines
  geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "darkgreen") +
  
  # horizontal line for your current n
  geom_hline(yintercept = 30, linetype = "dashed", color = "black") +
  
  # labels for vertical lines
  annotate("text", x = 5, y = 1900, label = "5%", color = "red", angle = 90, vjust = -0.5) +
  annotate("text", x = 10, y = 1900, label = "10%", color = "blue", angle = 90, vjust = -0.5) +
  annotate("text", x = 30, y = 1900, label = "30%", color = "darkgreen", angle = 90, vjust = -0.5) +
  
  
  # label for horizontal line
  annotate("text", x = 2, y = 30, label = "n = 30", color = "black",vjust = -0.5) +
  
  labs(
    title = "Number of flies required to detect difference in\n mean progeny at 80% power, 0.05 significance",
    x = "% change in mean progeny",
    y = "Required flies per group"
  ) +
  
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))
p2

grid.arrange(p1,p2,ncol=2)


