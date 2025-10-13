rm(list=ls(all=T))
library(ggplot2)
library(gridExtra)
library(tidyverse)

# Install FSA package if not already installed
#install.packages("FSA")
# Load the FSA package
library(FSA)
#*********************************************************************
#Upload data

data_homo <- read.table("ratio_WT_homo.txt", sep= "\t", header=T)
data_hetero <- read.table("ratio_WT_hetero.txt", sep= "\t", header=T)

#opposite direction
data_homo$haplotype <- factor(data_homo$haplotype, levels=c('VNR','ANR','VNK','VSR','ASR','VSK','ANK','ASK'), order=T)
data_hetero$haplotype <- factor(data_hetero$haplotype, levels=c('VNR','ANR','VNK','VSR','ASR','VSK','ANK','ASK'), order=T)

# Reshape to long format
long_data_homo <- data_homo %>%
  pivot_longer(cols = larva:adult, names_to = "stage", values_to = "value") %>%
  drop_na()
long_data_hetero <- data_hetero %>%
  pivot_longer(cols = larva:adult, names_to = "stage", values_to = "value") %>%
  drop_na()

long_data_homo$stage <- factor(long_data_homo$stage, levels=c("larva","pupa","adult"), order=T)
long_data_hetero$stage <- factor(long_data_hetero$stage, levels=c("larva","pupa","adult"), order=T)


p_homo <- ggplot(long_data_homo, aes(x = haplotype, y = value, fill = stage)) +
  #geom_hline(yintercept=1, linetype="dashed", size=0.8, color="darkgrey") +
  geom_bar(stat = "identity", color = "black", position = position_dodge(width = 0.8),width = 0.7) + 
  labs(x = "Haplotype", y = "ratio of observed vs expected") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),  # Keep horizontal lines
    panel.grid.minor.y = element_blank()
  ) + scale_fill_manual(values = c("whitesmoke","darkgrey","black"))
p_homo <- p_homo + ylim(0,1.8) + ggtitle("homozygous") +     theme(legend.position = "top")
p_homo


p_hetero <- ggplot(long_data_hetero, aes(x = haplotype, y = value, fill = stage)) +
  #geom_hline(yintercept=1, linetype="dashed", size=0.8, color="darkgrey") +
  geom_bar(stat = "identity",color = "black", position = position_dodge(width = 0.8),width = 0.7) + 
  labs(x = "Haplotype", y = "ratio of observed vs expected") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),  # Keep horizontal lines
    panel.grid.minor.y = element_blank()
  ) + scale_fill_manual(values = c("whitesmoke","darkgrey","black"))
p_hetero <- p_hetero + ylim(0,1.8) + ggtitle("heterozygous") +     theme(legend.position = "top")
p_hetero


grid.arrange(p_homo,p_hetero, ncol = 2)

