rm(list=ls(all=T))
library(ggplot2)
library(dplyr)
#setwd("C:/Users/borne/Documents/PostDoc/PostDoc_Project/trio_epistasis/clustering_softwares/MASS-PRF/Andrew_CDS_results/batch012/tables_only")
#setwd("C:/Users/borne/Documents/PostDoc/PostDoc_Project/trio_epistasis/clustering_softwares/MASS-PRF/Andrew_CDS_results/")
setwd("C:/Users/borne/Documents/PostDoc/PostDoc_Project/trio_epistasis/clustering_softwares/MASS-PRF/Andrew_CDS_results/Dmel_Dsim_Anc_v2/batch2/tables_only/several_ge4_some_lt20")

# List all files matching your pattern
file_list <- list.files(pattern = "*\\.txt$")
#Or use a txt file with files paths
#file_list <- readLines("adaptive_only_batch1.txt")

# Open a single PDF to save all plots
pdf("combined_plots_candidates_DmelDsimAnc_batch2.pdf", width = 8, height = 4)

# Loop through each file and plot
for (file_name in file_list) {
  
  # Read the data
  data <- read.table(file_name, header = TRUE)
  
  # Convert columns to numeric
  data$Gamma <- as.numeric(data$Gamma)
  data$Position <- as.numeric(data$Position)
  data$Upper_CI_Gamma <- as.numeric(data$Upper_CI_Gamma)
  data$Lower_CI_Gamma <- as.numeric(data$Lower_CI_Gamma)
  
  # Sort and group
  data <- data %>%
    arrange(Position) %>%
    mutate(
      color_group = factor(Lower_CI_Gamma >= 4, levels = c(FALSE, TRUE)),  # TRUE = highlight
      seg_id = cumsum(lag(color_group, default = first(color_group)) != color_group)
    )
  
  # Generate the plot
  p <- ggplot(data, aes(x = Position, y = Gamma)) +
    geom_errorbar(aes(ymax = Upper_CI_Gamma, ymin = Lower_CI_Gamma), color = "grey", alpha = 0.5) +
    geom_line(size = 0.8) +
    theme_minimal() +
    xlab("Codon position") +
    ylab("Selection intensity (γ)") +
    ggtitle(file_name) +
    geom_line(aes(group = seg_id, colour = color_group), size = 0.8) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    guides(colour = "none")
  
  # Print the plot to the PDF
  print(p)
}

# Close the PDF device
dev.off()
