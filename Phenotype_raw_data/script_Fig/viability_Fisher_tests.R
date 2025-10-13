rm(list=ls(all=T))

larva_homo <- read.table("larva_homo_Fisher_test.txt", sep= "\t", header=T)
#Can use the same script for: 
#larva_hetero <- read.table("larva_hetero_Fisher_test.txt", sep= "\t", header=T)
#pupa_homo <- read.table("pupa_homo_Fisher_test.txt", sep= "\t", header=T)
#pupa_hetero <- read.table("pupa_hetero_Fisher_test.txt", sep= "\t", header=T)
#adult_homo <- read.table("adult_homo_Fisher_test.txt", sep= "\t", header=T)
#adult_hetero <- read.table("adult_hetero_Fisher_test.txt", sep= "\t", header=T)


##Run p-value for each developmental stage
data <- larva_homo
data <- na.omit(data)


# Round up to integers
data$EYFP <- ceiling(data$EYFP)
data$WT <- ceiling(data$WT)

# Store p-values
pvals <- c()
comparisons <- c()

# Get ASK values
ask_eyfp <- data$EYFP[data$haplotype == "ASK"]
ask_wt <- data$WT[data$haplotype == "ASK"]

# Loop through other haplotypes
for (i in 1:nrow(data)) {
  hap <- data$haplotype[i]
  if (hap != "ASK") {
    hap_eyfp <- data$EYFP[i]
    hap_wt <- data$WT[i]
    
    # Create 2x2 table
    mat <- matrix(c(hap_eyfp, hap_wt, ask_eyfp, ask_wt),
                  nrow = 2,
                  byrow = TRUE,
                  dimnames = list(c(hap, "ASK"), c("EYFP", "WT")))
    
    # Fisher's exact test
    test <- fisher.test(mat)
    pvals <- c(pvals, test$p.value)
    comparisons <- c(comparisons, hap)
  }
}

# Adjusted p-values
pvals_adj <- p.adjust(pvals, method = "fdr")  # Or "bonferroni"

# Combine into a result data frame
results <- data.frame(
  Comparison = comparisons,
  Raw_P = pvals,
  Adjusted_P = pvals_adj
)

print(results)

