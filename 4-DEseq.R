################################################################################
# Loading dependencies
################################################################################

# Installing DESeq2

install.packages("BiocManager")

if (first == "install") {
  source("http://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  stop("Installation completed", call.=FALSE)
}

BiocManager::install("apeglm")

# Loading DESeq2

library("DESeq2")
library("apeglm")

################################################################################
# Processing featureCount output and setting up environment
################################################################################

# The output from featureCount is given as tabe separated files for each sample.
# Files need to be converted to a new format that merges all STM-related files
# and another for TCP4-related files



# Point to the working dir

setwd("C:/Users/c1831460/OneDrive - Cardiff University/Documents/DTP/Second Year/RNA-seq/featureCounts/")

# Create dataframe contianing output from 3-featureCounts.sh

STM <- list.files(
  pattern = "STM.*.csv")

TCP4 <- list.files(
  pattern = "TCP4.*.csv")

# Create tables to tell DESeq2 what the variables are
# The number of samplesheets (coldata) you need will depend on the number of
# timepoints/variables you want to compare

#STM
STM_coldata <- data.frame(
                      sample = c("STM_C3_S7",
                                 "STM_C4_S8",
                                 "STM_C5_S9",
                                 "STM_CD3_S10",
                                 "STM_CD4_S11",
                                 "STM_CD5_S12",
                                 "STM_D3_S4",
                                 "STM_D4_S5",
                                 "STM_D5_S6",
                                 "STM_M3_S1",
                                 "STM_M4_S2",
                                 "STM_M5_S3"
                                  ),
                      treatment = factor(c("C", "C", "C",
                                    "CD","CD","CD",
                                    "D","D","D",
                                    "M","M", "M"
                                    )),
                      replicate = c("1", "2", "3")
)

#TCP4
TCP4_coldata <- data.frame(
                      sample = c("TCP_C2_S19",
                                 "TCP_C4_S20",
                                 "TCP_C5_S21",
                                 "TCP_CO2_S22",
                                 "TCP_CO4_S23",
                                 "TCP_CO5_S24",
                                 "TCP_D2_S16",
                                 "TCP_D4_S17",
                                 "TCP_D5_S18",
                                 "TCP_M2_S13",
                                 "TCP_M4_S14",
                                 "TCP_M5_S15"
                                  ),
                      treatment = factor(c("C", "C", "C",
                                    "CD","CD","CD",
                                    "D","D","D",
                                    "M","M", "M"
                                   )),
                      replicate = c("1", "2", "3")
)

gene = list(STM = STM, TCP4 = TCP4)

################################################################################
# Main CMDs
################################################################################

# Create loop

colDatafield <- list(STM_coldata, TCP4_coldata)
namefield <- c("markdup", "rmdup")
genotypefield <- c("STM", "TCP4")
index1 = as.numeric(1)
index2 = as.numeric(1)


for (i in gene) {
  
  print(i)
  
  index1 = as.numeric(1)
  
  for(j in i) {
  
    print(j)
    # Read the data from standard input
  
    data <- as.matrix(read.csv(
      file = paste (j,
                    sep = ""),
      header = TRUE,
      row.names = "Geneid",
      sep = "\t"))
      
  # Use dds to combine the coldata and countdata matrix
  # State the variable/factor being analysed using the "design" flag
    coldata <- colDatafield[[index2]]
    dds <- DESeqDataSetFromMatrix(
      countData = data, 
      colData = coldata,
      design = ~ treatment)
    
    dds
    
  # As R will automatically choose the reference level for the variable/factor,
  # the control group needs to be defined and releveled
  
    dds$treatment <- relevel(
      dds$treatment, 
      ref = "M")
  
  # Run the DESeq command on the DESeq dataset
  
    dds <- DESeq(dds)
  
  # Look at the contrasts we can build
    resultsNames(dds)
    
# Generate a results table, and specify the contrast we want to build
  
    res <- results(
      dds,
      contrast=c(
        "treatment",
        "M",
        "D"))
  
  #OR...
    
#    res <- results(
#      dds,
#      name = "treatment_M_vs_D"
#    )

  # Log fold change shrinkage for visualisation can be done using different models
  
    resLFC <- lfcShrink(
      dds,
      coef= "treatment_M_vs_D",
      type= "apeglm")
  
    resLFC
  
  #OR...
  
#    resNorm <- lfcShrink(
#      dds, 
#      coef = "genotype__vs_", 
#      type = "normal")
    
#    resNorm
  
  # Plotting the normalised results with MA
  # Set the probability and log fold change thresholds
  
    xlim <- c(1,1e5); ylim <- c(-3,3)
    plotMA(
      resLFC, 
      xlim=xlim, 
      ylim=ylim, 
      main="apeglm")
  
  # The plot can be used to identify the rownumber of individual genes
  # interactively
  
  #  idx <- identify(
  #    res$baseMean, 
  #    res$log2FoldChange)
    
  #  rownames(res)[idx]
  
  # Extract all genes (independent of differential expression)
    print(file.path("DEGs", genotypefield[index2],
      paste( 
      genotypefield[index2], "_M_vs_D_", namefield[index1], "_DEGs.csv", 
      sep = "")))
    
    write.csv(resLFC, (
      file = file.path("DEGs", genotypefield[index2],
        paste( 
        genotypefield[index2], "_M_vs_D_", namefield[index1], "_DEGs.csv", 
        sep = ""))))
  
  # Extract significant upregulated and down regulated genes into separate
  # datasets
  
    upreg <- subset(
      resLFC, 
      log2FoldChange >=0.5 & padj <0.05)
    
    upreg <- upreg[ , c(-1,-3)]
  
    downreg <-subset(
      resLFC, 
      log2FoldChange <=-0.5 & padj <0.05)
    
    downreg <- downreg[ , c(-1,-3)]
  
  #write csv file of for up and down regualted gene IDs
  
    write.csv(upreg, (
      file = file.path("DEGs", genotypefield[index2],
                       paste( 
                         genotypefield[index2], "_M_vs_D_", namefield[index1], 
                         "_upreg.csv", 
                         sep = "")))
      )
      
    write.csv(upreg, (
      file = file.path("DEGs", genotypefield[index2],
                       paste( 
                         genotypefield[index2], "_M_vs_D_", namefield[index1], 
                         "_downreg.csv", 
                         sep = "")))
      )
    
    index1 <- index1+1
  
  }
  
  index2 <- index2+1
}
