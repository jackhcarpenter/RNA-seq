################################################################################
# Setup and loading dependencies
################################################################################

# Installing DESeq2

if (first == "install") {
  source("http://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  stop("Installation completed", call.=FALSE)
}

# Loading DESeq2

library("DESeq2")

# Point to the working dir

swd("path / to / your / dir")

# Create dataframe contianing output from 3-featureCounts.sh

rmdup <- list.files(
  pattern = " * rmdup.featurecount")


# Create tables to tell DESeq2 what the variables are
# The number of samplesheets (coldata) you need will depend on the number of
# timepoints/
# variables you want to compare

coldata <- data.frame(
  sample = c("sample_1", "sample_2", "sample_n"),
  genotype = c("genotype_1", "genotype_2", "genotype_n"),
  replicate = c("replicate_1", "replicate_2", "replicate_n")
)

################################################################################
# Main CMDs
################################################################################

# Create loop

for(i in rmdup) {

# Read the data from standard input

  data <- as.matrix(read.csv(
    file = paste (i, sep = ""),
    header = TRUE,
    row.names = "Geneid",
    sep = "\t"))
    
# Use dds to combine the coldata and countdata matrix
# State the variable/factor being analysed using the "design" flag

  dds <- DESeqDataSetFromMatrix(
    countData = data, 
    colData = coldata,
    design = ~ genotype)
  
  dds
  
# As R will automatically choose the reference level for the variable/factor,
# the control group needs to be defined and releveled

  dds$genotype <- relevel(
    dds$genotype, 
    ref = "WT")

# Run the DESeq command on the DESeq dataset

  dds <- DESeq(dds)

# Generate a results table, and specify the contrast we want to build

  res <- results(
    dds,
    contrast=c(
      "genotype",
      "WT",
      "TCP4"))

#OR...

  res <- results(
    dds,
    name = "genotype_TCP4_vs_WT")

# Log fold change shrinkage for visualisation can be done using different models

  resLFC <- lfcShrink(
    dds, 
    coef= "genotype_TCP4_vs_WT",
    type="apeglm")

  resLFC

#OR...

  resNorm <- lfcShrink(
    dds, 
    coef = "genotype_TCP4_vs_WT", 
    type = "normal")
  
  resNorm

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

  idx <- identify(
    res$baseMean, 
    res$log2FoldChange)
  
  rownames(res)[idx]

# Extract all genes (independent of differential expression)
  
  write.csv(resLFC, (
    file = paste( 
      i, 
      "_DEGs.csv", 
      sep = "")))

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
    file = paste( 
      i, 
      "_upreg.csv", 
      sep = "")))
    
  write.csv(downreg, (
    file = paste( 
      i, 
      "_downreg.csv", 
      sep = "")))
  
}