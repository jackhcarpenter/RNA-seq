################################################################################
# Loading dependencies
################################################################################

# Installing DESeq2

install.packages("BiocManager")
BiocManager::install("apeglm")
BiocManager::install("DESeq2")

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

setwd("C:/Users/c1831460/OneDrive - Cardiff University/Documents/DTP/Second Year/RNA-seq TCP4_STM/")

# Create dataframes containing output from 3.2-featureCounts_merged.sh

t3h <- list.files(
  pattern = "3h.*.csv")

t12h <- list.files(
  pattern = "12h.*.csv")

t24h <- list.files(
  pattern = "24h.*.csv")

t11d <- list.files(
  pattern = "11d.*.csv")

# Create tables to tell DESeq2 what the variables are
# The number of samplesheets (coldata) you need will depend on the number of
# timepoints/variables you want to compare

#3h
coldata_3h <- data.frame(
  sample = factor(c("Col-0_3h_1",
                    "Col-0_3h_2",
                    "Col-0_3h_3",
                    "STM_3h_1",
                    "STM_3h_2",
                    "STM_3h_3",
                    "TCP_3h_1",
                    "TCP_3h_2",
                    "TCP_3h_3"
  )),
  treatment = factor(c("Col", "Col", "Col",
                       "STM","STM","STM",
                       "TCP","TCP","TCP"
  )),
  replicate = factor(c("1", "2", "3",
                       "1", "2", "3",
                       "1", "2", "3")
  ))


#12h
coldata_12h <- data.frame(
  sample = factor(c("Col-0_12h_1",
                    "Col-0_12h_2",
                    "Col-0_12h_4",
                    "STM_12h_1",
                    "STM_12h_2",
                    "STM_12h_3",
                    "TCP_12h_1",
                    "TCP_12h_2",
                    "TCP_12h_4"
  )),
  treatment = factor(c("Col", "Col", "Col",
                       "STM","STM","STM",
                       "TCP","TCP","TCP"
  )),
  replicate = factor(c("1", "2", "3",
                       "1", "2", "3",
                       "1", "2", "3")
  ))


#24h
coldata_24h <- data.frame(
  sample = factor(c("Col-0_24h_1",
                    "Col-0_24h_3",
                    "Col-0_24h_4",
                    "STM_24h_2",
                    "STM_24h_3",
                    "STM_24h_5",
                    "TCP_24h_1",
                    "TCP_24h_3",
                    "TCP_24h_4"
  )),
  treatment = factor(c("Col", "Col", "Col",
                       "STM","STM","STM",
                       "TCP","TCP","TCP"
  )),
  replicate = factor(c("1", "2", "3",
                       "1", "2", "3",
                       "1", "2", "3")
  ))

#11d
coldata_11d <- data.frame(
  sample = factor(c("Col-0_11d_1",
                    "Col-0_11d_2",
                    "Col-0_11d_3",
                    "STM_11d_1",
                    "STM_11d_2",
                    "STM_11d_3",
                    "TCP_11d_1",
                    "TCP_11d_2",
                    "TCP4_11d_5"
  )),
  treatment = factor(c("Col", "Col", "Col",
                       "STM","STM","STM",
                       "TCP","TCP","TCP"
  )),
  replicate = factor(c("1", "2", "3",
                       "1", "2", "3",
                       "1", "2", "3")
  ))

################################################################################
# Main CMDs
################################################################################

# Create loop and var for loop

# List of timepoints to initiate the for loop
gene = list(t3h = t3h, t12h = t12h, t24h = t24h, t11d = t11d)

# For passing different timepoint data into coldata 
colDatafield <- list(coldata_3h, coldata_12h, coldata_24h, coldata_11d)

# Markduplicates and removeduplicates name list
namefield <- c("markdup", "rmdup")

# List of timepoint names as factor
timepointfield <- c("3h", "12h", "24h", "11d")

# List of genotypes as a factor
genotypefield <- c("STM","TCP")

# To switch between namefield variables
index1 = as.numeric(1)

# To switch between timepoint fields
index2 = as.numeric(1)

# To switch between genotype fields
index3 = as.numeric(1)


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
      row.names = 1,
      sep = ","))
    
    # Edit the col names of "data" to match the rownames of "coldata"
    
    colnames(data) <- sub("X.mnt.scratch.c1831460.RNA.seq.", "", colnames(data))
    colnames(data) <- sub(namefield[index1], "", colnames(data))
    colnames(data) <- sub(".merged.", "", colnames(data))
    colnames(data) <- sub(".bam", "", colnames(data))
      
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
      ref = "Col")
  
  # Run the DESeq command on the DESeq dataset
  
    dds <- DESeq(dds)
  
  # Look at the contrasts we can build
    resultsNames(dds)
    
  # Generate a results table, and specify the contrast we want to build
    
    index3 = as.numeric(1)
    
    for(x in genotypefield) {
  
      res <- results(
        dds,
        contrast=c(
          "treatment",
          genotypefield[index3],
          "Col"))
    
    #OR...
      
    #    res <- results(
    #      dds,
    #      name = "treatment_M_vs_D"
    #    )
  
    # Log fold change shrinkage for visualisation can be done using different models
      
      coef <- (paste0("treatment_", genotypefield[index3], "_vs_Col"))
    
      resLFC <- lfcShrink(
        dds,
        coef= coef,
        type= "apeglm")
    
      resLFC
    
    #OR...
    
    #    resNorm <- lfcShrink(
    #      dds, 
    #      coef = "timepoint__vs_", 
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
      print(file.path("DEGs", timepointfield[index2],
        paste( 
        timepointfield[index2], "_Col_vs_", genotypefield[index3], "_", namefield[index1], "_DEGs.csv", 
        sep = "")))
      
      write.csv(resLFC, (
        file = file.path("DEGs", timepointfield[index2],
          paste( 
          timepointfield[index2], "_Col_vs_", genotypefield[index3], "_", namefield[index1], "_DEGs.csv", 
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
        file = file.path("DEGs", timepointfield[index2],
                         paste( 
                           timepointfield[index2], "_Col_vs_", genotypefield[index3], "_", namefield[index1], 
                           "_upreg.csv", 
                           sep = "")))
        )
        
      write.csv(downreg, (
        file = file.path("DEGs", timepointfield[index2],
                         paste( 
                           timepointfield[index2], "_Col_vs_", genotypefield[index3], "_", namefield[index1], 
                           "_downreg.csv", 
                           sep = "")))
        )
      
      index3 <- index3+1
    
    }
    
    index1 <- index1+1
  
  }
  
  index2 <- index2+1
}
