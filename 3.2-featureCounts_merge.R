################################################################################
# Merging featureCount output files 
################################################################################

# Want to add all of the .featureCount.txt files together into one file

################################################################################
# Loading libraries and setting up environemnt
################################################################################

# Needed for dataframe manipulation
install.packages("dplyr")
library("dplyr")

# Workingdir
setwd("C:/Users/c1831460/OneDrive - Cardiff University/Documents/DTP/Second Year/RNA-seq/featureCounts/rawFC/")

# Making a list object containing files, separated by genotype
STMfiles <- list.files(pattern = "STM")

TCP4files <- list.files(pattern = "TCP")

# Making a list of file lists
genotypes <- list(STMfiles, TCP4files)

################################################################################
# Main CMDs
################################################################################

# Make a matrix that has the GeneID column so the columns from other
# featureCounts files can be added to it.
# Because the order of GeneIDs in feature counts is the same, any file can be
# used (can run example below).

#df1 <- read.table("STM_C3_S7.markdup.featurecounts.txt")
#df2 <-read.table("TCP_C2_S19.markdup.featurecounts.txt")
#print(df1$V1 %in% df2$V7)rm(df1, df2)

GeneID_col <- read.table("STM_C3_S7.markdup.featurecounts.txt", 
                         header = TRUE, 
                         sep = "\t")
GeneID_col <- pull(GeneID_col, 
                   var = +1, 
                   name = NULL)
GeneID_col <- as.matrix(GeneID_col)

# Initialise separate matricies for markdup and rmdup to be used in the if else loop
markdupGeneID <- GeneID_col
rmdupGeneID <- GeneID_col

# Character vector holding the genotype names and a counter to name files
Genotype_names <- c("STM", "TCP4")
Genotype_counter <-as.numeric(1)

# Outer for loop to switch between genotypes
for (x in genotypes) {
  
  # Inner for loop to iterate through file in the file list
  for (file in x) {
    
    # Inner if loop to delineate markdup and rmdup files in list with grep
    if (grepl("markdup", file)) {
      
      # Pass markdup file into markdup_merging_file and add the 7th coloum
      # representing counts into the final merged file markdupGeneID
      markdup_merging_file <- read.table(file, 
                                         header = TRUE, 
                                         sep = "\t")
      names(markdup_merging_file)[7] <- "counts"
      markdupGeneID <- bind_cols(markdupGeneID, markdup_merging_file$counts)
      #print(head(markdupGeneID))
      
  } else {
    
    # Pass rmdup file into rmdup_merging_file and add the 7th coloum
    # representing counts into the final merged file rmdupGeneID
    rmdup_merging_file <- read.table(file, 
                                     header = TRUE, 
                                     sep = "\t")
    #rename this column becuase it varies between files
    names(rmdup_merging_file)[7] <- "counts"
    rmdupGeneID <- bind_cols(rmdupGeneID, rmdup_merging_file$counts)
    #print(head(rmdupGeneID))
        
  }
    

  }
  
  # Write the markdup file
  write.csv(markdupGeneID, file = paste0(
    Genotype_names[Genotype_counter],"_markdup.csv"), row.names = FALSE)
  
  # Write the rmdup file
  write.csv(rmdupGeneID, file = paste0(
    Genotype_names[Genotype_counter],"_rmdup.csv"), row.names = FALSE)
  
  # Reset GeneID matricies
  markdupGeneID <- GeneID_col
  rmdupGeneID <- GeneID_col
  
  Genotype_counter <- Genotype_counter+1
    
}
