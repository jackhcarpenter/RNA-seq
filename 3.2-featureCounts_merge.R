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

#workingdir
setwd("C:/Users/jhami/OneDrive/Documents/DTP/Second Year/")

STMfiles <- list.files(pattern = "STM")

TCP4files <- list.files(pattern = "TCP")

genotypes <- list(STMfiles, TCP4files)

################################################################################
# Main CMDs
################################################################################

# Make a dataframe that has the GeneID column so the columns from other
# featureCounts files can be added to it.
# Because the order of GeneIDs in feature counts is the same, any file can be
# used (can run example below).

#df1 <- read.table("STM_C3_S7.markdup.featurecounts.txt")
#df2 <-read.table("TCP_C2_S19.markdup.featurecounts.txt")
#print(df1$V1 %in% df2$V7)rm(df1, df2)

GeneID_col <- read.table("STM_C3_S7.markdup.featurecounts.txt")
GeneID_col <- pull(GeneID_col, var = +1, name = NULL)

markdupGeneID <- GeneID_col
rmdupGeneID <- GeneID_col

# Oject holding the genotype names and a counter to name files
Genotype_names <- c("STM", "TCP4")
Genotype_counter <-as.numeric(1)

# inner for loop to switch between genotypes

for (x in genotypes) {
  
  #print(x)
  
  for (file in x) {
    
    if (grepl("markdup", file)) {
      
      #head(file)
    
      # Take each file and add to the GeneID_col
      markdup_merging_file <-read.table(file)
      markdupGeneID <- bind_cols(markdupGeneID, markdup_merging_file$V7)
      tail(markdupGeneID)
      
  } else {
      
       head(file)
      
        rmdup_merging_file <-read.table(file)
        rmdupGeneID <- bind_cols(rmdupGeneID, rmdup_merging_file$V7)
        
        tail(rmdupGeneID)
        
  }
    

  }
  
  markdupGeneID <- as.matrix(markdupGeneID)
  write.csv(markdupGeneID, file = paste(Genotype_names[Genotype_counter],"_markdup.csv", sep = "_", row.names = ))
  
  rmdupGeneID <- as.matrix(rmdupGeneID)
  write(rmdupGeneID, file = paste(Genotype_names[Genotype_counter],"_rmdup.csv", sep = "_", row.names = ))
  
  # Reset GeneID objects
  markdupGeneID <- GeneID_col
  rmdupGeneID <- GeneID_col
  
  Genotype_counter <- Genotype_counter+1
    
}








