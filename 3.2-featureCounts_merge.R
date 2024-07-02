################################################################################
# Merging featureCount output files 
################################################################################

# Want to add all of the .featureCount.txt files together into one file

################################################################################
# Loading libraries and setting up environement
################################################################################

# Needed for dataframe manipulation
#install.packages("dplyr")
library("dplyr")

# Workingdir
setwd("C:/Users/c1831460/OneDrive - Cardiff University/Documents/DTP/Second Year/RNA-seq TCP4_STM/featureCounts/")

# Making a list object containing files, separated by timepoint
files3h <- list.files(pattern = "3h")

files12h <- list.files(pattern = "12h")

files24h <- list.files(pattern = "24h")

files11d <- list.files(pattern = "11d")

# Making a list of file lists
timepoints <- list(files3h, files12h, files24h, files11d)

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

GeneID_col <- read.table("Col-0_3h_1_S10_markdup_featurecounts.txt", 
                         header = TRUE, 
                         sep = "\t")
GeneID_col <- pull(GeneID_col, 
                   var = +1, 
                   name = NULL)
GeneID_col <- as.matrix(GeneID_col)

# Initialise separate matricies for markdup and rmdup to be used in the if else loop
markdupGeneID <- GeneID_col
rmdupGeneID <- GeneID_col

# Character vector holding the timepoint names and a counter to name files
timepoint_names <- c("3h", "12h", "24h", "11d")
timepoint_counter <-as.numeric(1)

# Outer for loop to switch between timepoints
for (x in timepoints) {
  
  # Inner for loop to iterate through file in the file list
  for (file in x) {
    
    # Inner if loop to delineate markdup and rmdup files in list with grep
    if (grepl("markdup", file)) {
      
      # Pass markdup file into markdup_merging_file and add the 7th coloumn
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
  write.csv(markdupGeneID, file = paste0("../",
    timepoint_names[timepoint_counter],"_markdup.csv"), row.names = FALSE)
  
  # Write the rmdup file
  write.csv(rmdupGeneID, file = paste0("../",
    timepoint_names[timepoint_counter],"_rmdup.csv"), row.names = FALSE)
  
  # Reset GeneID matricies
  markdupGeneID <- GeneID_col
  rmdupGeneID <- GeneID_col
  
  timepoint_counter <- timepoint_counter+1
    
}
