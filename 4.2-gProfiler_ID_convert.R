################################################################################
# Loading dependencies
################################################################################

library("dplyr")
library("gprofiler2")

################################################################################
# Adding Gene Names to DEG files
################################################################################

# Point to the working dir

setwd("your/working/dir/DEGs")


# Read in all the files
files <- c(list.files("3h/", include.dirs = FALSE, full.names = TRUE),
           list.files("12h/", include.dirs = FALSE, full.names = TRUE),
           list.files("24h/", include.dirs = FALSE, full.names = TRUE),
           list.files("11d/", include.dirs = FALSE, full.names = TRUE))

# Creating subdirectories for output files - these lines should not be run before
# line 18

dir.create("3h/NamesAdded")
dir.create("12h/NamesAdded")
dir.create("24h/NamesAdded")
dir.create("11d/NamesAdded")

for (i in files) {

  # Read in each file as a csv  
  csvfile <- read.csv(file = paste0(i))
  # Change the name of first column to "geneID"
  names(csvfile)[1] <- "GeneID"
  # Extract the first/GeneID column and pass it into the query object to be read
  # by gconvert
  IDcolumn <- as.matrix(csvfile %>% select(1))
  query <- as.character(IDcolumn)

  GeneNames <- gconvert(query,
                        organism = "athaliana",
                        target = "ENSG",
                        numeric_ns = "",
                        mthreshold = 1,
                        filter_na = TRUE
                        )
  
  # Extract the name and description columns from the output
  GeneNameColumns <- GeneNames %>% 
    select(name, description)
  
  # Insert the name and description columns in between the first and second 
  # columns from the initial csv
  outputdf <- bind_cols(select(csvfile, 1), 
                        GeneNameColumns, 
                        select(csvfile, -1))
  outputdf <- as.data.frame(outputdf)
  
  # Remove parent directory (e.g. STM/) from the string held in the i variable
  i <- sub(".*/", "", i)

  # Write files into appropriate subdirectories depending on genotype
  if (grepl("3h", i)) {
      write.csv(outputdf, file = file.path("3h/NamesAdded", paste0(i, "")), 
                row.names = FALSE)
      } else if (grepl("12h", i)) {
        write.csv(outputdf, file = file.path("12h/NamesAdded", paste0(i, "")), 
                  row.names = FALSE)
      } else if (grepl("24h", i)) {
        write.csv(outputdf, file = file.path("24h/NamesAdded", paste0(i, "")), 
                  row.names = FALSE)
      } else if (grepl("11d", i)) {
        write.csv(outputdf, file = file.path("11d/NamesAdded", paste0(i, "")), 
                  row.names = FALSE)
    }
  }

