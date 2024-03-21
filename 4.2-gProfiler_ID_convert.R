###################################################################################
# Loading dependencies
################################################################################

library("dplyr")
library("gprofiler2")

################################################################################
# Adding Gene Names to DEG files
################################################################################

# Point to the working dir

setwd("C:/Users/c1831460/OneDrive - Cardiff University/Documents/DTP/Second Year/RNA-seq/featureCounts/DEGs/")


#files <- "STM/STM_C_vs_CD_markdup_DEGs.csv"
files <- c(list.files("TCP4/", include.dirs = FALSE, full.names = TRUE), list.files("STM/", include.dirs = FALSE, full.names = TRUE))

dir.create("STM/NamesAdded")
dir.create("TCP4/NamesAdded")

for (i in files) {
  
  csvfile <- read.csv(file = paste0(i))
  names(csvfile)[1] <- "GeneID"
  IDcolumn <- as.matrix(csvfile %>% select(1))
  query <- as.character(IDcolumn)

  GeneNames <- gconvert(query,
                        organism = "athaliana",
                        target = "ENSG",
                        numeric_ns = "",
                        mthreshold = 1,
                        filter_na = TRUE
                        )
  
  GeneNameColumns <- GeneNames %>% 
    select(name, description)
  
  outputdf <- bind_cols(select(csvfile, 1), 
                        GeneNameColumns, 
                        select(csvfile, -1))
  outputdf <- as.data.frame(outputdf)
  
  i <- sub(".*/", "", i)

  if (grepl("STM", i)) {
    write.csv(outputdf, file = file.path("STM/NamesAdded", paste0(i, "")), row.names = FALSE)
  } else {
    write.csv(outputdf, file = file.path("TCP4/NamesAdded", paste0(i, "")), row.names = FALSE)
  }
}
