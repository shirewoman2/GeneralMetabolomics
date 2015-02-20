# findfiles function


# This function searches for a given file and produces the full path for that 
# file.
# Input: 
#       1. File - a file name with its extension (character string)
#       2. Dir - the starting directory to search (character). This can be 
#          multiple levels in a tree or just the drive.
# Output: the path for that file

findfile <- function(File, StartDir) {
      require(stringr)
      
      OrigDir <- getwd()
      setwd(StartDir) 
      
      # Getting all the subdirectories
      AllDir <- sub("./", "", list.dirs(recursive = TRUE))
      
      # Removing git directories
      AllDir <- AllDir[!AllDir %in% AllDir[str_detect(AllDir, ".git")]]
      
      for (d in 1:length(AllDir)){
            setwd(AllDir[d])
            OutDir <- getwd()
            if (file.exists(File) == TRUE) { break }
            setwd(StartDir)
      }

      setwd(OrigDir)
      return(OutDir)
}

# # Example
# File <- "V2_predose_neg.xlsx"
# File <- "CYP2D6 DEX existent samples.xlsx"
# StartDir <- "D:/Users/Laura/Documents/Work/Lin Lab/CYP2D6 DEX"
# 
# findfile("Comparison of XCMS and MPP parameters.xlsx",
#          "D:/Users/Laura/Documents/Work/Lin Lab/CYP2D6 DEX")
# 
