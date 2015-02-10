# abslevelsMF function ----------------------------------------------------
# The point of running this function is to find which samples would be best to 
# fragment for a given MF by calculating their absolute levels.
# Input:
#    1. MassFeature = the name of the mass feature you're interested in 
#       fragmenting AS IT APPEARS in DF.
#    2. DF = a data.frame from data that have been processed with xcms but NOT
#       preprocessed. Must contain the following columns:
#             a. "MassFeature" = unique identifier for all mass features, 
#                including the mass feature listed in #1. 
#             b. "RT" = the retention time in minutes, 
#             c. "mz" = the m/z for that mass feature, and 
#             d. sample columns with peak areas. 
#    3. SampNames = names of the columns that contain absolute abundances
#    4. n (optional) = number of top mass features you want. Defaults to 3.
# Output: A data.frame with 2 columns: 1. the top n samples by abundance and 
# 2. the abundance values.

abslevelsMF <- function(MassFeature, DF, SampNames, n = 3){
      require(tidyr)
      
      MF <- DF[DF$MassFeature == MassFeature, SampNames]
      MF <- gather(MF, Sample, Abundance)
      MF <- arrange(MF, -Abundance)
      
      return(MF[1:n, ])
}


# # EXAMPLE
# setwd("C:/Users/Laura/Documents/SCOR project/84 SCOR EposU")
# EposU84.unproc <- read.csv("SCORMDZEposU84 peak table - RT above 2 min.csv")
# 
# SampNames <- names(EposU84.unproc)[9:96]
# 
# BestMF <- abslevelsMF("I114.0647R5.05", EposU84.unproc, SampNames, 5)
