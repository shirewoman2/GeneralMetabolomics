# rellevelsMF function ------------------------------------------------
# The point of running this function is to find the relative levels of any
# mass feature in preprocessed data.
# Input:
#    1. MassFeature = the name of the mass feature you're interested in 
#       fragmenting AS IT APPEARS in DF.
#    2. DF = a data.frame from data that have been preprocessed. Must contain
#       the following columns:
#             a. "MassFeature" = unique identifier for all mass features, 
#                including the mass feature listed in #1. 
#             b. "RT" = the retention time in minutes, 
#             c. "mz" = the m/z for that mass feature, and 
#             d. sample columns with peak areas. 
#    3. SampNames = names of the columns that contain relative abundance data
# Output: A long-format data.frame with 2 columns: 1. the sample, 2. the 
# relative abundance for the given mass feature.

rellevelsMF <- function(MassFeature, DF, SampNames){
      require(tidyr)
      
      MF <- DF[DF$MassFeature == MassFeature, SampNames]
      MF <- gather(MF, Sample, Abundance)
      
      return(MF)
}


# # EXAMPLE
# EposU84.proc <- read.csv("EposU84 - preprocessed.csv", skip = 1)
# SampleColumns <- names(EposU84.proc)[4:90]
# 
# I281.1498R5.72 <- rellevelsMF("I281.1498R5.72", EposU84.proc, SampleColumns)
