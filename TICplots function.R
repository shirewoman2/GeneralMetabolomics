# TICplots function

# This function plots TICs from mzdata files.
# Input:
#       a. Files - a data.frame with the following columns:
#             1. File - file names, including ".mzdata.xml"
#             2. SampleID (if you want to facet by sample ID)
# Output:
#       a. Plot - a ggplot object "Plot" with all TICs, facetted by File 
#       (default) or by SampleID when SampleID = TRUE
#       b. TICs - a data.frame with the TIC data when KeepTICs = TRUE

TICplots <- function(Files, SampleID = FALSE, KeepTICs = FALSE){
      require(xcms)
      require(plyr)
      require(ggplot2)
      
      TICs <- list()
      
      for (i in 1:nrow(Files)){
            Raw <- xcmsRaw(Files$File[i], profstep=0.01, profmethod="bin")
            TICs[[i]] <- data.frame(plotTIC(Raw))
            
            if (SampleID == TRUE){
                  TICs[[i]]$SampleID <- Files$SampleID[i]
            }
            
            TICs[[i]]$File <- Files$File[i]
                  
      }
      
      TICs <- rbind.fill(TICs)
      TICs <- plyr::rename(TICs, c("X1" = "Time", "X2" = "Intensity"))
      
      if (SampleID == FALSE) {
            Plot <<- ggplot(TICs, aes(x = Time, y = Intensity, color = File)) +
                  geom_line() +
                  facet_wrap(~ File)
      } else {
            Plot <<- ggplot(TICs, aes(x = Time, y = Intensity, color = File)) +
                  geom_line() +
                  facet_wrap(~ SampleID)
      }
      
      Plot
      
      if (KeepTICs == TRUE){
            TICs <<- TICs
      }
      
}