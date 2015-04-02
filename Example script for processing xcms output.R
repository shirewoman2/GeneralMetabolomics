# Example processing script

# This script processes the data from xcms.


# Housekeeping --------------------------------------------
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(lubridate)
library(stringr)

Dataset <- c("SFNEposP8", "SFNEnegP9", "SFNEposU10", "SFNEnegU11")
Datasets <- data.frame(Dataset = Dataset,
                       Mode = rep(c("Epos", "Eneg"), 2),
                       Matrix = rep(c("plasma", "urine"), each = 2))

Directory <- list("C:/Users/Laura/Documents/Sulforaphane project/8 SulfEposP",
                  "C:/Users/Laura/Documents/Sulforaphane project/9 SulfEnegP", 
                  "C:/Users/Laura/Documents/Sulforaphane project/SFNEposU10",
                  "C:/Users/Laura/Documents/Sulforaphane project/SFNEnegU11")
names(Directory) <- Dataset

MainDir <- c("C:/Users/Laura/Documents/Sulforaphane project")

# Loading data
RDataFiles <- paste(Dataset, "- filtered dataframe only.RData")
names(RDataFiles) <- Dataset

for (j in Dataset){
      setwd(Directory[j])
      load(RDataFiles[j])
}

Data.filtered <- c(SulfEposP8.filter, SulfEnegP9.filter, 
                      SulfEposU10.filter, SFNEnegU11.filter)
names(Data.filtered) <- Dataset

# Loading metadata
setwd(MainDir)
load("SFN metadata.RData")

# Selecting the files that were used with the xcms script
SampCol <- unlist(lapply(Data.filtered, function(x) names(x)))
GoodFiles <- Files[Files$SampCol %in% SampCol, ]

# Adding Dataset
GoodFiles <- join(GoodFiles, Datasets, by = c("Mode", "Matrix"))

# Breaking up by Dataset for matching elsewhere in script
GoodFiles <- dlply(GoodFiles, "Dataset")

GoodSamples <- list()
ClinSamples <- list()

# Selecting only the columns of interest in the data and renaming them sensibly
for (j in Dataset){
      
      GoodSamples[[j]] <- GoodFiles[[j]]$SampleID
      ClinSamples[[j]] <- GoodFiles[[j]]$SampleID[
            GoodFiles[[j]]$SampType == "clinical"]
      
      Data.filtered[[j]] <- Data.filtered[[j]][, c("MassFeature", "mz", "RT", 
                                                   GoodFiles[[j]]$SampCol)]
      names(Data.filtered[[j]]) <- c("MassFeature", "mz", "RT", 
                                     GoodFiles[[j]]$SampleID)
      
}



# Functions and themes ----------------------------------
# Theme for graphs made using ggplot2
ThemeLaura <- function (base_size = 12, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace% 
            theme(
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(colour = "black"),
                  axis.title.y = element_text(colour = "black", angle=0),
                  panel.background = element_rect(fill="white", color=NA),
                  panel.grid.minor.y = element_line(color="white"),
                  panel.grid.minor.x = element_line(color="white"),
                  panel.grid.major = element_line(colour = "white"),
                  plot.background = element_rect(fill="white", color=NA),
                  panel.border = element_rect(color="black", fill=NA),
                  strip.background = element_rect(color=NA, fill="white"),
                  legend.background = element_rect(color=NA, fill="white"),
                  legend.key = element_rect(color=NA, fill="white")
            )   
}

theme_set(ThemeLaura())


# Function for writing files with a note on the first line. 
my.write <- function(x, file, header, f = write.csv, ...){
      # create and open the file connection
      datafile <- file(file, open = 'wt')
      # close on exit
      on.exit(close(datafile))
      # if a header is defined, write it to the file
      if(!missing(header)) writeLines(header,con=datafile)
      # write the file using the defined function and required addition arguments  
      f(x, datafile,...)
}



# Processing the data ----------------------------------------------

# Starting from data that were already filtered by frequency to retain ions 
# detected initially in >= 25% of all samples. 

Hist <- list()
Data.TCCnorm <- list()
Data.log <- list()

for (j in Dataset){
      
      # Adding 1 to all sample columns
      Data.filtered[[j]][ , GoodSamples[[j]]] <- 
            Data.filtered[[j]][ , GoodSamples[[j]]] + 1
      
      # Calculate the TCC sum.
      TCCsum <- apply(Data.filtered[[j]][ , GoodSamples[[j]]], 
                      2, sum, na.rm=TRUE)
      
      
      # Normalize each sample by the TCC sum. 
      Data.TCCnorm[[j]] <- Data.filtered[[j]][ , GoodSamples[[j]]]
      for (s in 1:length(GoodSamples[[j]])){
            Data.TCCnorm[[j]][, GoodSamples[[j]][s]] <- 
                  (Data.filtered[[j]][, GoodSamples[[j]][s]]/ 
                         TCCsum[s])*1e6
            rm(s)
      }  
      
      
      # log10 transforming
      Data.log[[j]] <- data.frame(MassFeature = Data.filtered[[j]]$MassFeature, 
                                  mz = Data.filtered[[j]]$mz, 
                                  RT = Data.filtered[[j]]$RT,
                                  log10(Data.TCCnorm[[j]]))
      
      
      # Saving the files
      setwd(as.character(Directory[j]))
      FileNote <- paste("Dataset:", j,
                        "This file was generated on",Sys.Date(), 
                        "using the script SCOR MDZ preprocessing.R.")
      my.write(Data.log[[j]], paste(j,"- preprocessed.csv"), 
               header=FileNote, row.names=F)
      
      
      # Checking distributions. 
      Data.melt <- gather(Data.log[[j]], Sample, Abundance, -MassFeature, 
                          -mz, -RT)
      
      
      Hist[[j]] <- ggplot(Data.melt, aes(x=Abundance)) + 
            geom_density() +
            xlab("log10(Abundance)") + ggtitle(j)
      
      rm(Data.melt, TCCsum)
      
}


setwd(MainDir)
if (length(Dataset) == 2){
      png(paste(str_c(Dataset, collapse = " "), "Kernel density plots of MF abundances.png"), 
          height = 3, width = 8,
          units = "in", res=600)
      grid.arrange(Hist[[1]], Hist[[2]], ncol=2)
      dev.off()
      
} else {
      png(paste(str_c(Dataset, collapse = " "), "Kernel density plots of MF abundances.png"), 
          height = 6, width = 8,
          units = "in", res=600)
      grid.arrange(Hist[[1]], Hist[[2]], Hist[[3]], Hist[[4]], ncol=2)
      dev.off()
}

save(Data.log, Directory, Dataset, file="SFN main data.RData")

