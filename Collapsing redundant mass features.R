# Collapsing redundant mass features
# This script collapses features that should have been collapsed into one upon peak alignment. 
# Updated 2/9/13 -LS


# Calling up the data. 
# Set the directories where the original data are.
DirEposP <- "D:/My Documents/Work/Lin Lab/SCOR project/ESI+/Metabolomics plasma samples/SCOR plasma V"
DirEnegP <- "D:/My Documents/Work/Lin Lab/SCOR project/ESI-/ESI- XV"
DirAposP <- "D:/My Documents/Work/Lin Lab/SCOR project/APCI+/APCI+ VIII"
Directory <- c(DirEposP, DirEnegP, DirAposP)


FileIn <- c("20130209 SCOR EposP VA1.csv", "20130208 SCOR EnegP XV mass and RT.csv", "20130209 SCOR AposP VIIIA1.csv")
FileOut <- c("20130209 SCOR EposP VA1a collapsed.csv", "20130210 SCOR EnegP XVA collapsed - 15 ppm _2 min.csv", "20130209 SCOR AposP VIIIA1a collapsed.csv")

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


FileNote <- paste("This file was generated using the R script 'Collapsing redundant mass features.R' on", Sys.Date())


for (k in 1:3){
  
  # Opening file. Format: Samples in columns. 
  # 1st three columns should be 1. compound, 2. mass, 3. retention time. Subsequent columns should be samples.
  setwd(Directory[k])
  MassFeatures.A <- read.csv(FileIn[k], header=T)
  names(MassFeatures.A)[1:3] <- c("Compound", "Mass", "RT")
  MassFeatures.A$Compound <- as.character(MassFeatures.A$Compound)
  MassFeatures.B <- MassFeatures.A
  
  # Fill in the column numbers that contain samples. 
  SampleColumns <- c(4:length(MassFeatures.A))
  
  # Create a vector containing the samples. 
  Samples <- c(colnames(MassFeatures.A)[SampleColumns])
  
  # Set the mass accuracy desired in ppm. The final window will be double this.
  MassRange.ppm <- 15
  
  # Set the acceptable RT range in min. The final window will be double this.
  RTRange <- 0.2
  
  # Create a list where all the datasets will be stored. 
  MFmatch <- list()
  
  # For each entry in the MassFeatures.B dataset, check whether there are compounds
  # in the MassFeatures.A dataset that are within the RT and mass accuracy range
  # you set above. Make each of the dataframes into a list item in MFmatch.
  for (i in 1:nrow(MassFeatures.B)){
    B.Mass.i <- MassFeatures.B[i,2]
    B.RT.i <- MassFeatures.B[i,3]
    MFmatch[[i]] <- MassFeatures.A[MassFeatures.A$Mass>(B.Mass.i-(MassRange.ppm/1e6*B.Mass.i)) &
                                     MassFeatures.A$Mass<(B.Mass.i+(MassRange.ppm/1e6*B.Mass.i)) & 
                                     MassFeatures.A$RT>B.RT.i-RTRange & MassFeatures.A$RT<B.RT.i+RTRange,]
    
    MFmatch[[i]]$ppm <- c(1:nrow(MFmatch[[i]]))
    MFmatch[[i]]$RTdif <- c(1:nrow(MFmatch[[i]]))
    for (j in 1:nrow(MFmatch[[i]])){
      MFmatch[[i]]$ppm[j] <- abs((B.Mass.i-MFmatch[[i]]$Mass[j])/B.Mass.i*1e6)
      MFmatch[[i]]$RTdif[j] <- abs(B.RT.i-MFmatch[[i]]$RT[j])
    }
    
    rm(i, B.Mass.i, B.RT.i)
  }
  
  
  
  # Making a new data frame where the collapsed data will be stored.
  GoodData1 <- data.frame()
  
  
  for (i in 1:length(MFmatch)){
    GoodData1[i,1] <- MFmatch[[i]][1,1]
    
    for (j in 2:length(MassFeatures.A)){
      GoodData1[i,j] <- mean(MFmatch[[i]][,j], na.rm=T)  
    }
  }
  
  names(GoodData1) <- names(MassFeatures.A)
  GoodData2 <- unique(GoodData1)
  GoodData2$Compound2 <- GoodData2$Compound
  
  for (i in 1:nrow(GoodData2)){
    GoodData2$Compound2[i] <- paste("M", format(GoodData2$Mass[i], digits=8), "R", 
                                    format(GoodData2$RT[i], digits=8), sep="")
    
  }
  
  GoodData <- data.frame(GoodData2$Compound2, GoodData2[,2:(length(GoodData2)-1)])
  names(GoodData)[1] <- "Compound"
  
  
  
  
  my.write(GoodData, FileOut[k], row.names=F, header=FileNote)
  
  
}