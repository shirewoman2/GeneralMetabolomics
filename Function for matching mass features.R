library(plyr)

# This function will take two data.frames, X and Y (they can be called
# something other than that), and output a list of which mass features match 
# which in the other dataset. Data.frames must include the following columns:
#       1. MassFeature
#       2. mz
#       3. RT
# The default is that matches must be within 15 ppm and 0.2 min, but you can
# change the defaults when you call the function, eg. with defaults:
# MassFeatureMatch(Data1, Data2)
# with other settings for PPM and RTRange:
# MassFeatureMatch(Data1, Data2, PPM = 10, RTRange = 0.3)
# Output is a data.frame of all mass features in both original data.frames
# showing which ones match and what the m/z and RT differences are. The suffixes
# on the columns should make it clear which values came from which of the 
# original data.frames.


MassFeatureMatch <- function(X, Y, PPM = 15, RTRange = 0.2){
      DF.X <- X[, c("MassFeature", "mz", "RT")]
      DF.Y <- Y[, c("MassFeature", "mz", "RT")]
      
      # Renaming to keep track of which mass feature, m/z, and RT came from 
      # which dataset.
      names(DF.X) <- paste(names(DF.X), "X", sep = ".")
      names(DF.Y) <- paste(names(DF.Y), "Y", sep = ".")
      
      MFmatch <- list()
      Matched.Y <- c()
      MFname.X <- as.character(DF.X$MassFeature)
      mz.X <- DF.X$mz
      RT.X <- DF.X$RT
      
      # Checking each row in DF.X for any matches in DF.Y
      for (i in 1:nrow(DF.X)){
            
            MFmatch[[i]] <- DF.Y[DF.Y$mz.Y > (mz.X[i] - (PPM/1e6*mz.X[i])) 
                                 & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])) 
                                 & DF.Y$RT.Y > RT.X[i] - RTRange 
                                 & DF.Y$RT.Y < RT.X[i] + RTRange, ]
            
            Matched.Y[i] <- as.numeric(nrow(MFmatch[[i]]))
      }
      
      # Making a new data.frame to hold matched mass features
      Matches <- data.frame(MassFeature.X = DF.X$MassFeature.X,
                            MassFeature.Y = NA, 
                            mz.X = DF.X$mz.X, 
                            mz.Y = NA,
                            RT.X = DF.X$RT.X,
                            RT.Y = NA,
                            NumMatched = Matched.Y,
                            ppm = NA, 
                            RTdif = NA)
      
      # If there was more than 1 potential match, take the one that was closest
      # by m/z and then by RT. Calculate the difference in ppm and RT in min.
      for (i in which(Matched.Y > 0)){
            for (n in 1:nrow(MFmatch[[i]])){
                  MFmatch[[i]]$ppm[n] <- abs((mz.X[i]-MFmatch[[i]]$mz.Y[n])/mz.X[i]*1e6)
                  MFmatch[[i]]$RTdif[n] <- abs(RT.X[i]-MFmatch[[i]]$RT.Y[n])
                  MFmatch[[i]]$MassFeature.X <- MFname.X[i]
            }
            MFmatch[[i]] <- arrange(MFmatch[[i]], ppm, RTdif)
            
            Matches$MassFeature.Y[i] <- as.character(MFmatch[[i]]$MassFeature.Y[1])
            Matches$mz.Y[i] <- MFmatch[[i]]$mz.Y[1]
            Matches$RT.Y[i] <- MFmatch[[i]]$RT.Y[1]
            Matches$ppm[i] <- MFmatch[[i]]$ppm[1]
            Matches$RTdif[i] <- MFmatch[[i]]$RTdif[1]
            
      }
      
      # Add to the Matches data.frame any MFs in Y that weren't found in X.
      NoMatch <- DF.Y[DF.Y$MassFeature %in% 
                            setdiff(DF.Y$MassFeature, Matches$MassFeature.Y), ]
      
      Matches <- rbind.fill(Matches, NoMatch)
      Matches$mz.mean <- apply(Matches[, c("mz.X", "mz.Y")], 1, FUN = mean, 
                               na.rm=T)
      Matches <- arrange(Matches, mz.mean)
      
      names(Matches)
            
      ColSuffix.X <- as.character(deparse(substitute(X)))
      ColSuffix.Y <- as.character(deparse(substitute(Y)))
      
      names(Matches) <- sub("Y", ColSuffix.Y, names(Matches))
      names(Matches) <- sub("X", ColSuffix.X, names(Matches))
       
      return(Matches)
      
}


# Example of using the function MassFeatureMatch -----------------------------

# X
setwd("F:/Busulfan/Busulfan_Postdose_Data_10_2014/Busulfan EposP3")
BusulfCurrent <- read.csv("BusulfEposP3 preprocessed.csv", skip = 1)

# Y
setwd("G:/Data/Metabolomics/Laura/Busulfan/Initial pilot data/Busulf EposP/Busulf EposP V")
BusulfPilot <- read.csv("Busulf EposP VA.csv")
# This used "Mass" instead of m/z, so I'm switching back to m/z.
BusulfPilot$mz <- BusulfPilot$Mass + 1.0073


# Using the function to find matches.
BusulfMatches <- MassFeatureMatch(BusulfCurrent, BusulfPilot)



