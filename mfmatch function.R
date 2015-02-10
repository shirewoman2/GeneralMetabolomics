# mfmatch function --------------------------------------------------

# This function will take two data.frames, X and Y (they can be called
# something other than that), and output a list of which mass features match 
# which in the other dataset. Data.frames must include the following columns:
#       1. MassFeature
#       2. mz
#       3. RT
# The default is that matches must be within 15 ppm and 0.2 min, but you can
# change the defaults when you call the function, eg. with defaults:
# mfmatch(Data1, Data2)
# with other settings for PPM and RTRange:
# mfmatch(Data1, Data2, PPM = 10, RTRange = 0.3)
# Output is a data.frame of all mass features in both original data.frames
# showing which ones match and what the m/z and RT differences are. The suffixes
# on the columns should make it clear which values came from which of the 
# original data.frames.


mfmatch <- function(X, Y, PPM = 15, RTRange = 0.2){
      require(plyr)
      require(stringr)
      
      # Renaming data.frames and making sure that the mass features are unique.
      DF.X <- X[unique(X$MassFeature), c("MassFeature", "mz", "RT")]
      DF.Y <- Y[unique(Y$MassFeature), c("MassFeature", "mz", "RT")]
      
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
            
            if (is.na(DF.X$RT.X[i])){
                  MFmatch[[i]] <- DF.Y[
                        DF.Y$mz.Y > (mz.X[i] - (PPM/1e6*mz.X[i])) 
                        & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])), ]
                  
            } else {
                  
                  MFmatch[[i]] <- DF.Y[
                        DF.Y$mz.Y > (mz.X[i] - (PPM/1e6*mz.X[i])) 
                        & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])) 
                        & DF.Y$RT.Y > RT.X[i] - RTRange 
                        & DF.Y$RT.Y < RT.X[i] + RTRange, ]
                  MFmatch[[i]] <- MFmatch[[i]][
                        complete.cases(MFmatch[[i]]$mz.Y), ]
                  
                  if (nrow(MFmatch[[i]]) > 0) {
                        MFmatch[[i]] <- 
                              rbind(MFmatch[[i]], 
                                    DF.Y[DF.Y$mz.Y > (mz.X[i] - 
                                                            (PPM/1e6*mz.X[i])) 
                                         & DF.Y$mz.Y < (mz.X[i] + 
                                                              (PPM/1e6*mz.X[i])) 
                                         & is.na(DF.Y$RT.Y), ])
                        
                        
                  } else {
                        MFmatch[[i]] <- 
                              DF.Y[DF.Y$mz.Y > (mz.X[i] - 
                                                      (PPM/1e6*mz.X[i])) 
                                   & DF.Y$mz.Y < (mz.X[i] + 
                                                        (PPM/1e6*mz.X[i])) 
                                   & is.na(DF.Y$RT.Y), ]
                        
                  }
            }

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
      
      Matches <- dlply(Matches, "MassFeature.X")
      
      # If there was more than 1 potential match and the RT was listed, take 
      # the one that was closest by m/z and then by RT. If there was no RT 
      # listed, take all matches. 
      # Calculate the difference in m/z in ppm and RT in min.
      for (i in which(Matched.Y > 0)){
            
            if (is.na(Matches[[i]]$RT.X)){
                  Matches[[i]] <- MFmatch[[i]]
                  Matches[[i]]$ppm <- 
                        abs((mz.X[i] - Matches[[i]]$mz.Y)/mz.X[i]*1e6)
                  Matches[[i]]$MassFeature.X <- MFname.X[i]
                  
            } else {
                  
                  for (n in 1:nrow(MFmatch[[i]])){
                        MFmatch[[i]]$ppm[n] <- 
                              abs((mz.X[i] - MFmatch[[i]]$mz.Y[n])/mz.X[i]*1e6)
                        MFmatch[[i]]$RTdif[n] <- 
                              abs(RT.X[i] - MFmatch[[i]]$RT.Y[n])
                        MFmatch[[i]]$MassFeature.X <- MFname.X[i]
                  }
                  MFmatch[[i]] <- arrange(MFmatch[[i]], ppm, RTdif)
                  
                  Matches[[i]]$MassFeature.Y <- 
                        as.character(MFmatch[[i]]$MassFeature.Y[1])
                  Matches[[i]]$mz.Y <- MFmatch[[i]]$mz.Y[1]
                  Matches[[i]]$RT.Y <- MFmatch[[i]]$RT.Y[1]
                  Matches[[i]]$ppm <- MFmatch[[i]]$ppm[1]
                  Matches[[i]]$RTdif <- MFmatch[[i]]$RTdif[1]
                  
            }
      }
      
      Matches <- rbind.fill(Matches)
      
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
      
      OrigNames <- names(Matches)
      
      names(Matches)[str_detect(OrigNames, "Y$")] <- sub(
            "Y$", ColSuffix.Y, OrigNames)[str_detect(OrigNames, "Y$")]
      names(Matches)[str_detect(OrigNames, "X$")] <- sub(
            "X$", ColSuffix.X, OrigNames)[str_detect(OrigNames, "X$")]
      
      return(Matches)
      
}


# # EXAMPLE
# # Using mfmatch to compare MFs in SCOR MDZ EposU84 and in 
# # CYP2D6 DEX EposU2.
# 
# setwd("G:/Data/Metabolomics/Jessica/CYP2D6/mzData files 20131230/ESI+")
# CYP2D6DEX.EposU2 <- read.csv("CYP2D6 EposU 2 peak table - RT above 2 min.csv")
# CYP2D6DEX.EposU2 <- CYP2D6DEX.EposU2[, c("MassFeature", "mz", "RTmin")]
# CYP2D6DEX.EposU2 <- plyr::rename(CYP2D6DEX.EposU2, c("RTmin" = "RT"))
# 
# 
# setwd("C:/Users/Laura/Documents/SCOR project/84 SCOR EposU")
# SCORMDZ.EposU84 <- read.csv("SCORMDZEposU84 peak table - RT above 2 min.csv")
# SCORMDZ.EposU84 <- SCORMDZ.EposU84[, c("MassFeature", "mz", "RT")]
# 
# SCORMDZ.CYP2D6DEX.EposU <- mfmatch(SCORMDZ.EposU84, CYP2D6DEX.EposU2)
# 
