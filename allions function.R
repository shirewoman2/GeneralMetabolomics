# allions function --------------------------------------------------
# This function gathers into 1 data.frame the  EICs of the originally 
# identified mass feature plus all the possible ions CAMERA identified 
# as being associated with that mass feature. It uses the functions eic 
# and eicplot to do this.
# Input:
#       1. MF.df = a 1-row data.frame with the following columns:
#             a. the mass feature (MassFeature)
#             b. m/z (mz)
#             c. RT in minutes (RT)
#             d. ionization mode, as in "Epos", "Eneg", "Apos", or "Aneg" (Mode)
#             e. and matrix, as in "plasma", "urine", etc. (Matrix) 
#       2. Files = a data.frame with the following columns:
#             a. the names of the files, including ".mzdata.xml" (File)
#             b. the directory (Directory)
#             c. ionization mode (Mode)
#             d. and matrix (Matrix)
#       3. CameraList = the list output from running the camera function
# Output: a data.frame with the following columns:
#       1. Intensity - Intensity of the signal in counts
#       2. Dataset - the dataset as the mode and the matrix together
#       3. File - the file
#       4. MassFeature - the input mass feature
#       5. RT.original - the original RT of the mass feature in minutes
#       6. Matrix - the matrix
#       7. Mode - the mode
#       8. MassFeature.otherion - what the potentially matching ion is, i.e. 
#       "M if orig is M+K"
#       9. RT - the retention time for plotting the EIC
#       10. mz - the m/z of the potentially matching ion
#       11. IonType - the  type of ion the originally detected ion could be, 
#       i.e. "M+K"
#       12. Charge - the charge of the potentially matching ion
#       13. IsoGroup - if an isotope group was hypothesized (see CAMERA articles
#       for an explanation), what that group was


allions <- function(MF.df, Files, CameraList){
      
      require(plyr)
      require(dplyr)
      require(stringr)
      require(xcms)
      
      OrigDir <- getwd()
      setwd("I:/General LCMS scripts")
      source("eic function.R")
      setwd(OrigDir)
      
      OtherIonsIndex <- which(str_detect(names(CameraList), "otherions$"))
      Other.df <- CameraList[[OtherIonsIndex]]
      Other.df <- Other.df[Other.df$MassFeature == MF.df$MassFeature[1], ]
      for (m in 1:nrow(Other.df)){
            if(complete.cases(Other.df$mzOfM[m])){
                  Other.df$mz[m] <- as.numeric(Other.df$mzOfM[m])
            } else {
                  # If it's an adduct, look for the M+H or M-H ion since that's
                  # likely the most abundant peak.
                  if (str_detect(Other.df$Charge[m], "+")) {
                        Other.df$mz[m] <- Other.df$mzOfM[m] + 1.0073
                  } else {
                        Other.df$mz[m] <- Other.df$mzOfM[m] - 1.0073
                  }
                  
            }
      }
      
      # No need to plot ions if the main ion is just M. Removing those.
      Other.df <- Other.df[Other.df$IonType != "M", ]
      
      if (nrow(Other.df) > 0){
            
            for (m in 1:nrow(Other.df)){
                  if (Other.df$Charge[m] == "-" | Other.df$Charge[m] == "+") {
                        if (str_detect(Other.df$Charge[m], "+")){
                              Other.df$MassFeature.otherion[m] <- 
                                    paste("M+H if orig is", Other.df$IonType[m])
                        } else {
                              Other.df$MassFeature.otherion[m] <- 
                                    paste("M-H if orig is", Other.df$IonType[m])
                        }
                  } else {
                        Other.df$MassFeature.otherion[m] <- 
                              paste("ion with charge =", Other.df$Charge[m], 
                                    "\n if orig is", Other.df$IonType[m])
                  }
                  
            }
      }
      
      MF.df$MassFeature.otherion <- "originally detected ion"
      Other.df <- rbind.fill(Other.df, 
                             MF.df[, c("MassFeature", "mz", "RT", "MassFeature.otherion")])
      
      Other.df$Mode <- MF.df$Mode[1]
      Other.df$Matrix <- MF.df$Matrix[1]
      
      # Removing redundant rows
      Other.df <- subset(Other.df, MassFeature.otherion != "M+H if orig is M+H")
      Other.df <- subset(Other.df, MassFeature.otherion != "M-H if orig is M-H")
      
      EICs <- eic(Other.df, Files)
      
      EICs <- join(EICs, Other.df[, c("MassFeature", "mz", "IonType", 
                                      "mzOfM", "Charge", "IsoGroup",
                                      "MassFeature.otherion")], 
                   by = c("MassFeature", "MassFeature.otherion"))
      return(EICs)
}


# # EXAMPLE
# MF.df <- SCORMDZEnegP83.after2[SCORMDZEnegP83.after2$MassFeature == 
#                                      "I145.0677R5.63", 
#                                c("MassFeature", "mz", "RT")]
# MF.df$Mode <- "Eneg"
# MF.df$Matrix <- "plasma"
# 
# CameraList <- EnegP83.otherions
# 
# EICs.all <- allions(MF.df, Files, EnegP83.otherions)
