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
#             d. ionization mode (Mode)
#             e. and matrix (Matrix) 
#       2. Files = a data.frame with the following columns:
#             a. the files (File)
#             b. directory (Dir)
#             c. ionization mode (Mode)
#             d. and matrix (Matrix)
#       3. CameraList = the list output from running the camera function

allions <- function(MF.df, Files, CameraList){
      
      require(plyr)
      require(dplyr)
      
      OtherIonsIndex <- which(str_detect(names(CameraList), "otherions$"))
      Other.df <- CameraList[[OtherIonsIndex]]
      Other.df <- Other.df[Other.df$MassFeature == MF.df$MassFeature, ]
      Other.df$mz <- as.numeric(Other.df$Othermz)
      
      # No need to plot ions if the main ion is just M. Removing those.
      Other.df <- Other.df[Other.df$IonType != "M", ]
      
      if (nrow(Other.df) > 0){
            
            for (m in 1:nrow(Other.df)){
                  if (Other.df$Charge[m] == "-" | Other.df$Charge[m] == "+") {
                        Other.df$MassFeature.ion[m] <- paste("M if orig is", 
                                                             Other.df$IonType[m])
                  } else {
                        Other.df$MassFeature.ion[m] <- 
                              paste("ion with charge =", Other.df$Charge[m], 
                                    "\n if orig is", Other.df$IonType[m])
                  }
            }
      }
      
      MF.df$MassFeature.ion <- "originally detected ion"
      Other.df <- rbind.fill(Other.df, 
                             MF.df[, c("MassFeature", "mz", "RT", "MassFeature.ion")])
      
      Other.df$Mode <- MF.df$Mode[1]
      Other.df$Matrix <- MF.df$Matrix[1]
      
      EICs <- eic(Other.df, Files)
      EICs <- join(EICs, Other.df[, c("MassFeature", "mz", "IonType", 
                                      "Charge", "IsoGroup", "MassFeature.ion")], 
                   by = c("MassFeature", "MassFeature.ion"))
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
