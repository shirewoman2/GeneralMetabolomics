# eic function ----------------------------------------------------
# This function extracts EIC data for use in plotting EICs. Input:
#       1. MF.df = a data.frame with the following columns:
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
# for each file you want to extract the EIC data from.
# Column names must match. Output is a data.frame.


eic <- function(MF.df, Files, ppm = 15, PrintProgress = FALSE) {
      
      require(plyr)
      require(xcms)
      OrigDir <- getwd()
      MF.df$ModeMatrix <- paste(MF.df$Mode, MF.df$Matrix, sep = ".")
      Files$ModeMatrix <- paste(Files$Mode, Files$Matrix, sep = ".")
      
      # Only keeping files with the mode and matrix actually present in 
      # the list of mass features of interest.
      Files <- Files[Files$ModeMatrix %in% unique(MF.df$ModeMatrix), ]
      
      MF.list <- dlply(MF.df, c("Mode", "Matrix"))
      EIC <- list()
      Files$Directory <- as.character(Files$Directory)
      
      for (i in 1:nrow(Files)){
            
            setwd(Files$Dir[i])
            RawData <- xcmsRaw(paste0(Files$File[i], ".mzdata.xml"), 
                               profstep = 0, profmethod = "bin")
            
            Dataset <- paste(Files$Mode[i], Files$Matrix[i], sep = ".")
            
            EIC[[i]] <- list()
            
            for (m in 1:nrow(MF.list[[Dataset]])){
                  
                  MZRange <- cbind(MF.list[[Dataset]]$mz[m] - 
                                         ppm*MF.list[[Dataset]]$mz[m]/1e6,
                                   MF.list[[Dataset]]$mz[m] + 
                                         ppm*MF.list[[Dataset]]$mz[m]/1e6)
                  
                  RTRange <- t(cbind(MF.list[[Dataset]]$RT[m]*60 + c(-30, 30)))
                  EICdata <- getEIC(RawData, 
                                    mzrange = MZRange,
                                    rtrange = RTRange)
                  EIC[[i]][[m]] <- data.frame(EICdata@eic$xcmsRaw[[1]]) # names: "rt" and "intensity"
                  EIC[[i]][[m]]$Dataset <- Dataset
                  EIC[[i]][[m]]$File <- Files$File[i]
                  EIC[[i]][[m]]$MassFeature <- MF.list[[Dataset]]$MassFeature[m]
                  EIC[[i]][[m]]$RT.original <- MF.list[[Dataset]]$RT[m]
                  EIC[[i]][[m]]$Matrix <- MF.list[[Dataset]]$Matrix[m]
                  EIC[[i]][[m]]$Mode <- MF.list[[Dataset]]$Mode[m]
                  
                  if ("MassFeature.ion" %in% names(MF.df)) {
                        EIC[[i]][[m]]$MassFeature.ion <- MF.df$MassFeature.ion[m]
                  }
                  
            }
            
            EIC[[i]] <- rbind.fill(EIC[[i]])
            EIC[[i]]$Project <- Files$Project[i]
            
            if (PrintProgress == TRUE) { print(i) }
      }
      
      EIC <- rbind.fill(EIC)
      EIC$RT <- EIC$rt / 60
      EIC$rt <- NULL
      EIC <- plyr::rename(EIC, c("intensity" = "Intensity"))
      setwd(OrigDir)
      
      return(EIC)
      
}

# # EXAMPLE
# # Loading data 
# setwd("C:/Users/Laura/Documents/LCMS metabolomics")
# Files <- read.csv("List of metabolomics files.csv")
# 
# load("Significant MFs linear regression 3 - 20150125.RData")
# # Changing Sig.3$Mode to match what I use elsewhere. Probably should change 
# # this in the metadata and/or regression scripts. 
# Sig.3$Mode <- revalue(Sig.3$Mode, c("ESI+" = "Epos", 
#                                     "ESI-" = "Eneg"))
# 
# # Extracting EICs 
# 
# TopEICs <- eic(Sig.3, Files)
