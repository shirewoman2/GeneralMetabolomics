# General functions for checking mass features


# Housekeeping ---------------------------------------------------------
library(xcms)
library(ggplot2)
library(stringr)
library(XLConnect)
library(tidyr)
library(plyr)
library(dplyr)
library(lubridate)
library(XLConnect)
library(reshape2)

ThemeLaura <- function (base_size = 12, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace% 
            theme(
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(colour = "black"),
                  axis.title.y = element_text(colour = "black", angle=90),
                  panel.background = element_rect(fill="white", color=NA),
                  panel.grid.minor.y = element_line(color = NA),
                  panel.grid.minor.x = element_line(color = NA),
                  panel.grid.major = element_line(colour = NA),
                  plot.background = element_rect(fill="white", colour=NA),
                  panel.border = element_rect(color="black", fill=NA),
                  strip.background = element_rect(color=NA, fill="white"),
                  legend.background = element_rect(color=NA, fill=NA),
                  legend.key = element_rect(color=NA, fill=NA)
            )   
}

# Call up that theme before plotting graphs.
theme_set(ThemeLaura())


MainDir <- "C:/Users/Laura/Documents/SCOR project"

setwd(MainDir)


# eic function ============================================================
# Making a function to export EICs given a data.frame (MF) listing:
# 1. the mass feature (MassFeature)
# 2. m/z (mz)
# 3. RT in minutes (RT)
# 4. ionization mode (Mode)
# 5. and matrix (Matrix) 
# and also given a data.frame (Files) listing:
# 1. the files (File)
# 2. directory (Dir)
# 3. ionization mode (Mode)
# 4. and matrix (Matrix)
# for each file you want to extract the EIC data from.
# Column names must match. This creates a new list object called EIC.


eic <- function(MF, Files, ppm = 15) {
      
      require(plyr)
      OrigDir <- getwd()
      MF$ModeMatrix <- paste(MF$Mode, MF$Matrix, sep = ".")
      Files$ModeMatrix <- paste(Files$Mode, Files$Matrix, sep = ".")
      
      # Only keeping files with the mode and matrix actually present in 
      # the list of mass features of interest.
      Files <- Files[Files$ModeMatrix %in% unique(MF$ModeMatrix), ]
      
      MF <- dlply(MF, c("Mode", "Matrix"))
      EIC <- list()
      Files$Directory <- as.character(Files$Directory)
      
      for (i in 1:nrow(Files)){
            
            setwd(Files$Dir[i])
            RawData <- xcmsRaw(paste0(Files$File[i], ".mzdata.xml"), 
                               profstep = 0, profmethod = "bin")
            
            Dataset <- paste(Files$Mode[i], Files$Matrix[i], sep = ".")
            
            EIC[[i]] <- list()
            
            for (m in 1:nrow(MF[[Dataset]])){
                  
                  MZRange <- cbind(MF[[Dataset]]$mz[m] - 
                                         ppm*MF[[Dataset]]$mz[m]/1e6,
                                   MF[[Dataset]]$mz[m] + 
                                         ppm*MF[[Dataset]]$mz[m]/1e6)
                  
                  RTRange <- t(cbind(MF[[Dataset]]$RT[m]*60 + c(-30, 30)))
                  EICdata <- getEIC(RawData, 
                                    mzrange = MZRange,
                                    rtrange = RTRange)
                  EIC[[i]][[m]] <- data.frame(EICdata@eic$xcmsRaw[[1]]) # names: "rt" and "intensity"
                  EIC[[i]][[m]]$Dataset <- Dataset
                  EIC[[i]][[m]]$File <- Files$File[i]
                  EIC[[i]][[m]]$MassFeature <- MF[[Dataset]]$MassFeature[m]
                  EIC[[i]][[m]]$RT.original <- MF[[Dataset]]$RT[m]
                  EIC[[i]][[m]]$Matrix <- MF[[Dataset]]$Matrix[m]
                  EIC[[i]][[m]]$Mode <- MF[[Dataset]]$Mode[m]
                  
            }
            
            EIC[[i]] <- rbind.fill(EIC[[i]])
            EIC[[i]]$Project <- Files$Project[i]
            
      }
      
      EIC <- rbind.fill(EIC)
      EIC$RT <- EIC$rt / 60
      EIC$rt <- NULL
      EIC <- plyr::rename(EIC, c("intensity" = "Intensity"))
      setwd(OrigDir)
      
      return(EIC)
      
}

# Loading data ----------------------------------------------------
setwd(MainDir)
Files <- read.csv("List of metabolomics files.csv")

load("Significant MFs linear regression 3 - 20150125.RData")
# Changing Sig.3$Mode to match what I use elsewhere. Probably should change 
# this in the metadata and/or regression scripts. 
Sig.3$Mode <- revalue(Sig.3$Mode, c("ESI+" = "Epos", 
                                    "ESI-" = "Eneg"))

# Extracting EICs -------------------------------------------------

TopEICs <- eic(Sig.3, Files)


# MFplot function ----------------------------------------

# Input is 1 mass feature and the data.frame where all the EICs are.

MFplot <- function(MF, EICs) {
      MF.df <- EICs[EICs$MassFeature == MF, ]
      
      Plot <<- ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
            geom_line() + ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF)) +
            xlab("RT (min)") +
            geom_vline(data = MF.df, aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            theme(legend.position = "none") +
            facet_wrap(~ Project, scales = "free")
      Plot
      ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], paste0(MF, ".png")))
}


# Making graphs ---------------------------------------------------

setwd("Fragmentation results and data files/EnegP83 and EposU84 top MF plots")

for (m in 1:nrow(Sig.3)){
      
      MFplot(Sig.3$MassFeature[m], TopEICs)
}

setwd(MainDir)

# Matching to other MFs --------------------------------------------------

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
      require(plyr)
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

# Using MassFeatureMatch to compare MFs in SCOR MDZ EposU84 and in 
# CYP2D6 DEX EposU2.

setwd("G:/Data/Metabolomics/Jessica/CYP2D6/mzData files 20131230/ESI+")
CYP2D6DEX.EposU2 <- read.csv("CYP2D6 EposU 2 peak table - RT above 2 min.csv")
CYP2D6DEX.EposU2 <- CYP2D6DEX.EposU2[, c("MassFeature", "mz", "RTmin")]
CYP2D6DEX.EposU2 <- plyr::rename(CYP2D6DEX.EposU2, c("RTmin" = "RT"))


setwd("C:/Users/Laura/Documents/SCOR project/84 SCOR EposU")
SCORMDZ.EposU84 <- read.csv("SCORMDZEposU84 peak table - RT above 2 min.csv")
SCORMDZ.EposU84 <- SCORMDZ.EposU84[, c("MassFeature", "mz", "RT")]

SCORMDZ.CYP2D6DEX.EposU <- MassFeatureMatch(SCORMDZ.EposU84, CYP2D6DEX.EposU2)
