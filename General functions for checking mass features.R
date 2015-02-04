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


MainDir <- "C:/Users/Laura/Documents/LCMS metabolomics"

setwd(MainDir)


# eic function ============================================================
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
# setwd(MainDir)
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

# eicplot function ----------------------------------------
# This function plots EICs of the mass feature of interest, facetted by 
# project. It creates a new ggplot2 object, "Plot", that can be further modified.
# Input: 
#       1. MF = a character string listing the mass feature of interest
#       2. EICs = a data.frame with all the EIC info already extracted; this 
#          should be the output from the function "eic".

eicplot <- function(MF, EICs, Height = 8, Width = 8) {
      MF.df <- EICs[EICs$MassFeature == MF, ]
      
      Plot <<- ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
            geom_line() + ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF)) +
            xlab("RT (min)") +
            geom_vline(data = MF.df, aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            theme(legend.position = "none") +
            facet_wrap(~ Project, scales = "free")
      Plot
      ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], paste0(MF, ".png")),
             height = Height, width = Width)
}


# # EXAMPLE
# # Making graphs 
# setwd("Fragmentation results and data files/EnegP83 and EposU84 top MF plots")
# 
# for (m in 1:nrow(Sig.3)){
#       
#       eicplot(Sig.3$MassFeature[m], TopEICs)
# }
# 
# setwd(MainDir)

# Matching to other MFs --------------------------------------------------

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

# Absolute levels of MFs ----------------------------------------------------
# The point of running this function is to find which samples would be best to 
# fragment for a given MF.
# Input:
#    1. MassFeature = the name of the mass feature you're interested in 
#       fragmenting AS IT APPEARS in DF.
#    2. DF = a data.frame from data that have been processed with xcms but NOT
#       preprocessed. Must contain the following columns:
#             a. "MassFeature" = unique identifier for all mass features, 
#                including the mass feature listed in #1. 
#             b. "RT" = the retention time in minutes, 
#             c. "mz" = the m/z for that mass feature, and 
#             d. sample columns with peak areas. 
#    3. SampNames = names of the columns that contain absolute abundances
#    4. n (optional) = number of top mass features you want. Defaults to 3.
# Output: A data.frame with 2 columns: 1. the top n samples by abundance and 
# 2. the abundance values.

abslevelsMF <- function(MassFeature, DF, SampNames, n = 3){
      require(tidyr)
      
      MF <- DF[DF$MassFeature == MassFeature, SampNames]
      MF <- gather(MF, Sample, Abundance)
      MF <- arrange(MF, -Abundance)
      
      return(MF[1:n, ])
}


# # EXAMPLE
# setwd("C:/Users/Laura/Documents/SCOR project/84 SCOR EposU")
# EposU84.unproc <- read.csv("SCORMDZEposU84 peak table - RT above 2 min.csv")
# 
# SampNames <- names(EposU84.unproc)[9:96]
# 
# BestMF <- abslevelsMF("I114.0647R5.05", EposU84.unproc, SampNames, 5)


# Relative levels of MFs ------------------------------------------------
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


# Processing xcms data using CAMERA --------------------------
# This function processes an xcmsSet object using CAMERA to look for 
# other ions that could arise from the same mass feature. Input is an
# xcmsSet object (xset) and the ionization mode as "positive" or "negative".
# Output is a named list containing:
#       1. the CAMERA object (suffix is .xsa)
#       2. a data.frame of the output from CAMERA (.annot)
#       3. a data.frame of other potential ions arising from the same compound as
#          the original input mass feature (.otherions)

camera <- function(xset, Mode, PPM = 15, PVal = 0.0001) {
      require(xcms)
      require(CAMERA)
      require(plyr)
      require(stringr)
      
      # Create a CAMERA object
      xsa <- xsAnnotate(xset)
      
      # Find the compound groups based on RT
      xsa <- groupFWHM(xsa, perfwhm = 1.5)
      
      # Check that peaks in the same group correlate well enough to continue to be 
      # included. Generate pseudospectra.
      xsa <- groupCorr(xsa, cor_eic_th = 0.85, 
                       pval=PVal, calcIso=TRUE, calcCiS=TRUE)
      
      # Find the isotopes within each pseudospectrum.
      xsa <- findIsotopes(xsa, maxcharge = 3, maxiso = 8, 
                          ppm = PPM, mzabs = 0.015, intval = "maxo", 
                          minfrac = 0.25)
      
      # Find the adducts within each pseudospectrum.
      xsa <- findAdducts(xsa, ppm = PPM, mzabs = 0.015, 
                         multiplier = 3, polarity = Mode, 
                         rules = NULL, max_peaks = 100)
      
      # Make a table of the annotated peaks.
      xset.annot <- getPeaklist(xsa)
      
      # Once again, remove reference ions from the peak table.
      ifelse(Mode == "positive",
             RefIons <- c(121.0509, 922.0098), # ESI+       
             RefIons <- c(112.9856, 119.0363, 980.015)) # ESI-
      
      MF.refs <- list()
      
      for (i in 1:length(RefIons)){
            MF.refs[[i]] <- which(xset.annot$mz < 
                                        (RefIons[i] + PPM/1e6*RefIons[i]) &
                                        xset.annot$mz > 
                                        RefIons[i] - PPM/1e6*RefIons[i])
      }
      
      xset.annot <- xset.annot[-unlist(MF.refs), ]
      
      # Making a column with the mass feature name.
      xset.annot$MassFeature <- paste0("I", round((xset.annot$mz),
                                                  digits = 4), "R", 
                                       round(xset.annot$rt/60, 
                                             digits = 2))
      
      # Making a column with the RT in minutes. Note that this is different
      # from the column "rt", which is the RT in seconds. 
      xset.annot$RT <- round(xset.annot$rt/60, digits = 2)
      
      # Retaining only columns of interest in a logical order
      xset.annot <- xset.annot[, c("MassFeature", "mz", "RT", "isotopes",
                                   "adduct", "pcgroup")]
      
      # Using regex to extract meaningful info about isotopes and adducts.
      
      # Isotopes
      # Make the column "isotopes" be character data instead of the default, factor.
      # Replace empty cells with "NA".
      xset.annot$isotopes <- as.character(xset.annot$isotopes)
      xset.annot$isotopes[xset.annot$isotopes == ""] <- NA
      
      
      # Take a subset of the data that is only the mass features that appear
      # to have isotopes. 
      Iso <- xset.annot[complete.cases(xset.annot$isotopes), 
                        c("MassFeature", "mz", "RT", "isotopes", 
                          "pcgroup")]
      
      
      # Split the isotopes column into three pieces:
      # 1. xset.annot$IsoGroup (numeric isotope group)
      # 2. xset.annot$IonType (the isotope of that particular ion s/a "M" or "M+1")
      # 3. xset.annot$Charge (the charge of the isotope)
      IsoString <- str_split(Iso$isotopes, "\\]")
      IsoGroup <- sapply(IsoString, function(x) x[[1]])
      IsoGroup <- gsub("\\[", "", IsoGroup)
      
      IonType <- sapply(IsoString, function(x) x[[2]])
      IonType <- gsub("\\[", "", IonType)
      
      Charge <- sapply(IsoString, function(x) x[[3]])
      
      IsoList <- data.frame(MassFeature = Iso$MassFeature,
                            mz = Iso$mz,
                            RT = Iso$RT, 
                            pcgroup = Iso$pcgroup,
                            IsoGroup = IsoGroup, 
                            IonType = IonType, 
                            Charge = Charge)
      IsoList$IonType <- as.character(IsoList$IonType)
      
      IsoList$Othermz <- NA
      for (i in 1:nrow(IsoList)){
            
            if (IsoList$IonType[i] == "M") {
                  IsoList$Othermz[i] <- IsoList$mz[i]
            } else {
                  n <- as.numeric(str_sub(IsoList$IonType[i], 3, 3))
                  
                  IsoList$Othermz[i] <- IsoList$mz[i] - n * 1.00866
                  
            }
            
      }
      
      # Adducts
      # Make the column "adduct" be character data instead of the default, factor.
      # Replace empty cells with "NA".
      xset.annot$adduct <- as.character(xset.annot$adduct)
      xset.annot$adduct[xset.annot$adduct == ""] <- NA
      
      # Take a subset of the data that is only the mass features that appear
      # to have adduct. 
      Adduct <- xset.annot[complete.cases(xset.annot$adduct), 
                           c("MassFeature", "mz", "RT", "adduct", 
                             "pcgroup")]
      
      # Split the adduct column into one list for every possible adduct with 
      # each list having 4 pieces:
      # 1. IonType = type of adduct, eg. M+Cl
      # 2. mz.adduct = the m/z of that particular adduct
      # 3. MassFeature 
      # 4. pcgroup
      
      AdSplit <- str_split(Adduct$adduct, "\\[")
      AdList <- list()
      
      for (i in 1:nrow(Adduct)){
            
            IonType <- c()
            Charge <- c()
            Othermz <- c()
            
            for (m in 2:length(AdSplit[[i]])){
                  
                  IonType[m-1] <- unlist(str_extract(AdSplit[[i]][m], ".*\\]"))
                  IonType[m-1] <- gsub("\\]", "", IonType[m-1])
                  
                  Charge[m-1] <- unlist(str_extract(AdSplit[[i]][m], "\\].{1,2}"))
                  Charge[m-1] <- str_trim(gsub("\\]", "", Charge[m-1]))
                  
                  Othermz[m-1] <- str_trim(str_extract(AdSplit[[i]][m], " .*"))
            }
            
            AdList[[i]] <- data.frame(MassFeature = Adduct$MassFeature[i],
                                      mz = Adduct$mz[i],
                                      RT = Adduct$RT[i],
                                      pcgroup = Adduct$pcgroup[i],
                                      IonType = IonType,
                                      Charge = Charge,
                                      Othermz = as.numeric(Othermz))
      }
      
      AdList <- rbind.fill(AdList)
      
      
      OtherIons <- rbind.fill(AdList, IsoList)
      
      IonList <- list(xsa, xset.annot, OtherIons)
      names(IonList) <- c(
            paste0(as.character(deparse(substitute(xset))),".xsa"),
            paste0(as.character(deparse(substitute(xset))),".annot"),
            paste0(as.character(deparse(substitute(xset))),".otherions"))
      
      return(IonList)
      
}


# # EXAMPLE
# setwd("F:/SCOR/20120202 SCOR plasma and urine ESI neg")
# load("SCORMDZEnegP83 workspace.RData")
# 
# EnegP83.otherions <- camera(SCORMDZEnegP83, "negative")
# 
# SCORMDZEnegP83.annot <- EnegP83.otherions[["SCORMDZEnegP83.annot"]]
# OtherIons <- EnegP83.otherions[["SCORMDZEnegP83.otherions"]]


# EICs of all ions --------------------------------------------------
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


# allionplot function ----------------------------------------
# This function plots EICs of the mass feature of interest, facetted by 
# other ions. It creates a new ggplot2 object, "Plot.allion", that can be
# further modified.
# Input: the data.frame output from the function "allion".

allionplot <- function(allion.df, Height = 8, Width = 8, 
                       FileSuffix = "- all ions.png") {
      
      require(ggplot2)
      
      Plot.allion <<- ggplot(allion.df, aes(x = RT, y = Intensity, 
                                            color = File)) +
            geom_line() + ggtitle(paste(allion.df$Mode[1], 
                                        allion.df$Matrix[1], 
                                        allion.df$MassFeature[1])) +
            xlab("RT (min)") +
            geom_vline(data = allion.df, aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            theme(legend.position = "none") +
            facet_grid(MassFeature.ion ~ Project, scales = "free")
      Plot.allion
      ggsave(paste(allion.df$Mode[1], allion.df$Matrix[1], 
                   allion.df$MassFeature[1], FileSuffix),
             height = Height, width = Width)
}



# Specific ions ----------------------------------------
# Rather than have CAMERA search for likely related ions, look for specific 
# ones. Input is:
# 1. MF.df, a 1-row data.frame with:
#       a. MassFeature (character)
#       b. mz (numeric)
#       c. RT (numeric)
#       d. Mode (character or factor)
#       e. Matrix (character or factor)
# 2. Files = a data.frame with the following columns:
#       a. the files (File)
#       b. directory (Dir)
#       c. ionization mode (Mode)
#       d. and matrix (Matrix)
# 3. Ions, a list of which ions you want or the specific m/z to use. Specific 
#    ions to call by name are: "M+Na", "M+1", "M+2", "M+3", "M+Cl", "M-H2O".

specifions <- function(MF.df, Files, Ions) {
      
      PossibleIons <- data.frame(Ion = c("M+Na", 
                                            "M+1", 
                                            "M+2", 
                                            "M+3",
                                            "M+Cl", 
                                            "M-H2O"),
                                 MassAdd = c(22.98977- 1.0073, # positive mode: If molecule gains Na, wouldn't also gain H
                                             1.008665,
                                             2*1.008665,
                                             3*1.008665,
                                             34.96885+1.0073, # negative mode: If molecule gains Cl, wouldn't also lose H
                                             -18.01056))

      MF.df$MassFeature.ion <- "originally detected ion"
      
      for (i in 1:length(Ions)) {
            MF.df <- rbind(MF.df,
                           MF.df[1, ])
            
            if (Ions[[i]] %in% PossibleIons$Ion){
                  MF.df$mz[i+1] <- MF.df$mz[1] + 
                        PossibleIons$MassAdd[which(PossibleIons$Ion == 
                                                         Ions[[i]])]
            } else {
                  MF.df$mz[i+1] <- MF.df$mz[1] + Ions[[i]]
            }
            
            MF.df$MassFeature.ion[i+1] <- as.character(Ions[[i]])
            
      }
      
      EICs <- eic(MF.df, Files)
      EICs <- join(EICs, MF.df[, c("MassFeature", "mz", "MassFeature.ion")], 
                   by = c("MassFeature", "MassFeature.ion"))
      return(EICs)
      
}

