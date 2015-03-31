# CYP3A positive controls

# This script looks through datasets we already have and searches for 
# compounds thought to be associated with CYP3A activity. See notes in 
# Laura's notebook: SCOR MDZ tab, page "2/27/15 Consolidating positive 
# control info".

# Housekeeping -----------------------------------------------
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(xlsx)

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
# When there were any NA values for RT in Y, output also includes a data.frame 
# "Isomers" that shows which isomers were collapsed into one m/z for searching.
#
# !!! Important: If one of the data.frames includes compounds whose RT you 
# don't know, make that one Y. For somewhat arcane reasons, this works 
# much better.  !!! #


mfmatch <- function(X, Y, PPM = 15, RTRange = 0.2){
      require(plyr)
      require(stringr)
      
      DF.X <- X[, c("MassFeature", "mz", "RT")]
      DF.Y <- Y[, c("MassFeature", "mz", "RT")]
      
      DF.Y$MassFeature <- as.character(DF.Y$MassFeature)
      DF.X$MassFeature <- as.character(DF.X$MassFeature)
      
      DF.X <- arrange(DF.X, MassFeature)
      DF.Y <- arrange(DF.Y, MassFeature)
      
      # Need to condense compounds that are isomers without known RT since we
      # can't determine which would be the best match without RT
      if (anyNA(DF.Y$RT)) {
            
            DF.Y$mz.round <- round(DF.Y$mz, digits = 4)
            
            Isomers <- dlply(DF.Y[is.na(DF.Y$RT), ], "mz.round")
            
            for (r in 1:length(Isomers)){
                  Isomers[[r]]$MassFeature[1] <- str_c(Isomers[[r]]$MassFeature, 
                                                       collapse = " ")
                  Isomers[[r]] <- Isomers[[r]][1, ]
            }
            
            Isomers <- ldply(Isomers)
            Isomers$.id <- NULL
            
            Isomers <<- Isomers
            
            DF.Y <- rbind.fill(Isomers, DF.Y[complete.cases(DF.Y$RT), ])
            DF.Y$mz.round <- NULL
            
      }
      
      
      # Renaming to keep track of which mass feature, m/z, and RT came from 
      # which dataset.
      names(DF.X) <- paste(names(DF.X), "X", sep = ".")
      names(DF.Y) <- paste(names(DF.Y), "Y", sep = ".")
      
      MFmatch <- list()
      Matched.X <- c()
      MFname.X <- as.character(DF.X$MassFeature)
      mz.X <- DF.X$mz
      names(mz.X) <- DF.X$MassFeature.X
      RT.X <- DF.X$RT
      names(RT.X) <- DF.X$MassFeature.X
      
      # Checking each row in DF.X for any matches in DF.Y
      for (i in DF.X$MassFeature.X){
            
            MFmatch[[i]] <- DF.Y[
                  DF.Y$mz.Y > (mz.X[i] - (PPM/1e6*mz.X[i])) 
                  & DF.Y$mz.Y < (mz.X[i] + (PPM/1e6*mz.X[i])) 
                  & DF.Y$RT.Y > RT.X[i] - RTRange 
                  & DF.Y$RT.Y < RT.X[i] + RTRange, ]
            
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
            Matched.X[i] <- as.numeric(nrow(MFmatch[[i]]))            
      }
      
      # Making a new data.frame to hold matched mass features
      Matches <- data.frame(MassFeature.X = DF.X$MassFeature.X,
                            MassFeature.Y = NA, 
                            mz.X = DF.X$mz.X, 
                            mz.Y = NA,
                            RT.X = DF.X$RT.X,
                            RT.Y = NA,
                            NumMatched = Matched.X,
                            ppm = NA, 
                            RTdif = NA)
      
      Matches <- dlply(Matches, "MassFeature.X")
      
      # If there was more than 1 potential match and the RT was listed, take 
      # the one that was closest by m/z and then by RT. If there was no RT 
      # listed, take all matches. 
      # Calculate the difference in m/z in ppm and RT in min.
      for (i in names(Matched.X)[which(Matched.X > 0)]){
            
            for (n in 1:nrow(MFmatch[[i]])){
                  MFmatch[[i]]$ppm[n] <- 
                        abs((mz.X[i] - MFmatch[[i]]$mz.Y[n])/mz.X[i]*1e6)
                  MFmatch[[i]]$RTdif[n] <- 
                        abs(RT.X[i] - MFmatch[[i]]$RT.Y[n])
                  MFmatch[[i]]$MassFeature.X <- i
            }
            MFmatch[[i]] <- arrange(MFmatch[[i]], ppm, RTdif)
            
            Matches[[i]]$MassFeature.Y <- 
                  as.character(MFmatch[[i]]$MassFeature.Y[1])
            Matches[[i]]$mz.Y <- MFmatch[[i]]$mz.Y[1]
            Matches[[i]]$RT.Y <- MFmatch[[i]]$RT.Y[1]
            Matches[[i]]$ppm <- MFmatch[[i]]$ppm[1]
            Matches[[i]]$RTdif <- MFmatch[[i]]$RTdif[1]
            
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
      Files$Dir <- as.character(Files$Dir)
      
      for (i in 1:nrow(Files)){
            
            setwd(Files$Dir[i])
            RawData <- xcmsRaw(Files$File[i], profstep = 0, profmethod = "bin")
            
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


# eicplot function ----------------------------------------
# This function plots EICs of the mass feature of interest, facetted by 
# project. It creates a new ggplot2 object, "Plot", that can be further modified.
# Input: 
#       1. MF = a character string listing the mass feature of interest
#       2. EICs = a data.frame with all the EIC info already extracted; this 
#          should be the output from the function "eic".

eicplot <- function(MF, EICs, Height = 8, Width = 8) {
      
      require(ggplot2)
      
      MF.df <- EICs[EICs$MassFeature == MF, ]
      
      if("Project" %in% names(MF.df)){
            
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
      } else {
            Plot <<- ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
                  geom_line() + ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF)) +
                  xlab("RT (min)") +
                  geom_vline(data = MF.df, aes(xintercept = RT.original),
                             linetype = "dashed", size = 0.5, color = "gray50") +
                  theme(legend.position = "none")
            Plot
            ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], paste0(MF, ".png")),
                   height = Height, width = Width)
      }
}


# Positive controls ---------------------------------------
# Loading data
setwd("C:/Users/Laura/Documents/SCOR project")
PosCont <- read.xlsx("Positive controls.xlsx", sheetName = "unique",
                     colClasses = c(rep("character", 4), 
                                    "numeric",
                                    "character",
                                    rep("numeric", 3),
                                    rep("character", 2),
                                    "numeric"), 
                     colIndex = c(1:12), stringsAsFactors = FALSE)
# names(PosCont)
PosCont <- plyr::rename(PosCont, c("RT..min." = "RT",
                                   "Compound.name" = "MassFeature",
                                   "Matrix.in.which.compound.was.found" =
                                         "Matrix",
                                   "mode.used.to.detect.major.ion" = "Mode",
                                   "major.ion..m.z." = "mz"))

PosCont <- PosCont[str_detect(PosCont$Reason.for.interest, 
                              "CYP3A positive control") &
                         complete.cases(PosCont$RT), ]
PosCont$Mode[is.na(PosCont$Mode)] <- "either"

PosCont.Epos <- PosCont[PosCont$Mode %in% c("Epos", "either"), ]
PosCont.Epos$mz[is.na(PosCont.Epos$mz)] <- 
      PosCont.Epos$M.H.m.z[is.na(PosCont.Epos$mz)]

PosCont.Eneg <- PosCont[PosCont$Mode %in% c("Eneg", "either"), ]
PosCont.Eneg$mz[is.na(PosCont.Eneg$mz)] <- 
      PosCont.Eneg$M.H.m.z[is.na(PosCont.Eneg$mz)]


# Loading IVPO MDZ data -----------------------------------------------
setwd("G:/Data/Metabolomics/Laura/IV MDZ project")
load("IVPO MDZ metadata.RData")

# IVPOEposU3
setwd("G:/Data/Metabolomics/Laura/IV MDZ project/IVPOEposU3")
load("IVPO MDZ main data.RData")

IVPOEposU3 <- Data.log[["IVPOEposU3"]]

MissingFiles <- setdiff(Files$SampleID[Files$Mode == "Epos"],
                        names(IVPOEposU3))
Files[Files$SampleID %in% MissingFiles & Files$Mode == "Epos", ]
# Only missing files are "do not use" and master QC. Good.

Files <- Files[Files$Exists == TRUE & Files$Use == "use", ]


# MDZEnegU2
setwd("G:/Data/Metabolomics/Laura/IV MDZ project/MDZEnegU2")
IVPOEnegU2 <- read.csv("MDZEnegU2a.csv", skip = 1)

# Looking for positive controls ----------------------------------------
# IVPOEposU3
PosCont.IVPO.Epos <- mfmatch(PosCont.Epos, IVPOEposU3, RTRange = 0.5)
PosCont.IVPO.Epos <- PosCont.IVPO.Epos[
      complete.cases(PosCont.IVPO.Epos$MassFeature.PosCont.Epos) &
            complete.cases(PosCont.IVPO.Epos$MassFeature.IVPOEposU3), ]
PosCont.IVPO.Epos <- plyr::rename(PosCont.IVPO.Epos, 
                                  c("MassFeature.IVPOEposU3" = "MassFeature"))
PosCont.IVPO.Epos <- join(PosCont.IVPO.Epos, IVPOEposU3, by = "MassFeature", 
                          type = "left")

PosCont.IVPO.Epos <- gather(PosCont.IVPO.Epos, SampleID, Abundance, 
                            -MassFeature.PosCont.Epos, -MassFeature,
                            -mz.PosCont.Epos, -mz.IVPOEposU3, 
                            -RT.PosCont.Epos, -RT.IVPOEposU3,
                            -NumMatched, -ppm, -RTdif, -mz.mean,
                            -mz, -RT)
PosCont.IVPO.Epos <- join(PosCont.IVPO.Epos, Samples, by = "SampleID", 
                          type = "left")

setwd("C:/Users/Laura/Documents/SCOR project")
write.xlsx(PosCont.IVPO.Epos, file = "IVPO MDZ positive controls.xlsx",
           sheetName = "Epos")

# LEFT OFF HERE !!!! ####
IVPO.Epos <- eic()




# IVPOEnegU2
PosCont.IVPO.Eneg <- mfmatch(PosCont.Eneg, IVPOEnegU2, RTRange = 0.5)
PosCont.IVPO.Eneg <- PosCont.IVPO.Eneg[
      complete.cases(PosCont.IVPO.Eneg$MassFeature.PosCont.Eneg) &
            complete.cases(PosCont.IVPO.Eneg$MassFeature.IVPOEnegU2), ]
PosCont.IVPO.Eneg <- plyr::rename(PosCont.IVPO.Eneg, 
                                  c("MassFeature.IVPOEnegU2" = "MassFeature"))
PosCont.IVPO.Eneg <- join(PosCont.IVPO.Eneg, IVPOEnegU2, by = "MassFeature", 
                          type = "left")

PosCont.IVPO.Eneg <- gather(PosCont.IVPO.Eneg, SampleID, Abundance, 
                            -MassFeature.PosCont.Eneg, -MassFeature,
                            -mz.PosCont.Eneg, -mz.IVPOEnegU2, 
                            -RT.PosCont.Eneg, -RT.IVPOEnegU2,
                            -NumMatched, -ppm, -RTdif, -mz.mean,
                            -mz, -RT)
PosCont.IVPO.Eneg <- join(PosCont.IVPO.Eneg, Samples, by = "SampleID", 
                          type = "left")









# SFN data -------------------------------------------------------
setwd("C:/Users/Laura/Documents/Sulforaphane project")
load("SFN metadata.RData")
load("SCOR MDZ main data.RData")

SFNFiles <- Files
SFNSamples <- Samples
SFN.log <- Data.log

SFNFiles$SampCol <- make.names(paste0(SFNFiles$File, ".mzdata"))
SFNFiles$Dataset[SFNFiles$Mode == "Epos" & SFNFiles$Matrix == "plasma"] <- 
      "SFNEposP8"
SFNFiles$Dataset[SFNFiles$Mode == "Eneg" & SFNFiles$Matrix == "plasma"] <- 
      "SFNEnegP9"
SFNFiles <- dlply(SFNFiles, "Dataset")

load("SFN metadata.RData")
SFNSamples2 <- Samples


# SFNEposP8 ----------------------------------------------------
SFN.EposP <- SFN.log[["SFNEposP8"]]

PosCont.SFN.EposP <- mfmatch(PosCont.Epos, SFN.EposP, RTRange = 0.5)
PosCont.SFN.EposP <- PosCont.SFN.EposP[
      complete.cases(PosCont.SFN.EposP$MassFeature.PosCont.Epos) &
            complete.cases(PosCont.SFN.EposP$MassFeature.SFN.EposP), ]

PosCont.SFN.EposP <- plyr::rename(PosCont.SFN.EposP, 
                                 c("MassFeature.SFN.EposP" = "MassFeature"))
PosCont.SFN.EposP <- join(PosCont.SFN.EposP, SFN.EposP, by = "MassFeature", 
                         type = "left")

PosCont.SFN.EposP.long <- gather(PosCont.SFN.EposP, SampleID, Abundance, 
                                -MassFeature.PosCont.Epos, -MassFeature,
                                -mz.PosCont.Epos, -mz.SFN.EposP, 
                                -RT.PosCont.Epos, -RT.SFN.EposP,
                                -NumMatched, -ppm, -RTdif, -mz.mean,
                                -mz, -RT)
PosCont.SFN.EposP.long <- join(PosCont.SFN.EposP.long, 
                              SFNFiles[["SFNEposP8"]], 
                              by = "SampleID", 
                              type = "left")
PosCont.SFN.EposP.long <- join(PosCont.SFN.EposP.long, SFNSamples2, 
                              type = "left")
PosCont.SFN.EposP.long <- subset(PosCont.SFN.EposP.long, Use == "use" &
                                       SampType == "clinical" &
                                       complete.cases(Effector))
PosCont.SFN.EposP.long$Tx <- paste(PosCont.SFN.EposP.long$Effector,
                                   PosCont.SFN.EposP.long$Day)
PosCont.SFN.EposP.long[which(is.na(PosCont.SFN.EposP.long$Effector)), ]

setwd("C:/Users/Laura/Documents/SCOR project")
ggplot(PosCont.SFN.EposP.long, aes(x = Tx, y = Abundance, fill = Tx)) +
      geom_boxplot() +
      facet_wrap(~ MassFeature.PosCont.Epos, scales = "free")

ggsave("SFN plasma positive control compounds - boxplot.png")

write.xlsx(PosCont.SFN.EposP.long, file = "SFN positive controls.xlsx",
           sheetName = "ESI pos plasma")

# Checking EICs
PosCont.SFN.Epos$Mode <- "Epos"
PosCont.SFN.Epos$Matrix <- "plasma"

SFNFiles$File <- paste0(SFNFiles$File, ".mzdata.xml")
SFNFiles$Dir <- "F:/Sulforaphane/All SFN mzdata files"
SFNFiles.l <- dlply(SFNFiles, c("Mode", "Matrix"), function(x) 
      x[sample(10), ])
names(SFNFiles.l)

SFN.Epos.EIC <- eic(PosCont.SFN.Epos, SFNFiles.l[["Epos.plasma"]])

for (i in 1:nrow(PosCont.SFN.Epos)){
      MF.df <- subset(SFN.Epos.EIC, MassFeature == PosCont.SFN.Epos$MassFeature[i])
      
      ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
            geom_line() + 
            ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF.df$MassFeature[1],
                          PosCont.SFN.Epos$MassFeature.PosCont.Epos[i])) +
            xlab("RT (min)") +
            geom_vline(data = MF.df, aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            theme(legend.position = "none")
      
      ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], MF.df$MassFeature[1],
                   paste0(PosCont.SFN.Epos$MassFeature.PosCont.Epos[i], ".png")),
             height = 8, width = 8)
      
}


write.xlsx(PosCont.SFN.EposP.long, file = "SFN positive controls.xlsx",
           sheetName = "EICs for ESIpos positive controls", append = TRUE)



# SFNEnegP9 ----------------------------------------------------
SFN.EnegP <- SFN.log[["SFNEnegP9"]]

PosCont.SFN.EnegP <- mfmatch(PosCont.Eneg, SFN.EnegP, RTRange = 0.5)
PosCont.SFN.EnegP <- PosCont.SFN.EnegP[
      complete.cases(PosCont.SFN.EnegP$MassFeature.PosCont.Eneg) &
            complete.cases(PosCont.SFN.EnegP$MassFeature.SFN.EnegP), ]

PosCont.SFN.EnegP <- plyr::rename(PosCont.SFN.EnegP, 
                                  c("MassFeature.SFN.EnegP" = "MassFeature"))
PosCont.SFN.EnegP <- join(PosCont.SFN.EnegP, SFN.EnegP, by = "MassFeature", 
                          type = "left")

PosCont.SFN.EnegP.long <- gather(PosCont.SFN.EnegP, SampleID, Abundance, 
                                 -MassFeature.PosCont.Eneg, -MassFeature,
                                 -mz.PosCont.Eneg, -mz.SFN.EnegP, 
                                 -RT.PosCont.Eneg, -RT.SFN.EnegP,
                                 -NumMatched, -ppm, -RTdif, -mz.mean,
                                 -mz, -RT)
PosCont.SFN.EnegP.long <- join(PosCont.SFN.EnegP.long, 
                               SFNFiles[["SFNEnegP9"]], 
                               by = "SampleID", 
                               type = "left")
PosCont.SFN.EnegP.long <- join(PosCont.SFN.EnegP.long, SFNSamples2, 
                               type = "left")
PosCont.SFN.EnegP.long <- subset(PosCont.SFN.EnegP.long, Use == "use" &
                                       SampType == "clinical" &
                                       complete.cases(Effector))
PosCont.SFN.EnegP.long$Tx <- paste(PosCont.SFN.EnegP.long$Effector,
                                   PosCont.SFN.EnegP.long$Day)
PosCont.SFN.EnegP.long[which(is.na(PosCont.SFN.EnegP.long$Effector)), ]

setwd("C:/Users/Laura/Documents/SCOR project")
ggplot(PosCont.SFN.EnegP.long, aes(x = Tx, y = Abundance, fill = Tx)) +
      geom_boxplot() +
      facet_wrap(~ MassFeature.PosCont.Eneg, scales = "free")

ggsave("SFN plasma positive control compounds - boxplot.png")

write.xlsx(PosCont.SFN.EnegP.long, file = "SFN positive controls.xlsx",
           sheetName = "ESI neg plasma", append = TRUE)

# Checking EICs
PosCont.SFN.Eneg$Mode <- "Eneg"
PosCont.SFN.Eneg$Matrix <- "plasma"

SFNFiles$File <- paste0(SFNFiles$File, ".mzdata.xml")
SFNFiles$Dir <- "F:/Sulforaphane/All SFN mzdata files"
SFNFiles.l <- dlply(SFNFiles, c("Mode", "Matrix"), function(x) 
      x[sample(10), ])
names(SFNFiles.l)

SFN.Eneg.EIC <- eic(PosCont.SFN.Eneg, SFNFiles.l[["Eneg.plasma"]])

for (i in 1:nrow(PosCont.SFN.Eneg)){
      MF.df <- subset(SFN.Eneg.EIC, MassFeature == PosCont.SFN.Eneg$MassFeature[i])
      
      ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
            geom_line() + 
            ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF.df$MassFeature[1],
                          PosCont.SFN.Eneg$MassFeature.PosCont.Eneg[i])) +
            xlab("RT (min)") +
            geom_vline(data = MF.df, aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            theme(legend.position = "none")
      
      ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], MF.df$MassFeature[1],
                   paste0(PosCont.SFN.Eneg$MassFeature.PosCont.Eneg[i], ".png")),
             height = 8, width = 8)
      
}


write.xlsx(SFN.Eneg.EIC, file = "SFN positive controls.xlsx",
           sheetName = "EICs for ESIneg positive controls", append = TRUE)



# Saving --------------------------------------
save.image("CYP3A positive controls workspace.RData")
