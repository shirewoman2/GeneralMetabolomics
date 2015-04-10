# ToyEposP dataset
# 4/9/15 LS

# This script generates some toy data to play with.

# Notes to the users of this script ---------------------------------


###       ITEMS YOU MAY WANT TO ADJUST EACH TIME YOU RUN THIS SCRIPT:       ###
# 1. Set the directories for:
#       a. all your mzdata files (RawDataDir),
#       b. where you want your final data to be saved (MainDir)
# 
# 2. Find and replace "ToyEposP" with some single-word name that has meaning to
# you and doesn't start with a number or symbol (those are the rules for naming
# objects in R). Now, your objects will be clearly named and the output files
# will be similarly named.
# 
# 3. Load the specific metadata file(s) you need for this project. It's much 
# easier to have a separate script where you figure out all your file names, 
# which ones you want to use, all your samples, all your subjects, etc. and
# then just save those objects as .RData in that script and then load those
# objects here. That way, you have only one place -- and an obvious, trackable,
# and reproducible place -- where you note which files and samples you're 
# using and why.
# 
# 4. If you have more than 25 or so samples, split up your samples when you run
# xcmsSet() so that you save your progress occasionally. For example, set up 
# the code like this:
#       ToyEposP.xs1 <- xcmsSet(GoodSamples[1:25], method = "centWave",  ppm=15, 
#                                 peakwidth=c(4,12), snthresh = 20, 
#                                 mzCenterFun="apex", prefilter=c(10, 10000),
#                                 integrate = 1, fitgauss= TRUE)
#       save(ToyEposP.xs1, file="ToyEposP xs1.RData")
# 
#       ToyEposP.xs2 <- xcmsSet(GoodSamples[26:50], method = "centWave",  ppm=15, 
#                                 peakwidth=c(4,12), snthresh = 20, 
#                                 mzCenterFun="apex", prefilter=c(10, 10000),
#                                 integrate = 1, fitgauss= TRUE)
#       save(ToyEposP.xs2, file="ToyEposP xs2.RData")
# etc. until you've done peak picking on all of your samples. That way, if the
# computer shuts down or crashes, you've saved your progress and don't have to 
# start over. To load those files into memory, here is the code:
#     load("ToyEposP xs1.RData")
# A note here: This seems to be an obvious place for a for loop rather than 
# having to type out each set of xcmsSet commands for each subset of samples. 
# However, xcsmSet objects are S4 objects (the regular objects you typically use
# in R are S3) and require different treatment. S4 objects CAN do some cool 
# stuff, though, and that's why, when we do ToyEposP <- group(...) and 
# then, subsequently, ToyEposP <- retcor(...) and so forth, we're actually 
# NOT overwriting; we are filling different "slots" in the S4 object that 
# is ToyEposP.
# 
# 5. If you don't have much signal for some samples, you may want to adjust the
# xcmsSet.prefilter setting to be less stringent. If you are getting a lot of
# compounds you think are noise, make it more stringent. 
# 
# 6. I added a requirement for the quality checking step that any mass feature
# we look at must be present in at least 25% of all samples. You can change that
# number to whatever you want, though. See the "Quality control" section. 
# 
# 7. Set the ionization mode (see the section "Housekeeping") as appropriate 
# for your data for to remove the appropriate set of reference ions.
#
# Email me (Laura) at shirewoman2@gmail.com if you have questions.
#



# Housekeeping -------------------------------------------------
# This script uses the package xcms to generate a list of ions that have been
# picked, aligned, RT corrected, realigned, recursively filled, and filtered. 
# It then performs a check on the quality of the data extraction.

# Dataset: ToyEposP

IonizationMode <- "Epos" # Change to "Epos", "Eneg", "Apos", or "Aneg" as 
# appropriate. Capital letter is for ESI or APCI and pos or neg is for 
# positive or negative.

library(plyr)
library(ggplot2)
library(gridExtra)
library(xcms)
library(xlsx)
library(dplyr)
library(seqinr)
library(lubridate)
library(stringr)
library(tidyr)

# Set path of the working directory, i.e. the folder where your files are. 
# Note that you need forward slashes!
MainDir <- "F:/Example data/ToyEposP"
RawDataDir <- "F:/Example data"

# Loading metadata
# These data came from the Progestin project, so loading those metadata.
setwd("F:/Progestin")
load("Progestin metadata.RData")
setwd(RawDataDir)
Files <- subset(Files, File %in% sub(".mzdata.xml", "", 
                                     list.files(pattern = ".mzdata.xml")))

# Samples to use
GoodSamples <- paste0(Files$File, ".mzdata.xml")

# Getting the names of the output data.frames' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(GoodSamples))


# General functions ----------------------------------------------
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

# Peak picking ------------------------------------------
setwd(RawDataDir)

# Making a note of start time for this step
ToyEposP.tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
SNthresh <- 15
Prefilter <- c(10, 5000)

ToyEposP.xs1 <- xcmsSet(GoodSamples, method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(ToyEposP.xs1, file="ToyEposP xs1.RData")


# Tracking how long this step takes
ToyEposP.tpick.final <- Sys.time()
ToyEposP.tpick <- as.numeric(difftime(ToyEposP.tpick.final, 
                                      ToyEposP.tpick.init, units = "mins"))
write.csv(ToyEposP.tpick, "ToyEposP tpick.csv")



# Initial peak alignment -------------------------------------------------------
ToyEposP.tgroup.init <- Sys.time()

# Peak alignment
ToyEposP <- group(ToyEposP.xs1, 
                  method = "density", bw = 4, minfrac = 0,
                  minsamp = 1, mzwid = 0.007, max = 100)
# If you split your samples into multiple xcmsSet steps, change the first bit
# of the above command to:
# group (c(ToyEposP.xs1, ToyEposP.xs2, .. ToyEposP.xsn), method = ...)

ToyEposP.tgroup.final <- Sys.time()
ToyEposP.tgroup <- as.numeric(difftime(ToyEposP.tgroup.final, 
                                       ToyEposP.tgroup.init, units="mins"))
write.csv(ToyEposP.tgroup, "ToyEposP tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(ToyEposP.xs1)


# If you need to reload the xcmsSet objects b/c of a glitch or to remove 
# specific files, uncomment and run the lines below:
# setwd(RawDataDir)
# load("ToyEposP xs1.RData")


# RT correction -------------------------------------------------------
ToyEposP.tretcor.init <- Sys.time() 

# Performing retention time correction
ToyEposP <- retcor(ToyEposP, method = "peakgroups", 
                   missing = 0.1*length(GoodSamples), 
                   extra = 2*length(GoodSamples), smooth = "loess", 
                   family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("ToyEposP RTcor plot.png", width = 4, height = 5, units = "in", res = 300)
plotrt(ToyEposP, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

ToyEposP.tretcor.final <- Sys.time()
ToyEposP.tretcor <- as.numeric(difftime(ToyEposP.tretcor.final, ToyEposP.tretcor.init, units="mins"))
write.csv(ToyEposP.tretcor, "ToyEposP tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
ToyEposP.tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
ToyEposP <- group(ToyEposP, method = "density", 
                  minsamp = 1, minfrac = 0, mzwid = 0.007, 
                  bw = 2, max=100)

# Making a data.frame of the data before the recursive peak filling step
ToyEposP.unfilled <- peakTable(ToyEposP)
ToyEposP.unfilled$MassFeature <- paste("I", round((
      ToyEposP.unfilled$mz),digits=4), "R", round(
            ToyEposP.unfilled$rt/60, digits=2), sep="")

ToyEposP.tgroup2.final <- Sys.time()
ToyEposP.tgroup2 <- as.numeric(difftime(ToyEposP.tgroup2.final, 
                                        ToyEposP.tgroup2.init, units="mins"))
write.csv(ToyEposP.tgroup2, "ToyEposP tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
ToyEposP.tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
ToyEposP <- fillPeaks(ToyEposP)

ToyEposP.tfillPeaks.final <- Sys.time()
ToyEposP.tfillPeaks <- as.numeric(difftime(ToyEposP.tfillPeaks.final, 
                                           ToyEposP.tfillPeaks.init, 
                                           units="mins"))
write.csv(ToyEposP.tfillPeaks, "ToyEposP tfillPeaks.csv")

setwd(MainDir)
save(ToyEposP, file = "ToyEposP xcmsSet object.RData")

# Generate a data.frame with all the peaks ----------------------------------
ToyEposP.tPeakTable.init <- Sys.time()

# Making a data.frame of the recursively filled data
ToyEposP.allpeaks <- peakTable(ToyEposP)

# Setting up a function to count length of columns with 0 or NA in the 
# data before recursive peak filling step. Will filter by this later.
Detected <- function(x) {length(x) - length(which(is.na(x) | x == 0))}
# Counting
ToyEposP.allpeaks$Count <- apply(ToyEposP.unfilled[, SampCol], 
                                 MARGIN = 1, FUN = Detected)

# Making a column with the mass feature name
ToyEposP.allpeaks$MassFeature <- paste("I", round((
      ToyEposP.allpeaks$mz),digits=4), "R", round(
            ToyEposP.allpeaks$rt/60, digits=2), sep="")

# Making a column with the mass feature name as xcms sets it
ToyEposP.allpeaks$groupname <- groupnames(ToyEposP)

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
ToyEposP.allpeaks$RT <- ToyEposP.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "Epos" | IonizationMode == "Apos",
       RefIons <- c(121.0509, 922.0098), # Epos and Apos       
       RefIons <- c(112.9856, 119.0363, 980.015)) # Eneg and Aneg

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(ToyEposP.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 ToyEposP.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
ToyEposP.noref <- ToyEposP.allpeaks[-unlist(RefMFs), ]

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
ToyEposP.after2 <- subset(ToyEposP.noref, ToyEposP.noref$rt 
                          > 120)

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
ToyEposP.filter <- subset(ToyEposP.after2, Count >= 0.25*length(SampCol))
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
setwd(MainDir)
write.csv(ToyEposP.filter, 
          "ToyEposP peak table - mass features must be present in at least 25% of all samples.csv")

save(ToyEposP.after2, ToyEposP.allpeaks, ToyEposP.filter, ToyEposP.noref,
     ToyEposP.unfilled, file = "ToyEposP - all main dataframes.RData")
save(ToyEposP.filter, file = "ToyEposP - filtered dataframe only.RData")

setwd(RawDataDir)
ToyEposP.tPeakTable.final <- Sys.time()
ToyEposP.tPeakTable <- as.numeric(difftime(ToyEposP.tPeakTable.final, 
                                           ToyEposP.tPeakTable.init, 
                                           units="mins"))
write.csv(ToyEposP.tPeakTable, "ToyEposP tPeakTable.csv")

# Quality control ----------------------------------------
ToyEposP.tQC.init <- Sys.time()

setwd(MainDir)

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
set.seed(253)

MFs <- ToyEposP.filter$groupname[runif(30, min = 1, 
                                       max = length(
                                             ToyEposP.filter$groupname))]
RandSamp <- round(runif(10, min = 1, max = length(Samples)))
write.csv(MFs, paste(Sys.Date(), 
                     "ToyEposP randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "ToyEposP randomly selected samples.csv"))

LastSamp <- GoodSamples[RandSamp[length(RandSamp)]]

EIC.uncorrected <- list()
EIC.corrected <- list()

# This next step will take some time to process, so don't expect instant results. 
for (i in 1:length(MFs)){
      EIC.uncorrected[[i]] <- getEIC(ToyEposP, rt="raw", 
                                     groupidx=MFs[i], sampleidx=RandSamp)
      EIC.corrected[[i]] <- getEIC(ToyEposP, rt="corrected", 
                                   groupidx=MFs[i], sampleidx=RandSamp)  
}

ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp)-1), "red")

setwd(RawDataDir)
xset.raw <- xcmsRaw(LastSamp, profstep=0.01, profmethod="bin")

setwd(MainDir)
pdf("ToyEposP quality check.pdf", 8.5,11)

# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the 1st sample for that compound with a 
# dashed horizontal line where the calculated m/z is.

par(mfrow=c(4,3), mar=c(3,3,3,0.5))
for(i in 1:length(MFs)){
      
      plot(EIC.uncorrected[[i]], ToyEposP, groupidx=1, rtrange=60, 
           col=MyColors, main=MFs[i])
      mtext(paste(i, ToyEposP.filter$MassFeature[
            ToyEposP.filter$groupname == MFs[i]]), side=3, line=-1, 
            adj=0, padj=0, cex=0.8)
      plot(EIC.corrected[[i]], ToyEposP, groupidx=1, rtrange=60, 
           col=MyColors)
      
      RT <- ToyEposP.filter$rt[ToyEposP.filter$groupname == MFs[i]]
      RTRange <- c(RT-30, RT+30)
      
      mz <- ToyEposP.filter$mz[ToyEposP.filter$groupname == MFs[i]]
      mzRange <- c(mz-0.02, mz+0.02)
      mzRange.poly.low <- mz- mz*(0.5*PPM)/1e6
      mzRange.poly.up <- mz*(0.5*PPM)/1e6 + mz
      
      plotRaw(xset.raw, mzrange=mzRange, rtrange=RTRange, log=FALSE)
      abline(h=mz, lty=2, col="gray35")
      mtext(paste("abund =", round(ToyEposP.filter[
            ToyEposP.filter$groupname == MFs[i], 
            SampCol[RandSamp[length(RandSamp)]]], digits=0)), 
            side=3, line=-1, adj=0, padj=0, cex=0.8)
      polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
              c(mzRange.poly.up, mzRange.poly.up, 
                mzRange.poly.low, mzRange.poly.low), 
              col=col2alpha("blue", alpha=0.1), border=NA)
      abline(v=RT, lty=2, col="gray35")
      
}

dev.off()

save(EIC.corrected, EIC.uncorrected, xset.raw, 
     file = "ToyEposP QC data.RData")

setwd(RawDataDir)
ToyEposP.tQC.final <- Sys.time()
ToyEposP.tQC <- as.numeric(difftime(ToyEposP.tQC.final, 
                                    ToyEposP.tQC.init, units = "mins"))
write.csv(ToyEposP.tQC, "ToyEposP tQC.csv")



# Calculating processing times for each step --------------------------------
setwd(RawDataDir)

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes. If any of the objects required to 
# construct "Times" don't exist, read all of them from file.
if(any(exists(c("ToyEposP.tpick", "ToyEposP.tgroup",
                "ToyEposP.tretcor", "ToyEposP.tretcor",
                "ToyEposP.tgroup2", "ToyEposP.tfillPeaks",
                "ToyEposP.tPeakTable", "ToyEposP.tQC")) == FALSE)) {
      ToyEposP.tpick <- read.csv("ToyEposP tpick.csv")[1,2]
      ToyEposP.tgroup <- read.csv("ToyEposP tgroup.csv")[1,2]
      ToyEposP.tretcor <- read.csv("ToyEposP tretcor.csv")[1,2]
      ToyEposP.tgroup2 <- read.csv("ToyEposP tgroup2.csv")[1,2]
      ToyEposP.tfillPeaks <- read.csv("ToyEposP tfillPeaks.csv")[1,2]
      ToyEposP.tPeakTable <- read.csv("ToyEposP tPeakTable.csv")[1,2]
      ToyEposP.tQC <- read.csv("ToyEposP tQC.csv")[1,2]
}

Times <- data.frame(rbind(ToyEposP.tpick, ToyEposP.tgroup, 
                          ToyEposP.tretcor, ToyEposP.tgroup2, 
                          ToyEposP.tfillPeaks, 
                          ToyEposP.tPeakTable, ToyEposP.tQC))
Times$Step <- c("pick peaks", # tpick
                "group 1", # tgroup
                "retcor", # tretcor
                "group 2", # tgroup2
                "fill peaks", # tfillPeaks
                "peak table", # tPeakTable
                "QC") # tQC
names(Times) <- c("Duration", "Step")
row.names(Times) <- 1:nrow(Times)

# Calculate the total time in hours. 
TotalTime <- c("Total time in hours", format(sum(Times$Duration)/60, digits=2))
Times$Duration <- round(Times$Duration, digits=2)
Times <- Times[, c("Step", "Duration")]

setwd(MainDir)
write.csv(rbind(Times, TotalTime), "ToyEposP processing times.csv", row.names=F)
write.csv(data.frame(NumSamp = length(GoodSamples),
                     PeakPickPerSamp = ToyEposP.tpick/length(GoodSamples)),
          "ToyEposP peak picking time per sample.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for ToyEposP") +
      guides(fill=FALSE)
ggsave("ToyEposP processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for ToyEposP") +
      guides(fill=FALSE)
ggsave("ToyEposP processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                      "Count" = c(nrow(ToyEposP.unfilled),
                                  nrow(ToyEposP.noref),
                                  nrow(ToyEposP.after2),
                                  nrow(ToyEposP.filter)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "ToyEposP counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("ToyEposP bar chart of numbers of mass features at each step.png")


# Processing the data ----------------------------------------------
# Selecting only the columns of interest in the data and renaming them sensibly
ToyEposP.filter <- ToyEposP.filter[, c("MassFeature", "mz", "RT", 
                                       Files$SampCol)]
names(ToyEposP.filter) <- c("MassFeature", "mz", "RT", 
                            Files$SampleID)

# Starting from data that were already filtered by frequency to retain ions 
# detected initially in >= 25% of all samples. 
# Adding 1 to all sample columns
ToyEposP.filter[ , Files$SampleID] <- ToyEposP.filter[ , Files$SampleID] + 1

# Calculate the TCC sum.
TCCsum <- apply(ToyEposP.filter[ , Files$SampleID], 2, sum, na.rm=TRUE)


# Normalize each sample by the TCC sum. 
ToyEposP.TCCnorm <- ToyEposP.filter[ , Files$SampleID]
for (s in 1:length(Files$SampleID)){
      ToyEposP.TCCnorm[, Files$SampleID[s]] <- 
            (ToyEposP.filter[, Files$SampleID[s]]/TCCsum[s])*1e6
      rm(s)
}  


# log10 transforming
ToyEposP.log <- data.frame(MassFeature = ToyEposP.filter$MassFeature, 
                           mz = ToyEposP.filter$mz, 
                           RT = ToyEposP.filter$RT,
                           log10(ToyEposP.TCCnorm))


# Saving the files
setwd(MainDir)
FileNote <- paste("Dataset: ToyEposP --",
                  "This file was generated on", Sys.Date(), 
                  "using the script 'Example data.R.'")
my.write(ToyEposP.log, "ToyEposP - preprocessed.csv", 
         header=FileNote, row.names=F)


# Checking distributions. 
ToyEposP.melt <- gather(ToyEposP.log, Sample, Abundance, -MassFeature, 
                        -mz, -RT)


ggplot(ToyEposP.melt, aes(x=Abundance)) + 
      geom_density() +
      xlab("log10(Abundance)") + ggtitle("ToyEposP")
ggsave("Kernel density plot of ToyEposP abundance data.png")

rm(ToyEposP.melt, TCCsum)

save(ToyEposP.log, file="ToyEposP main data.RData")


# Applying various functions ------------------------------------------
# Files data.frame will require mode, matrix, and directory for some functions.
# Adding those.
Files$Mode <- "Epos"
Files$Matrix <- "plasma"
Files$Directory <- RawDataDir

# Files will also need the full name, i.e. XXXX.mzdata.xml. Adding the file extension.
Files$File <- paste0(Files$File, ".mzdata.xml")

# CAMERA and allions functions -------------------------------------------------------
setwd("F:/General LCMS scripts")
source("camera function.R")
source("allions function.R")
source("eic function.R")
setwd(MainDir)

ToyEposP.cam <- camera(ToyEposP, "Epos")
save(ToyEposP.cam, file = "ToyEposP camera object.RData")

# Randomly selecting some ions with possible adducts or isotopes to see what
# their EICs look like. 
ToyEposP.otherions <- ToyEposP.cam[["ToyEposP.otherions"]]
set.seed(206)
AdductsOfInterest <- subset(ToyEposP.otherions, complete.cases(NeutralMassOfM) &
                                  MassFeature %in% ToyEposP.filter$MassFeature)
IsoOfInterest <- subset(ToyEposP.otherions, complete.cases(IsoGroup) &
                              MassFeature %in% ToyEposP.filter$MassFeature)

IonsOfInterest <- rbind(AdductsOfInterest[
      runif(5, min = 1, max = nrow(AdductsOfInterest)), ],
      IsoOfInterest[runif(5, min = 1, max = nrow(IsoOfInterest)), ])
IonsOfInterest$Matrix <- "plasma"
IonsOfInterest$Mode <- "Epos"

EICs <- list()

for (i in 1:nrow(IonsOfInterest)){
      DF <- IonsOfInterest[i, ]
      EICs[[i]] <- allions(DF, Files, ToyEposP.cam)
      ggplot(EICs[[i]], aes(x = RT, y = Intensity, color = File)) +
            geom_line() + 
            xlab("RT (min)") + ggtitle(IonsOfInterest$MassFeature[i]) +
            geom_vline(data = EICs[[i]], aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            facet_wrap(~ MassFeature.otherion) +
            theme(legend.position = "none")
      ggsave(paste(IonsOfInterest$MassFeature[i], "EICs of CAMERA-detected ions.png"))
}


# mfmatch function -------------------------------------------------
# Matching CYP3A positive controls to these data.
setwd("F:/General LCMS scripts")
source("mfmatch function.R")
setwd(MainDir)

setwd("C:/Users/Laura/Documents/SCOR project")
CYP3A <- read.xlsx("Positive controls.xlsx", sheetName = "unique")
CYP3A <- plyr::rename(CYP3A, c("mode.used.to.detect.major.ion" = "Mode",
                               "RT..min." = "RT",
                               "major.ion..m.z." = "mz",
                               "Compound.name" = "MassFeature"))
CYP3A$RT <- as.numeric(as.character(CYP3A$RT))
CYP3A$mz <- as.numeric(as.character(CYP3A$mz))
CYP3A <- subset(CYP3A, Mode == "Epos" & complete.cases(RT))

Match <- mfmatch(CYP3A, ToyEposP.after2, PPM = 20, RTRange = 0.5)
Match <- subset(Match, complete.cases(MassFeature.CYP3A))
# No hits. Bummer.

# # Positive control to make sure the function is matching: Added a compound
# # that was already in the ToyEposP data.
# CYP3A <- rbind(CYP3A[, c("MassFeature", "mz", "RT")],
#                data.frame(MassFeature = "I251.1802R8.53",
#                           mz = 251.18015,
#                           RT = 8.529446))
# Match <- mfmatch(CYP3A, ToyEposP.filter, PPM = 20, RTRange = 0.5)
# Match <- subset(Match, complete.cases(MassFeature.CYP3A))
# # Good. That's working. 


# eic and eicplot functions ------------------------------------
setwd("F:/General LCMS scripts")
source("eic function.R")
source("eicplot function.R")
setwd(MainDir)

# Picking an ion to look at. Making it one with a really good signal.
IonSignal <- apply(ToyEposP.filter[, Files$SampleID], 1, mean)
hist(IonSignal)
AbundIons <- ToyEposP.filter[which(rank(-IonSignal) %in% 1:10), 
                            c("MassFeature", "mz", "RT")]
AbundIons$Mode <- "Epos"
AbundIons$Matrix <- "plasma"

AbundEICs <- eic(AbundIons, Files)

for (i in 1:nrow(AbundIons)){
      eicplot(AbundIons$MassFeature[i], AbundEICs)
}


# specific ions function --------------------------------------------------
setwd("F:/General LCMS scripts")
source("specifions function.R")
setwd(MainDir)

AbundEICs2 <- specifions(AbundIons[1, ], Files, c("M+Na", 
                                                  "M+1", 
                                                  "M+2",
                                                  "M-H2O"))
AbundEICs2$Group <- paste(AbundEICs2$File, AbundEICs2$MassFeature)

ggplot(AbundEICs2, aes(x = RT, y = Intensity, color = File, 
                       group = Group)) +
      geom_line() + ggtitle(AbundEICs2$MassFeature[1]) +
      facet_wrap(~ MassFeature.ion)

# This doesn't look quite how it should. Having some trouble with grouping.


# worklist function ------------------------------------------------
# This example of how the worklist function works is being made with data 
# from the Progestin project. 
setwd("F:/General LCMS scripts")
source("Generate worklist function.R")
setwd(MainDir)

# # Load this if it's not already in your workspace:
# setwd("F:/Progestin")
# load("Progestin metadata.RData")

# NB: If you need to prep over several days (prepping 75 samples is a long day,
# for example), it's a good idea to randomize the treatment groups so that you 
# don't prep all the samples from one treatment group on one day and all
# the samples from another a different day. That will assuredly create
# artifacts. This example is for the entire set of samples in this project,
# which I actually prepped over the course of 4 days.

Run <- subset(Samples, Subject %in% unique(Samples$Subject)[1:15])
Run$FileLabel <- gsub("\\.", " ", Run$SampleID) # Need a label

ToyWorklist <- worklist(Run, Date = 20150410,
                        Project = "Toy data", Matrix = "plasma",
                        FilePath = "D:\\MassHunter\\Data\\Laura\\Toy data\\")





