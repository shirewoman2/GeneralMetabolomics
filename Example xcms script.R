# ProgEnegU2 xcms
# 3/31/15 LS

# Notes to the users of this script ---------------------------------


###       ITEMS YOU MAY WANT TO ADJUST EACH TIME YOU RUN THIS SCRIPT:       ###
# 1. Set the directories for:
#       a. all your mzdata files (RawDataDir),
#       b. where you want your final data to be saved (MainDir)
# 
# 2. Find and replace "ProgEnegU2" with some single-word name that has meaning to
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
#       ProgEnegU2.xs1 <- xcmsSet(GoodSamples[1:25], method = "centWave",  ppm=15, 
#                                 peakwidth=c(4,12), snthresh = 20, 
#                                 mzCenterFun="apex", prefilter=c(10, 10000),
#                                 integrate = 1, fitgauss= TRUE)
#       save(ProgEnegU2.xs1, file="ProgEnegU2 xs1.RData")
# 
#       ProgEnegU2.xs2 <- xcmsSet(GoodSamples[26:50], method = "centWave",  ppm=15, 
#                                 peakwidth=c(4,12), snthresh = 20, 
#                                 mzCenterFun="apex", prefilter=c(10, 10000),
#                                 integrate = 1, fitgauss= TRUE)
#       save(ProgEnegU2.xs2, file="ProgEnegU2 xs2.RData")
# etc. until you've done peak picking on all of your samples. That way, if the
# computer shuts down or crashes, you've saved your progress and don't have to 
# start over. To load those files into memory, here is the code:
#     load("ProgEnegU2 xs1.RData")
# A note here: This seems to be an obvious place for a for loop rather than 
# having to type out each set of xcmsSet commands for each subset of samples. 
# However, xcsmSet objects are S4 objects (the regular objects you typically use
# in R are S3) and require different treatment. S4 objects CAN do some cool 
# stuff, though, and that's why, when we do ProgEnegU2 <- group(...) and 
# then, subsequently, ProgEnegU2 <- retcor(...) and so forth, we're actually 
# NOT overwriting; we are filling different "slots" in the S4 object that 
# is ProgEnegU2.
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

# Dataset: ProgEnegU2

IonizationMode <- "Eneg" # Change to "Epos", "Eneg", "Apos", or "Aneg" as 
# appropriate. Capital letter is for ESI or APCI and pos or neg is for 
# positive or negative.

library(plyr)
library(ggplot2)
library(gridExtra)
library(xcms)
library(XLConnect)
library(dplyr)
library(seqinr)
library(lubridate)

# Set path of the working directory, i.e. the folder where your files are. 
# Note that you need forward slashes!
MainDir <- "F:/Progestin/Progestin raw data/ProgEnegU2"
RawDataDir <- "F:/Progestin/Progestin raw data"

# Loading metadata
setwd(RawDataDir)
load("Progestin metadata.RData")

# Samples to use
GoodSamples <- paste0(Files$File[Files$Mode == IonizationMode & 
                                       Files$Matrix == "plasma" &
                                       Files$Use == "use" &
                                       Files$SampType == "clinical"],
                      ".mzdata.xml")

# Getting the names of the output data.frames' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(GoodSamples))

# Getting the name of the column that will contain the number of times
# an ion was detected.
CountCol <- make.names(basename(RawDataDir))



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


# Peak picking ------------------------------------------
setwd(RawDataDir)

# Making a note of start time for this step
ProgEnegU2.tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
SNthresh <- 15
Prefilter <- c(10, 5000)

ProgEnegU2.xs1 <- xcmsSet(GoodSamples[1:25], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs1, file="ProgEnegU2 xs1.RData")

ProgEnegU2.xs2 <- xcmsSet(GoodSamples[26:50], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs2, file="ProgEnegU2 xs2.RData")

ProgEnegU2.xs3 <- xcmsSet(GoodSamples[51:75], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs3, file="ProgEnegU2 xs3.RData")

ProgEnegU2.xs4 <- xcmsSet(GoodSamples[76:100], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs4, file="ProgEnegU2 xs4.RData")

ProgEnegU2.xs5 <- xcmsSet(GoodSamples[101:125], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs5, file="ProgEnegU2 xs5.RData")

ProgEnegU2.xs6 <- xcmsSet(GoodSamples[126:150], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs6, file="ProgEnegU2 xs6.RData")

ProgEnegU2.xs7 <- xcmsSet(GoodSamples[151:175], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs7, file="ProgEnegU2 xs7.RData")

ProgEnegU2.xs8 <- xcmsSet(GoodSamples[176:200], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs8, file="ProgEnegU2 xs8.RData")

ProgEnegU2.xs9 <- xcmsSet(GoodSamples[201:225], method = "centWave",  ppm=PPM, 
                          peakwidth=c(4,12), 
                          snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                Prefilter,
                          integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs9, file="ProgEnegU2 xs9.RData")

ProgEnegU2.xs10 <- xcmsSet(GoodSamples[226:250], method = "centWave",  ppm=PPM, 
                           peakwidth=c(4,12), 
                           snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                 Prefilter,
                           integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs10, file="ProgEnegU2 xs10.RData")

ProgEnegU2.xs11 <- xcmsSet(GoodSamples[251:257], method = "centWave",  ppm=PPM, 
                           peakwidth=c(4,12), 
                           snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                                 Prefilter,
                           integrate = 1, fitgauss= TRUE)
save(ProgEnegU2.xs11, file="ProgEnegU2 xs11.RData")


# Tracking how long this step takes
ProgEnegU2.tpick.final <- Sys.time()
ProgEnegU2.tpick <- as.numeric(difftime(ProgEnegU2.tpick.final, 
                                        ProgEnegU2.tpick.init, units = "mins"))
write.csv(ProgEnegU2.tpick, "ProgEnegU2 tpick.csv")



# Initial peak alignment -------------------------------------------------------
ProgEnegU2.tgroup.init <- Sys.time()

# Peak alignment
ProgEnegU2 <- group(c(ProgEnegU2.xs1, ProgEnegU2.xs2, ProgEnegU2.xs3,
                      ProgEnegU2.xs4, ProgEnegU2.xs5, ProgEnegU2.xs6,
                      ProgEnegU2.xs7, ProgEnegU2.xs8, ProgEnegU2.xs9,
                      ProgEnegU2.xs10, ProgEnegU2.xs11), 
                    method = "density", bw = 4, minfrac = 0,
                    minsamp = 1, mzwid = 0.007, max = 100)
# If you split your samples into multiple xcmsSet steps, change the first bit
# of the above command to:
# group (c(ProgEnegU2.xs1, ProgEnegU2.xs2, .. ProgEnegU2.xsn), method = ...)

ProgEnegU2.tgroup.final <- Sys.time()
ProgEnegU2.tgroup <- as.numeric(difftime(ProgEnegU2.tgroup.final, 
                                         ProgEnegU2.tgroup.init, units="mins"))
write.csv(ProgEnegU2.tgroup, "ProgEnegU2 tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(ProgEnegU2.xs1, ProgEnegU2.xs2, ProgEnegU2.xs3,
   ProgEnegU2.xs4, ProgEnegU2.xs5, ProgEnegU2.xs6,
   ProgEnegU2.xs7, ProgEnegU2.xs8, ProgEnegU2.xs9,
   ProgEnegU2.xs10, ProgEnegU2.xs11)


# RT correction -------------------------------------------------------
ProgEnegU2.tretcor.init <- Sys.time() 

# Performing retention time correction
ProgEnegU2 <- retcor(ProgEnegU2, method = "peakgroups", 
                     missing = 0.1*length(GoodSamples), 
                     extra = 2*length(GoodSamples), smooth = "loess", 
                     family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("ProgEnegU2 RTcor plot.png", width = 4, height = 5, units = "in", res = 300)
plotrt(ProgEnegU2, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

ProgEnegU2.tretcor.final <- Sys.time()
ProgEnegU2.tretcor <- as.numeric(difftime(ProgEnegU2.tretcor.final, ProgEnegU2.tretcor.init, units="mins"))
write.csv(ProgEnegU2.tretcor, "ProgEnegU2 tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
ProgEnegU2.tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
ProgEnegU2 <- group(ProgEnegU2, method = "density", 
                    minsamp = 1, minfrac = 0, mzwid = 0.007, 
                    bw = 2, max=100)

# Making a data.frame of the data before the recursive peak filling step
ProgEnegU2.unfilled <- peakTable(ProgEnegU2)
ProgEnegU2.unfilled$MassFeature <- paste("I", round((
      ProgEnegU2.unfilled$mz),digits=4), "R", round(
            ProgEnegU2.unfilled$rt/60, digits=2), sep="")

ProgEnegU2.tgroup2.final <- Sys.time()
ProgEnegU2.tgroup2 <- as.numeric(difftime(ProgEnegU2.tgroup2.final, 
                                          ProgEnegU2.tgroup2.init, units="mins"))
write.csv(ProgEnegU2.tgroup2, "ProgEnegU2 tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
ProgEnegU2.tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
ProgEnegU2 <- fillPeaks(ProgEnegU2)

ProgEnegU2.tfillPeaks.final <- Sys.time()
ProgEnegU2.tfillPeaks <- as.numeric(difftime(ProgEnegU2.tfillPeaks.final, 
                                             ProgEnegU2.tfillPeaks.init, 
                                             units="mins"))
write.csv(ProgEnegU2.tfillPeaks, "ProgEnegU2 tfillPeaks.csv")

# Generate a data.frame with all the peaks ----------------------------------
ProgEnegU2.tPeakTable.init <- Sys.time()

# Making a data.frame of the recursively filled data
ProgEnegU2.allpeaks <- peakTable(ProgEnegU2)

# Changing the name of the count column
names(ProgEnegU2.allpeaks)[names(ProgEnegU2.allpeaks) == CountCol] <- "Count"

# Making a column with the mass feature name
ProgEnegU2.allpeaks$MassFeature <- paste("I", round((
      ProgEnegU2.allpeaks$mz),digits=4), "R", round(
            ProgEnegU2.allpeaks$rt/60, digits=2), sep="")

# Making a column with the mass feature name as xcms sets it
ProgEnegU2.allpeaks$groupname <- groupnames(ProgEnegU2)

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
ProgEnegU2.allpeaks$RT <- ProgEnegU2.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "Epos" | IonizationMode == "Apos",
       RefIons <- c(121.0509, 922.0098), # Epos and Apos       
       RefIons <- c(112.9856, 119.0363, 980.015)) # Eneg and Aneg

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(ProgEnegU2.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 ProgEnegU2.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
ProgEnegU2.noref <- ProgEnegU2.allpeaks[-unlist(RefMFs), ]

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
ProgEnegU2.after2 <- subset(ProgEnegU2.noref, ProgEnegU2.noref$rt 
                            > 120)

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
ProgEnegU2.filter <- subset(ProgEnegU2.after2, Count >= 0.25*length(SampCol))
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
setwd(MainDir)
write.csv(ProgEnegU2.filter, 
          "ProgEnegU2 peak table - mass features must be present in at least 25% of all samples.csv")

save(ProgEnegU2.after2, ProgEnegU2.allpeaks, ProgEnegU2.filter, ProgEnegU2.noref,
     ProgEnegU2.unfilled, file = "ProgEnegU2 - all main dataframes.RData")
save(ProgEnegU2.filter, file = "ProgEnegU2 - filtered dataframe only.RData")

setwd(RawDataDir)
ProgEnegU2.tPeakTable.final <- Sys.time()
ProgEnegU2.tPeakTable <- as.numeric(difftime(ProgEnegU2.tPeakTable.final, 
                                             ProgEnegU2.tPeakTable.init, 
                                             units="mins"))
write.csv(ProgEnegU2.tPeakTable, "ProgEnegU2 tPeakTable.csv")

# Quality control ----------------------------------------
ProgEnegU2.tQC.init <- Sys.time()

setwd(MainDir)

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
set.seed(253)

MFs <- ProgEnegU2.filter$groupname[runif(30, min = 1, 
                                         max = length(
                                               ProgEnegU2.filter$groupname))]
RandSamp <- GoodSamples[runif(10, min = 1, max = length(GoodSamples))]
write.csv(MFs, paste(Sys.Date(), 
                     "ProgEnegU2 randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "ProgEnegU2 randomly selected samples.csv"))

LastSamp <- RandSamp[length(RandSamp)]

EIC.uncorrected <- list()
EIC.corrected <- list()

# This next step will take some time to process, so don't expect instant results. 
for (i in 1:length(MFs)){
      EIC.uncorrected[[i]] <- getEIC(ProgEnegU2, rt="raw", 
                                     groupidx=MFs[i], sampleidx=RandSamp)
      EIC.corrected[[i]] <- getEIC(ProgEnegU2, rt="corrected", 
                                   groupidx=MFs[i], sampleidx=RandSamp)  
}

ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp)-1), "red")

setwd(RawDataDir)
xset.raw <- xcmsRaw(GoodSamples[LastSamp], profstep=0.01, profmethod="bin")

setwd(MainDir)
pdf("ProgEnegU2 quality check.pdf", 8.5,11)

# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the 1st sample for that compound with a 
# dashed horizontal line where the calculated m/z is.

par(mfrow=c(4,3), mar=c(3,3,3,0.5))
for(i in 1:length(MFs)){
      
      plot(EIC.uncorrected[[i]], ProgEnegU2, groupidx=1, rtrange=60, 
           col=MyColors, main=MFs[i])
      mtext(paste(i, ProgEnegU2.filter$MassFeature[
            ProgEnegU2.filter$groupname == MFs[i]]), side=3, line=-1, 
            adj=0, padj=0, cex=0.8)
      plot(EIC.corrected[[i]], ProgEnegU2, groupidx=1, rtrange=60, 
           col=MyColors)
      
      RT <- ProgEnegU2.filter$rt[ProgEnegU2.filter$groupname == MFs[i]]
      RTRange <- c(RT-30, RT+30)
      
      mz <- ProgEnegU2.filter$mz[ProgEnegU2.filter$groupname == MFs[i]]
      mzRange <- c(mz-0.02, mz+0.02)
      mzRange.poly.low <- mz- mz*(0.5*PPM)/1e6
      mzRange.poly.up <- mz*(0.5*PPM)/1e6 + mz
      
      plotRaw(xset.raw, mzrange=mzRange, rtrange=RTRange, log=FALSE)
      abline(h=mz, lty=2, col="gray35")
      mtext(paste("abund =", round(ProgEnegU2.filter[
            ProgEnegU2.filter$groupname == MFs[i], SampCol[LastSamp]], digits=0)), 
            side=3, line=-1, adj=0, padj=0, cex=0.8)
      polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
              c(mzRange.poly.up, mzRange.poly.up, 
                mzRange.poly.low, mzRange.poly.low), 
              col=col2alpha("blue", alpha=0.1), border=NA)
      abline(v=RT, lty=2, col="gray35")
      
}

dev.off()

save(EIC.corrected, EIC.uncorrected, xset.raw, 
     file = "ProgEnegU2 QC data.RData")

setwd(RawDataDir)
ProgEnegU2.tQC.final <- Sys.time()
ProgEnegU2.tQC <- as.numeric(difftime(ProgEnegU2.tQC.final, 
                                      ProgEnegU2.tQC.init, units = "mins"))
write.csv(ProgEnegU2.tQC, "ProgEnegU2 tQC.csv")



# Calculating processing times for each step --------------------------------
setwd(RawDataDir)

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes.
Times <- data.frame(rbind(ProgEnegU2.tpick, ProgEnegU2.tgroup, 
                          ProgEnegU2.tretcor, ProgEnegU2.tgroup2, 
                          ProgEnegU2.tfillPeaks, 
                          ProgEnegU2.tPeakTable, ProgEnegU2.tQC))
Times$Step <- c("pick peaks", # ProgEnegU2.tpick
                "group 1", # ProgEnegU2.tgroup
                "retcor", # ProgEnegU2.tretcor
                "group 2", # ProgEnegU2.tgroup2
                "fill peaks", # ProgEnegU2.tfillPeaks
                "peak table", # ProgEnegU2.tPeakTable
                "QC") # ProgEnegU2.tQC
names(Times) <- c("Duration", "Step")
row.names(Times) <- 1:nrow(Times)

# Calculate the total time in hours. 
TotalTime <- c("Total time in hours", format(sum(Times$Duration)/60, digits=2))
Times$Duration <- round(Times$Duration, digits=2)
Times <- Times[, c("Step", "Duration")]

setwd(MainDir)
write.csv(rbind(Times, TotalTime), "ProgEnegU2 processing times.csv", 
          row.names=F)
write.csv(data.frame(NumSamp = length(GoodSamples),
                     PeakPickPerSamp = ProgEnegU2.tpick/length(GoodSamples)),
          "ProgEnegU2 peak picking time per sample.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for ProgEnegU2") +
      guides(fill=FALSE)
ggsave("ProgEnegU2 processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for ProgEnegU2") +
      guides(fill=FALSE)
ggsave("ProgEnegU2 processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                      "Count" = c(nrow(ProgEnegU2.unfilled),
                                  nrow(ProgEnegU2.noref),
                                  nrow(ProgEnegU2.after2),
                                  nrow(ProgEnegU2.filter)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "ProgEnegU2 counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("ProgEnegU2 bar chart of numbers of mass features at each step.png")


# Saving final workspace ------------------------------
setwd(MainDir)
save(ProgEnegU2, file = "ProgEnegU2 xcmsSet object.RData")
# save.image("ProgEnegU2 workspace.RData") # This saves EVERYTHING that is 
# currently in your workspace, which is a pretty huge file. Skip this step if 
# you don't think you'll need that. 
