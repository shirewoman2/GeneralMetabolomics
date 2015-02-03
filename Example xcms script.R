# Example xcms script
# 2/3/15 LS

# Notes to the users of this script ---------------------------------


###       ITEMS YOU WILL WANT TO ADJUST EACH TIME YOU RUN THIS SCRIPT:       ###
# 1. Set the working directory (RawDataDir).

# 2. If you have more than 25 or so samples, split up your samples when you run
# xcmsSet() so that you save your progress occasionally. For example, set up 
# the code like this:
#     xs1 <- xcmsSet(Samples[1:25], method = "centWave",  ppm=15, peakwidth=c(4,12), 
#     snthresh = 20, mzCenterFun="apex", prefilter=c(10, 10000),
#     integrate = 1, fitgauss= TRUE)
#     save(xs1, file="xcmsSet objects.RData")
#
#     xs2 <- xcmsSet(Samples[26:50], method = "centWave",  ppm=15, peakwidth=c(4,12), 
#     snthresh = 20, mzCenterFun="apex", prefilter=c(10, 10000),
#     integrate = 1, fitgauss= TRUE)
#     save(c(xs1, xs2), file="xcmsSet objects.RData")
# etc. until you've done peak picking on all of your samples. That way, if the
# computer shuts down or crashes, you've saved your progress and don't have to 
# start over. To load those files into memory, here is the code:
#     load("xcmsSet objects.RData")

# 3. Parameters I commonly tweak to get things to look right:
#     xcmsSet() -- snthresh and prefilter
#     group() -- bw

# 4. Find and replace "SCORMDZEnegP83" with some single-word name that has meaning to
# you and doesn't start with a number or symbol (those are the rules for naming
# objects in R). Now, your objects will be clearly named and the output files
# will be similarly named.

# 5. I added a requirement for the quality checking step that any mass feature
# we look at must be present in at least 25% of all samples. You can change that
# number to whatever you want, though. See the "Quality control" section. 

# 6. Set the ionization mode (see the section "Housekeeping") as appropriate 
# for your data for CAMERA to work properly and for the appropriate set of 
# reference ions to be removed from your data.



###                          GENERAL TIPS                                   ###
# Do yourself a favor and make a metadata file that includes information on
# each sample. Also include any info you might later want to sort or filter by
# such as matrix, date prepped, date run, what kind of sample is it (clinical,
# QC, Master QC, water), treatment group, whether that injection was a good one
# if some of them were not, etc. Then, use that file throughout your analyses, 
# including here, where you can use it to pick which samples you want to 
# work with. 



# Housekeeping -------------------------------------------------
# This script uses the package xcms to generate a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled. It then performs
# a check on the quality of the data extraction, collapses ions that are likely
# isotopic peaks, and checks on the quality of the selection of isotopes.

# Dataset: SCORMDZEnegP83

IonizationMode <- "negative" # Change to "positive" or "negative" as appropriate.

library(plyr)
library(ggplot2)
library(stringr)
library(xcms)


# Set path of the working directory, i.e. the folder where your files are. 
# Note that you need forward slashes!
RawDataDir <- "F:/SCOR/20120202 SCOR plasma and urine ESI neg"

MainDir <- "C:/Users/Laura/Documents/SCOR project/79 SCOR EnegP"

# Getting samples and peak picking ------------------------------------------

setwd("F:/SCOR")
Meta <- read.csv("SCOR MDZ metadata.csv", skip = 1, 
                 na.strings = c("has not run", ""))

# Removing files that have not run.
Meta <- Meta[complete.cases(Meta$File), ]
Meta$File <- paste0(Meta$File, ".mzdata.xml")

# Changing the ionization mode to match what I use elsewhere
Meta$Mode <- revalue(Meta$Mode, c("APCI+" = "Apos",
                                  "ESI+" = "Epos",
                                  "ESI-" = "Eneg"))

# Looking for the appropriate files
setwd(RawDataDir)

Samples <- Meta$File[which(Meta$Mode == "Eneg" & Meta$Matrix == "plasma" &
                                 Meta$Subject != "Master QC")]

SampleFileCheck <- file.exists(Samples)
# OK.

# Getting the names of the output tables' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(Samples))

# Making a note of start time for this step
tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
SNthresh <- 30
Prefilter <- c(10, 5000)

SCORMDZEnegP83.xs1 <- xcmsSet(Samples[1:15], method = "centWave",  ppm=PPM, 
                              peakwidth=c(4,12), 
                              snthresh = SNthresh, mzCenterFun="apex", 
                              prefilter = Prefilter,
                              integrate = 1, fitgauss= TRUE)
save(SCORMDZEnegP83.xs1, file="SCORMDZEnegP83 xs1.RData")


SCORMDZEnegP83.xs2 <- xcmsSet(Samples[16:32], method = "centWave",  ppm=PPM, 
                              peakwidth=c(4,12), 
                              snthresh = SNthresh, mzCenterFun="apex", 
                              prefilter = Prefilter,
                              integrate = 1, fitgauss= TRUE)
save(SCORMDZEnegP83.xs2, file="SCORMDZEnegP83 xs2.RData")

# Tracking how long this step takes
tpick.final <- Sys.time()
tpick <- as.numeric(difftime(tpick.final, tpick.init, units="mins")) 

write.csv(tpick, "tpick.csv")


# Initial peak alignment -------------------------------------------------------
tgroup.init <- Sys.time()

# Peak alignment
SCORMDZEnegP83 <- group(c(SCORMDZEnegP83.xs1, SCORMDZEnegP83.xs2), 
                        method = "density", bw = 4, minfrac = 0,
                        minsamp = 1, mzwid = 0.007, max = 100)

tgroup.final <- Sys.time()
tgroup <- as.numeric(difftime(tgroup.final, tgroup.init, units="mins"))
write.csv(tgroup, "tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(xs1, xs2, xs3)


# RT correction -------------------------------------------------------
tretcor.init <- Sys.time() 

# Performing retention time correction
SCORMDZEnegP83 <- retcor(SCORMDZEnegP83, method = "peakgroups", 
                         missing = 0.1*length(Samples), 
                         extra = 2*length(Samples), smooth = "loess", 
                         family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("SCORMDZEnegP83 RTcor plot.png", width = 6, height = 8, units = "in", res = 300)
plotrt(SCORMDZEnegP83, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

tretcor.final <- Sys.time()
tretcor <- as.numeric(difftime(tretcor.final, tretcor.init, units="mins"))
write.csv(tretcor, "tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
SCORMDZEnegP83 <- group(SCORMDZEnegP83, method = "density", 
                        minsamp = 1, minfrac = 0, mzwid = 0.007, 
                        bw = 2, max=100)

# Making a table of the data before the recursive peak filling step
SCORMDZEnegP83.unfilled <- peakTable(SCORMDZEnegP83)

setwd(MainDir)
write.csv(SCORMDZEnegP83.unfilled, 
          "SCORMDZEnegP83 peak table - unfilled.csv")
setwd(RawDataDir)

tgroup2.final <- Sys.time()
tgroup2 <- as.numeric(difftime(tgroup2.final, tgroup2.init, units="mins"))
write.csv(tgroup2, "tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
SCORMDZEnegP83 <- fillPeaks(SCORMDZEnegP83)

tfillPeaks.final <- Sys.time()
tfillPeaks <- as.numeric(difftime(tfillPeaks.final, tfillPeaks.init, 
                                  units="mins"))
write.csv(tfillPeaks, "tfillPeaks.csv")

# Generate a table with all the peaks ----------------------------------------
tPeakTable.init <- Sys.time()

# Making a table of the recursively filled data
SCORMDZEnegP83.allpeaks <- peakTable(SCORMDZEnegP83)

# Making a column with the mass feature name
SCORMDZEnegP83.allpeaks$MassFeature <- paste("I", round((
      SCORMDZEnegP83.allpeaks$mz),digits=4), "R", round(
            SCORMDZEnegP83.allpeaks$rt/60, digits=2), sep="")

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
SCORMDZEnegP83.allpeaks$RT <- SCORMDZEnegP83.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "positive",
       RefIons <- c(121.0509, 922.0098), # ESI+       
       RefIons <- c(112.9856, 119.0363, 980.015)) # ESI-

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(SCORMDZEnegP83.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 SCORMDZEnegP83.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
SCORMDZEnegP83.noref <- SCORMDZEnegP83.allpeaks[-unlist(RefMFs), ]

setwd(MainDir)
write.csv(SCORMDZEnegP83.noref, "SCORMDZEnegP83 peak table without ref ions.csv")

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
SCORMDZEnegP83.after2 <- subset(SCORMDZEnegP83.noref, SCORMDZEnegP83.noref$rt 
                                > 120)
write.csv(SCORMDZEnegP83.after2, 
          "SCORMDZEnegP83 peak table - RT above 2 min.csv")

# Counting how many times a mass feature was detected BEFORE recursive peak
# filling step.
Count.df <- SCORMDZEnegP83.unfilled[as.numeric(row.names(
      SCORMDZEnegP83.after2)), SampCol]

SCORMDZEnegP83.after2$Count <- rep(NA, nrow(SCORMDZEnegP83.after2))

for (m in 1:nrow(Count.df)){
      SCORMDZEnegP83.after2$Count[m] <- length(which(!is.na(as.numeric(
            Count.df[m, ]))))
}

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
SCORMDZEnegP83.filter <- SCORMDZEnegP83.after2[SCORMDZEnegP83.after2$Count >= 
                                                     0.25*length(SampCol), ]
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
rm(Count.df)
write.csv(SCORMDZEnegP83.filter, 
          "SCORMDZEnegP83 peak table - mass features must be present in at least 25% of all samples.csv")

setwd(RawDataDir)
# Saving the workspace up to here. 
save.image("SCORMDZEnegP83 workspace.RData")

tPeakTable.final <- Sys.time()
tPeakTable <- as.numeric(difftime(tPeakTable.final, tPeakTable.init, 
                                  units="mins"))
write.csv(tPeakTable, "tPeakTable.csv")

# Quality control ----------------------------------------
tQC.init <- Sys.time()

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
MFs <- as.numeric(sample(row.names(SCORMDZEnegP83.filter), 30))
RandSamp <- as.numeric(sample(length(Samples), 10))
write.csv(MFs, paste(Sys.Date(), 
                     "SCORMDZEnegP83 randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "SCORMDZEnegP83 randomly selected samples.csv"))

EIC.uncorrected <- list()
EIC.corrected <- list()

# Looking at the extracted ion chromatograms of each sample and mass feature.
# This will take some time to process, so don't expect instant results. 
for (m in MFs){
      EIC.uncorrected[[m]] <- getEIC(SCORMDZEnegP83, rt = "raw", 
                                     groupidx = m, sampleidx = RandSamp)
      EIC.corrected[[m]] <- getEIC(SCORMDZEnegP83, rt = "corrected", 
                                   groupidx = m, sampleidx = RandSamp)  
}

# Setting the palette for drawing the EICs. The last sample in "RandSamp"
# will be colored red and will be the one used to generate the raw mass 
# spectral data.
ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp) - 1), "red")

xset.raw <- xcmsRaw(Samples[RandSamp[10]], profstep = 0.01, profmethod = "bin")

setwd(MainDir)
pdf("SCORMDZEnegP83 peak picking and alignment quality check.pdf", 8.5, 11)

# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the last sample in RandSamp. A dashed 
# horizontal line shows where the calculated m/z is; a dashed vertical line 
# shows the median RT for that compound. A semi-transparent gray rectangle shows 
# the compound's m/z within a window equal to the value of PPM.

par(mfrow=c(4,3), mar=c(3,3,3,0.5))
for(i in 1:30){
      m <- MFs[i]
      plot(EIC.uncorrected[[m]], SCORMDZEnegP83, groupidx = 1, 
           rtrange = 60, col = MyColors, main = MFs[m])
      mtext(paste(i, SCORMDZEnegP83.allpeaks$MassFeature[m]), side = 3, line = -1, 
            adj = 0, padj = 0, cex = 0.8)
      plot(EIC.corrected[[m]], SCORMDZEnegP83, groupidx = 1, 
           rtrange = 60, col = MyColors)
      
      RT <- SCORMDZEnegP83.allpeaks$rt[m]
      RTRange <- c(RT - 30, RT + 30)
      
      mz <- SCORMDZEnegP83.allpeaks$mz[m]
      mzRange <- c(mz - 0.02, mz + 0.02)
      mzRange.poly.low <- mz - mz * 7.5/1e6
      mzRange.poly.up <- mz * 7.5/1e6 + mz
      
      plotRaw(xset.raw, mzrange = mzRange, rtrange = RTRange, log = FALSE)
      abline(h = mz, lty = 2, col = "gray35")
      mtext(paste("abund =", round(SCORMDZEnegP83.allpeaks[m, (length(RandSamp))], 
                                   digits = 0)), 
            side = 3, line = -1, adj = 0, padj = 0, cex = 0.8)
      polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
              c(mzRange.poly.up, mzRange.poly.up, mzRange.poly.low, 
                mzRange.poly.low), 
              col = col2alpha("blue", alpha = 0.1), border=NA)
      abline(v = RT, lty = 2, col = "gray35")
      
}

dev.off()

setwd(RawDataDir)

tQC.final <- Sys.time()
tQC <- as.numeric(difftime(tQC.final, tQC.init, units = "mins"))
write.csv(tQC, "tQC.csv")

# Calculating processing times for each step --------------------------------------------
tpick <- read.csv("tpick.csv")
tgroup <- read.csv("tgroup.csv")
tretcor <- read.csv("tretcor.csv")
tgroup2 <- read.csv("tgroup2.csv")
tfillPeaks <- read.csv("tfillPeaks.csv")
tPeakTable <- read.csv("tPeakTable.csv")
tQC <- read.csv("tQC.csv")

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes.
Times <- rbind(tpick, tgroup, tretcor, tgroup2, tfillPeaks, tPeakTable, tQC)
Times$Step <- c("pick peaks", # tpick
                "group 1", # tgroup
                "retcor", # tretcor
                "group 2", # tgroup2
                "fill peaks", # tfillPeaks
                "peak table", # tPeakTable
                "QC") # tQC

Times <- rename(Times, c("x" = "Duration"))

# Calculate the total time in hours. 
TotalTime <- c("Total time in hours", format(sum(Times$Duration)/60, digits=2))
Times$Duration <- round(Times$Duration, digits=2)
Times <- Times[, c("Step", "Duration")]

setwd(MainDir)
write.csv(rbind(Times, TotalTime), "SCORMDZEnegP83 processing times.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SCORMDZEnegP83") +
      guides(fill=FALSE)
ggsave("SCORMDZEnegP83 processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SCORMDZEnegP83") +
      guides(fill=FALSE)
ggsave("SCORMDZEnegP83 processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                      "Count" = c(nrow(SCORMDZEnegP83.unfilled),
                                  nrow(SCORMDZEnegP83.noref),
                                  nrow(SCORMDZEnegP83.after2),
                                  nrow(SCORMDZEnegP83.filter)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "SCORMDZEnegP83 counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("SCORMDZEnegP83 bar chart of numbers of mass features at each step.png")


# Saving final workspace ------------------------------
setwd(RawDataDir)
save(SCORMDZEnegP83, file = "SCORMDZEnegP83 xcmsSet object.RData")
save.image("SCORMDZEnegP83 workspace.RData") # This saves EVERYTHING that is currently
# in your workspace, which is a pretty huge file. Skip this step if you don't 
# think you'll need that. 
