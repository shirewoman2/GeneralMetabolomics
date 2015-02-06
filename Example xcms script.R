# DataEnegP xcms
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

# 4. Find and replace "DataEnegP" with some single-word name that has meaning to
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

# Dataset: DataEnegP

IonizationMode <- "negative" # Change to "positive" or "negative" as appropriate.

library(plyr)
library(ggplot2)
library(gridExtra)
library(xcms)
library(CAMERA)
library(XLConnect)

# Set path of the working directory, i.e. the folder where your files are. 
# Note that you need forward slashes!
RawDataDir <- "D:/Users/Laura/Documents/Work/Lin Lab/LCMS metabolomics/Whole blood"

MainDir <- paste0(RawDataDir, "/WBEnegP4")

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

# Multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      require(grid)
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
            # Make the panel
            # ncol: Number of columns of plots
            # nrow: Number of rows needed, calculated from # of cols
            layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                             ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
            print(plots[[1]])
            
      } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            
            # Make each plot, in the correct location
            for (i in 1:numPlots) {
                  # Get the i,j matrix positions of the regions that contain this subplot
                  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                  
                  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
            }
      }
}

# Getting samples and peak picking ------------------------------------------
# Using an Excel file with metadata to select which samples I want to process.
setwd(RawDataDir)
wb <- loadWorkbook("20141014 Whole blood workup worklist.xlsx")
Meta <- readWorksheet(wb, "Sheet2")
Meta$File <- paste0(Meta$File, ".mzdata.xml")
Samples <- Meta$File[file.exists(Meta$File) & Meta$Mode == "ESI-"]


# Getting the names of the output data.frames' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(Samples))

# Making a note of start time for this step
tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
SNthresh <- 15
Prefilter <- c(10, 5000)

DataEnegP.xs1 <- xcmsSet(Samples, method = "centWave",  ppm=PPM, 
                         peakwidth=c(4,12), 
                         snthresh = SNthresh, mzCenterFun="apex", prefilter = Prefilter,
                         integrate = 1, fitgauss= TRUE)
save(DataEnegP.xs1, file="DataEnegP xs1.RData")


# Tracking how long this step takes
tpick.final <- Sys.time()
tpick <- as.numeric(difftime(tpick.final, tpick.init, units="mins")) 

write.csv(tpick, "tpick.csv")


# Initial peak alignment -------------------------------------------------------
tgroup.init <- Sys.time()

# Peak alignment
DataEnegP <- group(c(DataEnegP.xs1), method = "density", bw = 4, minfrac = 0,
                   minsamp = 1, mzwid = 0.007, max = 100)
# If you split your samples into multiple xcmsSet steps, change the first bit
# of the above command to:
# group (c(DataEnegP.xs1, DataEnegP.xs2, .. DataEnegP.xsn), method = ...)

tgroup.final <- Sys.time()
tgroup <- as.numeric(difftime(tgroup.final, tgroup.init, units="mins"))
write.csv(tgroup, "tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(DataEnegP.xs1)


# RT correction -------------------------------------------------------
tretcor.init <- Sys.time() 

# Performing retention time correction
DataEnegP <- retcor(DataEnegP, method = "peakgroups", 
                    missing = 0.1*length(Samples), 
                    extra = 2*length(Samples), smooth = "loess", 
                    family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("DataEnegP RTcor plot.png", width = 4, height = 5, units = "in", res = 300)
plotrt(DataEnegP, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

tretcor.final <- Sys.time()
tretcor <- as.numeric(difftime(tretcor.final, tretcor.init, units="mins"))
write.csv(tretcor, "tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
DataEnegP <- group(DataEnegP, method = "density", 
                   minsamp = 1, minfrac = 0, mzwid = 0.007, 
                   bw = 2, max=100)

# Making a data.frame of the data before the recursive peak filling step
DataEnegP.unfilled <- peakTable(DataEnegP)

setwd(MainDir)
write.csv(DataEnegP.unfilled, 
          "DataEnegP peak table - unfilled.csv")
setwd(RawDataDir)

tgroup2.final <- Sys.time()
tgroup2 <- as.numeric(difftime(tgroup2.final, tgroup2.init, units="mins"))
write.csv(tgroup2, "tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
DataEnegP <- fillPeaks(DataEnegP)

tfillPeaks.final <- Sys.time()
tfillPeaks <- as.numeric(difftime(tfillPeaks.final, tfillPeaks.init, 
                                  units="mins"))
write.csv(tfillPeaks, "tfillPeaks.csv")

# Generate a data.frame with all the peaks ----------------------------------
tPeakTable.init <- Sys.time()

# Making a data.frame of the recursively filled data
DataEnegP.allpeaks <- peakTable(DataEnegP)

# Making a column with the mass feature name
DataEnegP.allpeaks$MassFeature <- paste("I", round((
      DataEnegP.allpeaks$mz),digits=4), "R", round(
            DataEnegP.allpeaks$rt/60, digits=2), sep="")

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
DataEnegP.allpeaks$RT <- DataEnegP.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "positive",
       RefIons <- c(121.0509, 922.0098), # ESI+       
       RefIons <- c(112.9856, 119.0363, 980.015)) # ESI-

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(DataEnegP.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 DataEnegP.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
DataEnegP.noref <- DataEnegP.allpeaks[-unlist(RefMFs), ]

# Counting how many times a mass feature was detected BEFORE recursive peak
# filling step.
nonzero <- function(x) {length(x[complete.cases(as.numeric(x))])}
DataEnegP.noref$Count <- apply(DataEnegP.noref[, SampCol], 1, nonzero)

# Changing order of columns
DataEnegP.noref <- DataEnegP.noref[, c("MassFeature", "mz", "RT", 
                                       "mzmax", "mzmin", "rt", "rtmax", 
                                       "rtmin", "Count", SampCol)]

setwd(MainDir)
write.csv(DataEnegP.noref, "DataEnegP peak table without ref ions.csv")

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
DataEnegP.after2 <- subset(DataEnegP.noref, DataEnegP.noref$rt 
                           > 120)
write.csv(DataEnegP.after2, 
          "DataEnegP peak table - RT above 2 min.csv")

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
DataEnegP.filter <- DataEnegP.after2[DataEnegP.after2$Count >= 
                                           0.25*length(SampCol), ]
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
write.csv(DataEnegP.filter, 
          "DataEnegP peak table - mass features must be present in at least 25% of all samples.csv")

save(DataEnegP.after2, DataEnegP.allpeaks, DataEnegP.filter, DataEnegP.noref,
     DataEnegP.unfilled, file = "DataEnegP - all main dataframes.RData")
save(DataEnegP.filter, file = "DataEnegP - filtered dataframe only.RData")

setwd(RawDataDir)
tPeakTable.final <- Sys.time()
tPeakTable <- as.numeric(difftime(tPeakTable.final, tPeakTable.init, 
                                  units="mins"))
write.csv(tPeakTable, "tPeakTable.csv")

# Quality control ----------------------------------------
tQC.init <- Sys.time()

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
MFs <- as.numeric(sample(1:nrow(DataEnegP.filter), 30))
RandSamp <- as.numeric(sample(length(Samples), 10))
write.csv(MFs, paste(Sys.Date(), 
                     "DataEnegP randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "DataEnegP randomly selected samples.csv"))

EIC.uncorrected <- list()
EIC.corrected <- list()

# Looking at the extracted ion chromatograms of each sample and mass feature.
# This will take some time to process, so don't expect instant results. 
for (i in 1:length(MFs)){
      EIC.uncorrected[[i]] <- getEIC(DataEnegP, rt = "raw", 
                                     groupidx = MFs[i], sampleidx = RandSamp)
      EIC.corrected[[i]] <- getEIC(DataEnegP, rt = "corrected", 
                                   groupidx = MFs[i], sampleidx = RandSamp)  
}

# Setting the palette for drawing the EICs. The last sample in "RandSamp"
# will be colored red and will be the one used to generate the raw mass 
# spectral data.
ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp) - 1), "red")

xset.raw <- xcmsRaw(Samples[RandSamp[10]], profstep = 0.01, profmethod = "bin")

setwd(MainDir)
# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the last sample in RandSamp. A dashed 
# horizontal line shows where the calculated m/z is; a dashed vertical line 
# shows the median RT for that compound. A semi-transparent gray rectangle shows 
# the compound's m/z within a window equal to the value of PPM.

Plots <- list()

EICuncor <- list()
EICcor <- list()
MSraw <- list()
MFparam <- list()

EICplot <- list()

for(i in 1:30){
      m <- MFs[i]
      
      EICuncor[[i]] <- list()
      EICcor[[i]] <- list()
      
      for (s in 1:length(RandSamp)){
            EICuncor[[i]][[s]] <- data.frame(EIC.uncorrected[[i]]@eic[[s]][[1]])
            EICuncor[[i]][[s]]$Sample <- RandSamp[s]
            EICuncor[[i]][[s]]$SpecType <- "uncorrected"
            EICcor[[i]][[s]] <- data.frame(EIC.corrected[[i]]@eic[[s]][[1]])
            EICcor[[i]][[s]]$Sample <- as.factor(RandSamp[s])
            EICcor[[i]][[s]]$SpecType <- "corrected"
      }
      
      EICuncor[[i]] <- rbind.fill(EICuncor[[i]])
      EICuncor[[i]]$Sample <- as.factor(EICuncor[[i]]$Sample)
      EICcor[[i]] <- rbind.fill(EICcor[[i]])
      EICcor[[i]]$Sample <- as.factor(EICcor[[i]]$Sample)
      
      RT <- DataEnegP.allpeaks$rt[m]
      RTRange <- c(RT - 30, RT + 30)
      
      mz <- DataEnegP.allpeaks$mz[m]
      mzRange <- c(mz - 0.02, mz + 0.02)
      
      MSraw[[i]] <- data.frame(plotRaw(xset.raw, mzrange = mzRange, 
                                       rtrange = RTRange, log = FALSE))
      MFparam[[i]] <- data.frame(MassFeature = DataEnegP.allpeaks$MassFeature[m],
                                 mz = mz,
                                 RT = RT/60,
                                 mzRange.rect.low = mz - mz * 7.5/1e6,
                                 mzRange.rect.high = mz * 7.5/1e6 + mz)

      EICplot[[3*i-2]] <- ggplot(EICuncor[[i]], aes(x = rt, y = intensity, 
                                                    color = Sample)) +
            geom_line() + 
            scale_color_manual(values = MyColors) +
            ggtitle("Uncorrected") +
            ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 1.5,
                              hjust = -0.05, size = 3,
                              label = DataEnegP.filter$MassFeature[MFs[i]]) +
            theme(legend.position = "none", 
                  text = element_text(size = 8),
                  axis.text.y = element_text(angle = 90, hjust = 0.5))
      
      EICplot[[3*i-1]] <- ggplot(EICcor[[i]], aes(x = rt, y = intensity, 
                                                  color = Sample)) +
            geom_line() + 
            scale_color_manual(values = MyColors) +
            ggtitle("Corrected") +
            ggplot2::annotate("text", x = -Inf, y = Inf, vjust = 1.5,
                              hjust = -0.05, size = 3,
                              label = DataEnegP.filter$MassFeature[MFs[i]]) +
            theme(legend.position = "none", 
                  text = element_text(size = 8),
                  axis.text.y = element_text(angle = 90, hjust = 0.5))
      
      
      EICplot[[3*i]]<- ggplot(MSraw[[i]], aes(x = time/60, y = mz, color = intensity)) +
            geom_point() + xlab("RT (min)") + ylab("m/z") +
            ggplot2::annotate("rect", xmin = min(MSraw[[i]]$time/60), 
                              xmax = max(MSraw[[i]]$time/60), 
                              ymin = MFparam[[i]]$mzRange.rect.low, 
                              ymax = MFparam[[i]]$mzRange.rect.high,
                              alpha = 0.1) +
            geom_vline(xintercept = MFparam[[i]]$RT, linetype="dashed") +
            geom_hline(yintercept = MFparam[[i]]$mz, linetype = "dashed") +
            theme(legend.position = "none", 
                  text = element_text(size = 8),
                  axis.text.y = element_text(angle = 90, hjust = 0.5, 
                                             size = 6))
      
      
      
      
      
}

p <- do.call(marrangeGrob, c(EICplot, ncol = 3, nrow = 4))
ggsave("DataEnegP quality check.pdf", p, width = 8.5, height = 11)

setwd(RawDataDir)
tQC.final <- Sys.time()
tQC <- as.numeric(difftime(tQC.final, tQC.init, units = "mins"))
write.csv(tQC, "tQC.csv")

# Calculating processing times for each step --------------------------------------------
setwd(RawDataDir)

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes.
Times <- data.frame(rbind(tpick, tgroup, tretcor, tgroup2, tfillPeaks, 
                          tPeakTable, tQC))
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
write.csv(rbind(Times, TotalTime), "DataEnegP processing times.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for DataEnegP") +
      guides(fill=FALSE)
ggsave("DataEnegP processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for DataEnegP") +
      guides(fill=FALSE)
ggsave("DataEnegP processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                      "Count" = c(nrow(DataEnegP.unfilled),
                                  nrow(DataEnegP.noref),
                                  nrow(DataEnegP.after2),
                                  nrow(DataEnegP.filter)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "DataEnegP counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("DataEnegP bar chart of numbers of mass features at each step.png")


# Saving final workspace ------------------------------
setwd(RawDataDir)
save(DataEnegP, file = "DataEnegP xcmsSet object.RData")
# save.image("DataEnegP workspace.RData") # This saves EVERYTHING that is currently
# in your workspace, which is a pretty huge file. Skip this step if you don't 
# think you'll need that. 
