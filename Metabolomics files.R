# Metabolomics files

# This script puts together a list of metabolomics files that I've generated
# or used.

# Housekeeping ====================================================
library(stringr)
library(XLConnect)
library(tidyr)
library(plyr)
library(dplyr)
library(lubridate)
library(reshape2)


MainDir <- "C:/Users/Laura/Documents/LCMS metabolomics"

# Loading data ==================================================
# Busulf --------------------------------------------------
setwd("F:/Busulfan/Busulfan_Postdose_Data_10_2014")
load("Busulfan retrospective tidy metadata.RData")

BusulfFiles <- Meta
rm(Meta, Meta.all, Meta.all.clin, Meta.clin, Clin)

BusulfFiles <- plyr::rename(BusulfFiles, c("drugs.present.in.sample" = "Special",
                                           "UPN" = "Subject"))
BusulfFiles$Matrix <- "plasma"
BusulfFiles$Date <- ymd(str_sub(BusulfFiles$File, 1, 8))
BusulfFiles$Project <- "Busulf"
BusulfFiles$Directory <- "F:/Busulfan/Busulfan_Postdose_Data_10_2014/ESI_neg"
BusulfFiles$Directory[BusulfFiles$Mode == "Epos"] <- 
      "F:/Busulfan/Busulfan_Postdose_Data_10_2014/ESI_pos"

BusulfFiles <- BusulfFiles[, c("SampleID", "Subject", "Special", "File", "Mode",
                               "Matrix", "SampType", "Date", "Project", 
                               "Directory")]

# CHC2 ------------------------------------------------------------
setwd("C:/Users/Laura/Documents/Farm Worker Families (CHC2)/CHC2 statistical analyses for DEOHS")
library(xlsx)
# I couldn't get R to read from the file "CHC2 metadata v3.xlsx" for some
# reason, so this file is just the copied and "paste values"ed verion 
# of the "metadata - long" tab of the former.
CHC2Files <- read.xlsx("CHC2 metadata v3 long.xlsx", sheetName = "Sheet1",
                       colClasses = rep("character", 17))
CHC2Files <- plyr::rename(CHC2Files, c("Filename" = "File", 
                                       "Use.this." = "Use", 
                                       "Sample.type" = "SampType",
                                       "Sample" = "SampleID",
                                       "Date.run" = "Date"))
CHC2Files$Mode[CHC2Files$Dataset == "EposU"] <- "Epos"
CHC2Files$Mode[CHC2Files$Dataset == "EnegU"] <- "Eneg"

CHC2Files$Matrix <- "urine"

CHC2Files$Directory[CHC2Files$Mode == "Epos"] <- "F:/CHC2/EposU"
CHC2Files$Directory[CHC2Files$Mode == "Eneg"] <- "F:/CHC2/EnegU"

CHC2Files$SampType <- revalue(CHC2Files$SampType, c("M" = "Master QC",
                                                    "Q" = "QC", 
                                                    "U" = "clinical"))

CHC2Files$Date <- mdy(CHC2Files$Date)
CHC2Files$Project <- "CHC2"

# CYP2D6 DEX ------------------------------------------------
setwd("F:/CYP2D6 DEX")
load("CYP2D6 metadata.RData")

DEXFiles <- Meta.all
rm(Meta, Meta.all, Clin, NonexistentInjections, Samples.long, Dates, 
   ExistentSamples)

DEXFiles$Matrix <- "urine"
DEXFiles$Date <- mdy(str_sub(DEXFiles$File, 1, 8))
DEXFiles$Project <- "CYP2D6 DEX"
DEXFiles$Directory <- paste0("F:/CYP2D6 DEX/", DEXFiles$Dir)


# CYP2D6 Metoprolol -----------------------------------------
setwd("C:/Users/Laura/Documents/CYP2D6 metoprolol")
MetopFiles <- read.csv("CYP2D6 metoprolol metadata.csv", skip=1, 
                       na.strings=c("","NA", "#N/A", "Discontinued study", 
                                    "#VALUE!", "has not run"))
MetopFiles$SampleID <- as.character(MetopFiles$SampleID)

MetopFiles <- plyr::rename(MetopFiles, c("PregStage" = "Special",
                                       "SampleType" = "SampType"))

MetopFiles <- subset(MetopFiles, , c("SampleID", "Subject", "Special",
                                   "Matrix", "EposFile", "EnegFile", 
                                   "SampType"))

MetopFiles <- gather(MetopFiles, Mode, File, -SampleID, -Subject, -Special, 
                    -Matrix, -SampType)
MetopFiles$Mode <- revalue(MetopFiles$Mode, c("EposFile" = "Epos", 
                                            "EnegFile" = "Eneg"))

MetopFiles <- subset(MetopFiles, complete.cases(File))
MetopFiles$Date <- ymd(str_sub(MetopFiles$File, 1, 8))
MetopFiles$Project <- "CYP2D6 Metoprolol"
MetopFiles$Directory[MetopFiles$Matrix == "urine" &
                          MetopFiles$Mode == "Epos"] <- 
      "F:/Metoprolol/Metoprolol EposU"
MetopFiles$Directory[MetopFiles$Matrix == "urine" &
                          MetopFiles$Mode == "Eneg"] <- 
      "F:/Metoprolol/Metoprolol EnegU"
MetopFiles$Directory[MetopFiles$Matrix == "plasma" &
                          MetopFiles$Mode == "Epos"] <- 
      "F:/Metoprolol/Metoprolol EposP"
MetopFiles$Directory[MetopFiles$Matrix == "plasma" &
                          MetopFiles$Mode == "Eneg"] <- 
      "F:/Metoprolol/Metoprolol EnegP"


# IVPO MDZ files --------------------------------------------------
setwd("G:/Data/Metabolomics/Laura/IV MDZ project")
load("IVPO MDZ metadata.RData")

IVPOFiles <- Files
rm(Samples, Subjects, Visits, Files)
IVPOFiles <- plyr::rename(IVPOFiles, c("SubjectID" = "Subject",
                                       "DateTime" = "Date",
                                       "DoseType" = "Special"))
IVPOFiles$Project <- "IVPO MDZ"
IVPOFiles$Directory <- "G:/Data/Metabolomics/Laura/IV MDZ project/IVPO MDZ raw data"


# MnPS files ---------------------------------------------
setwd("F:/Mn exposure/Mn exposure Puget Sound workers")
load("Mn Puget Sound metadata.RData")
MnPSFiles <- Meta
rm(Meta)

MnPSFiles <- plyr::rename(MnPSFiles, c("Group" = "Special", 
                                       "SubjID" = "Subject"))
MnPSFiles$SampleID <- as.character(MnPSFiles$SampleID)
MnPSFiles$Date <- ymd(str_sub(MnPSFiles$File, 1,8))
MnPSFiles$Project <- "MnPS"
MnPSFiles$Directory <- "F:/Mn exposure/Mn exposure Puget Sound workers/Mn exposure Puget Sound workers raw data"


# MnWI --------------------------------------------------
setwd("F:/Mn exposure/Mn WI")
load("MnWI metadata.RData")

MnWIFiles <- Files
rm(Files, Samples, Subjects)

MnWIFiles <- plyr::rename(MnWIFiles, c("DateTime" = "Date"))
MnWIFiles$Project <- "MnWI"
MnWIFiles$Directory <- "F:/Mn exposure/Mn WI/MnWI raw data"

# Progestin files --------------------------------------------------------
setwd("F:/Progestin")
load("Progestin metadata.RData")

ProgFiles <- Files
ProgFiles <- join(ProgFiles, Samples[, c("SampleID", "Subject", "Visit")], 
                  by = "SampleID", type = "left")
ProgFiles <- plyr::rename(ProgFiles, c("DateTime" = "Date",
                                       "Visit" = "Special"))

rm(Files, Samples, Subjects)

ProgFiles$Project <- "Progestin"
ProgFiles$Directory <- "F:/Progestin/Progestin raw data"


# SCOR DIG files ---------------------------------------------------
setwd("F:/SCOR DIG")
load("SCOR DIG metadata.RData")

DIGFiles <- Meta.all
DIGFiles <- plyr::rename(DIGFiles, c("DateTime" = "Date",
                                     "Preg" = "Special"))

DIGFiles$Project <- "SCOR DIG"
DIGFiles$Directory <- "F:/SCOR DIG/SCOR DIG raw data"

rm(Meta.all, Meta.clin, PK, PKType, Samples)


# SCOR MDZ files -----------------------------------

SCORDir <- c("G:/Data/Metabolomics/Laura/SCOR project/20110620 SCOR plasma/20130212 SCOR EposP VI", # EposP
             "G:/Data/Metabolomics/Laura/SCOR project/20120202 SCOR plasma and urine ESI neg", # EnegP
             "G:/Data/Metabolomics/Laura/SCOR project/20110211 SCOR subjects 1-6", # EposU #1
             "G:/Data/Metabolomics/Laura/SCOR project/20110224 SCOR subjects 7-15", # EposU #2
             "G:/Data/Metabolomics/Laura/SCOR project/20120202 SCOR plasma and urine ESI neg") # EnegU

setwd("C:/Users/Laura/Documents/SCOR project")

SCORmeta <- read.xlsx("all SCOR samples.xlsx", sheetName = "ClinicalSampleRuns")
SCORmeta <- SCORmeta[, c("Sample.name", "Subject", "Matrix", 
                         "X3rd.trimester.or.post.partum",
                         "time.point", "ESI..file.name",
                         "ESI..date.sample.ran", "APCI.date.sample.ran", 
                         "APCI.file.name", "ESI..date.sample.ran.1", 
                         "ESI..file.name.1")]

SCORmeta <- plyr::rename(SCORmeta, c("Sample.name" = "SampleID",
                                     "X3rd.trimester.or.post.partum" = "Special",
                                     "time.point" = "TimePoint",
                                     "ESI..date.sample.ran" = "Date.Epos",
                                     "ESI..file.name" = "File.Epos",
                                     "APCI.date.sample.ran" = "Date.Apos",
                                     "APCI.file.name" = "File.Apos",
                                     "ESI..date.sample.ran.1" = "Date.Eneg",
                                     "ESI..file.name.1" = "File.Eneg"))
names(SCORmeta)

# For subsequent steps, can't have identical SampleID on same day, mode, and 
# matrix. Fixing that for QC samples and just removing Master QC and blanks.
SCORmeta <- SCORmeta[!(SCORmeta$Subject %in%  c("blank", "Master QC", 
                                                "master QC", "water")), ]
SCORmeta$SampleID <- as.character(SCORmeta$SampleID)
SCORmeta$SampleID[SCORmeta$Subject == "QC"] <- 
      paste("QC", 1:length(SCORmeta$SampleID[SCORmeta$Subject == "QC"]))

# Getting data into best format. Making dates character for easier back conversion
SCORmeta$Date.Epos <- as.character(SCORmeta$Date.Epos)
SCORmeta$Date.Eneg <- as.character(SCORmeta$Date.Eneg)
SCORmeta$Date.Apos <- as.character(SCORmeta$Date.Apos)

SCORFiles <- gather(data = SCORmeta, Key, Value, -Subject, -Matrix, 
                    -Special, -SampleID, -TimePoint) %>%
      separate(Key, c("Type", "Mode"), "\\.") 
SCORFiles <- SCORFiles[complete.cases(SCORFiles$Value), ]
SCORFiles <- spread(SCORFiles, Type, Value)
SCORFiles$Date <- ymd(SCORFiles$Date)

# Adding directory. Note that I haven't bothered adding the Apos directory
# since it seems unlikely we'll use those data any time soon. 
SCORFiles$Directory[SCORFiles$Matrix == "plasma" & 
                          SCORFiles$Mode == "Epos"] <- SCORDir[1]

SCORFiles$Directory[SCORFiles$Matrix == "plasma" & 
                          SCORFiles$Mode == "Eneg"] <- SCORDir[2]

SCORFiles$Directory[SCORFiles$Matrix == "urine" & 
                          SCORFiles$Mode == "Epos" &
                          SCORFiles$Date == ymd("2011-02-11")] <- SCORDir[3]

SCORFiles$Directory[SCORFiles$Matrix == "urine" & 
                          SCORFiles$Mode == "Epos" &
                          SCORFiles$Date == ymd("2011-02-23")] <- SCORDir[4]

SCORFiles$Directory[SCORFiles$Matrix == "urine" & 
                          SCORFiles$Mode == "Eneg"] <- SCORDir[5]

SCORFiles$Project <- "SCOR MDZ"
SCORFiles$SampType <- "clinical"
SCORFiles$SampType[SCORFiles$Subject == "QC"] <- "QC"

SCORFiles <- subset(SCORFiles, complete.cases(File) & File != "has not run")

rm(SCORmeta, SCORDir)


# SCOR MDZ fragmentation files ----------------------------------------------
# 10/27/14
SCOR20141027 <- data.frame(Project = "SCOR MDZ fragment", 
                           Directory = "F:/SCOR/20141027 SCOR MDZ fragmentation",
                           File = c("20141027 SCOR MDZ EposP74 S01 PP a",
                                    "20141027 SCOR MDZ EposP74 S01 PP b",
                                    "20141027 SCOR MDZ EposP74 pooled A",
                                    "20141027 SCOR MDZ EposU76 QC urine A",
                                    "20141027 SCOR MDZ EposU76 QC urine B",
                                    "20141027 SCOR MDZ EposU76 S07 PP 6to24",
                                    "20141027 SCOR MDZ EnegU S10 PP 6to24"), 
                           SampleID = c(rep("S01 PP", 2),
                                        "pooled",
                                        rep("QC", 2),
                                        "S07 PP C", 
                                        "S10 PP C"),
                           Matrix = c(rep("plasma", 3), rep("urine", 4)),
                           Special = c("PP", "PP", "PP", "QC", "QC", "PP", "PP"),
                           Date = ymd("2014-10-27"),
                           Mode = c(rep("Epos", 6), "Eneg"))
# 1/3/15
setwd("C:/Users/Laura/Documents/SCOR project")
SCOR20150103 <- read.xlsx("Fragmentation plans - edited.xlsx", 
                          sheetName = "worklist")
SCOR20150103$Project <- "SCOR MDZ fragment"
SCOR20150103$Directory <- "F:/SCOR/20150103 SCOR MDZ fragmentation"
SCOR20150103$Date <- ymd("20150103")
SCOR20150103 <- SCOR20150103[SCOR20150103$MSn == "MS", ]

# Removing bad injection
SCOR20150103 <- SCOR20150103[SCOR20150103$File != 
                                   "20150103 SCOR MDZ EposP SFN 2214 310 2", ]


# 1/26/15
setwd("C:/Users/Laura/Documents/SCOR project/Fragmentation results and data files/20150126 SCOR MDZ fragmentation")
SCOR20150126 <- read.xlsx("20150126 SCOR MDZ fragmentation worklist.xlsx", 
                          sheetName = "worklist")
SCOR20150126 <- SCOR20150126[which(SCOR20150126$MSn == "MS"), ]
SCOR20150126$Project <- "SCOR MDZ fragment"
SCOR20150126$Directory <- "F:/SCOR/20150126 SCOR MDZ fragmentation"
SCOR20150126$Date <- ymd("20150126")
SCOR20150126$SampleID <- SCOR20150126$Sample

SCORFragmentFiles <- rbind.fill(SCOR20141027, 
                                SCOR20150103,
                                SCOR20150126)

rm(SCOR20141027, SCOR20150103, SCOR20150126)



# SFN files ---------------------------------------------
setwd("F:/Sulforaphane/All SFN mzdata files")
load("SFN metadata.RData")

SFNFiles <- subset(Files, Exists = TRUE)
SFNFiles$Directory <- "F:/Sulforaphane/All SFN mzdata files"
SFNFiles$Project <- "SFN"
SFNFiles <- join(SFNFiles, Samples, by = "SampleID", type = "left")
SFNFiles$Special <- SFNFiles$InductStat
SFNFiles$Date <- mdy(SFNFiles$Date)

rm(Files, Samples)



# Whole blood vs plasma data ----------------------------------------
setwd("F:/Whole blood vs plasma")
WBFiles <- read.csv("20140929 ESI+ whole blood IS peaks.csv", skip = 1)

names(WBFiles)
WBFiles <- plyr::rename(WBFiles, c("Sample" = "SampleID"))

# Keeping only the columns I want
WBFiles <- WBFiles[, c("SampleID", "Matrix", "File", "Date")]

# Only keeping single replicate of each sample and plasma only. Adding column
# for ionization mode.
WBFiles <- WBFiles[WBFiles$Matrix == "plasma", ]
WBFiles <- WBFiles[!(duplicated(WBFiles$SampleID)), ]
WBFiles$Mode <- "Epos"

WBFiles$Directory <- "F:/Whole blood vs plasma/Whole blood ESI+"

# Loading the Eneg data for the same day
WBFiles.Eneg <- read.csv("20140930 ESI- whole blood IS peaks.csv", skip = 1)

names(WBFiles.Eneg)
WBFiles.Eneg <- plyr::rename(WBFiles.Eneg, c("Sample" = "SampleID"))

# Keeping only the columns I want
WBFiles.Eneg <- WBFiles.Eneg[, c("SampleID", "Matrix", "File", "Date")]

# Only keeping single replicate of each sample and plasma only. Adding column
# for ionization mode.
WBFiles.Eneg <- WBFiles.Eneg[WBFiles.Eneg$Matrix == "plasma", ]
WBFiles.Eneg <- WBFiles.Eneg[!(duplicated(WBFiles.Eneg$SampleID)), ]
WBFiles.Eneg$Mode <- "Eneg"
WBFiles.Eneg$Directory <- "F:/Whole blood vs plasma/Whole blood ESI-"

WBFiles <- rbind(WBFiles, WBFiles.Eneg)
WBFiles$Project <- "Whole blood vs plasma"
WBFiles$Date <- mdy_hm(WBFiles$Date)

# Removing .d from file name.
WBFiles$File <- sub(".d$", "", as.character(WBFiles$File))

rm(WBFiles.Eneg)


# Putting all the files together ---------------------------------------
Files <- rbind.fill(BusulfFiles, CHC2Files, DEXFiles, DIGFiles, IVPOFiles, 
                    MetopFiles, MnPSFiles, MnWIFiles, ProgFiles, SCORFiles, 
                    SCORFragmentFiles, SFNFiles, WBFiles)

Files <- subset(Files, , c("Project", "SampleID", "Subject", "Special", 
                           "SampType", "File", "Directory", 
                           "Mode", "Matrix",
                           "Date"))

# Removing any missing files or directories. 
Files <- Files[complete.cases(Files$File) & complete.cases(Files$Directory), ]

# Removing trailing or leading white spaces from the file names. They sometimes
# got in there b/c of the way I had set up my Excel files to make the file name.
Files$File <- str_trim(Files$File)

# Adding dates for files that are missing them.
Files$Date[is.na(Files$Date)] <- ymd(str_sub(Files$File[is.na(Files$Date)], 
                                             1, 8))

# Checking numbers
ddply(Files, c("Mode", "Matrix"), function(x) nrow(x))
ddply(Files, c("Project", "Mode"), function(x) nrow(x))

# Saving the data ------------------------------------------
setwd(MainDir)
write.csv(Files, "Metabolomics files.csv", row.names = FALSE)
save(Files, file = "Metabolomics files.RData")

