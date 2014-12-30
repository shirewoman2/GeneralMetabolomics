# Master QC tracking

# This script looks at how internal standard peaks in Master QC samples have
# changed over time. 

# Housekeeping ===================================
library(XLConnect)
library(plyr)
library(reshape2)
library(tidyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(stringr)

ThemeLaura <- function (base_size = 12, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace% 
            theme(
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(colour = "black"),
                  axis.title.y = element_text(colour = "black", angle=0),
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


# Mn Puget Sound =======================================================

# 20141117 Epos run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/Mn exposure project")
R20141117 <- loadWorkbook("20141117 ESI pos CEEH.xlsx")
setMissingValue(R20141117, value = "")
R20141117 <- readWorksheet(R20141117, sheet = "Sheet1", startRow = 2)

names(R20141117)

R20141117 <- plyr::rename(R20141117, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

R20141117 <- R20141117[, -c(4, 8, 12, 16)]
IS <- c("MT", "RSG", "MDZ", "Prog")

names(R20141117)[4:15] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20141117$Matrix <- "urine"

# Reshaping data to long format
R20141117 <- R20141117 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20141117$SampType <- "clinical"
R20141117$SampType[str_detect(R20141117$Sample, "*QC*")] <- "QC"
R20141117$SampType[str_detect(R20141117$Sample, "*MQC*")] <- "Master QC"

# Adding a column for the run
R20141117$Run <- "Mn Puget Sound 20141117"

# Adding a column for the mode
R20141117$Mode <- "Epos"


# 20141118 EnegU run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/Mn exposure project")
R20141118 <- loadWorkbook("20141118 EnegU IS quant report.xlsx")
setMissingValue(R20141118, value = "")
R20141118 <- readWorksheet(R20141118, sheet = "Sheet1", startRow = 2)

names(R20141118)

R20141118 <- plyr::rename(R20141118, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("Salicylic", "Pred", "Stearic")

R20141118 <- R20141118[, -c(4, 8, 12)]

names(R20141118)[4:12] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20141118$Matrix <- "urine"

# Reshaping data to long format
R20141118 <- R20141118 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20141118$SampType <- "clinical"
R20141118$SampType[str_detect(R20141118$Sample, "*QC*")] <- "QC"
R20141118$SampType[str_detect(R20141118$Sample, "*MQC*")] <- "Master QC"

# Adding a column for the run
R20141118$Run <- "Mn Puget Sound 20141118"

# Adding a column for the mode
R20141118$Mode <- "Eneg"



# Whole blood =========================================

# 20140930 Eneg run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/LCMS metabolomics/Whole blood")
R20140930 <- loadWorkbook("2014030 ESI- whole blood quant report.xlsx")
setMissingValue(R20140930, value = "")
R20140930 <- readWorksheet(R20140930, sheet = "Sheet1", startRow = 1)

names(R20140930)

R20140930 <- plyr::rename(R20140930, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("Salicylic", "Pred", "Stearic")

names(R20140930)[5:13] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20140930 <- R20140930 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140930$SampType <- "clinical"
R20140930$SampType[str_detect(R20140930$Sample, "*QC*")] <- "QC"
R20140930$SampType[str_detect(R20140930$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20140930$Run <- "Whole blood 20140930"

# Adding a column for the mode
R20140930$Mode <- "Eneg"


# 20141014 Epos run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/LCMS metabolomics/Whole blood")
R20141014 <- loadWorkbook("20141014 ESI+ whole blood quant report.xlsx")
setMissingValue(R20141014, value = "")
R20141014 <- readWorksheet(R20141014, sheet = "Sheet1", startRow = 2)

names(R20141014)

R20141014 <- plyr::rename(R20141014, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

R20141014 <- R20141014[, -c(4, 8, 12, 16)]
IS <- c("MT", "RSG", "MDZ", "Prog")

names(R20141014)[4:15] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20141014$Matrix <- "plasma"
R20141014$Matrix[str_detect(R20141014$Sample, "*blood*")] <- "whole blood"

# Reshaping data to long format
R20141014 <- R20141014 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20141014$SampType <- "clinical"
R20141014$SampType[str_detect(R20141014$Sample, "*QC*")] <- "QC"
R20141014$SampType[str_detect(R20141014$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20141014$Run <- "Whole blood 20141014"

# Adding a column for the mode
R20141014$Mode <- "Epos"


# Busulfan ========================================================

# 20141017 Eneg run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/Busulfan project/Busulfan retrospective study")
R20141017 <- loadWorkbook("Busulfan EnegP IS peaks.xlsx")
setMissingValue(R20141017, value = "")
R20141017 <- readWorksheet(R20141017, sheet = "Sheet1", startRow = 2)

names(R20141017)

R20141017 <- plyr::rename(R20141017, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

R20141017 <- R20141017[, -c(4, 8, 12)]

IS <- c("Salicylic", "Pred", "Stearic")

names(R20141017)[4:12] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20141017$Matrix <- "plasma"
R20141017$Matrix[str_sub(R20141017$File, 9, 9) %in% c("M", "U")] <- "urine"
R20141017$Matrix[str_sub(R20141017$File, 9, 9) == "B"] <- "water"

# Reshaping data to long format
R20141017 <- R20141017 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20141017$SampType <- "clinical"
R20141017$SampType[str_detect(R20141017$Sample, "*QC*")] <- "QC"
R20141017$SampType[str_detect(R20141017$Sample, "*Master*")] <- "Master QC"

# Adding a column for the run
R20141017$Run <- "Busulf 20141017"

# Adding a column for the mode
R20141017$Mode <- "Eneg"


# 20141015 Epos run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/Busulfan project/Busulfan retrospective study")
R20141015 <- loadWorkbook("Busulfan EposP IS peaks.xlsx")
setMissingValue(R20141015, value = "")
R20141015 <- readWorksheet(R20141015, sheet = "Sheet1", startRow = 2)

names(R20141015)

R20141015 <- plyr::rename(R20141015, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

R20141015 <- R20141015[, -c(4, 8, 12, 16)]

IS <- c("MT", "RSG", "MDZ", "Prog")

names(R20141015)[4:15] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20141015$Matrix <- "plasma"
R20141015$Matrix[str_sub(R20141015$File, 9, 9) %in% c("M", "U")] <- "urine"
R20141015$Matrix[str_sub(R20141015$File, 9, 9) == "B"] <- "water"

# Reshaping data to long format
R20141015 <- R20141015 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20141015$SampType <- "clinical"
R20141015$SampType[str_detect(R20141015$Sample, "*QC*")] <- "QC"
R20141015$SampType[str_detect(R20141015$Sample, "*Master*")] <- "Master QC"

# Adding a column for the run
R20141015$Run <- "Busulf 20141015"

# Adding a column for the mode
R20141015$Mode <- "Epos"


# SCOR MDZ ===================================================

# 20110211 EposU run --------------------------------------------
setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR MDZ/ESI+/Metabolomics urine samples/20110211 SCOR subjects 1-6")
R20110211 <- loadWorkbook("20110211 SCOR subjects 1-6.xlsx")
setMissingValue(R20110211, value = "")
R20110211 <- readWorksheet(R20110211, sheet = "data", startRow = 1,
                           endRow = 45)
names(R20110211)
R20110211 <- plyr::rename(R20110211, 
                          c("Sample.type" = "SampType",
                            "file.name" = "File",
                            "X5.MT.RT" = "MT.RT",
                            "X5.MT.Area" = "MT.Area",
                            "PROG.Area" = "Prog.Area",
                            "PROG.RT" = "Prog.RT"))

R20110211$DateTime <- ymd("20110211")
R20110211$Sample <- paste(R20110211$Subject, 
                          R20110211$X3rd.trimester.or.post.partum,
                          R20110211$time.point)

IS <- c("MT", "RSG", "MDZ", "Prog")

# Switching the values to drop "urine" from QC and Master QC sample types.
R20110211$SampType[R20110211$SampType == "QC urine"] <- "QC"
R20110211$SampType[R20110211$SampType == "Master QC urine"] <- "Master QC"

# Changing sample type of "urine" to "clinical".
R20110211$SampType[R20110211$SampType == "urine"] <- "clinical"

# Removing blank and water injections.
R20110211 <- R20110211[R20110211$Matrix == "urine", ]

# Keeping only required columns.
R20110211 <- R20110211[, c("Sample", "File", "DateTime", "Matrix", "SampType",
                           paste(rep(IS, each = 2), sep=".", 
                                 rep(c("RT", "Area"), length(IS))))]

# Reshaping data to long format
R20110211 <- R20110211 %>% gather(key, Value, -Sample, -File, -DateTime,
                                  -Matrix, -SampType) %>%
      separate(key, c("IS", "Type"), "\\.") %>% spread(Type, Value)

# Adding a column for the run
R20110211$Run <- "SCOR MDZ 20110211"

# Adding a column for the mode
R20110211$Mode <- "Epos"


# 20110223 EposU run --------------------------------------------
setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR MDZ/ESI+/Metabolomics urine samples/20110223 SCOR subjects 7-15")
R20110223 <- loadWorkbook("20110223 SCOR subjects 7-15.xlsx")
setMissingValue(R20110223, value = "")
R20110223 <- readWorksheet(R20110223, sheet = "20110223", startRow = 1,
                           endRow = 63)
names(R20110223)
R20110223 <- plyr::rename(R20110223, 
                          c("file.name" = "File",
                            "X5.MT.RT" = "MT.RT",
                            "X5.MT.Area" = "MT.Area",
                            "ROS.Area" = "RSG.Area",
                            "ROS.RT" = "RSG.RT"))

R20110223$DateTime <- ymd("20110223")
R20110223$Sample <- paste(R20110223$Subject, 
                          R20110223$X3rd.trimester.or.post.partum,
                          R20110223$time.point)

IS <- c("MT", "RSG", "MDZ", "Prog")

# Adding sample type column
R20110223$SampType <- "clinical"
R20110223$SampType[R20110223$Subject == "Master QC"] <- "Master QC"
R20110223$SampType[R20110223$Subject == "QC"] <- "QC"

# Removing blank and water injections.
R20110223 <- R20110223[!(R20110223$Subject %in% c("blank", "water")), ]

# Fixing the matrix when it was "urine-blank"
R20110223$Matrix <- "urine"

# Keeping only required columns.
R20110223 <- R20110223[, c("Sample", "File", "DateTime", "Matrix", "SampType",
                           paste(rep(IS, each = 2), sep=".", 
                                 rep(c("RT", "Area"), length(IS))))]

# Reshaping data to long format
R20110223 <- R20110223 %>% gather(key, Value, -Sample, -File, -DateTime,
                                  -Matrix, -SampType) %>%
      separate(key, c("IS", "Type"), "\\.") %>% spread(Type, Value)

# Adding a column for the run
R20110223$Run <- "SCOR MDZ 20110223"

# Adding a column for the mode
R20110223$Mode <- "Epos"


# 20110620 EposP run --------------------------------------------
setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR MDZ/ESI+/Metabolomics plasma samples/20110620 SCOR plasma samples run #3")
R20110620 <- loadWorkbook("20110620 IS peak areas.xlsx")
setMissingValue(R20110620, value = "")
R20110620 <- readWorksheet(R20110620, sheet = "IS peak areas", startRow = 1,
                           endRow = 40)
names(R20110620)
R20110620 <- plyr::rename(R20110620, 
                          c("Data.File" = "File",
                            "Acq..Date.Time" = "DateTime",
                            "Sample.Type" = "SampType",
                            "X5.MT.RT" = "MT.RT",
                            "X5.MT.Height" = "MT.Height",
                            "X5.MT.Area" = "MT.Area"))

R20110620$Sample <- paste(R20110620$Subject, 
                          R20110620$Gestation.state)

# Fixing the sample type when it was "QC plasma" and when it was just "plasma"
R20110620$SampType[R20110620$SampType == "QC plasma"] <- "QC"
R20110620$SampType[R20110620$SampType == "plasma"] <- "clinical"

# Removing blank and water injections.
R20110620 <- R20110620[R20110620$SampType != "water", ]

# Removing samples that were for adjusting for matrix
R20110620 <- R20110620[! str_detect(R20110620$Name, "adjust"), ]

# Adding a column for the matrix
R20110620$Matrix <- "plasma"
R20110620$Matrix[R20110620$SampType == "Master QC"] <- "urine"

IS <- c("MT", "RSG", "MDZ", "Prog")


# Keeping only required columns.
R20110620 <- R20110620[, c("Sample", "File", "DateTime", "Matrix", "SampType",
                           paste(rep(IS, each = 2), sep=".", 
                                 rep(c("RT", "Area"), length(IS))))]

# Reshaping data to long format
R20110620 <- R20110620 %>% gather(key, Value, -Sample, -File, -DateTime,
                                  -Matrix, -SampType) %>%
      separate(key, c("IS", "Type"), "\\.") %>% spread(Type, Value)

# Adding a column for the run
R20110620$Run <- "SCOR MDZ 20110620"

# Adding a column for the mode
R20110620$Mode <- "Epos"




# 20130201 EnegP run =======================================
setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR MDZ/ESI-")
R20130201 <- loadWorkbook("20130201 SCOR EnegP IS peak areas and RTs.xlsx")
setMissingValue(R20130201, value = "")
R20130201 <- readWorksheet(R20130201, sheet = "Quant report", startRow = 2)

R20130201$Sample <- NULL

R20130201 <- R20130201[, -c(1,2,5,6)]
IS <- c("Salicylic", "Pred", "Stearic")

names(R20130201) <- c("Sample", "File", "DateTime", 
                      paste(rep(c("RT", "Area", "Height"), length(IS)),
                            rep(IS, each = 3), sep="."))

# Reshaping data to long format
R20130201 <- R20130201 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20130201$SampType <- "clinical"
R20130201$SampType[str_detect(R20130201$Sample, "*QC*")] <- "QC"
R20130201$SampType[str_detect(R20130201$Sample, "*aster*")] <- "Master QC"

# Adding a column for the matrix
R20130201$Matrix <- "plasma"

# Adding a column for the run
R20130201$Run <- "SCOR MDZ 20130201"

# Adding a column for the mode
R20130201$Mode <- "Eneg"

# 20120202 EnegU and EnegP run ======================================
setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR MDZ/ESI-/20120202 SCOR urine and plasma ESI-")
R20120202 <- loadWorkbook("ESI- IS quant results with smoothing.xlsx")
setMissingValue(R20120202, value = "")
R20120202 <- readWorksheet(R20120202, sheet = "peak areas", startRow = 1)

names(R20120202)

R20120202$Subject <- NULL
R20120202$Sample.type <- NULL

R20120202 <- plyr::rename(R20120202, c("Name" = "Sample",
                                       "Data.File" = "File", 
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("Salicylic", "Pred", "Stearic")

names(R20120202) <- c("Sample", "Matrix", "File", "DateTime", 
                      paste(rep(c("RT", "Area"), length(IS)),
                            rep(IS, each = 2), sep="."))

# Reshaping data to long format
R20120202 <- R20120202 %>% gather(key, Value, -Sample, -Matrix, -File, 
                                  -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding column for height
R20120202$Height <- NA

# Adding a column for the type of sample
R20120202$SampType <- "clinical"
R20120202$SampType[str_detect(R20120202$Sample, "*QC*")] <- "QC"
R20120202$SampType[str_detect(R20120202$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20120202$Run <- "SCOR MDZ 20120202"

# Adding a column for the mode
R20120202$Mode <- "Eneg"


# CHC2 =========================================================

# 20121008 EnegU run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/CHC2 project")
R20121008 <- loadWorkbook("EnegU IS peaks quant report.xlsx")
setMissingValue(R20121008, value = "")
R20121008 <- readWorksheet(R20121008, sheet = "Sheet1", startRow = 1)

names(R20121008)

R20121008 <- plyr::rename(R20121008, c("Date" = "DateTime"))

IS <- c("Salicylic", "Pred", "Stearic")

names(R20121008)[4:12] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20121008 <- R20121008 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20121008$SampType <- "clinical"
R20121008$SampType[str_detect(R20121008$Sample, "*QC*")] <- "QC"
R20121008$SampType[str_detect(R20121008$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20121008$Run <- "CHC2 20121008"

# Adding a column for the mode
R20121008$Mode <- "Eneg"

# Adding a column for the matrix
R20121008$Matrix <- "urine"


# 20110504 EposU run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/CHC2 project")
R20110504 <- loadWorkbook("EposU IS peaks quant report.xlsx")
setMissingValue(R20110504, value = "")
R20110504 <- readWorksheet(R20110504, sheet = "quant report", startRow = 1)

names(R20110504)

R20110504 <- plyr::rename(R20110504, c("Date" = "DateTime"))

IS <- c("MT", "RSG", "MDZ", "Prog")

names(R20110504)[4:15] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20110504 <- R20110504 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20110504$SampType <- "clinical"
R20110504$SampType[str_detect(R20110504$Sample, "*QC*")] <- "QC"
R20110504$SampType[str_detect(R20110504$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20110504$Run <- "CHC2 20110504"

# Adding a column for the mode
R20110504$Mode <- "Epos"

# Adding a column for the matrix
R20110504$Matrix <- "urine"



# CYP2D6 Metoprolol ========================================================
# 20140630 EposP run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/CYP2D6 metoprolol")

R20140630 <- loadWorkbook("Metoprolol EposP IS peaks.xlsx")
setMissingValue(R20140630, value = "")
R20140630 <- readWorksheet(R20140630, sheet = "Sheet1", startRow = 2)

names(R20140630)

R20140630 <- plyr::rename(R20140630, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("MT", "RSG", "MDZ", "Prog")

names(R20140630)[4:15] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20140630 <- R20140630 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140630$SampType <- "clinical"
R20140630$SampType[str_detect(R20140630$Sample, "*QC*")] <- "QC"
R20140630$SampType[str_detect(R20140630$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20140630$Run <- "CYP2D6 Metoprolol 20140630"

# Adding a column for the mode
R20140630$Mode <- "Epos"

# Adding a column for the matrix
R20140630$Matrix <- "plasma"
R20140630$Matrix[R20140630$SampType == "Master QC"] <- "urine"


# 20140714 EposU run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/CYP2D6 metoprolol")
R20140714 <- loadWorkbook("Metoprolol EposU IS peaks.xlsx")
setMissingValue(R20140714, value = "")
R20140714 <- readWorksheet(R20140714, sheet = "Sheet1", startRow = 2)

names(R20140714)

R20140714 <- plyr::rename(R20140714, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("MT", "RSG", "MDZ", "Prog")

names(R20140714)[4:15] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20140714 <- R20140714 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140714$SampType <- "clinical"
R20140714$SampType[str_detect(R20140714$Sample, "*QC*")] <- "QC"
R20140714$SampType[str_detect(R20140714$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20140714$Run <- "CYP2D6 Metoprolol 20140714"

# Adding a column for the mode
R20140714$Mode <- "Epos"

# Adding a column for the matrix
R20140714$Matrix <- "urine"
# R20140714$Matrix[R20140714$SampType == "Master QC"] <- "urine"


# 20140716 EnegU run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/CYP2D6 metoprolol")
R20140716 <- loadWorkbook("Metoprolol EnegU IS peaks.xlsx")
setMissingValue(R20140716, value = "")
R20140716 <- readWorksheet(R20140716, sheet = "Sheet1", startRow = 2)

names(R20140716)

R20140716 <- plyr::rename(R20140716, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("Salicylic", "Pred", "Stearic")

names(R20140716)[4:12] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20140716 <- R20140716 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140716$SampType <- "clinical"
R20140716$SampType[str_detect(R20140716$Sample, "*QC*")] <- "QC"
R20140716$SampType[str_detect(R20140716$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20140716$Run <- "CYP2D6 Metoprolol 20140716"

# Adding a column for the mode
R20140716$Mode <- "Eneg"

# Adding a column for the matrix
R20140716$Matrix <- "urine"
# R20140716$Matrix[R20140716$SampType == "Master QC"] <- "urine"


# 20140701 EnegP run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/CYP2D6 metoprolol")
R20140701 <- loadWorkbook("Metoprolol EnegP IS peaks.xlsx")
setMissingValue(R20140701, value = "")
R20140701 <- readWorksheet(R20140701, sheet = "Sheet1", startRow = 2)

names(R20140701)

R20140701 <- plyr::rename(R20140701, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

IS <- c("Salicylic", "Pred", "Stearic")

names(R20140701)[4:12] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Reshaping data to long format
R20140701 <- R20140701 %>% gather(key, Value, -Sample, -File, -DateTime) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140701$SampType <- "clinical"
R20140701$SampType[str_detect(R20140701$Sample, "*QC*")] <- "QC"
R20140701$SampType[str_detect(R20140701$Sample, "*aster*")] <- "Master QC"

# Adding a column for the run
R20140701$Run <- "CYP2D6 Metoprolol 20140701"

# Adding a column for the mode
R20140701$Mode <- "Eneg"

# Adding a column for the matrix
R20140701$Matrix <- "plasma"
R20140701$Matrix[R20140701$SampType == "Master QC"] <- "urine"


# SCOR DIG ===================================================
# 20140327 Epos run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR DIG")
R20140327 <- loadWorkbook("20140327 SCOR DIG Epos IS peaks.xlsx")
setMissingValue(R20140327, value = "")
R20140327 <- readWorksheet(R20140327, sheet = "Sheet1", startRow = 2)

names(R20140327)

R20140327 <- plyr::rename(R20140327, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

R20140327 <- R20140327[, -c(4, 8, 12, 16, 20)]
IS <- c("MT", "RSG", "MDZ", "DIG", "Prog")

names(R20140327)[4:18] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20140327$Matrix <- "plasma"
R20140327$Matrix[str_sub(R20140327$File, 9, 9) %in% c("M", "U")] <- "urine"
R20140327$Matrix[str_sub(R20140327$File, 9, 9) == "B"] <- "water"

# Reshaping data to long format
R20140327 <- R20140327 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140327$SampType <- "clinical"
R20140327$SampType[str_detect(R20140327$Sample, "*QC*")] <- "QC"
R20140327$SampType[str_detect(R20140327$Sample, "*Master*")] <- "Master QC"

# Adding a column for the run
R20140327$Run <- "SCOR DIG 20140327"

# Adding a column for the mode
R20140327$Mode <- "Epos"

# Removing DIG b/c not useful for comparisons
R20140327 <- R20140327[R20140327$IS != "DIG", ]

# 20140222 Epos run ==============================================

setwd("D:/Users/Laura/Documents/Work/Lin Lab/SCOR DIG")
R20140222 <- loadWorkbook("20140327 SCOR DIG Eneg IS peaks.xlsx")
setMissingValue(R20140222, value = "")
R20140222 <- readWorksheet(R20140222, sheet = "Sheet1", startRow = 2)

names(R20140222)

R20140222 <- plyr::rename(R20140222, c("Name" = "Sample",
                                       "Data.File" = "File",
                                       "Acq..Date.Time" = "DateTime"))

R20140222 <- R20140222[, -c(4, 8, 12)]

IS <- c("Salicylic", "Pred", "Stearic")

names(R20140222)[4:12] <- paste(rep(c("RT", "Area", "Height"), length(IS)),
                                rep(IS, each = 3), sep=".")

# Adding a column for the matrix
R20140222$Matrix <- "plasma"
R20140222$Matrix[str_sub(R20140222$File, 9, 9) %in% c("M", "U")] <- "urine"
R20140222$Matrix[str_sub(R20140222$File, 9, 9) == "B"] <- "water"

# Reshaping data to long format
R20140222 <- R20140222 %>% gather(key, Value, -Sample, -File, -DateTime, 
                                  -Matrix) %>%
      separate(key, c("Type", "IS"), "\\.") %>% spread(Type, Value)

# Adding a column for the type of sample
R20140222$SampType <- "clinical"
R20140222$SampType[str_detect(R20140222$Sample, "*QC*")] <- "QC"
R20140222$SampType[str_detect(R20140222$Sample, "*Master*")] <- "Master QC"

# Adding a column for the run
R20140222$Run <- "SCOR DIG 20140222"

# Adding a column for the mode
R20140222$Mode <- "Eneg"




# Putting all the data.frames together ==================================
setwd("D:/Users/Laura/Documents/Work/Lin Lab/LCMS metabolomics")

DF <- list(R20110211, R20110223, R20110620, R20141017, 
           R20110504, R20120202, R20121008, R20130201,
           R20140222, R20140327, R20140630, R20140701,
           R20140714, R20140716, R20140930, R20141014,
           R20141015, R20141117, R20141118)

# Some of my data.frames had the class wrong for DateTime. Checking that.
llply(DF, function(x) class(x$DateTime))
# Good. Problem fixed.

Data <- rbind.fill(DF)
Data <- Data[Data$Matrix != "water" & Data$Matrix != "whole blood", ]

Runs <- data.frame(Run = Data$Run,
                   Date = str_sub(Data$Run, nchar(Data$Run)-7, nchar(Data$Run)))
Runs$Date <- as.numeric(as.character(Runs$Date))
Runs <- arrange(Runs, Date)
Runs <- unique(Runs)

Runs$PreppedBy <- "LS"
Runs$PreppedBy[Runs$Run %in% c("CHC2 20110504", "CHC2 20121008",
                               "Busulf 20141015")] <- "TS"

Data$Run <- factor(Data$Run, levels = Runs$Run)

# Replacing "NA" in peak area column with 0 since those peaks were not 
# detected because their peak areas and heights were 0, not just missing.
Data$Area[is.na(Data$Area)] <- 0
# Data$Height[is.na(Data$Height)] <- 0  # Sometimes, I didn't collect data on 
#the height, only the area, so I don't want to set these values as zero.

Data <- join(Data, Runs, by = "Run") 

# Final, complete data
Data <- Data[, c("Sample", "File", "DateTime", "Matrix", "SampType", 
                 "IS", "Area", "Height", "RT", "Run", "Mode", "Date", 
                 "PreppedBy")]
Data <- arrange(Data, Date, DateTime, File)
Data$Date <- as.factor(Data$Date)

save(Data, Runs, DF, file = "IS peak areas for several runs.RData")


# Epos IS peak areas ===================================================
Epos <- Data[Data$Mode == "Epos" & Data$SampType == "clinical", ]

MQC.pos <- Data[Data$Mode == "Epos" & Data$SampType == "Master QC" &
                      Data$DateTime < "2014-09-01", ]

windows()
ggplot(MQC.pos, aes(x = Date, y = Area, fill = Run)) + geom_boxplot() + 
      facet_wrap(~ IS, scales = "free") + ggtitle("Master QC ESI+") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("Master QC Epos boxplot of IS peak areas.png", height = 8, width = 14)

windows()
ggplot(Epos, aes(x = Run, y = Area, fill = Run)) + geom_boxplot() + 
      facet_grid(IS ~ Matrix, scale = "free") + 
      ggtitle("All Samples ESI+") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("All Epos clinical samples - boxplot of IS peak areas by matrix.png", 
       height = 12, width = 8)


# Eneg IS peak areas ===================================================
Eneg <- Data[Data$Mode == "Eneg" & Data$SampType == "clinical", ]

MQC.neg <- Data[Data$Mode == "Eneg" & Data$SampType == "Master QC" &
                      Data$DateTime < "2014-09-01", ]


ggplot(MQC.neg, aes(x = Date, y = Area, fill = Run)) + geom_boxplot() + 
      facet_wrap(~ IS, scales = "free") + ggtitle("Master QC ESI-") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("Master QC Eneg boxplot of IS peak areas.png", height = 8, width = 14)

windows()
ggplot(Eneg, aes(x = Run, y = Area, fill = Run)) + geom_boxplot() + 
      facet_grid(IS ~ Matrix, scale = "free") + 
      ggtitle("All clinical samples ESI-") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("All Eneg clinical samples - boxplot of IS peak areas by matrix.png", 
       height = 12, width = 8)


# Just checking recent runs that I prepped. Removing 5-MT b/c that's not even
# in the ACN sln any more. 
LSrecent <- Data[Data$PreppedBy == "LS" & 
                       Data$DateTime > "2014-01-01" &
                       Data$SampType != "Master QC" &
                       Data$IS != "MT", ]

LS.Epos <- LSrecent[LSrecent$Mode == "Epos", ]

ggplot(LS.Epos, aes(x = Run, y = Area, fill = Run)) + 
      geom_boxplot() +
      facet_grid(IS ~ Matrix, scale = "free") + 
      ggtitle("All recent LS samples ESI+") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("All recent Epos LS samples - boxplot of IS peak areas.png", 
       height = 8, width = 8)


LS.Eneg <- LSrecent[LSrecent$Mode == "Eneg", ]

ggplot(LS.Eneg, aes(x = Run, y = Area, fill = Run)) + 
      geom_boxplot() +
      facet_grid(IS ~ Matrix, scale = "free") + 
      ggtitle("All recent LS samples ESI-") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("All recent Eneg LS samples - boxplot of IS peak areas.png", 
       height = 8, width = 8)



# Checking RTs
MeanRT <- ddply(Data, c("Run", "IS", "Matrix"), function(x)
                Ave = mean(x$RT, na.rm = T))
MeanRT <- plyr::rename(MeanRT, c("V1" = "MeanRT"))

RT.Epos <- MeanRT[MeanRT$IS %in% c("MT", "MDZ", "RSG", "Prog"), ]
RT.Eneg <- MeanRT[MeanRT$IS %in% c("Salicylic", "Stearic", "Pred"), ]

windows()
ggplot(RT.Epos, aes(x = Run, y = MeanRT, fill = Matrix, color = Matrix)) +
      geom_boxplot() + 
      facet_wrap( ~ IS, scale = "free", ncol = 2) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("Boxplots comparing RTs for all Epos IS.png", 
       height = 8, width = 10)


ggplot(RT.Eneg, aes(x = Run, y = MeanRT, fill = Matrix, color = Matrix)) +
      geom_boxplot() + 
      facet_wrap( ~ IS, scale = "free", ncol = 2) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("Boxplots comparing RTs for all Eneg IS.png", 
       height = 8, width = 10)
