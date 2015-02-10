# This function generates a worklist for metabolomics samples where the sample
# order is randomized, there is a QC injection every 10 injections, and the 
# worklist starts and ends with a Master QC injection. 
# Input: 
#       1. Samples -- a data.frame containing columns for:
#             a. SampleID - a vector of the sample IDs as you would like them 
#             to show up in the file names. These must be unique.
#             b. Matrix -- the matrix of the sample
#       2. Date -- a number for how you'd like to represent today's 
#       date, which will be at the front of each file. No symbols allowed.
#       3. Matrix -- the sample matrix (this is set up to make a worklist for
#                     one matrix at a time. For multiple matrices,
#                     make multiple worklists.)
#      3. Project -- a brief abbreviation or code word for which project this
#       is for. Example: "MnPS" or "SCOR MDZ".
#       4. Qnum.start -- number of QC injections you'd like at the start of the
#       run. (Defaults to 3)
#       5. Qnum.end -- number of QC injections you'd like at the end of the run. 
#       (Defaults to 1.)
#       6. Mode -- a character string of which ionization modes you'd like to run. 
#       Defaults to c("Epos", "Eneg") but other options are "Apos" and "Aneg".

worklist <- function(Samples, Date, Matrix, FilePath,
                     Qnum.start = 3, Qnum.end = 1, 
                     MQnum.start = 1, MQnum.end = 1,
                     Mode = c("Epos", "Eneg")) {
      
      require(plyr)
      require(tidyr)
      require(dplyr)
      require(stringr)
      
      set.seed(98032)
      
      Samples <- data.frame(SampleID = SampleID, 
                            RandNum = rnorm(length(SampleID)))
      Samples <- arrange(Samples, RandNum)
      
      # Setting up QC injections
      Qnum.start <- Qnum.start-1
      Qnum <- Qnum.start + Qnum.end + round(nrow(Samples)/10)
      QCsamp <- data.frame(SampleID = paste0("QC", 1:Qnum))
      QCsamp$VialPos <- "P1-A2"
      
      # Setting up Master QC injections
      MQCsamp <- data.frame(SampleID = paste0("MQC", 1:(MQnum.start +
                                                              MQnum.end)))
      MQCsamp$VialPos <- "P1-A1"
      
      if ("Epos" %in% Mode) {
            Samples$File.Epos <- paste(
                  Date, Project, paste0("Epos", toupper(str_sub(Matrix, 1, 1))),
                                       Samples$SampleID)
            QCsamp$File.Epos <- paste(
                  Date, Project, paste0("Epos", toupper(str_sub(Matrix, 1, 1))), 
                                        QCsamp$SampleID)
            MQCsamp$File.Epos <- paste(
                  Date, Project, paste0("Epos", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$SampleID)
            
      }
      
      if ("Eneg" %in% Mode){
            Samples$File.Eneg <- paste(
                  Date, Project, paste0("Eneg", toupper(str_sub(Matrix, 1, 1))),
                                       Samples$SampleID)
            QCsamp$File.Eneg <- paste(
                  Date, Project, paste0("Eneg", toupper(str_sub(Matrix, 1, 1))),
                  QCsamp$SampleID)
            MQCsamp$File.Eneg <- paste(
                  Date, Project, paste0("Eneg", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$SampleID)
      }
      
      if ("Apos" %in% Mode){
            Samples$File.Apos <- paste(
                  Date, Project, paste0("Apos", toupper(str_sub(Matrix, 1, 1))),
                  Samples$SampleID)
            QCsamp$File.Apos <- paste(
                  Date, Project, paste0("Apos", toupper(str_sub(Matrix, 1, 1))),
                  QCsamp$SampleID)
            MQCsamp$File.Apos <- paste(
                  Date, Project, paste0("Apos", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$SampleID)
      }
      
      if ("Aneg" %in% Mode){
            Samples$File.Aneg <- paste(
                  Date, Project, paste0("Aneg", toupper(str_sub(Matrix, 1, 1))),
                  Samples$SampleID)
            QCsamp$File.Aneg <- paste(
                  Date, Project, paste0("Aneg", toupper(str_sub(Matrix, 1, 1))),
                  QCsamp$SampleID)
            MQCsamp$File.Aneg <- paste(
                  Date, Project, paste0("Aneg", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$SampleID)
      }
           
      # Setting the vial positions
      VialPos <- paste0(rep(c("P1-", "P2-"), each = 54), # 6 rows, 9 columns 
                        rep(paste0(rep(LETTERS[1:6], each = 9), c(1:9)), 2))
      
      Samples$VialPos <- VialPos[3:(nrow(Samples) + 2)]
      
      # Putting together the worklist
      Worklist <- list()
      
      for (i in 1:ceiling(nrow(Samples)/10)){
            Worklist[[i]] <- rbind.fill(QCsamp[i+Qnum.start, ],
                                   Samples[
                                         ((1+10*(i-1)):(10*i)), ])
      }
      
      if (nrow(Samples) %% 10 > 0) {
            Worklist[[ceiling(nrow(Samples)/10)]] <- rbind.fill(
                  QCsamp[Qnum-1, ],
                  Samples[((nrow(Samples) %/% 10)*10+1):nrow(Samples), ])
      }
      
      Worklist[[1]] <- rbind.fill(MQCsamp[1:MQnum.start, ], 
                                  QCsamp[1:Qnum.start, ],
                                  Worklist[[1]])

      Worklist[[ceiling(nrow(Samples)/10)]] <- 
            rbind.fill(Worklist[[ceiling(nrow(Samples)/10)]], 
                       QCsamp[(Qnum-Qnum.end+1):Qnum, ],
                  MQCsamp[(MQnum.start+1):(MQnum.start+MQnum.end), ])
      
      Worklist <- rbind.fill(Worklist)
      Worklist <- Worklist[, c("SampleID", "VialPos", 
                               names(Worklist)[names(Worklist) %in% 
                                                     paste0("File.", Mode)])]
      
      Worklist$InjectionOrder <- 1:nrow(Worklist)
      
      # Converting the worklist to long format
      Worklist <- Worklist %>% gather(key, File, -SampleID, -VialPos, 
                                      -InjectionOrder) %>% 
            separate(key, c("type", "Mode"), "\\.")
      
      Worklist$type <- NULL
      
      # Adding the method to use
      Methods <- data.frame(Mode = c("Epos", "Eneg", "Apos", "Aneg"),
                            Method = c("metabolomics+ESI.m",
                                       "metabolomicsnegESI.m",
                                       "metabolomicsAPCI+.m",
                                       "metabolomicsAPCIneg.m")) # !!! NEED TO CHECK APCI METHODS!!!!
      Worklist <- join(Worklist, Methods, by = "Mode")
      
      # Adding the file path
      Worklist$FilePath <- paste0(FilePath, Worklist$File, ".d")
      
      # Noting sample type
      Worklist$SampType <- "clinical"
      Worklist$SampType[str_detect(Worklist$SampleID, "QC")] <- "QC"
      Worklist$SampType[str_detect(Worklist$SampleID, "MQC")] <- "Master QC"
      Worklist <- Worklist[, c("SampleID", "VialPos", "Method", "FilePath",
                               "InjectionOrder", "Mode", "SampType", "File")]
      
      return(Worklist)
            
}


# Example
library(xlsx)
setwd("D:/Users/Laura/Documents/Work/Lin Lab/Mn exposure project")
Meta <- read.xlsx("subject information north star metabolomics.xlsx", 
                   sheetName = 4)
SampleID <- as.character(Meta$LabID)
Project <- "MnPS"
Matrix <- "urine"
FilePath <- "D:\\MassHunter\\Data\\Laura\\Mn exposure\\Mn exposure Puget Sound workers\\"
Date <- 20150209

MyWorklist <- worklist(SampleID, Date, Matrix, FilePath)

