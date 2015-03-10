# This function generates a worklist for metabolomics samples where the sample
# order is randomized, there is a QC injection every 10 injections, and the 
# worklist starts and ends with a Master QC injection. 
# Input: 
#       1. Samples -- a data.frame with the following columns:
#             a. FileLabel - a vector of the sample IDs as you would like them 
#             to show up in the file names. These must be unique and should NOT
#             contain any special characters, including periods because it
#             makes things more challenging downstream. They CAN contain spaces.
#             b. SampleID -- a unique identifier that ties everything back to the 
#             samples, something that will allow you to connect the worklist 
#             with the samples in your data.frame. These should follow 
#             standard R requirements for column names and CAN contain 
#             periods but NO spaces.
#       2. Date -- a number for today's date in the format YYYYMMDD, which 
#       will be at the front of each file. No symbols allowed.
#       3. Project -- a brief abbreviation or code word for which project this
#       is for. Example: "MnPS" or "SCOR MDZ".
#       4. Matrix -- the sample matrix (this is set up to make a worklist for
#                    one matrix at a time. For multiple matrices,
#                    make multiple worklists.)
#       5. Mode -- a character string of which ionization modes you'd like to run. 
#       Defaults to c("Epos", "Eneg") but other options are "Apos" and "Aneg".
#       6. FilePath -- a character string with the path where you'd like to
#       save the files. Format is tricky; it must be like this:
#             "D:\\MassHunter\\Data\\Laura\\Mn exposure\\Mn exposure Puget Sound workers\\"
#       7. Column -- a character string listing which column was used. Options
#       are "SB-Aq" and "TOSOH". Default is "SB-Aq".
#       8. Qnum.start -- number of QC injections you'd like at the start of the
#       run. (Defaults to 3)
#       9. Qnum.end -- number of QC injections you'd like at the end of the run. 
#       (Defaults to 1.)
#       10. MQnum.start -- number of Master QC injections you'd like at the start of the
#       run. (Defaults to 1)
#       11. MQnum.end -- number of Master QC injections you'd like at the end of the run. 
#       (Defaults to 1.)
#       12. Initials -- defaults to "LS" and doesn't do anything unless 
#       "Labels" is set to TRUE.
#       13. Labels -- if set to TRUE, will create a data.frame "Labels" that 
#       is set up to have 21 rows, just like the ToughTag layout, for 
#       printing vial labels for samples. The format will be: 
#       "Project SampleID dd/mm/yyyy Initials". Copy and paste this into Word 
#       with the ToughTag template and print.
##    !!!   Set the seed if you want to create the same worklist every time you
###   !!!   run this script. i.e. set.seed(1234)

worklist <- function(Samples, Date, Project, Matrix, 
                     Mode = c("Epos", "Eneg"),
                     FilePath, Column = "SB-Aq",
                     Qnum.start = 3, Qnum.end = 1, 
                     MQnum.start = 1, MQnum.end = 1,
                     Initials = "LS", Labels = FALSE) {
      
      require(plyr)
      require(tidyr)
      require(dplyr)
      require(stringr)
      require(lubridate)
      
      Samples$RandNum <- rnorm(nrow(Samples))
      Samples <- Samples[, c("SampleID", "FileLabel", "RandNum")]
      
      # Setting up QC injections
      Qnum.start.2 <- Qnum.start-1
      Qnum <- Qnum.start + Qnum.end + nrow(Samples) %/%10
      QCsamp <- data.frame(FileLabel = paste0("QC", 1:Qnum))
      QCsamp$VialPos <- "P1-A2"
      
      # Setting up Master QC injections
      MQCsamp <- data.frame(FileLabel = paste0("MQC", 1:(MQnum.start +
                                                               MQnum.end)))
      MQCsamp$VialPos <- "P1-A1"
      
      if ("Epos" %in% Mode) {
            Samples$File.Epos <- paste(
                  Date, Project, paste0("Epos", toupper(str_sub(Matrix, 1, 1))),
                  Samples$FileLabel)
            QCsamp$File.Epos <- paste(
                  Date, Project, paste0("Epos", toupper(str_sub(Matrix, 1, 1))), 
                  QCsamp$FileLabel)
            MQCsamp$File.Epos <- paste(
                  Date, Project, paste0("Epos", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$FileLabel)
            
      }
      
      if ("Eneg" %in% Mode){
            Samples$File.Eneg <- paste(
                  Date, Project, paste0("Eneg", toupper(str_sub(Matrix, 1, 1))),
                  Samples$FileLabel)
            QCsamp$File.Eneg <- paste(
                  Date, Project, paste0("Eneg", toupper(str_sub(Matrix, 1, 1))),
                  QCsamp$FileLabel)
            MQCsamp$File.Eneg <- paste(
                  Date, Project, paste0("Eneg", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$FileLabel)
      }
      
      if ("Apos" %in% Mode){
            Samples$File.Apos <- paste(
                  Date, Project, paste0("Apos", toupper(str_sub(Matrix, 1, 1))),
                  Samples$FileLabel)
            QCsamp$File.Apos <- paste(
                  Date, Project, paste0("Apos", toupper(str_sub(Matrix, 1, 1))),
                  QCsamp$FileLabel)
            MQCsamp$File.Apos <- paste(
                  Date, Project, paste0("Apos", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$FileLabel)
      }
      
      if ("Aneg" %in% Mode){
            Samples$File.Aneg <- paste(
                  Date, Project, paste0("Aneg", toupper(str_sub(Matrix, 1, 1))),
                  Samples$FileLabel)
            QCsamp$File.Aneg <- paste(
                  Date, Project, paste0("Aneg", toupper(str_sub(Matrix, 1, 1))),
                  QCsamp$FileLabel)
            MQCsamp$File.Aneg <- paste(
                  Date, Project, paste0("Aneg", toupper(str_sub(Matrix, 1, 1))),
                  MQCsamp$FileLabel)
      }
      
      # Setting the vial positions
      VialPos <- paste0(rep(c("P1-", "P2-"), each = 54), # 6 rows, 9 columns 
                        rep(paste0(rep(LETTERS[1:6], each = 9), c(1:9)), 2))
      
      Samples$VialPos <- VialPos[3:(nrow(Samples) + 2)]
      Samples <- arrange(Samples, RandNum)
      
      # Putting together the worklist
      Worklist <- list()
      
      for (i in 1:(nrow(Samples) %/% 10)){
            Worklist[[i]] <- rbind.fill(QCsamp[i+Qnum.start.2, ],
                                        Samples[
                                              ((1+10*(i-1)):(10*i)), ])
      }
      
      if (nrow(Samples) %% 10 > 0) {
            Worklist[[nrow(Samples) %/% 10 + 1]] <- rbind.fill(
                  QCsamp[nrow(Samples) %/% 10 + Qnum.start, ],
                  Samples[((nrow(Samples) %/% 10)*10+1):nrow(Samples), ])
      }
      
      Worklist[[1]] <- rbind.fill(MQCsamp[1:MQnum.start, ], 
                                  QCsamp[1:Qnum.start.2, ],
                                  Worklist[[1]])
      
      Worklist[[nrow(Samples) %/% 10 + 1]] <- 
            rbind.fill(Worklist[[nrow(Samples) %/% 10 + 1]], 
                       QCsamp[(Qnum-Qnum.end+1):Qnum, ],
                       MQCsamp[(MQnum.start+1):(MQnum.start+MQnum.end), ])
      
      Worklist <- rbind.fill(Worklist)
      Worklist <- Worklist[, c("SampleID", "FileLabel", "VialPos", 
                               names(Worklist)[names(Worklist) %in% 
                                                     paste0("File.", Mode)])]
      
      Worklist$InjectionOrder <- 1:nrow(Worklist)
      
      # Converting the worklist to long format
      Worklist <- Worklist %>% gather(key, File, -SampleID, -FileLabel,
                                      -VialPos, 
                                      -InjectionOrder) %>% 
            separate(key, c("type", "Mode"), "\\.")
      
      Worklist$type <- NULL
      
      Worklist$Column <- Column
      
      # Adding the method to use
      Methods <- data.frame(Mode = rep(c("Epos", "Eneg", "Apos", "Aneg"), 2),
                            Column = rep(c("SB-Aq", "TOSOH"), each = 4),
                            Method = c("metabolomics+ESI.m",
                                       "metabolomicsnegESI.m",
                                       "metabolomicsAPCI.m",
                                       "metabolomicsAPCIneg.m", # Not sure of APCI- method.
                                       "polar_TOSOH_metabolomics+ESI.m", 
                                       "polar_TOSOH_metabolomicsnegESI.m", 
                                       "polar_TOSOH_metabolomics+APCI.m",
                                       NA)) 
      Worklist <- join(Worklist, Methods, by = c("Mode", "Column"))
      
      # Adding the file path
      Worklist$FilePath <- paste0(FilePath, Worklist$File, ".d")
      
      # Noting sample type
      Worklist$SampType <- "clinical"
      Worklist$SampType[str_detect(Worklist$FileLabel, "QC")] <- "QC"
      Worklist$SampType[str_detect(Worklist$FileLabel, "MQC")] <- "Master QC"
      
      Worklist$SampleID <- as.character(Worklist$SampleID)
      
      Worklist$SampleID[Worklist$SampType == "QC"] <- 
            Worklist$FileLabel[Worklist$SampType == "QC"]
      
      Worklist$SampleID[Worklist$SampType == "Master QC"] <- 
            Worklist$FileLabel[Worklist$SampType == "Master QC"]
      
      Worklist <- Worklist[, c("SampleID", "VialPos", "Method", 
                               "FilePath", "FileLabel", 
                               "InjectionOrder", "Mode", "SampType", "File",
                               "Column")]
      
      if (Labels == TRUE) {
            Samples <- arrange(Samples, VialPos)
            
            Lab.string <- paste(Project, Samples$FileLabel, 
                                format(ymd(Date), format="%m/%d/%Y"),
                                Initials)
            Lab.string <- c(Lab.string, 
                            paste(Project, "QC", 
                                  format(ymd(Date), format="%m/%d/%Y"),
                                  Initials))
            NCol <- length(Lab.string) %/% 21
            LastCol <- length(Lab.string ) %% 21
            
            Col <- list()
            for (n in 1:NCol) {
                  Col[[n]] <- Lab.string[((n-1)*21+1):(21*n)]
            }
            
            Col[[NCol+1]] <- Lab.string[(NCol*21)+(1:LastCol)]
            Col[[NCol+1]] <- c(Col[[NCol+1]], rep("", 21-LastCol))
            
            names(Col) <- LETTERS[1:(NCol+1)]
            
            template <- data.frame(A = Col[[1]])
            for (n in 2:length(Col)) {
                  template <- cbind(template, Col[[n]])
            }
            
            names(template) <- names(Col)
            
            template <<- template
      }
      
      return(Worklist)
      
}

# Example
library(xlsx)
setwd("D:/Users/Laura/Documents/Work/Lin Lab/Mn exposure project")
Samples <- read.xlsx("subject information north star metabolomics.xlsx", 
                   sheetName = 4)
Samples$FileLabel <- gsub("\\.", " ", as.character(Meta$LabID))
Samples <- plyr::rename(Samples, c("LabID" = "SampleID"))
Samples <- arrange(Samples, SampleID)
Date <- 20150209
Project <- "MnPS"
Matrix <- "urine"
Mode = c("Epos", "Eneg")
FilePath <- "D:\\MassHunter\\Data\\Laura\\Mn exposure\\Mn exposure Puget Sound workers\\"
Column = "SB-Aq"
Qnum.start = 3
Qnum.end = 1
MQnum.start = 1
MQnum.end = 1
Initials = "LS"
Labels = TRUE

MyWorklist <- worklist(Samples = Samples, Date = Date, Project = Project, 
                       Matrix = Matrix, FilePath = FilePath, Labels = TRUE)

# Checking that all the samples are in the worklist.
Samples$SampleID %in% MyWorklist$SampleID

# Checking that I don't have any replicate file names.
anyDuplicated(Worklist$File)

