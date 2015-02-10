# camera function --------------------------
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
