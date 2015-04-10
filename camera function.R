# camera function --------------------------
# This function processes an xcmsSet object using CAMERA to look for 
# other ions that could arise from the same mass feature. Input is an
# xcmsSet object (xset) and the ionization mode as "positive" or "negative".
# Input is:
#       1. a grouped xcmsSet object (xset)
#       2. the ionization mode as "Epos", "Eneg", "Apos", or "Aneg" (Mode)
#       3. the resolution in ppm to be used for matching adducts
#       and isotopes (PPM); defaults to 15
#       4. Since CAMERA looks at the correlation of chromatograms between ions
#       as part of how it determines which ions are, in fact, adducts or isotopes,
#       the p value cutoff required. True adducts and ions should correlate 
#       quite well in the raw data, so the default is 0.0001. (PVal)
# Output is a named list containing:
#       1. the CAMERA object (root part of the name is whatever the input 
#       xcmsSet object was called, suffix is .xsa)
#       2. a data.frame of the output from CAMERA (.annot). There will be one row for
#       every ion in the original data after the reference ions were removed. 
#       Columns:
#             a. the mass feature name (MassFeature), 
#             b. the m/z (mz), 
#             c. the retention time in minutes (RT), 
#             d. any isotopes CAMERA thinks it has found (isotopes), 
#             e. any adducts it thinks it has found (adduct), and 
#             f. the pseudospectrum group (pcgroup). "Pseudospectra" are the 
#             chromatograms that CAMERA hypothesizes could be related, i.e. they
#             follow the same pattern in intensity. Pseudospectra do not necessarily
#             have m/z that would actually make them the same compound, though, so
#             think of grouping the ions into pseudospectra as simply the first 
#             step in making a hypothesis about which ions might be related. 
#       3. a data.frame (suffix is .otherions) where each row is an ion that had 
#       at least one other potentially related ion, so >= 1 other isotope or 
#       adduct. The columns are:
#             a. MassFeature, 
#             b. mz, 
#             c. RT, 
#             d. pseudospectrum group (pcgroup), 
#             e. the type of ion CAMERA thinks this ion might be (IonType), 
#             f. the possible charge of the ion (Charge), 
#             g. the m/z of the other, potentially related ion (mzOfM), 
#             h. and a group number where other ions with the same number
#             are potentially isotopes of this one (IsoGroup). 


# Nota bene: The xcmsSet object must have already had the group() function
# applied to it or this will not work. That's something inherent in the 
# CAMERA package, not something I implemented. The problem with that is that
# you then have to run this on your ENTIRE dataset, which will be time 
# consuming and RAM intensive if you've got very many samples. You can
# always subset your xcmsSet object to contain only a handful of samples, but
# then you have to redo the xcms::group() function and then, once you've
# run that function, you'll have to figure out which ions correspond to which
# in the full dataset.   - Laura

camera <- function(xset, Mode, PPM = 15, PVal = 0.0001) {
      require(xcms)
      require(CAMERA)
      require(plyr)
      require(stringr)
      
      if(Mode %in% c("Epos", "Apos")){
            Mode2 <- "positive"
      } else {
            Mode2 <- "negative"
      }
      
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
                         multiplier = 3, polarity = Mode2, 
                         rules = NULL, max_peaks = 100)
      
      # Make a table of the annotated peaks.
      xset.annot <- getPeaklist(xsa)
      
      # Once again, remove reference ions from the peak table.
      ifelse(Mode %in% c("Epos", "Apos"),
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
                            pcgroup = as.numeric(as.character(Iso$pcgroup)),
                            IsoGroup = as.numeric(as.character(IsoGroup)), 
                            IonType = as.character(IonType), 
                            Charge = as.character(Charge),
                            stringsAsFactors = FALSE)
      
      IsoList$mzOfM <- NA
      for (i in 1:nrow(IsoList)){
            
            if (IsoList$IonType[i] == "M") {
                  IsoList$mzOfM[i] <- IsoList$mz[i]
            } else {
                  n <- as.numeric(str_sub(IsoList$IonType[i], 3, 3))
                  
                  IsoList$mzOfM[i] <- IsoList$mz[i] - n * 1.00866
                  
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
      # 2. NeutralMassOfM = the neutral mass of that particular adduct
      # 3. MassFeature 
      # 4. pcgroup
      
      AdSplit <- str_split(Adduct$adduct, "\\[")
      AdList <- list()
      
      for (i in 1:nrow(Adduct)){
            
            IonType <- c()
            Charge <- c()
            NeutralMassOfM <- c()
            
            for (m in 2:length(AdSplit[[i]])){
                  
                  IonType[m-1] <- unlist(str_extract(AdSplit[[i]][m], ".*\\]"))
                  IonType[m-1] <- gsub("\\]", "", IonType[m-1])
                  
                  Charge[m-1] <- unlist(str_extract(AdSplit[[i]][m], "\\].{1,2}"))
                  Charge[m-1] <- str_trim(gsub("\\]", "", Charge[m-1]))
                  
                  NeutralMassOfM[m-1] <- str_trim(str_extract(AdSplit[[i]][m], " .*"))
            }
            
            AdList[[i]] <- data.frame(MassFeature = Adduct$MassFeature[i],
                                      mz = Adduct$mz[i],
                                      RT = Adduct$RT[i],
                                      pcgroup = as.numeric(as.character(
                                            Adduct$pcgroup[i])),
                                      IonType = IonType,
                                      Charge = Charge,
                                      NeutralMassOfM = 
                                            as.numeric(NeutralMassOfM),
                                      stringsAsFactors = FALSE)
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

