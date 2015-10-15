# specifions function ----------------------------------------
# Rather than have CAMERA search for likely related ions, look for specific 
# ones. Input is:
# 1. MF.df, a 1-row data.frame with:
#       a. MassFeature (character)
#       b. mz (numeric)
#       c. RT (numeric)
#       d. Mode (character or factor)
#       e. Matrix (character or factor)
# 2. Files = a data.frame with the following columns:
#       a. the files (File)
#       b. directory (Directory)
#       c. ionization mode (Mode)
#       d. and matrix (Matrix)
# 3. Ions, a list of which ions you want or the specific m/z to use. Specific 
#    ions to call by name are: "M+Na", "M+1", "M+2", "M+3", "M+Cl", "M-H2O", 
#    "M-1", "M-2", "M-3". Default is all of those.
# Output is data.frame with all the EICs.



specifions <- function(MF.df, Files, 
                       Ions = c("M+Na", "M-Na", "M+1", "M+2", "M-1", "M-2", 
                                "M-H2O", "M+H2O")) {
      
      require(plyr)
      
      OrigDir <- getwd()
      setwd("I:/General LCMS scripts")
      source("eic function.R")
      setwd(OrigDir)
      
      PossibleIons <- data.frame(Ion = c("M+Na",
                                         "M-Na",
                                         "M+K",
                                         "M-K",
                                         "M+Cl", 
                                         "M-Cl",
                                         "M+H2O",
                                         "M-H2O",
                                         "M+1", 
                                         "M+2", 
                                         "M+3",
                                         "M-1",
                                         "M-2", 
                                         "M-3"),
                                 MassAdd = c(22.98977- 1.0073, 
                                             -22.98977+ 1.0073,
                                             39.0983 - 1.0073,
                                             -39.0983 + 1.0073,
                                             34.96885+1.0073, 
                                             -34.96885-1.0073, 
                                             15.999405 + 2*1.0073,
                                             -(15.999405 + 2*1.0073),
                                             1.008665,
                                             2*1.008665,
                                             3*1.008665,
                                             -1.008665,
                                             -2*1.008665,
                                             -3*1.008665),
                                 Mode = c(rep("pos", 4), "neg", "neg", 
                                          rep(NA, 8)))
      
      MF.df$MassFeature.ion <- "originally detected ion"
      
      for (i in 1:length(Ions)) {
            MF.df <- rbind(MF.df,
                           MF.df[1, ])
            
            if (Ions[[i]] %in% PossibleIons$Ion){
                  MF.df$mz[i+1] <- MF.df$mz[1] + 
                        PossibleIons$MassAdd[which(PossibleIons$Ion == 
                                                         Ions[[i]])]
            } else {
                  MF.df$mz[i+1] <- MF.df$mz[1] + Ions[[i]]
            }
            
            MF.df$MassFeature.ion[i+1] <- as.character(Ions[[i]])
            
      }
      
      EICs <- eic(MF.df, Files)
      if ("MassFeature.ion" %in% names(EICs) & 
                "MassFeature.ion" %in% names(MF.df)){
            EICs <- join(EICs, MF.df[, c("MassFeature", "mz", "MassFeature.ion")], 
                         by = c("MassFeature", "MassFeature.ion"))
      } else {
            EICs <- join(EICs, MF.df[, c("MassFeature", "mz", "MassFeature.ion")], 
                         by = c("MassFeature"))
      }
      
      
      return(EICs)
      
}
