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
#       b. directory (Dir)
#       c. ionization mode (Mode)
#       d. and matrix (Matrix)
# 3. Ions, a list of which ions you want or the specific m/z to use. Specific 
#    ions to call by name are: "M+Na", "M+1", "M+2", "M+3", "M+Cl", "M-H2O".
# !!!!! You must have sourced the function "eic" for this to run. !!!!!

specifions <- function(MF.df, Files, Ions) {
      
      require(plyr)
      
      PossibleIons <- data.frame(Ion = c("M+Na", 
                                         "M+1", 
                                         "M+2", 
                                         "M+3",
                                         "M+Cl", 
                                         "M-H2O",
                                         "M-1",
                                         "M-2", 
                                         "M-3"),
                                 MassAdd = c(22.98977- 1.0073, # positive mode: If molecule gains Na, wouldn't also gain H
                                             1.008665,
                                             2*1.008665,
                                             3*1.008665,
                                             34.96885+1.0073, # negative mode: If molecule gains Cl, wouldn't also lose H
                                             -18.01056,
                                             -1.008665,
                                             -2*1.008665,
                                             -3*1.008665))
      
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
      EICs <- join(EICs, MF.df[, c("MassFeature", "mz", "MassFeature.ion")], 
                   by = c("MassFeature", "MassFeature.ion"))
      return(EICs)
      
}
