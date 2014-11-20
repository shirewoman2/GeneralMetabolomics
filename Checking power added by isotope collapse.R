# Checking power added by isotope collapse
# 10/9/14


# Housekeeping =================================
library(ggplot2)
library(plyr)
library(reshape2)
library(xlsx)


setwd("D:/Users/Laura/Documents/Work/Lin Lab/LCMS metabolomics")

Numbers <- read.xlsx2("Checking value of collapsing isotopes in terms of statistical power.xlsx",
                     sheetName = "NumMFs", colClasses = c("character", rep("numeric", 3)))
names(Numbers) <- c("Dataset", "NumWithout", "NumWith", "PercDecrease")
str(Numbers)

PowerGain <- function(NumWithout, NumWith, p) {
      Without <- p.adjust(p = p, method = "BH", n = NumWithout)
      With <- p.adjust(p = p, method = "BH", n = NumWith)
      PercFewer <- (NumWithout - NumWith)/NumWithout
      
      print(paste0(round(PercFewer*100, digits = 0),"% fewer compounds with collapsing"))
      print(paste("Adjusted p value without collapsing =", 
                  signif(Without, digits = 2)))
      print(paste("Adjusted p value with collapsing =",
                  signif(With, digits = 2)))
      print(paste0(round((Without - With)/Without*100, digits = 0),
                   "% lower adjusted p value"))
}

i=5
PowerGain(Numbers$NumWithout[i], Numbers$NumWith[i], p = 5e-06)

