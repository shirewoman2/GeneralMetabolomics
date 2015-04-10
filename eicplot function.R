# eicplot function ----------------------------------------
# This function plots EICs of the mass feature of interest, facetted by 
# project if you input a project if Project is one of the columns in the
# input data.frame. It creates a new ggplot2 object, "Plot", that can be 
# further modified.
# Input: 
#       1. MF = a character string listing the mass feature of interest
#       2. EICs = a data.frame with all the EIC info already extracted; this 
#          should be the output from the function "eic".

eicplot <- function(MF, EICs, Height = 8, Width = 8) {
      
      require(ggplot2)
      
      MF.df <- EICs[EICs$MassFeature == MF, ]
      
      if("Project" %in% names(MF.df)){
            
            Plot <<- ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
                  geom_line() + ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF)) +
                  xlab("RT (min)") +
                  geom_vline(data = MF.df, aes(xintercept = RT.original),
                             linetype = "dashed", size = 0.5, color = "gray50") +
                  theme(legend.position = "none") +
                  facet_wrap(~ Project, scales = "free")
            Plot
            ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], paste0(MF, ".png")),
                   height = Height, width = Width)
      } else {
            Plot <<- ggplot(MF.df, aes(x = RT, y = Intensity, color = File)) +
                  geom_line() + ggtitle(paste(MF.df$Mode[1], MF.df$Matrix[1], MF)) +
                  xlab("RT (min)") +
                  geom_vline(data = MF.df, aes(xintercept = RT.original),
                             linetype = "dashed", size = 0.5, color = "gray50") +
                  theme(legend.position = "none")
            Plot
            ggsave(paste(MF.df$Mode[1], MF.df$Matrix[1], paste0(MF, ".png")),
                   height = Height, width = Width)
      }
}


# # EXAMPLE
# # Making graphs 
# setwd("Fragmentation results and data files/EnegP83 and EposU84 top MF plots")
# 
# for (m in 1:nrow(Sig.3)){
#       
#       eicplot(Sig.3$MassFeature[m], TopEICs)
# }
# 
# setwd("C:/Users/Laura/Documents/LCMS metabolomics")
