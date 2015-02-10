# allionplot function ----------------------------------------
# This function plots EICs of the mass feature of interest, facetted by 
# other ions. It creates a new ggplot2 object, "Plot.allion", that can be
# further modified.
# Input: the data.frame output from the function "allion".

allionplot <- function(allion.df, Height = 8, Width = 8, 
                       FileSuffix = "- all ions.png") {
      
      require(ggplot2)
      
      Plot.allion <<- ggplot(allion.df, aes(x = RT, y = Intensity, 
                                            color = File)) +
            geom_line() + ggtitle(paste(allion.df$Mode[1], 
                                        allion.df$Matrix[1], 
                                        allion.df$MassFeature[1])) +
            xlab("RT (min)") +
            geom_vline(data = allion.df, aes(xintercept = RT.original),
                       linetype = "dashed", size = 0.5, color = "gray50") +
            theme(legend.position = "none") +
            facet_grid(MassFeature.ion ~ Project, scales = "free")
      Plot.allion
      ggsave(paste(allion.df$Mode[1], allion.df$Matrix[1], 
                   allion.df$MassFeature[1], FileSuffix),
             height = Height, width = Width)
}
