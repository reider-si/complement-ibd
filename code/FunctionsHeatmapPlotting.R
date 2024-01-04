prepare_df_hmp <- function(datafr, cutoffval) {
  group <- paste(datafr$DISEASE, datafr$LESIONAL)           # paste Disease and Lesional group together
  datafr<-datafr[3:26]                                      # remove DISEASE and LESIONAL column, not longer needed
  datafr <- t(datafr)                                       # transpose the dataframe
  colnames(datafr)<-group                                   # Rename column headers 
  
  datafr <- rownames_to_column(as.data.frame(unclass(datafr)), var="Target")        # Target as a column, not as row names
  rm(group)
  datafr.norm <- datafr[,3:6]/datafr[,2]                          # Dividing everything by the healthy control (i.e. first row)
  datafr.norm <- cbind(datafr$Target, datafr.norm)                            # constructs normalized dataframe
  names(datafr.norm)<- c("Target", "CD NL", "CD L", "UC NL", "UC L")          # Rename columns
  rm(datafr)
  
  datafr.norm.melt <- melt(datafr.norm, id.vars = "Target")                               # Melting df
  
  datafr.norm.melt$logval <- log10(datafr.norm.melt$value)                                # Transforming value column to log10
  datafr.norm.melt$logvalr <- round(datafr.norm.melt$logval, 2)                           # Round to 2 digits
  datafr.norm.melt$logvalp = datafr.norm.melt$logval
  datafr.norm.melt$logvalp[datafr.norm.melt$logvalr > cutoffval] <- cutoffval             # set logval 3 values above cutoffval to cutoffval
  datafr.norm.melt$logvalp[datafr.norm.melt$logvalr < -cutoffval] <- -cutoffval           # set logval 3 values below cutoffval to cutoffval
  datafr.norm.melt                                                                        # return the dataframe
}

plot_heatmap <- function(datafr, caption_string, dpival, file_suffix) {
  pal <- diverge_hcl(9, h = c(260, 8), c=100, l = c(45,100), power = c(0.1, 0.9), fixup = T)
  print(ggplot(datafr, aes(x=variable, y = Target))+
          geom_tile(aes(fill = logvalp))+
          scale_fill_gradientn(colours = pal, limits = c(-cutoffval, cutoffval))+
          labs(x="", fill = "", title = caption_string)+
          scale_x_discrete(position = "top")+
          scale_y_discrete(limits = rev(levels(factor(datafr$Target))))+
          theme_classic()+
          theme(text = element_text(size=14))
  )
  ggsave(file = paste0("results/qpcr/heatmap-log ", file_suffix, ".pdf"), height = 7, width = 5, dev = "pdf")
  ggsave(file = paste0("results/qpcr/heatmap-log ", file_suffix, ".png"), dpi = dpival, height = 7, width = 5)
}