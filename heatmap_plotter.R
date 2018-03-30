heatmap_plotter <- function(temp_table, allele, suffix, group_colors, out_path){
  pdf(paste0(out_path,"heatmaps/",allele,"/", allele, suffix,".pdf"), width = 10, height = 10)
  heatmap.2(x = as.matrix(temp_table),
            trace = "none",
            dendrogram = "column", 
            labRow = F, 
            keysize = 0.5, 
            key = F,
            distfun = distfun,
            density.info = 'none',
            hclustfun = hclust.ave,
            col = mypalette,
            ColSideColors = group_colors)
  dev.off()
}