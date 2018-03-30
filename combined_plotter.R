# combined_plotter


group_color_mapper <- function(x, out_path){
  # group colors
  MM <-as.vector(read.csv(paste0(out_path, "MMs"), header = F)[,1])
  MCL <-as.vector(read.csv(paste0(out_path,"MCLs"), header = F)[,1])
  CLL <-as.vector(read.csv(paste0(out_path,"CLLs"), header = F)[,1])
  AML <-as.vector(read.csv(paste0(out_path,"AMLs"), header = F)[,1])
  CML <-as.vector(read.csv(paste0(out_path,"CMLs"), header = F)[,1])
  if(x %in% MM){
    return("orange")
  }else if(x %in% MCL){
    return("yellow")
  }else if(x %in% CLL){
    return("red")
  }else if(x %in% AML){
    return("green")
  }else if(x %in% CML){
    return("darkgreen")
  }else{
    print(paste0("no_mapping for: ", x))
  }
}

combined_plotter <- function(temp_table, allele, suffix, out_path){
  print(suffix)
  
  temp_table <- temp_table[which(rowSums(temp_table)>0),]
  temp_table <- temp_table[,which(colSums(temp_table)>0)]
  group_colors <- sapply(colnames(temp_table), function(x) group_color_mapper(x, out_path))
  
  heatmap_plotter(temp_table, allele = allele, suffix = suffix, group_colors = group_colors, out_path)
  # venny
  peptide_groups <- peptide_group_mapper(temp_table, out_path)
  group_frame <- data.frame(matrix("",ncol=4, nrow= max(sapply(peptide_groups, length))))
  colnames(group_frame) <- names(peptide_groups)
  for(item in names(peptide_groups)){
    group_frame[,item] <- c(peptide_groups[[item]], rep("", (nrow(group_frame)-length(peptide_groups[[item]]))))
  }
  write.csv(group_frame,paste0(out_path,"vennies/", allele,"/", allele,suffix, "_lineage_proteins.csv"), row.names=T)
  
  overlap <- venn(peptide_groups)
  
  out_frame <- data.frame(matrix("",ncol=15, nrow= max(sapply(attributes(overlap)$intersections, length))))
  colnames(out_frame) <- names(attributes(overlap)$intersections)
  for(item in names(attributes(overlap)$intersections)){
    out_frame[,item] <- c(attributes(overlap)$intersections[[item]], rep("", (nrow(out_frame)-length(attributes(overlap)$intersections[[item]]))))
  }
  write.csv(out_frame,paste0(out_path,"vennies/", allele,"/", allele,suffix, ".csv"), row.names=T)
  venn.diagram(x= peptide_groups, filename = paste0(out_path,"vennies/",allele,"/",allele,suffix,".png"), main = paste0(allele," Peptides") , imagetype = "png")
  # venn statistics
  inner_overlap <- attributes(overlap)$intersections$"MM:CLL:AML:CML"
  inner_overlap_matrix <- data.frame(matrix(0, nrow = length(inner_overlap), ncol = 4))
  row.names(inner_overlap_matrix) <- inner_overlap
  colnames(inner_overlap_matrix) <- c("MM/MCL", "CLL", "AML", "CML")
  MM <-as.vector(read.csv(paste0(out_path, "MMs"), header = F)[,1])
  MCL <-as.vector(read.csv(paste0(out_path,"MCLs"), header = F)[,1])
  CLL <-as.vector(read.csv(paste0(out_path,"CLLs"), header = F)[,1])
  AML <-as.vector(read.csv(paste0(out_path,"AMLs"), header = F)[,1])
  CML <-as.vector(read.csv(paste0(out_path,"CMLs"), header = F)[,1])
  
  for(pep in row.names(inner_overlap_matrix)){
    samples <- colnames(temp_table)[which(temp_table[pep,]==1)]
    for(s in samples){
      if(s %in% MM | s %in% MCL){
        inner_overlap_matrix[pep, "MM/MCL"] <- inner_overlap_matrix[pep, "MM/MCL"] +1
      }else if(s %in% CLL){
        inner_overlap_matrix[pep, "CLL"] <- inner_overlap_matrix[pep, "CLL"] +1
      }else if(s %in% AML){
        inner_overlap_matrix[pep, "AML"] <- inner_overlap_matrix[pep, "AML"] +1
      }else if(s %in% CML){
        inner_overlap_matrix[pep, "CML"] <- inner_overlap_matrix[pep, "CML"] +1
      }
    }
  }
  write.csv(inner_overlap_matrix,paste0(out_path,"vennies/", allele,"/", allele,suffix, "_inner_overlap_statistics.csv"), row.names=T)
  
  
  # graph
  cutoffs <- c(0, 0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  distance_method <- c(T,F)
  for(cut in cutoffs){
    print(cut)
    dir.create(file.path(paste0(out_path,"overlap_graph/",cut,"/", allele)), showWarnings = FALSE)
    for(jaccard in distance_method){
      overlap_matrix_protein <- calculate_overlap_matrix(temp_table, cut, jaccard)
      if(jaccard){
        write.csv(overlap_matrix_protein,paste0(out_path,"overlap_graph/",cut,"/", allele,"/", allele,suffix, "_jaccard.csv"), row.names=T)
        group_colors <- sapply(colnames(overlap_matrix_protein), function(x) group_color_mapper(x, out_path))
        graph_plotter(overlap_matrix = overlap_matrix_protein, allele = allele, group_colors = group_colors, suffix = paste0(suffix, "_jaccard"), cut = cut, out_path)
        
      }else{
        write.csv(overlap_matrix_protein,paste0(out_path,"overlap_graph/",cut,"/", allele,"/", allele,suffix, ".csv"), row.names=T)
        group_colors <- sapply(colnames(overlap_matrix_protein), function(x) group_color_mapper(x, out_path))
        graph_plotter(overlap_matrix = overlap_matrix_protein, allele = allele, group_colors = group_colors, suffix = suffix, cut = cut, out_path)
        
        
      }
      
    }}
}