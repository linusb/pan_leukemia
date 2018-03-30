# calculate_overlap_matrix

calculate_overlap_matrix <- function(temp_table, cut, jaccard){
  overlap_matrix_protein <- data.frame(matrix(ncol=ncol(temp_table), nrow=ncol(temp_table)))
  colnames(overlap_matrix_protein) <- colnames(temp_table)
  row.names(overlap_matrix_protein) <- colnames(temp_table)
  # calculating overlap
  for(i in colnames(temp_table)){
    for(j in colnames(temp_table)){
      if(j == i){
        overlap_matrix_protein[i,j] <- 0
      }else{
        if(jaccard){
          overlap <- 1-dist(rbind(temp_table[,i], temp_table[,j]),method = "binary")
        }else{
          i_peptides <- row.names(temp_table)[which(temp_table[,i]==1)]
          j_peptides <- row.names(temp_table)[which(temp_table[,j]==1)]
          overlap <- length(intersect(i_peptides, j_peptides))/ 
            ((length(i_peptides)+ length(j_peptides))/2)
        }
        if(overlap<cut){
          overlap_matrix_protein[i,j] <- 0
        }else{
          overlap_matrix_protein[i,j] <- overlap
        }
      }
    }
  }
  return(overlap_matrix_protein)
}