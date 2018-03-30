allotype_binder_table_calculator <- function(binders, allele){
  allotype_peptides <- binders[which(!is.na(binders[, paste0(allele, "_rank")])),]
  allotype_binders <- allotype_peptides[which(apply(cbind(
    allotype_peptides[, paste0(allele, "_aff") ]<=500,
    allotype_peptides[, paste0(allele, "_rank") ]<=2),
    1, any)),]
  
  # create 0-1-datatable
  allotype_binder_table <- data.frame(matrix(0, nrow = length(unique(unlist(allotype_binders[,"Seq"]))), 
                                             ncol = length(unique(unlist(allotype_binders[,"Patient"])))))
  row.names(allotype_binder_table) <- unique(unlist(allotype_binders[,"Seq"]))
  colnames(allotype_binder_table) <- unique(unlist(allotype_binders[,"Patient"]))
  
  for(patient in colnames(allotype_binder_table)){
    temp_data <- unlist(allotype_binders[which(allotype_binders[,"Patient"]==patient), "Seq"])
    for(pep in temp_data){
      allotype_binder_table[pep,patient] <- 1
    }
  }
  return(allotype_binder_table)
}