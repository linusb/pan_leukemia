protein_table_calculator <- function(tumor_specific_allotype_binder_table, peptide_protein_list){
  # get all possible proteins
  proteins <- c()
  for(pep in row.names(tumor_specific_allotype_binder_table)){
    proteins <- unique(c(proteins,peptide_protein_list[[pep]]))
  }
  # create protein table
  protein_table <- data.frame(matrix(0, nrow = length(proteins), ncol = ncol(tumor_specific_allotype_binder_table)))
  row.names(protein_table) <- proteins
  colnames(protein_table) <- colnames(tumor_specific_allotype_binder_table)
  
  for(patient in colnames(tumor_specific_allotype_binder_table)){
    temp_patient <- tumor_specific_allotype_binder_table[,patient]
    # remove all 0 rows
    peptides <- row.names(tumor_specific_allotype_binder_table)[which(temp_patient==1)]
    # extract patient proteins
    patient_proteins <-c()
    for(pep in peptides){
      patient_proteins <- unique(c(patient_proteins,peptide_protein_list[[pep]]))
    }
    # set protein to 1
    for(patient_prot in patient_proteins){
      protein_table[patient_prot,patient] <- 1
    }
  }
  return(protein_table)
}