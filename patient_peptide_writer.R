#patient_peptide_writer

patient_peptide_writer <- function(temp_table, allele, suffix, out_path){
  # create folders
  dir.create(file.path(paste0(out_path,"listen/",allele)), showWarnings = FALSE)
  dir.create(file.path(paste0(out_path,"listen/",allele,"/",suffix)), showWarnings = FALSE)
  
  # write tumor specific peptides for each allele an patient
  for(p in colnames(temp_table)){
    p_specific_tumor_binders <- temp_table[,p]
    write.table(x = row.names(temp_table)[which(p_specific_tumor_binders==1)], file = paste0(out_path,"listen/",allele,"/", suffix,"/", p,".csv"), col.names = F, row.names = F)
  }
  
}