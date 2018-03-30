# peptide_group_mapper
peptide_group_mapper <- function(temp_table, out_path){
  MM <-as.vector(read.csv(paste0(out_path, "MMs"), header = F)[,1])
  MCL <-as.vector(read.csv(paste0(out_path,"MCLs"), header = F)[,1])
  CLL <-as.vector(read.csv(paste0(out_path,"CLLs"), header = F)[,1])
  AML <-as.vector(read.csv(paste0(out_path,"AMLs"), header = F)[,1])
  CML <-as.vector(read.csv(paste0(out_path,"CMLs"), header = F)[,1])
  
  
  peptide_groups <- list()
  peptide_groups[["MM"]] <- temp_table[,which(colnames(temp_table)%in% MM | colnames(temp_table)%in% MCL)]
  if(typeof(peptide_groups[["MM"]])=="list"){
    peptide_groups[["MM"]]  <- row.names(peptide_groups[["MM"]][which(rowSums(peptide_groups[["MM"]])>0),]) 
  }else{
    peptide_groups[["MM"]]  <- row.names(temp_table[which(peptide_groups[["MM"]]==1),]) 
  }
  peptide_groups[["CLL"]] <- temp_table[,which(colnames(temp_table)%in% CLL)]
  if(typeof(peptide_groups[["CLL"]])=="list"){
    peptide_groups[["CLL"]]  <- row.names(peptide_groups[["CLL"]][which(rowSums(peptide_groups[["CLL"]])>0),]) 
  }else{
    peptide_groups[["CLL"]]  <- row.names(temp_table[which(peptide_groups[["CLL"]]==1),]) 
  }
  peptide_groups[["AML"]] <- temp_table[,which(colnames(temp_table)%in% AML)]
  if(typeof(peptide_groups[["AML"]])=="list"){
    peptide_groups[["AML"]]  <- row.names(peptide_groups[["AML"]][which(rowSums(peptide_groups[["AML"]])>0),]) 
  }else{
    peptide_groups[["AML"]]  <- row.names(temp_table[which(peptide_groups[["AML"]]==1),]) 
  }
  peptide_groups[["CML"]] <- temp_table[,which(colnames(temp_table)%in% CML)]
  if(typeof(peptide_groups[["CML"]])=="list"){
    peptide_groups[["CML"]]  <- row.names(peptide_groups[["CML"]][which(rowSums(peptide_groups[["CML"]])>0),]) 
  }else{
    peptide_groups[["CML"]]  <- row.names(temp_table[which(peptide_groups[["CML"]]==1),]) 
  }

  
  return(peptide_groups)
}