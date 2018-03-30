# reads for each patient the peptides
raw_peptide_reader <- function(path){
  #path <- "../../rohdaten/csv/individual_patient_peptides/"
  files <- list.files(path = path, pattern = "classI.csv")
  result <- list()
  for(f in files){
    result[[strsplit(f, split = "\\_")[[1]][1]]] <- as.vector(read.csv(paste0(path,f))[,1])
  }
  return(result)
}