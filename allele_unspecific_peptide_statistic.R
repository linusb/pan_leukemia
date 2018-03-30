# count the number of peptides of each sample

files <- list.files(path = "../../rohdaten/csv/individual_patient_peptides/", pattern = "classI.csv")

result <- data.frame(matrix(ncol= length(files), nrow = 1))
colnames(result) <- sapply(files, function(x){
  return(strsplit(x, split = "_")[[1]][1])
})
row.names(result) <- "Number of peptides"
for(f in files){
  result[1, strsplit(f, split = "_")[[1]][1]] <- nrow(read.csv(paste0("../../rohdaten/csv/individual_patient_peptides/",f), header = F) )
}

write.table(result, file = "allele_unspecific_peptide_statistic.csv", sep = ";", dec = ",")
