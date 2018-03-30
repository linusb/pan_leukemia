#main_complete heatmap

library("seqinr")
library(dplyr)
#library(data.table)
library(gplots)
library("VennDiagram")
library(igraph)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

source("prediction_reader.R")
source("allotype_binder_table_calculator.R")
source("protein_mapper.R")
source("protein_table_calculator.R")
source("heatmap_plotter.R")
source("peptide_group_mapper.R")
source("graph_plotter.R")
source("calculate_overlap_matrix.R")
source("combined_plotter.R")
source("patient_peptide_writer.R")
source("raw_peptide_reader.R")



# reading predicted peptides
path <- "../../predictions/"
out_path <- "global_heatmap/"
dir.create(file.path(paste0(out_path)), showWarnings = FALSE)


# extract patient data and read tumor and benign data
#all_patients <- unique(binders[,"Patient"])
#write.csv(all_patients, file=paste0(out_path,"patients.csv"), row.names = F)
#tumor <- as.vector(read.csv(paste0(out_path,"tumor.csv"), header = F)[,1])

patient_peptides <- raw_peptide_reader("../../rohdaten/csv/individual_patient_peptides/")
patient_peptides <- lapply(patient_peptides, function(x){
  return(sapply(x, toupper))
})
write.csv(names(patient_peptides), file=paste0(out_path,"all_samples.csv"), row.names = F)

benign <- as.vector(read.csv(paste0(out_path,"benign.csv"), header = F)[,1])
benign_peptides <- unique(do.call(c,patient_peptides[which(names(patient_peptides) %in% benign)]))


tumor <- as.vector(read.csv(paste0(out_path,"tumor.csv"), header = F)[,1])
tumor_samples <- patient_peptides[which(names(patient_peptides) %in% tumor)]
tumor_peptides <- unique(do.call(c,patient_peptides[which(names(patient_peptides) %in% tumor)]))


# create 0-1 matrix
tumor_matrix <- data.frame(matrix(0, ncol=length(tumor_samples), nrow = length(tumor_peptides)))
row.names(tumor_matrix) <- tumor_peptides
colnames(tumor_matrix) <- names(tumor_samples)

for(sample in names(tumor_samples)){
  print(sample)
  for(pep in tumor_samples[[sample]]){
    tumor_matrix[pep,sample] <- 1
  }
}
temp_matrix <- cbind(row.names(tumor_matrix),tumor_matrix)
colnames(temp_matrix)[1] <- "Seq"
temp_matrix[,1] <- sapply(temp_matrix[,1], as.character)
#peptide_protein_list <- protein_mapper(temp_matrix, "swissprotHUMANwoi_130927.fasta")
#save(x= peptide_protein_list, file = paste0(out_path,"peptide_protein_list.RData"))
load(file = paste0(out_path, "peptide_protein_list.RData"))

protein_matrix <- protein_table_calculator(tumor_specific_allotype_binder_table = tumor_matrix, peptide_protein_list = peptide_protein_list)

protein_matrix <- protein_matrix[which(rowSums(protein_matrix)>1),]

# Heatmap properties
hclust.ave <- function(x) hclust(x, method = "complete")
mypalette <- colorRampPalette(c('white','black'))(n=2)
distfun <- function(x) dist(x, method = "binary")
# group colors
MM <-as.vector(read.csv("MMs", header = F)[,1])
MCL <-as.vector(read.csv("MCLs", header = F)[,1])
CLL <-as.vector(read.csv("CLLs", header = F)[,1])
AML <-as.vector(read.csv("AMLs", header = F)[,1])
CML <-as.vector(read.csv("CMLs", header = F)[,1])

dir.create(file.path(paste0(out_path,"listen/")), showWarnings = FALSE)

dir.create(file.path(paste0(out_path,"overlap_graph/")), showWarnings = FALSE)
dir.create(file.path(paste0(out_path,"heatmaps/")), showWarnings = FALSE)
dir.create(file.path(paste0(out_path,"vennies/")), showWarnings = FALSE)

dir.create(file.path(paste0(out_path,"heatmaps/", "all")), showWarnings = FALSE)
dir.create(file.path(paste0(out_path,"vennies/", "all")), showWarnings = FALSE)
cutoffs <- c(0, 0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for(cut in cutoffs){
  print(cut)
  dir.create(file.path(paste0(out_path,"overlap_graph/",cut)), showWarnings = FALSE)
}




combined_plotter(protein_matrix, "all", "_all_proteins", out_path = out_path)


