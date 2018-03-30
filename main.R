# creates a heatmap based on binary peptide data

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

# reading predicted peptides
path <- "../../predictions/"
out_path <- "scored_substraction/"
dir.create(file.path(paste0("scored_substraction/")), showWarnings = FALSE)


binders <- prediction_reader(path = path)
write.csv(file = paste0(out_path,"all_binders.csv"), x= binders, row.names = F)
save(x= binders, file= paste0(out_path,"all_binders.RData"))



# extract patient data and read tumor and benign data
all_patients <- unique(binders[,"Patient"])
write.csv(all_patients, file=paste0(out_path,"patients.csv"), row.names = F)
tumor <- as.vector(read.csv(paste0(out_path,"tumor.csv"), header = F)[,1])
benign <- as.vector(read.csv(paste0(out_path,"benign.csv"), header = F)[,1])

# get all alleles
alleles <- unique(sapply(colnames(binders)[4:ncol(binders)], function(x){
  return(strsplit(x, split = "_")[[1]][1])
}))

#####################
# Protein mapping   #
#####################
# mapping peptides to proteins

peptide_protein_list <- protein_mapper(binders, "swissprotHUMANwoi_130927.fasta")
save(x= peptide_protein_list, file = paste0(out_path,"peptide_protein_list.RData"))
#load(file = "peptide_protein_list.RData")




# write statistic tables
statistic_table <- data.frame(matrix(ncol= length(all_patients), nrow = length(alleles)))
row.names(statistic_table) <- alleles
colnames(statistic_table) <- all_patients
for(p in all_patients){
  patient_binders <- binders[which(binders$Patient == p),]
  for(a in alleles){
    statistic_table[a,p] <- nrow(patient_binders[which(patient_binders[,paste0(a,"_aff")]<=500 | patient_binders[,paste0(a,"_rank")]<=2 ),])
  }
}
write.table(statistic_table, file = paste0(out_path,"number_of_binders_per_patient.csv"), sep = ';', dec = ",", row.names = T, col.names = NA)

# Heatmap properties
hclust.ave <- function(x) hclust(x, method = "complete")
mypalette <- colorRampPalette(c('white','black'))(n=2)
distfun <- function(x) dist(x, method = "binary")

dir.create(file.path(paste0(out_path,"listen/")), showWarnings = FALSE)

dir.create(file.path(paste0(out_path,"overlap_graph/")), showWarnings = FALSE)
dir.create(file.path(paste0(out_path,"heatmaps/")), showWarnings = FALSE)
dir.create(file.path(paste0(out_path,"vennies/")), showWarnings = FALSE)
cutoffs <- c(0, 0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for(cut in cutoffs){
  print(cut)
  dir.create(file.path(paste0(out_path,"overlap_graph/",cut)), showWarnings = FALSE)
}

for(allele in alleles){
  print(allele)
  # get all allotype typings
  #allele <- "A.02.01"
  dir.create(file.path(paste0(out_path,"heatmaps/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0(out_path,"vennies/", allele)), showWarnings = FALSE)
  
  
  allotype_binder_table <- allotype_binder_table_calculator(binders = binders, allele = allele)
  
  tumor_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% tumor)]
  benign_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% benign)]
  
  ######################
  # TUMOR SPECIFIC     #
  ######################
  # extracting only tumor specific peptides
  tumor_specific_allotype_binder_table <- tumor_allotype_binder_table[0,]
  for(pep in row.names(tumor_allotype_binder_table)){
    if(sum(benign_allotype_binder_table[pep,])==0){
      tumor_specific_allotype_binder_table <- rbind(tumor_specific_allotype_binder_table, tumor_allotype_binder_table[pep,])
    } 
  }
  
  
  combined_plotter(tumor_specific_allotype_binder_table, allele, "_peptide_tumor_specific", out_path)
  # write tumor specific peptides for each allele an patient
  patient_peptide_writer(tumor_specific_allotype_binder_table, allele = allele, suffix = "peptide_tumor_specific", out_path)
  
  # remove one hit wonders
  tumor_specific_allotype_binder_table <- tumor_specific_allotype_binder_table[which(rowSums(tumor_specific_allotype_binder_table)>1),]
  
  combined_plotter(tumor_specific_allotype_binder_table, allele, "_peptide_tumor_specific_no_one_hit", out_path)
  patient_peptide_writer(tumor_specific_allotype_binder_table, allele = allele, suffix = "peptide_tumor_specific_no_one_hit", out_path)
  
  ######################
  # TUMOR SPECIFIC     #
  # peptide_Treshold   #
  ######################
  
  # remove one hit wonders
  tumor_specific_allotype_binder_table <- tumor_specific_allotype_binder_table[which(rowSums(tumor_specific_allotype_binder_table)>1),]
  # 5 threshold
  tumor_specific_allotype_binder_table_cut <- tumor_specific_allotype_binder_table[,which(colSums(tumor_specific_allotype_binder_table)>5)]
  
  combined_plotter(tumor_specific_allotype_binder_table_cut, allele, "_peptide_tumor_specific_5_cut", out_path)
  patient_peptide_writer(tumor_specific_allotype_binder_table_cut, allele = allele, suffix = "peptide_tumor_specific_5_cut", out_path)
  
  ######################
  # NOT TUMOR SPECIFIC #
  ######################
  # remove one hit wonders
  tumor_allotype_binder_table <- tumor_allotype_binder_table[which(rowSums(tumor_allotype_binder_table)>1),]
  
  
  combined_plotter(tumor_allotype_binder_table, allele, "_peptide_tumor_UNspecific", out_path)
  patient_peptide_writer(tumor_allotype_binder_table, allele = allele, suffix = "peptide_tumor_UNspecific", out_path)
  
  ######################
  # NOT TUMOR SPECIFIC #
  # peptide_cutoff     #
  ######################
  # remove one hit wonders
  tumor_allotype_binder_table <- tumor_allotype_binder_table[which(rowSums(tumor_allotype_binder_table)>1),]
  tumor_allotype_binder_table_cut <- tumor_allotype_binder_table[,which(colSums(tumor_allotype_binder_table)>100)]
  
  combined_plotter(tumor_allotype_binder_table_cut, allele, "_peptide_tumor_UNspecific_100cut", out_path)
  patient_peptide_writer(tumor_allotype_binder_table_cut, allele = allele, suffix = "peptide_tumor_UNspecific_100cut", out_path)

  ######################
  # PROTEINS           #
  ######################

  
  #############################
  # tumor UNspecific proteins #
  #############################
  # get all allotype typings
  allotype_binder_table <- allotype_binder_table_calculator(binders, allele)
  tumor_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% tumor)]
  
  protein_table <- protein_table_calculator(tumor_allotype_binder_table, peptide_protein_list)
  combined_plotter(protein_table, allele, "_protein_tumor_UNspecific", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_UNspecific", out_path)
  # NO one hit wonders
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  combined_plotter(protein_table, allele, "_protein_tumor_UNspecific_no_one_hit", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_UNspecific_no_one_hit", out_path)
  
  ####################################
  # only peptides with one protein   #
  ####################################
  allotype_binder_table <- allotype_binder_table_calculator(binders = binders, allele = allele)
  tumor_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% tumor)]
  
  # tumor unspecific
  tumor_unspecific_allotype_binder_table <- tumor_allotype_binder_table
  
  # get all possible proteins
  proteins <- c()
  for(pep in row.names(tumor_unspecific_allotype_binder_table)){
    if(length(peptide_protein_list[[pep]])==1){
      proteins <- unique(c(proteins,peptide_protein_list[[pep]]))
    }
  }
  # create protein table
  protein_table <- data.frame(matrix(0, nrow = length(proteins), ncol = ncol(tumor_unspecific_allotype_binder_table)))
  row.names(protein_table) <- proteins
  colnames(protein_table) <- colnames(tumor_unspecific_allotype_binder_table)
  
  for(patient in colnames(tumor_unspecific_allotype_binder_table)){
    temp_patient <- tumor_unspecific_allotype_binder_table[,patient]
    # remove all 0 rows
    peptides <- row.names(tumor_unspecific_allotype_binder_table)[which(temp_patient==1)]
    # extract patient proteins
    patient_proteins <-c()
    for(pep in peptides){
      if(length(peptide_protein_list[[pep]])==1){
        patient_proteins <- unique(c(patient_proteins,peptide_protein_list[[pep]]))
      }
    }
    # set protein to 1
    for(patient_prot in patient_proteins){
      protein_table[patient_prot,patient] <- 1
    }
  }
  
  combined_plotter(protein_table, allele, "_protein_tumor_UNspecific_single", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_UNspecific_single", out_path)
  # remove one hit wonders
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  
  combined_plotter(protein_table, allele, "_protein_tumor_UNspecific_no_one_hit_single", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_UNspecific_no_one_hit_single", out_path)
  ###########################
  # tumor specific proteins #
  ###########################
  # get all allotype typings
  allotype_binder_table <- allotype_binder_table_calculator(binders, allele)
  
  tumor_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% tumor)]
  benign_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% benign)]
  # tumor specific
  tumor_specific_allotype_binder_table <- tumor_allotype_binder_table[0,]
  for(pep in row.names(tumor_allotype_binder_table)){
    if(sum(benign_allotype_binder_table[pep,])==0){
      tumor_specific_allotype_binder_table <- rbind(tumor_specific_allotype_binder_table, tumor_allotype_binder_table[pep,])
    } 
  }
  
  protein_table <- protein_table_calculator(tumor_specific_allotype_binder_table, peptide_protein_list)
  
  
  combined_plotter(protein_table, allele, "_protein_tumor_specific", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_specific", out_path)
  
  # remove one hit wonders
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  
  combined_plotter(protein_table, allele, "_protein_tumor_specific_no_one_hit", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_specific_no_one_hit", out_path)
  
  ####################################
  # only peptides with one protein   #
  ####################################
  allotype_binder_table <- allotype_binder_table_calculator(binders = binders, allele = allele)
  
  tumor_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% tumor)]
  benign_allotype_binder_table <- allotype_binder_table[,which(colnames(allotype_binder_table) %in% benign)]
  
  # tumor specific
  tumor_specific_allotype_binder_table <- tumor_allotype_binder_table[0,]
  for(pep in row.names(tumor_allotype_binder_table)){
    if(sum(benign_allotype_binder_table[pep,])==0){
      tumor_specific_allotype_binder_table <- rbind(tumor_specific_allotype_binder_table, tumor_allotype_binder_table[pep,])
    } 
  }
  
  # get all possible proteins
  proteins <- c()
  for(pep in row.names(tumor_specific_allotype_binder_table)){
    if(length(peptide_protein_list[[pep]])==1){
      proteins <- unique(c(proteins,peptide_protein_list[[pep]]))
    }
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
      if(length(peptide_protein_list[[pep]])==1){
        patient_proteins <- unique(c(patient_proteins,peptide_protein_list[[pep]]))
      }
    }
    # set protein to 1
    for(patient_prot in patient_proteins){
      protein_table[patient_prot,patient] <- 1
    }
  }
  
  combined_plotter(protein_table, allele, "_protein_tumor_specific_single", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_specific_single", out_path)
  #remove one hit wonders
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  combined_plotter(protein_table, allele, "_protein_tumor_specific_no_one_hit_single", out_path)
  patient_peptide_writer(protein_table, allele = allele, suffix = "protein_tumor_specific_no_one_hit_single", out_path)
}
