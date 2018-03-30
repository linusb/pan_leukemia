# creates a heatmap based on binary peptide data

library(dplyr)
#library(data.table)
library(gplots)

source("prediction_reader.R")
source("allotype_binder_table_calculator.R")
source("protein_mapper.R")
source("protein_table_calculator.R")
source("heatmap_plotter.R")
source("peptide_group_mapper.R")
source("graph_plotter.R")
source("calculate_overlap_matrix.R")

# reading predicted peptides
path <- "../../predictions/"

binders <- prediction_reader(path = path)
write.csv(file = "all_binders.csv", x= binders, row.names = F)
save(x= binders, file= "all_binders.RData")



# extract patient data and read tumor and benign data
all_patients <- unique(binders[,"Patient"])
write.csv(all_patients, file="patients.csv", row.names = F)
tumor <- as.vector(read.csv("tumor.csv", header = F)[,1])
benign <- as.vector(read.csv("benign.csv", header = F)[,1])

# get all alleles
alleles <- unique(sapply(colnames(binders)[4:ncol(binders)], function(x){
  return(strsplit(x, split = "_")[[1]][1])
}))


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
write.table(statistic_table, file = "number_of_binders_per_patient.csv", sep = ';', dec = ",", row.names = T, col.names = NA)

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

group_color_mapper <- function(x){
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
  }
}
cutoffs <- c(0, 0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for(cut in cutoffs){
  print(cut)
  dir.create(file.path(paste0("overlap_graph/", cut)), showWarnings = FALSE)
}
dir.create(file.path(paste0("listen/")), showWarnings = FALSE)
dir.create(file.path(paste0("listen/tumor_specific_peptides/")), showWarnings = FALSE)
dir.create(file.path(paste0("listen/tumor_specific_peptides_no_one_hit_wonders/")), showWarnings = FALSE)
dir.create(file.path(paste0("listen/not_tumor_specific_peptides_no_one_hit_wonders/")), showWarnings = FALSE)

for(allele in alleles){
  # get all allotype typings
  #allele <- "A.02.01"
  dir.create(file.path(paste0("heatmaps/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0("listen/tumor_specific_peptides/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0("listen/tumor_specific_peptides_no_one_hit_wonders/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0("listen/not_tumor_specific_peptides_no_one_hit_wonders/", allele)), showWarnings = FALSE)
  
  
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
  
  # write tumor specific peptides for each allele an patient
  for(p in colnames(tumor_specific_allotype_binder_table)){
    p_specific_tumor_binders <- tumor_specific_allotype_binder_table[,p]
    write.table(x = row.names(tumor_specific_allotype_binder_table)[which(p_specific_tumor_binders==1)], file = paste0("listen/tumor_specific_peptides/", allele,"/", p,".csv"), col.names = F, row.names = F)
  }
  
  # remove one hit wonders
  tumor_specific_allotype_binder_table <- tumor_specific_allotype_binder_table[which(rowSums(tumor_specific_allotype_binder_table)>1),]
  
  # write tumor specific peptides for each allele an patient
  for(p in colnames(tumor_specific_allotype_binder_table)){
    p_specific_tumor_binders <- tumor_specific_allotype_binder_table[,p]
    write.table(x = row.names(tumor_specific_allotype_binder_table)[which(p_specific_tumor_binders==1)], file = paste0("listen/tumor_specific_peptides_no_one_hit_wonders/", allele,"/", p,".csv"), col.names = F, row.names = F)
  }
  
  group_colors <- sapply(colnames(tumor_specific_allotype_binder_table), group_color_mapper)
  
  heatmap_plotter(tumor_specific_allotype_binder_table, allele = allele, suffix = "_peptide_tumor_specific", group_colors = group_colors)
  
  
  ######################
  # TUMOR SPECIFIC     #
  # peptide_Treshold   #
  ######################
  
  # remove one hit wonders
  tumor_specific_allotype_binder_table <- tumor_specific_allotype_binder_table[which(rowSums(tumor_specific_allotype_binder_table)>1),]
  # 5 threshold
  tumor_specific_allotype_binder_table_cut <- tumor_specific_allotype_binder_table[,which(colSums(tumor_specific_allotype_binder_table)>5)]
  
  group_colors <- sapply(colnames(tumor_specific_allotype_binder_table_cut), group_color_mapper)
  heatmap_plotter(tumor_specific_allotype_binder_table_cut, allele = allele, suffix = "_peptide_tumor_specific_5_cut", group_colors = group_colors)
  
  ######################
  # NOT TUMOR SPECIFIC #
  ######################
  # remove one hit wonders
  tumor_allotype_binder_table <- tumor_allotype_binder_table[which(rowSums(tumor_allotype_binder_table)>1),]
  # write tumor specific peptides for each allele an patient
  for(p in colnames(tumor_allotype_binder_table)){
    p_specific_tumor_binders <- tumor_allotype_binder_table[,p]
    write.table(x = row.names(tumor_allotype_binder_table)[which(p_specific_tumor_binders==1)], file = paste0("listen/not_tumor_specific_peptides_no_one_hit_wonders/", allele,"/", p,".csv"), col.names = F, row.names = F)
  }
  
  group_colors <- sapply(colnames(tumor_allotype_binder_table), group_color_mapper)
  heatmap_plotter(tumor_allotype_binder_table, allele = allele, suffix = "_peptide", group_colors = group_colors)
  
  
  ######################
  # NOT TUMOR SPECIFIC #
  # peptide_cutoff     #
  ######################
  # remove one hit wonders
  tumor_allotype_binder_table <- tumor_allotype_binder_table[which(rowSums(tumor_allotype_binder_table)>1),]
  tumor_allotype_binder_table_cut <- tumor_allotype_binder_table[,which(colSums(tumor_allotype_binder_table)>100)]
  
  group_colors <- sapply(colnames(tumor_allotype_binder_table_cut), group_color_mapper)
  heatmap_plotter(tumor_allotype_binder_table_cut, allele = allele, suffix = "_peptide_100cut", group_colors = group_colors)
  
}


dir.create(file.path(paste0("listen/tumor_specific_proteins/")), showWarnings = FALSE)
dir.create(file.path(paste0("listen/tumor_specific_proteins_no_one_hit_wonders/")), showWarnings = FALSE)
dir.create(file.path(paste0("listen/single_tumor_specific_protein/")), showWarnings = FALSE)

for(allele in alleles){
  #allele <- "A.02.01"
  dir.create(file.path(paste0("heatmaps/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0("listen/tumor_specific_proteins/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0("listen/tumor_specific_proteins_no_one_hit_wonders/", allele)), showWarnings = FALSE)
  dir.create(file.path(paste0("listen/single_tumor_specific_protein/", allele)), showWarnings = FALSE)
  
  
  #############################
  # tumor UNspecific proteins #
  #############################
  # get all allotype typings
  allotype_binder_table <- allotype_binder_table_calculator(binders, allele)
  
  protein_table <- protein_table_calculator(allotype_binder_table, peptide_protein_list)
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_tumor_UNspecific_protein", group_colors = group_colors)
  
  # NO one hit wonders
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_no_one_hit_tumor_UNspecific_protein", group_colors = group_colors)
  ####################################
  # only peptides with one protein   #
  ####################################
  allotype_binder_table <- allotype_binder_table_calculator(binders = binders, allele = allele)
  
  # tumor specific
  tumor_unspecific_allotype_binder_table <- allotype_binder_table[0,]
  for(pep in row.names(allotype_binder_table)){
    if(sum(benign_allotype_binder_table[pep,])==0){
      tumor_unspecific_allotype_binder_table <- rbind(tumor_unspecific_allotype_binder_table, allotype_binder_table[pep,])
    } 
  }
  
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
  
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_single_tumor_UNspecific_protein", group_colors = group_colors)
  
  
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_no_one_hit_single_tumor_UNspecific_protein", group_colors = group_colors)
  
  
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
  
  
  # write tumor specific peptides for each allele an patient
  for(p in colnames(protein_table)){
    p_specific_tumor_binders <- protein_table[,p]
    write.table(x = row.names(protein_table)[which(p_specific_tumor_binders==1)], file = paste0("listen/tumor_specific_proteins/", allele,"/", p,".csv"), col.names = F, row.names = F)
  }
  
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_tumor_specific_protein", group_colors = group_colors)
  
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  
  # write tumor specific peptides for each allele an patient
  for(p in colnames(protein_table)){
    p_specific_tumor_binders <- protein_table[,p]
    write.table(x = row.names(protein_table)[which(p_specific_tumor_binders==1)], file = paste0("listen/tumor_specific_proteins_no_one_hit_wonders/", allele,"/", p,".csv"), col.names = F, row.names = F)
  }
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_no_one_hit_tumor_specific_protein", group_colors = group_colors)
  
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
  # write tumor specific peptides for each allele an patient
  for(p in colnames(protein_table)){
    p_specific_tumor_binders <- protein_table[,p]
    write.table(x = row.names(protein_table)[which(p_specific_tumor_binders==1)], file = paste0("listen/single_tumor_specific_protein/", allele,"/", p,".csv"), col.names = F, row.names = F)
  }
  
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_single_tumor_specific_protein", group_colors = group_colors)
  
  
  protein_table <- protein_table[which(rowSums(protein_table)>1),]
  
  group_colors <- sapply(colnames(protein_table), group_color_mapper)
  heatmap_plotter(protein_table, allele = allele, suffix = "_no_one_hit_single_tumor_specific_protein", group_colors = group_colors)
  
}
