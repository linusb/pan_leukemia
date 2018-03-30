protein_mapper <- function(binders, swissprotfile){

  library("seqinr")
  swissprot <- read.fasta(swissprotfile, seqtype = "AA", as.string = T) 
  number_of_peptides <- length(unique(binders[,"Seq"]))
  counter <- 0
  peptide_protein_list <- sapply(unique(binders[,"Seq"]),function(x){
    if(counter %% 20 == 0){
      print(paste0(counter,"/",number_of_peptides))
    }
    #proteins <- c()
    proteins <- sapply(swissprot, function(prot){
      if(grepl(prot[[1]], pattern = x, fixed = T)){
        return(strsplit(attr(prot, "name"), split = "\\|")[[1]][2])
      }else{
        return(NA)
      }
    })
    counter <<- counter +1
    return(proteins[!is.na(proteins)])
  })
  return(peptide_protein_list)
}