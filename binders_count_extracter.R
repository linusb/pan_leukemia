# extracts the number of binders per sample

binders <- read.csv("complete_substraction/all_binders.csv")


alleles <- unique(sapply(colnames(binders)[4:ncol(binders)],function(x){
  strsplit(x, split="_")[[1]][1]
}))
patients <- unique(binders$Patient)
result <- data.frame(matrix(ncol=length(alleles)+1))
colnames(result) <- c("Patient",alleles)
for(p in patients){
  p_binders <- binders[which(binders$Patient == p),]
  number_of_binders <- NA
  for(a in alleles){
    a_binders <- p_binders[which(p_binders[,paste0(a,"_aff")] <= 500 | p_binders[,paste0(a,"_rank")]<= 2),]
    if(is.na(number_of_binders)[1]){
      number_of_binders <- nrow(a_binders)
    }else{
    number_of_binders <- c(number_of_binders, nrow(a_binders))
    }
  }
  result <- rbind(result, c(p,number_of_binders))
}

write.csv(file = "binders_count.csv", x= result[2:nrow(result),], row.names = F)
