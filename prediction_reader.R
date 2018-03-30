
library(dplyr)
prediction_reader <- function(path){
  files <- list.files(path = path, pattern = "*.txt")
  
  data <- list()
  for(file in files){
    print(file)
    temp_file <- read.csv(paste0(path,file), sep = "\t", stringsAsFactors = F)
    temp_file_split_data <- temp_file
    temp_file_split_data <- cbind(strsplit(x = file, split = "_")[[1]][1],temp_file_split_data)
    colnames(temp_file_split_data)[1] <- "Patient"
    temp_file_split_data[, "Patient"] <- data.frame(as.character(temp_file_split_data[, "Patient"]), stringsAsFactors=FALSE)
    
    data[[file]] <- temp_file_split_data
  }
  
  #all_predictions <- rbindlist(l = data, fill = T, use.names = T)
  all_predictions <- as.data.frame(bind_rows(data))
  
  remove(data)
  
  all_predictions_split <- data.frame(matrix(ncol = ((ncol(all_predictions)*2)-3), nrow = nrow(all_predictions)))
  headers <- c()
  for(i in seq(from=4, to = ncol(all_predictions), by=1)){
    headers <- c(headers, paste0(colnames(all_predictions)[i],"_aff"))
    headers<- c(headers, paste0(colnames(all_predictions)[i],"_rank"))
  }
  colnames(all_predictions_split) <- c(colnames(all_predictions)[1:3],headers)
  all_predictions_split[,1:3] <- all_predictions[,1:3]
  
  for(i in seq(from=4, to = ncol(all_predictions), by=1)){
    print(colnames(all_predictions)[i])
    print("aff")
    all_predictions_split[, paste0(colnames(all_predictions)[i], "_aff")] <- sapply(all_predictions[,i], 
                                                                                    function(x){
                                                                                      if(!is.na(x)){
                                                                                        return(as.numeric(strsplit(gsub(x=x, pattern = "[()]", replacement = ""), split = ',')[[1]][1]))
                                                                                      }else{
                                                                                        return(NA)
                                                                                      }
                                                                                    })
    print("rank")
    all_predictions_split[, paste0(colnames(all_predictions)[i], "_rank")] <- sapply(all_predictions[,i], function(x){
      if(!is.na(x)){
        return(as.numeric(strsplit(gsub(x=x, pattern = "[()]", replacement = ""), split = ',')[[1]][2]))
      }else{
        return(NA)
      }
    })
  }
  
  
  remove(all_predictions)
  save(x = all_predictions_split, file = paste0(out_path,"all_predictions.RData"))
  write.csv(x = all_predictions_split, file = paste0(out_path,"all_predictions.csv"), row.names = F)
  
  binders <- all_predictions_split[which(apply(cbind(
    apply(all_predictions_split[, which(grepl(colnames(all_predictions_split), pattern = "aff"))
                                ]<=500, 1, any)
    ,apply(all_predictions_split[, which(grepl(colnames(all_predictions_split), pattern = "rank"))
                                 ]<=2, 1, any)
  ),1, any)),]
  remove(all_predictions_split)
  return(binders)
}