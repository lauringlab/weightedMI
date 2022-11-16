#Utility functions

readMSA <- function(file_path){
  t <- readLines(file_path)
  id_index <- which(str_detect(t, ">"))
  ids <- gsub(">", "", t[id_index])
  nseqs <- length(id_index)
  
  seqs <- numeric(nseqs)
  for(i in 1:(nseqs-1)){
    seqs[i] <- paste(t[(id_index[i]+1):(id_index[i+1]-1)], collapse = "")
  }
  seqs[length(id_index)] <- paste(t[(id_index[nseqs]+1):length(t)], collapse = "")
  
  seq_mat <- matrix(unlist(strsplit(seqs, "")), ncol = nchar(seqs[1]), byrow = TRUE)
  rownames(seq_mat) <- ids
  
  return(seq_mat)
}

getNetworkInput <- function(mi_output){
  return(data.frame(
    a = mi_output$V1,
    b = mi_output$V2,
    strength = mi_output$mip
  ))
}











