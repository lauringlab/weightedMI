#Utility functions

#' Read in an MSA in fasta format
#'
#' This function reads in an MSA in fasta format and generates an MSA matrix.
#' The MSA matrix can then be used as input for calculateMI
#'
#'
#' @param file_path path to fasta file
#' @return The MSA formatted as a matrix
#' @export
readMSA <- function(file_path){
  t <- readLines(file_path)
  id_index <- which(grepl(">", t))
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

#' Generate input format for associationsubgraphs network visualization
#'
#' This function takes the output from calculateMI and generates a data.frame
#' that can be used as input for the associationsubgraphs package
#'
#'
#' @param mi_output output from calculateMI
#' @return A dataframe containing columns 'a', 'b', and 'strength' (MIp)
#' @export
getNetworkInput <- function(mi_output){
  return(data.frame(
    a = mi_output$V1,
    b = mi_output$V2,
    strength = mi_output$mip
  ))
}








