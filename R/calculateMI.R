#Functions for calculating MI and weighted MI

library(data.table)
library(foreach)
library(doSNOW)

#Main function
calculateMI <- function(msa_matrix, weighted = TRUE, groups = NULL,
                        weights = "equal", ncores = 2){
  
  entropies <- vector(length = ncol(msa_matrix))
  for(i in 1:ncol(msa_matrix)){
    mat_i <- as.matrix(msa_matrix[,i])
    rownames(mat_i) <- rownames(msa_matrix)
    if(weighted){
      entropies[i] = get_ent_weighted(mat_i, groups, weights)
    }
    else{
      entropies[i] = get_ent_unweighted(mat_i)
    }
  }
  
  non_zero_e <- which(entropies!=0)
  joint_cols <- as.data.frame(t(combn(non_zero_e, 2)))
  
  cl <- makeCluster(ncores, type="SOCK") # for 4 cores machine
  clusterExport(cl, c("get_ent_unweighted", "get_ent_weighted","calculateMI"))
  registerDoSNOW (cl)
  
  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # parallelization with vectorization
  joint_ents <- foreach(i = 1:nrow(joint_cols), .combine="c",
                        .packages=c('data.table'), .options.snow = opts) %dopar%
    {
      joint_mat <- as.matrix(paste(msa_matrix[,joint_cols[i, "V1"]],
                         msa_matrix[,joint_cols[i, "V2"]]))
      rownames(joint_mat) <- rownames(msa_matrix)
      
      if(weighted){
        get_ent_weighted(joint_mat, groups, weights = weights)
      }
      else{
        get_ent_unweighted(joint_mat)
      }
    }
  close(pb)
  stopCluster(cl)
  
  return(format_mi(unlist(entropies), joint_ents, joint_cols))
}

#Get Shannon entropy
get_ent_unweighted <- function(mat){
  t <- prop.table(table(mat))
  return(sum(-t*log2(t)))
}

#Get weighted Shannon entropy
get_ent_weighted <- function(mat, df, weights){
  groups <- unique(df$date)
  
  group_freqs <- list()
  for(i in 1:length(groups)){
    mat_i <- mat[df[df$date == groups[i],]$rowname,]
    group_freqs[[i]] <- as.list(prop.table(table(mat_i)))
  }
  
  group_freqs_bind <- rbindlist(group_freqs, fill = TRUE)
  group_freqs_bind[is.na(group_freqs_bind)] <- 0
  
  if(weights == "equal"){
    mean_freqs <- colMeans(group_freqs_bind)
  }
  else{
    mean_freqs <- colSums(group_freqs_bind*weights)
  }
  
  return(sum(-log2(mean_freqs)*mean_freqs))
}

#Format MI results, regardless of weighting
format_mi <- function(entropies, joint_entropies, joint_cols){
  joint_cols$joint_entropy <- joint_entropies
  joint_cols$v1_entropy <- entropies[joint_cols$V1]
  joint_cols$v2_entropy <- entropies[joint_cols$V2]
  joint_cols$mi <- (joint_cols$v1_entropy + joint_cols$v2_entropy) - joint_cols$joint_entropy
  
  meanMI <- mean(joint_cols$mi)
  
  non_zero_e <- which(entropies!=0)
  mean_MIs <- numeric(length(non_zero_e))
  for (i in non_zero_e){
    mean_MIs[i] <- mean(joint_cols[joint_cols$V1 == i | joint_cols$V2 == i,]$mi)
  }
  
  joint_cols$v1_meanMI <- mean_MIs[joint_cols$V1]
  joint_cols$v2_meanMI <- mean_MIs[joint_cols$V2]
  joint_cols$apc <- (joint_cols$v1_meanMI*joint_cols$v2_meanMI)/meanMI
  joint_cols$mip <- joint_cols$mi - joint_cols$apc
  joint_cols$Group <- paste("[", joint_cols$V1, ";", joint_cols$V2, "]", sep = "")
  
  return(joint_cols)
}
