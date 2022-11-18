#Functions for calculating MI and weighted MI

#' Calculate MI, raw or weighted
#'
#' This function calculates the MI, APC, and MIp from an MSA loaded with readMSA.
#' If weighted is set to FALSE, then only msa_matrix is needed as input. If
#' weighted is set to TRUE and weights is set to "equal", then a "groups"
#' dataframe is required that contains an 'ID' column that matches the msa_matrix
#' and a 'group' column. Otherwise, a weights dataframe is required that
#' contains a 'group' column and a 'weight' column.
#'
#'
#' @param msa_matrix Matrix format of MSA, from readMSA
#' @param weighted either TRUE, for weighted MI , or FALSE, for raw MI
#' @param groups a dataframe with 'ID' and 'group'
#' @param weights either "equal" or a dataframe with 'group' and 'weight'
#' @param ncores Number of cores for parallelization
#' @return A data.frame containing entropy, MI, APC, MIp, and z-score (MIp)
#' @importFrom foreach "%dopar%"
#' @export
calculateMI <- function(msa_matrix, weighted = FALSE, groups = NULL,
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
  
  cl <- snow::makeCluster(ncores, type="SOCK") # for 4 cores machine
  snow::clusterExport(cl, list = c("get_ent_unweighted",
                            "get_ent_weighted",
                            "calculateMI"),
                      )
  doSNOW::registerDoSNOW (cl)
  
  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # parallelization with vectorization
  joint_ents <- foreach::foreach(i = 1:nrow(joint_cols), .combine="c",
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
  snow::stopCluster(cl)
  
  return(format_mi(unlist(entropies), joint_ents, joint_cols))
}

#Putting entropy functions here because of weird cluster issues
#Get Shannon entropy
#' @export
get_ent_unweighted <- function(mat){
  t <- prop.table(table(mat))
  return(sum(-t*log2(t)))
}

#Get weighted Shannon entropy
#' @export
get_ent_weighted <- function(mat, df, weights){
  groups <- sort(unique(df$group))
  
  group_freqs <- list()
  for(i in 1:length(groups)){
    mat_i <- mat[df[df$group == groups[i],]$ID,]
    group_freqs[[i]] <- as.list(prop.table(table(mat_i)))
  }
  
  group_freqs_bind <- data.table::rbindlist(group_freqs, fill = TRUE)
  group_freqs_bind[is.na(group_freqs_bind)] <- 0
  
  if(is.null(dim(weights))){
    mean_freqs <- colMeans(group_freqs_bind)
  }
  else{
    mean_freqs <- colSums(group_freqs_bind*weights$weight)
  }
  
  return(sum(-log2(mean_freqs)*mean_freqs))
}


#Format MI results, regardless of weighting
format_mi <- function(entropies, joint_entropies, joint_cols){
  joint_cols$joint_entropy <- joint_entropies
  joint_cols$v1_entropy <- entropies[joint_cols$V1]
  joint_cols$v2_entropy <- entropies[joint_cols$V2]
  joint_cols$mi <- (joint_cols$v1_entropy+joint_cols$v2_entropy)-joint_cols$joint_entropy
  
  meanMI <- mean(joint_cols$mi)
  
  non_zero_e <- which(entropies!=0)
  mean_MIs <- numeric(length(non_zero_e))
  for (i in non_zero_e){
    mean_MIs[i] <- mean(joint_cols[joint_cols$V1 == i | joint_cols$V2 == i,]$mi)
  }
  
  joint_cols$Group <- paste("[",joint_cols$V1,";",joint_cols$V2, "]",sep = "")
  joint_cols$v1_meanMI <- mean_MIs[joint_cols$V1]
  joint_cols$v2_meanMI <- mean_MIs[joint_cols$V2]
  joint_cols$apc <- (joint_cols$v1_meanMI*joint_cols$v2_meanMI)/meanMI
  joint_cols$mip <- joint_cols$mi - joint_cols$apc
  mip_mean <- mean(joint_cols$mip)
  mip_sd <- sd(joint_cols$mip)
  joint_cols$z_score <- (joint_cols$mip-mip_mean)/mip_sd
  
  return(joint_cols)
}
