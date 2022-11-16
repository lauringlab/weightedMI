#Functions for calculating MI and weighted MI

library(data.table)
library(foreach)
library(doSNOW)

#Main function
calculateMI <- function(msa_matrix,
                        weighted = TRUE,
                        groups = NULL,
                        weights = "equal",
                        ncores = 2){
  
  entropies <- get_entropies(msa_matrix, weighted = weighted, groups = groups)
  non_zero_e <- which(entropies!=0)
  joint_cols <- as.data.frame(t(combn(non_zero_e, 2)))
  
  cl <- makeCluster(ncores, type="SOCK") # for 4 cores machine
  clusterExport(cl, c("get_ent_meaned", "get_ent_unweighted", "calculateMI"))
  registerDoSNOW (cl)
  
  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # parallelization with vectorization
  joint_ents <- foreach(i = 1:nrow(joint_cols), .combine="c",
                        .packages=c('data.table'),
                        .options.snow = opts) %dopar%
    {
      joint_mat <- as.matrix(paste(msa_matrix[,joint_cols[i, "V1"]],
                         msa_matrix[,joint_cols[i, "V2"]]))
      rownames(joint_mat) <- rownames(msa_matrix)
      
      if(weighted){
        get_ent_meaned(joint_mat, groups)
      }
      else{
        get_ent_unweighted(joint_mat)
      }
    }
  
  close(pb)
  stopCluster(cl)
  
  return(format_mi(unlist(entropies), joint_ents, joint_cols))
}

get_entropies <- function(mat, weighted = TRUE, groups = NULL, weights = NULL){
  entropies <- vector(length = ncol(mat))
  for(i in 1:ncol(mat)){
    mat_i <- as.matrix(mat[,i])
    rownames(mat_i) <- rownames(mat)
    
    if(weighted){
      entropies[i] = get_ent_meaned(mat_i, groups)
    }
    else{
      entropies[i] = get_ent_unweighted(mat_i)
    }
  }
  return(entropies)
}

#Get Shannon entropy
get_ent_unweighted <- function(mat){
  t <- prop.table(table(mat))
  
  return(sum(-t*log2(t)))
}

#Get equal-weighted Shannon entropy
get_ent_meaned <- function(mat, df){
  groups <- unique(df$date)
  
  t1f <- list()
  for(i in 1:length(groups)){
    mat_i <- mat[df[df$date == groups[i],]$rowname,]
    t1f[[i]] <- as.list(prop.table(table(mat_i)))
  }
  
  t1frbind <- rbindlist(t1f, fill = TRUE)
  t1frbind[is.na(t1frbind)] <- 0
  t1frw <- colMeans(t1frbind, na.rm = TRUE)
  
  return(sum(-log2(t1frw)*t1frw))
}

#Format MI results, regardless of weighting
format_mi <- function(entropies, joint_entropies, joint_cols){
  joint_cols$joint_entropy <- joint_entropies
  non_zero_e <- which(entropies!=0)
  
  v1_entropy <- numeric(nrow(joint_cols))
  v2_entropy <- numeric(nrow(joint_cols))
  for (i in 1:nrow(joint_cols)) {
    v1_entropy[i] <- entropies[joint_cols$V1[i]]
    v2_entropy[i] <- entropies[joint_cols$V2[i]]
  }
  
  joint_cols$v1_entropy <- v1_entropy
  joint_cols$v2_entropy <- v2_entropy
  joint_cols$mi <- (v1_entropy + v2_entropy) - joint_entropies
  
  meanMI <- mean(joint_cols$mi)
  
  mean_MIs <- numeric(length(non_zero_e))
  for (i in non_zero_e){
    mean_MIs[i] <- mean(pull(joint_cols[joint_cols$V1 == i | joint_cols$V2 == i,], mi))
  }
  
  v1_meanMI <- numeric(nrow(joint_cols))
  v2_meanMI <- numeric(nrow(joint_cols))
  for (i in 1:nrow(joint_cols)) {
    v1_meanMI[i] <- mean_MIs[joint_cols$V1[i]]
    v2_meanMI[i] <- mean_MIs[joint_cols$V2[i]]
  }
  
  joint_cols$v1_meanMI <- v1_meanMI
  joint_cols$v2_meanMI <- v2_meanMI
  joint_cols$apc <- (v1_meanMI*v2_meanMI)/meanMI
  joint_cols$mip <- joint_cols$mi - joint_cols$apc
  joint_cols$Group <- paste("[", joint_cols$V1, ";", joint_cols$V2, "]", sep = "")
  
  return(joint_cols)
}
