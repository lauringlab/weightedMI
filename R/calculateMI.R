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
  if(weighted){
    mi <- calculate_weighted_mi(msa_matrix, groups, weights, ncores)
  }
  else{
    mi <- calculate_raw_mi(msa_matrix, ncores)
  }
  
  return(mi)
}

get_entropies <- function(mat, weighted = TRUE, groups = NULL, weights = NULL){
  entropies <- vector(length = ncol(mat))
  if(weighted){
    for(i in 1:ncol(mat)){
      entropies[i] = get_ent_meaned(i, groups, mat)
    }
  }
  else{
    for(i in 1:ncol(mat)){
      entropies[i] = get_ent_unweighted(i, mat)
    }
  }
  return(entropies)
}

#Calculate MI without any weighting
calculate_raw_mi <- function(mat, ncores){
  entropies <- get_entropies(mat, weighted = FALSE)
  non_zero_e <- which(entropies!=0)
  
  joint_cols <- as.data.frame(t(combn(non_zero_e, 2)))
  
  cl <- makeCluster(ncores, type="SOCK") # for 4 cores machine
  clusterExport(cl, c("get_jointent_unweighted", "calculate_raw_mi"))
  registerDoSNOW (cl)
  
  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # parallelization with vectorization
  joint_ents <- foreach(i = 1:nrow(joint_cols), .combine="c", .packages=c('data.table'),
                        .options.snow = opts) %dopar%
    {
      get_jointent_unweighted(joint_cols[i, "V1"], joint_cols[i, "V2"], mat)
    }
  
  close(pb)
  stopCluster(cl)
  
  return(format_mi(entropies, joint_ents, joint_cols))
}

#Get Shannon entropy
get_ent_unweighted <- function(residue, mat){
  t <- prop.table(table(mat[,residue]))
  
  return(sum(-t*log2(t)))
}

#Get joint entropy
get_jointent_unweighted <- function(residue_a,residue_b, mat){
  jointcol <- as.matrix(paste(mat[,residue_a],mat[,residue_b],sep=""))
  
  t <- prop.table(table(jointcol))
  
  return(sum(-t*log2(t)))
}

#Caculate MI with weights (equal or defined)
calculate_weighted_mi <- function(mat, groups, weights, ncores){
  ents <- get_entropies(mat, groups = groups)
  
  non_zero_e <- which(ents!=0)
  
  print("doing combn step")
  joint_cols <- as.data.frame(t(combn(non_zero_e, 2)))
  
  print("starting clusters")
  cl <- makeCluster(ncores, type="SOCK") # for 4 cores machine
  clusterExport(cl, c("get_jointent_meaned", "calculate_weighted_mi"))
  registerDoSNOW (cl)
  
  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # parallelization with vectorization
  joint_ents <- foreach(i = 1:nrow(joint_cols), .combine="c", .packages=c('data.table'),
                        .options.snow = opts) %dopar%
    {
      get_jointent_meaned(joint_cols[i, "V1"], joint_cols[i, "V2"], groups, mat)
    }
  
  close(pb)
  stopCluster(cl)
  
  ents <- unlist(ents)
  return(format_mi(ents, joint_ents, joint_cols))
}

#Get equal-weighted Shannon entropy
get_ent_meaned <- function(pos, df, mat){
  groups <- unique(df$date)
  
  t1f <- list()
  for(i in 1:length(groups)){
    mat_i <- mat[df[df$date == groups[i],]$rowname,]
    if(is.null(dim(mat_i))){
      t1f[[i]] = as.list(1)
      names(t1f[[i]]) <- mat_i[pos]
    }
    else{
      t1f[[i]] <- as.list(prop.table(table(mat_i[,pos])))
    }
  }
  
  t1frbind <- rbindlist(t1f, fill = TRUE)
  t1frbind[is.na(t1frbind)] <- 0
  t1frw <- colMeans(t1frbind, na.rm = TRUE)
  
  return(sum(-log2(t1frw)*t1frw))
}

#Get equal-weighted joint entropy
get_jointent_meaned <- function(a, b, df, mat){
  groups <- unique(df$date)

  t12f <- list()
  for(i in 1:length(groups)){
    mat_i <- mat[df[df$date == groups[i],]$rowname,]
    if(is.null(dim(mat_i))){
      t12f[[i]] = as.list(1)
      names(t12f[[i]]) <- paste(mat_i[a], mat_i[b])
    }
    else{
      t12f[[i]] <- as.list(prop.table(table(paste(mat_i[,a], mat_i[,b]))))
    }
  }
  
  t12frbind <- rbindlist(t12f, fill = TRUE)
  t12frbind[is.na(t12frbind)] <- 0
  t12frw <- colMeans(t12frbind, na.rm = TRUE)
  
  return(sum(-log2(t12frw)*t12frw))
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
