#Tests for calculateMI functions
# seq_times <- readRDS("../timeMI/data/seq_times_091122.rda")
# seq_times_rndate <- seq_times[c('rowname', 'date')]
# 
# h3n2_pol <- readMSA("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Processing/pol_aa.fasta")
# 
# h3n2_pol_small <- h3n2_pol[,1:100]
# 
# h3n2_mi_noHA <- calculateMI(h3n2_pol_small, weighted = FALSE)
# 
# h3n2_meaned_mi_noHA <- calculateMI(h3n2_pol_small, groups = seq_times_rndate)
# h3n2_weighted_mi_noHA <- calculateMI(h3n2_pol_small, groups = seq_times_rndate)
# 
# 
# #Compare
# mi_weighted_res <- readRDS("../timeMI/data/MI_results_091122/h3n2_meaned_mi_noHA.rds")
# mi_weighted_res_small <- mi_weighted_res[mi_weighted_res$V1 <= 100 &
#                                            mi_weighted_res$V2 <= 100,]
# 
# sum(h3n2_weighted_mi_noHA$mi == mi_weighted_res_small$mi)
# sum(h3n2_weighted_mi_noHA$v1_entropy == mi_weighted_res_small$v1_entropy)
# 
# net <- getNetworkInput(h3n2_weighted_mi_noHA)





