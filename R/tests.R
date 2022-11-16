#Tests for calculateMI functions


library(tidyverse)
library(purrr)
library(data.table)
library(foreach)
library(doSNOW)
library(readxl)
library(lubridate)
library(Biostrings)

seq_times <- readRDS("../timeMI/data/seq_times_091122.rda")
seq_times_rndate <- seq_times %>%
  select(rowname, date)

h3n2_pol <- readAAStringSet("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Processing/pol_aa.fasta") %>%
  as.matrix()

h3n2_pol_small <- h3n2_pol[,1:100]

h3n2_mi_noHA <- calculateMI(h3n2_pol_small, weighted = FALSE)

h3n2_meaned_mi_noHA <- calculateMI(h3n2_pol_small, groups = seq_times_rndate)


#Compare
mi_weighted_res <- readRDS("../timeMI/data/MI_results_091122/h3n2_meaned_mi_noHA.rds") %>%
  filter(V1 <= 100 &
           V2 <= 100)

sum(h3n2_meaned_mi_noHA$mi == mi_weighted_res$mi)
sum(h3n2_meaned_mi_noHA$v1_entropy == mi_weighted_res$v1_entropy)


t <- as.matrix(paste(h3n2_pol[,4],h3n2_pol[,8]))
t2 <- paste(h3n2_pol[,4],h3n2_pol[,8])

rownames(t) <- rownames(h3n2_pol)
names(t2) <- rownames(h3n2_pol)
t3 <- as.matrix(t2)

t4 <- t3[c("EPI_ISL_200691", "EPI_ISL_200685"),]

t5 <- t3[c("EPI_ISL_200691"),]

names(t5) <- t4

prop.table(table(t4))


