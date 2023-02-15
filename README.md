# weightedMI

A package for computing weighted MI and entropy from multiple sequence alignments

## Installation

Install weightedMI from github using devtools:

    
    # install.packages("devtools")
    devtools::install_github("lauringlab/weightedMI")
    

## Usage

We will use an example subset of SARS-CoV-2 Spike sequences. We can use the `readMSA` function to read in a fasta file:

    
    library(weightedMI)
    pb2_fasta <- system.file("extdata", "pb2_sample.fasta", package = "weightedMI")
    pb2 <- readMSA(pb2_fasta)
    

### Calculating the raw MI

We can use the `calculateMI` function to calculate the raw MI

    
    rawMI <- calculateMI(pb2)
    

### Calculating the equal-weighted MI

The `calculateMI` function to can also calculated the equal-weighted MI when a "groups" data frame is also provided

    
    data(pb2_meta)
    equalweightedMI <- calculateMI(pb2, weighted = TRUE, groups = pb2_meta)
    
If weights are known for each group, these can be provided as a "weights" data frame

    

### Network output

Lastly, the output from `calculateMI` can be re-formatted to work with the [associationsubgraphs](https://github.com/nstrayer/associationsubgraphs) package

    
    netMI <- getNetworkInput(equalweightedMI)
    
