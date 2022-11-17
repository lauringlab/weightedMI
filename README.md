# weightedMI

A package for computing weighted MI and entropy from multiple sequence alignments

## Installation

Install weightedMI from github using devtools:

    
    # install.packages("devtools")
    devtools::install_github("lauringlab/weightedMI")
    

## Usage

We will use an example subset of SARS-CoV-2 Spike sequences. We can use the `readMSA` function to read in a fasta file:

    
    library(weightedMI)
    # spike_msa <- readMSA()
    

### Calculating the raw MI

We can use the `calculateMI` function to calculate the raw MI

    
    # rawMI <- calculateMI(spike_msa)
    # head(rawMI)
    

### Calculating the equal-weighted MI

The `calculateMI` function to can also calculated the equal-weighted MI when a "groups" data frame is also provided

    
    # head(groups)
    # equalweightedMI <- calculateMI(spike_msa, weighted = TRUE, groups = groups)
    # head(equalweightedMI)
    

### Calculating the weighted MI

If weights are known for each group, these can be provided as a "weights" data frame

    
    # head(weights)
    # weightedMI <- calculateMI(spike_msa, weighted = TRUE, groups = groups,
    #                           weights = weights)
    # head(weightedMI)
    

### Network output

Lastly, the output from `calculateMI` can be re-formatted to work with the [associationsubgraphs](https://github.com/nstrayer/associationsubgraphs) package

    
    # netMI <- getNetworkInput(weightedMI)
    # head(netMI)
    
