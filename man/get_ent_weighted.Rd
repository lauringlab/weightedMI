% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateMI.R
\name{get_ent_weighted}
\alias{get_ent_weighted}
\title{Get Shannon entropy, weighted}
\usage{
get_ent_weighted(mat, df, weights)
}
\arguments{
\item{mat}{Matrix format of MSA, from readMSA}

\item{df}{dataframe containing ID and group, where ID matches rownames of mat}

\item{weights}{a dataframe containing group and weight, specifying the weight for each group in df}
}
\value{
A table containing the AA frequencies from mat
}
\description{
This function calculates the weighted Shannon entropy of a column or pair of
columns from an MSA
}
