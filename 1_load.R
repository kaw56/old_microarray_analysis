# load Brian's RMA_data

## @knitr load
library(ggplot2)
library(plyr)
library(reshape2)
library(knitr)
library(xtable)

# dataframe of array data
array_data <- read.table("RMA_DATA.txt", header = TRUE)

# dataframe of contigs used for the array that have mitochondrial sequence
mito_contigs <- read.table("454_contigs_for_arrays.fastaRESULTS_summary.txt",
                           col.names =c("contig_name", "gene", "evalue"))

