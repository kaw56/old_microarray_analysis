# load Brian's RMA_data

## @knitr load
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(knitr)

# dataframe of array data
array_data <- read.table("eurydice_arrays/For R/RMA_DATA.txt", header = TRUE)

# dataframe of contigs used for the array that have mitochondrial sequence
mito_contigs <- read.table("Array_Analysis/Using_Ep_blast_v2.2.1/Briancontigs/454_contigs_for_arrays.fastaRESULTS_summary.txt",
                           col.names =c("contig_name", "gene", "evalue"))

