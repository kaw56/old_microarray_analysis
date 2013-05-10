# load Brian's RMA_data

setwd("/Users/oj/Documents/4-5-4 Arrays/eurydice_arrays/") 
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)

# dataframe of array data
array.data <- read.table("For R/RMA_DATA.txt", header = TRUE)

# dataframe of contigs used for the array that have mitochondrial sequence
mito_contigs <- read.table("454_contigs_for_arrays.fasta-RESULTS.txt", 
                           col.names =c("contig_name", "gene", "start", "stop"))

