# load Brian's RMA_data

## @knitr load
library(ggplot2)
library(plyr)
library(reshape2)
library(knitr)
library(xtable)

# dataframe of array data
array_data <- read.table("Data/RMA_DATA.txt", header = TRUE)

# dataframe of key between probes and contigs
key <- read.table("Data/array_probeset_key.txt", header = TRUE)

# dataframe of contigs used for the array that have mitochondrial sequence
mito_contigs <- read.table("Data/Array_probeset.txtRESULTS_summary.txt",
                           col.names =c("probe", "gene", "score", "length", "start", "end"))

