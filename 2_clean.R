## @knitr clean
library(plyr)
library(reshape2)

# clean up and reshape RMA_data so that ddply and ggplot can be used on it

# get the perfect hits
perfect_hits <- subset(mito_contigs, evalue == 0)

#drop the rev com duplicates
perfect_hits <- perfect_hits[!grepl("rc_contig\\w\\w\\w\\w\\w", perfect_hits$contig_name, perl = TRUE),]
perf_hits_no_rrnL <- subset(perfect_hits, gene != "rrnL")

# select the perfect hit mitochondrial contigs
mito_arrays <- array_data[array_data$Probeset.ID %in% perfect_hits$contig_name,]

# add a gene name column
mito_arrays$gene_name <- perfect_hits[mito_arrays$Probeset.ID %in% perfect_hits$contig_name, "gene"]

without_rrnl <- subset(mito_arrays, gene_name != "rrnL")

##################
# data reshaping #
##################



# seperate out time points 
HW1 <- mito_arrays[, c(1, 2:4, 14)]
HW1 <- renaming(HW1)

LW1 <- mito_arrays[, c(1, 5:7, 14)]
LW1 <-renaming(LW1) 

HW2 <- mito_arrays[, c(1, 8:10, 14)]
HW2 <- renaming(HW2) 

LW2 <- mito_arrays[, c(1, 11:13, 14)]
LW2 <- renaming(LW2) 

# long replicate data
HW1 <- make_long(HW1)
LW1 <- make_long(LW1)
HW2 <- make_long(HW2)
LW2 <- make_long(LW2)

# add time column
HW1$time <- rep("HW1", length(HW1$Probeset.ID)) 
LW1$time <- rep("LW1", length(HW1$Probeset.ID)) 
HW2$time <- rep("HW2", length(HW1$Probeset.ID)) 
LW2$time <- rep("LW2", length(HW1$Probeset.ID)) 

#combine back to 1 dataframe
tidal_arrays.long <- rbind(HW1, LW1, HW2, LW2)
rm(HW1, LW1, HW2, LW2)


# rearrange columns
arrays_long <- tidal_arrays.long[c("Probeset.ID", "gene", "time", "replicate", "expression")]

# make time a factor
arrays_long$time <- factor(tidal_arrays.long$time, 
                           levels = c("HW1", "LW2", "HW2", "LW1"))

