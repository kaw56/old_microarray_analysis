## @knitr clean
library(plyr)
library(reshape2)

# clean up and reshape RMA_data so that ddply and ggplot can be used on it

# get the perfect hits
perfect_hits <- subset(mito_contigs, evalue == 0)

#drop the rev com duplicates
perfect_hits_forward <- perfect_hits[!grepl("rc_contig\\w\\w\\w\\w\\w", perfect_hits$contig_name, perl = TRUE),]
perf_hits_forward_no_rrnL <- subset(perfect_hits_forward, gene != "rrnL")

perfect_hits_reverse<- perfect_hits[grepl("rc_contig\\w\\w\\w\\w\\w", perfect_hits$contig_name, perl = TRUE),]
perf_hits_reverse_no_rrnL <- subset(perfect_hits_reverse, gene != "rrnL")

# select the perfect hit mitochondrial contigs
mito_arrays_forward <- array_data[array_data$Probeset.ID %in% perfect_hits_forward$contig_name,]
mito_arrays_reverse <- array_data[array_data$Probeset.ID %in% perfect_hits_reverse$contig_name,]
# add a gene name column
mito_arrays_forward$gene_name <- perfect_hits[mito_arrays$Probeset.ID %in% perfect_hits_forward$contig_name, "gene"]
mito_arrays_reverse$gene_name <- perfect_hits[mito_arrays$Probeset.ID %in% perfect_hits_reverse$contig_name, "gene"]

without_rrnl <- subset(mito_arrays, gene_name != "rrnL")

##################
# data reshaping #
##################



# seperate out time points
HW1_forward <- mito_arrays_forward[, c(1, 2:4, 14)]
HW1_forward <- renaming(HW1_forward)

HW1_reverse <- mito_arrays_reverse[, c(1, 2:4, 14)]
HW1_reverse <- renaming(HW1_reverse)

LW1_forward <- mito_arrays_forward[, c(1, 5:7, 14)]
LW1_forward <-renaming(LW1_forward)

LW1_reverse <- mito_arrays_reverse[, c(1, 5:7, 14)]
LW1_reverse <-renaming(LW1_reverse)

HW2_forward <- mito_arrays_forward[, c(1, 8:10, 14)]
HW2_forward <- renaming(HW2_forward)

HW2_reverse <- mito_arrays_reverse[, c(1, 8:10, 14)]
HW2_reverse <- renaming(HW2_reverse)

LW2_forward <- mito_arrays_forward[, c(1, 11:13, 14)]
LW2_forward <- renaming(LW2_forward)

LW2_reverse <- mito_arrays_reverse[, c(1, 11:13, 14)]
LW2_reverse <- renaming(LW2_reverse)

# long replicate data
HW1_forward <- make_long(HW1_forward)
LW1_forward <- make_long(LW1_forward)
HW2_forward <- make_long(HW2_forward)
LW2_forward <- make_long(LW2_forward)

HW1_reverse <- make_long(HW1_reverse)
LW1_reverse <- make_long(LW1_reverse)
HW2_reverse <- make_long(HW2_reverse)
LW2_reverse <- make_long(LW2_reverse)

# add time column
HW1_forward$time <- rep("HW1", length(HW1_forward$Probeset.ID))
LW1_forward$time <- rep("LW1", length(HW1_forward$Probeset.ID))
HW2_forward$time <- rep("HW2", length(HW1_forward$Probeset.ID))
LW2_forward$time <- rep("LW2", length(HW1_forward$Probeset.ID))

HW1_reverse$time <- rep("HW1", length(HW1_reverse$Probeset.ID))
LW1_reverse$time <- rep("LW1", length(HW1_reverse$Probeset.ID))
HW2_reverse$time <- rep("HW2", length(HW1_reverse$Probeset.ID))
LW2_reverse$time <- rep("LW2", length(HW1_reverse$Probeset.ID))

#combine back to 1 dataframe
tidal_arrays_forward.long <- rbind(HW1_forward, LW1_forward, HW2_forward, LW2_forward)

tidal_arrays_reverse.long <- rbind(HW1_reverse, LW1_reverse, HW2_reverse, LW2_reverse)

rm(HW1_forward, LW1_forward, HW2_forward, LW2_forward, HW1_reverse, LW1_reverse, HW2_reverse, LW2_reverse)


# rearrange columns
arrays_long_forward <- tidal_arrays_forward.long[c("Probeset.ID", "gene", "time", "replicate", "expression")]
arrays_long_reverse <- tidal_arrays_reverse.long[c("Probeset.ID", "gene", "time", "replicate", "expression")]

# make time a factor
arrays_long_forward$time <- factor(tidal_arrays_forward.long$time,
                           levels = c("HW1", "LW2", "HW2", "LW1"))
arrays_long_reverse$time <- factor(tidal_arrays_reverse.long$time,
                           levels = c("HW1", "LW2", "HW2", "LW1"))

