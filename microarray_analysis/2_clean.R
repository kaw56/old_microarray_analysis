## @knitr clean
library(plyr)
library(reshape2)

# clean up and reshape RMA_data so that ddply and ggplot can be used on it

# selecting representative probes
# finding probes that are not duplicated
mito_contigs <- mito_contigs[!duplicated(mito_contigs$probe),]

# specificity
strong <- mito_contigs[mito_contigs$score >= 48,]
moderate <- mito_contigs[mito_contigs$score >= 32 & mito_contigs$score <= 47,]
weak <- mito_contigs[mito_contigs$score < 32,]

# Add contig name 
strong$contig_name <- key[key$ProbeID %in% strong$probe, "PrimaryAccession"]

mito_arrays <- array_data[array_data$Probeset.ID %in% strong$contig_name,]

# add gene name
mito_arrays$gene <- strong[strong$contig_name %in% mito_arrays$Probeset.ID, "gene"]

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


