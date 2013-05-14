library(plyr)
library(reshape2)

# clean up and reshape RMA_data so that ddply and ggplot can be used on it

# get the perfect hits
perfect_hits <- subset(mito_contigs, evalue == 0)

# select the perfect hit mitochondrial contigs
mito_arrays <- array_data[array_data$Probeset.ID %in% perfect_hits$contig_name,]

# add a gene name column
mito_arrays$gene_name <- perfect_hits[mito_arrays$Probeset.ID %in% perfect_hits$contig_name, "gene"]

without_rrnl <- subset(mito_arrays, gene_name != "rrnL")

##################
# data reshaping #
##################

# (seperate out each timepoint to average expression) (there is probably a loopy way to do this but r confuses and frightens me when it comes to loops)
HW1 <- mito_arrays[, c(1, 2:4, 14)]
names(HW1) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")
LW1 <- mito_arrays[, c(1, 5:7, 14)]
names(LW1) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")
HW2 <- mito_arrays[, c(1, 8:10, 14)]
names(HW2) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")
LW2 <- mito_arrays[, c(1, 11:13, 14)]
names(LW2) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")

# long replicate data
HW1 <- melt(HW1, id.vars=c("Probeset.ID", "gene"), 
            variable.name="replicate", 
            value.name="expression")
LW1 <- melt(LW1, id.vars=c("Probeset.ID", "gene"), 
            variable.name="replicate", 
            value.name="expression")
HW2 <- melt(HW2, id.vars=c("Probeset.ID", "gene"), 
            variable.name="replicate", 
            value.name="expression")
LW2 <- melt(LW2, id.vars=c("Probeset.ID", "gene"), 
            variable.name="replicate", 
            value.name="expression")

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

