library(plyr)
library(reshape2)

# clean up and reshape RMA_data so that ddply and ggplot can be used on it

# remove duplicate contigs (some contigs have multiple hits to mitochondrial
# sequence)
de_duplicate <- subset(mito_contigs, !duplicated(mito_contigs$contig_name))

# remove the over represented rrnL ribosomal RNA 
clean_mito <- subset(mito_contigs, gene!="rrnL")
clean_mito <- droplevels(clean_mito)

# select the mitochondrial contigs
mito_arrays <- array_data[array_data$Probeset.ID %in% mito_contigs$contig_name,]

# add a gene name column
mito_arrays$gene_name <- de_duplicate[mito_arrays$Probeset.ID %in% de_duplicate$contig_name, "gene"]

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
arrays.long <- tidal_arrays.long[c("Probeset.ID", "gene", "time", "replicate", "expression")]

# make time a factor
arrays.long$time <- factor(tidal_arrays.long$time, 
                           levels = c("HW1", "LW1", "HW2", "LW2"))

#############################
# summarising and filtering #
#############################

# collapse levels to give tidal data
tidal <- arrays.long
levels(tidal$time) <- c("HW", "LW", "HW", "LW")
tidal_average <- ddply(tidal, c("Probeset.ID", "gene", "time"), 
                       summarise, 
                       mean_expression = mean(expression), 
                       sd = sd(expression), 
                       n = length(expression), 
                       se = sd/sqrt(n))

# collapse to give circadian data
circadian <- arrays.long
levels(circadian$time) <- c("day", "day", "night", "night")
circadian_average <- ddply(circadian, c("Probeset.ID", "time", "gene"), 
                           summarise, 
                           mean_expression = mean(expression), 
                           sd = sd(expression), 
                           n = length(expression), 
                           se = sd/sqrt(n))

# average biological replicates at time point (plus standard deviation and standard error)
arrays_average <- ddply(arrays.long, c("Probeset.ID", "gene", "time"), 
                        summarise, 
                        mean_expression = mean(expression), 
                        sd = sd(expression), 
                        n = length(expression), 
                        se = sd/sqrt(n))

# filter out unchanging sequences tidal
tidal_average_wide <-dcast(tidal_average, Probeset.ID + gene ~ time, value.var="mean_expression")
tidal_change <- ddply(tidal_average_wide, c("Probeset.ID", "gene"), summarise, change = abs(HW - LW))
tidal_filter_set <- tidal_change[ tidal_change$change > 0.5 ,]
tidal_average_filtered <- tidal_average[tidal_average$Probeset.ID %in% tidal_filter_set$Probeset.ID,]

# filter out unchanging circadian
circadian_average_wide <-dcast(circadian_average, Probeset.ID + gene ~ time, value.var="mean_expression")
circadian_change <- ddply(circadian_average_wide, c("Probeset.ID", "gene"), summarise, change = abs(day - night))
circadian_filter_set <-circadian_change[ circadian_change$change > 0.5 ,]
circadian_average_filtered <- circadian_average[circadian_average$Probeset.ID %in% circadian_filter_set$Probeset.ID,]