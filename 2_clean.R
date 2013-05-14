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
arrays.long <- tidal_arrays.long[c("Probeset.ID", "gene", "time", "replicate", "expression")]

# make time a factor
arrays.long$time <- factor(tidal_arrays.long$time, 
                           levels = c("HW1", "LW2", "HW2", "LW1"))

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

circadian_ttest <- ddply(circadian, c("Probeset.ID", "gene"),
                         summarise,
                         pvalue = t.test(expression ~ time)$p.value)




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
tidal_filter_set <- tidal_change[ tidal_change$change > 0.3 ,]
tidal_average_filtered <- tidal_average[tidal_average$Probeset.ID %in% tidal_filter_set$Probeset.ID,]

# filter out unchanging circadian
circadian_average_wide <-dcast(circadian_average, Probeset.ID + gene ~ time, value.var="mean_expression")
circadian_change <- ddply(circadian_average_wide, c("Probeset.ID", "gene"), summarise, change = abs(day - night))
circadian_filter_set <-circadian_change[ circadian_change$change > 0.5 ,]
circadian_average_filtered <- circadian_average[circadian_average$Probeset.ID %in% circadian_filter_set$Probeset.ID,]