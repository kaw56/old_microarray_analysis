#############################
# summarising and filtering #
#############################

## @knitr filter

# collapse levels to give tidal and circdian data
tidal <- arrays_long
levels(tidal$time) <- c("HW", "LW", "HW", "LW")

circadian <- arrays_long
levels(circadian$time) <- c("day", "day", "night", "night")

# average of the biological replicates at each time point 
# (plus standard deviation and standard error)
arrays_average <- summary_stats(arrays_long)
circadian_average <- summary_stats(circadian)
tidal_average <- summary_stats(tidal)

# t-test to find significant changes
circadian_t_test <- perform_t_test(circadian)
tidal_t_test <- perform_t_test(tidal)

# filter out unchanging circadian
circadian_filter_set <-circadian_t_test[circadian_t_test$pvalue < 0.1 ,]
circadian_average_filtered <- circadian_average[circadian_average$Probeset.ID %in% circadian_filter_set$Probeset.ID,]

# filter out unchanging sequences tidal
tidal_filter_set <- tidal_t_test[ tidal_t_test$change < 0.2 ,]
tidal_average_filtered <- tidal_average[tidal_average$Probeset.ID %in% tidal_filter_set$Probeset.ID,]

# number of rrnL hits
rrnl_num <- nrow(perfect_hits) - nrow(perf_hits_no_rrnL)

##########
# tables #
##########

## @knitr contig_table


#perfect hits table
perf_table <- perf_hits_no_rrnL
# drop the evalues
perf_table$evalue <- NULL
# reorder columns
perf_table <- perf_table[c("gene", "contig_name")]
# sort by gene
perf_table <- perf_table[order(perf_table$gene),]

perf_table_print <- xtable(perf_table,
                            caption = "List of represented mitochondrial genes and their associated contigs")

print.xtable(perf_table_print, include.rownames = FALSE)

##########
# graphs #
##########

## @knitr graphs

# line graph of each contig average over biological replicates at each timepoint
raw_timecourse <- BaseLineGraph(arrays_average) 
raw_timecourse + scale_color_brewer("gene", palette = "Set3") +
    ggtitle("Mitochondrial gene expression: all time points")

# line graph for each contig: tidal
tidal <- BaseLineGraph(tidal_average)    
tidal + scale_color_brewer("gene", palette = "Set3") +
    ggtitle("Mitochondrial gene expression: tidal time")


# filtered graph: tidal (plus error bars)
tidal_filtered <- BaseLineGraph(tidal_average_filtered)                          
tidal_filtered + 
    ggtitle("Mitochondrial gene expression: tidal, significant change (p < 0.05)") +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

# line graph for each contig: circadian
circa <- BaseLineGraph(circadian_average)
circa + scale_color_brewer("gene", palette = "Set3") + 
    ggtitle("Mitochondrial gene expression: circadian time") 
    

# filtered graph: circadian
circa_filtered <- BaseLineGraph(circadian_average_filtered)
circa_filtered + 
    ggtitle("Mitochondrial gene expression: circadian, significant change (p < 0.05)") + 
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

