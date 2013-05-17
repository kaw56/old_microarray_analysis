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
circadian_filter_set <-circadian_t_test[circadian_t_test$pvalue < 0.05 ,]
circadian_average_filtered <- circadian_average[circadian_average$Probeset.ID %in% circadian_filter_set$Probeset.ID,]

# filter out unchanging sequences tidal
tidal_filter_set <- tidal_t_test[ tidal_t_test$pvalue < 0.05 ,]
tidal_average_filtered <- tidal_average[tidal_average$Probeset.ID %in% tidal_filter_set$Probeset.ID,]

# number of rrnL hits
rrnl_num <- nrow(high_hits) - nrow(hi_hits_no_rrnL)

##########
# tables #
##########

## @knitr gene_table

hi_table <- droplevels(mito_arrays$gene_name)
hi_table_print <- summary(hi_table)
hi_table_print <- as.data.frame(hi_table_print)
hi_table_print <- cbind(levels(hi_table), hi_table_print)
hi_table_print <- xtable(hi_table_print,
                         caption = "Number of contigs that will give information about gene expression for each mitochondrial gene")

print.xtable(hi_table_print, include.rownames = FALSE)

## @knitr circa_t_test
circa_t_table <- circadian_t_test
circa_t_table <- make_t_table(circa_t_table)

circa_t_table_print <- xtable(circa_t_table,
                              caption = "P values of t-test comparing day time and night time mitochondrial gene expression")

print.xtable(circa_t_table_print, include.rownames = FALSE) 


## @knitr tidal_t_test
tidal_t_table <- tidal_t_test
tidal_t_table <- make_t_table(tidal_t_table)

tidal_t_table_print <- xtable(tidal_t_table,
                              caption = "P values of t-test comparing high tide and low tide mitochondrial gene expression")

print.xtable(tidal_t_table_print, include.rownames = FALSE)

##########
# graphs #
##########

## @knitr whole_time_course

# line graph of each contig average over biological replicates at each timepoint
raw_timecourse <- BaseLineGraph(arrays_average) 
raw_timecourse + scale_color_brewer("gene", palette = "Set3") 

## @knitr tidal_contig_graph

# line graph for each contig: tidal
tidal <- BaseLineGraph(tidal_average)    
tidal + scale_color_brewer("gene", palette = "Set3") 


## @knitr tidal_filter_graph

# filtered graph: tidal (plus error bars)
tidal_filtered <- BaseLineGraph(tidal_average_filtered)                          
tidal_filtered +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)


## @knitr circadian_contig_graph

# line graph for each contig: circadian
circa <- BaseLineGraph(circadian_average)
circa + scale_color_brewer("gene", palette = "Set3")  
    
## @knitr circa_filter_graph

# filtered graph: circadian
circa_filtered <- BaseLineGraph(circadian_average_filtered)
circa_filtered + 
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

## @knitr
