#############################
# summarising and filtering #
#############################

## @knitr filter_forward

# collapse levels to give tidal and circdian data
tidal_forward <- arrays_long_forward
levels(tidal_forward$time) <- c("HW", "LW", "HW", "LW")

circadian_forward <- arrays_long_forward
levels(circadian_forward$time) <- c("day", "day", "night", "night")

# average of the biological replicates at each time point
# (plus standard deviation and standard error)
arrays_average_forward <- summary_stats(arrays_long_forward)
circadian_average_forward <- summary_stats(circadian_forward)
tidal_average_forward <- summary_stats(tidal_forward)

# t-test to find significant changes
circadian_t_test_forward <- perform_t_test(circadian_forward)
tidal_t_test_forward <- perform_t_test(tidal_forward)

# filter out unchanging circadian
circadian_filter_set_forward <-circadian_t_test_forward[circadian_t_test_forward$pvalue < 0.05 ,]
circadian_average_filtered_forward <- circadian_average_forward[circadian_average_forward$Probeset.ID %in% circadian_filter_set_forward$Probeset.ID,]

# filter out unchanging sequences tidal
tidal_filter_set_forward <- tidal_t_test_forward[ tidal_t_test_forward$pvalue < 0.05 ,]
tidal_average_filtered_forward <- tidal_average_forward[tidal_average_forward$Probeset.ID %in% tidal_filter_set_forward$Probeset.ID,]

# number of rrnL hits
rrnl_num_forward <- nrow(perfect_hits_forward) - nrow(perf_hits_forward_no_rrnL)

##########
# tables #
##########

# looking at the represented genes
probe_per_gene <- ddply(strong, "gene", summarise, n = length(score))

## @knitr contig_table_forward
perf_table_forward <- perf_hits_forward_no_rrnL
# drop the evalues
perf_table_forward$evalue <- NULL
# reorder columns
perf_table_forward <- perf_table_forward[c("gene", "contig_name")]
# sort by gene
perf_table_forward <- perf_table_forward[order(perf_table_forward$gene),]

perf_table_print_forward <- xtable(perf_table_forward,
                            caption = "List of represented mitochondrial genes and their associated contigs")

print.xtable(perf_table_print_forward, include.rownames = FALSE)

## @knitr circa_t_test_forward
circa_t_table_forward <- circadian_t_test_forward
circa_t_table_forward <- make_t_table(circa_t_table_forward)

circa_t_table_print_forward <- xtable(circa_t_table_forward,
                              caption = "P values of t-test comparing day time and night time mitochondrial gene expression")

print.xtable(circa_t_table_print_forward, include.rownames = FALSE)


## @knitr tidal_t_test_forward
tidal_t_table_forward <- tidal_t_test_forward
tidal_t_table_forward <- make_t_table(tidal_t_table_forward)

tidal_t_table_print_forward <- xtable(tidal_t_table_forward,
                              caption = "P values of t-test comparing high tide and low tide mitochondrial gene expression")

print.xtable(tidal_t_table_print_forward, include.rownames = FALSE)

##########
# graphs #
##########

## @knitr whole_time_course_forward

# line graph of each contig average over biological replicates at each timepoint
raw_timecourse_forward <- BaseLineGraph(arrays_average_forward)
raw_timecourse_forward + scale_color_brewer("gene", palette = "Set3")

## @knitr tidal_contig_graph_forward

# line graph for each contig: tidal
tidal_forward <- BaseLineGraph(tidal_average_forward)
tidal_forward + scale_color_brewer("gene", palette = "Set3")


## @knitr tidal_filter_graph_forward

# filtered graph: tidal (plus error bars)
tidal_filtered_forward <- BaseLineGraph(tidal_average_filtered_forward)
tidal_filtered_forward +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)


## @knitr circadian_contig_graph_forward

# line graph for each contig: circadian
circa_forward <- BaseLineGraph(circadian_average_forward)
circa_forward + scale_color_brewer("gene", palette = "Set3")

## @knitr circa_filter_graph_forward

# filtered graph: circadian
circa_filtered_forward <- BaseLineGraph(circadian_average_filtered_forward)
circa_filtered_forward +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

## @knitr


#############################
# summarising and filtering #
#############################

## @knitr filter_reverse

# collapse levels to give tidal and circdian data
tidal_reverse <- arrays_long_reverse
levels(tidal_reverse$time) <- c("HW", "LW", "HW", "LW")

circadian_reverse <- arrays_long_reverse
levels(circadian_reverse$time) <- c("day", "day", "night", "night")

# average of the biological replicates at each time point
# (plus standard deviation and standard error)
arrays_average_reverse <- summary_stats(arrays_long_reverse)
circadian_average_reverse <- summary_stats(circadian_reverse)
tidal_average_reverse <- summary_stats(tidal_reverse)

# t-test to find significant changes
circadian_t_test_reverse <- perform_t_test(circadian_reverse)
tidal_t_test_reverse <- perform_t_test(tidal_reverse)

# filter out unchanging circadian
circadian_filter_set_reverse <-circadian_t_test_reverse[circadian_t_test_reverse$pvalue < 0.05 ,]
circadian_average_filtered_reverse <- circadian_average_reverse[circadian_average_reverse$Probeset.ID %in% circadian_filter_set_reverse$Probeset.ID,]

# filter out unchanging sequences tidal
tidal_filter_set_reverse <- tidal_t_test_reverse[ tidal_t_test_reverse$pvalue < 0.05 ,]
tidal_average_filtered_reverse <- tidal_average_reverse[tidal_average_reverse$Probeset.ID %in% tidal_filter_set_reverse$Probeset.ID,]

# number of rrnL hits
rrnl_num_reverse <- nrow(perfect_hits_reverse) - nrow(perf_hits_reverse_no_rrnL)

##########
# tables #
##########

## @knitr contig_table_reverse
perf_table_reverse <- perf_hits_reverse_no_rrnL
# drop the evalues
perf_table_reverse$evalue <- NULL
# reorder columns
perf_table_reverse<- perf_table_reverse[c("gene", "contig_name")]
# sort by gene
perf_table_reverse <- perf_table_reverse[order(perf_table_reverse$gene),]

perf_table_print_reverse <- xtable(perf_table_reverse,
                           caption = "List of represented mitochondrial genes and their associated contigs")

print.xtable(perf_table_print_reverse, include.rownames = FALSE)

## @knitr circa_t_test_reverse
circa_t_table_reverse <- circadian_t_test_reverse
circa_t_table_reverse <- make_t_table(circa_t_table_reverse)

circa_t_table_print_reverse <- xtable(circa_t_table_reverse,
                              caption = "P values of t-test comparing day time and night time mitochondrial gene expression")

print.xtable(circa_t_table_print_reverse, include.rownames = FALSE)


## @knitr tidal_t_test
tidal_t_table_reverse <- tidal_t_test_reverse
tidal_t_table_reverse <- make_t_table(tidal_t_table_reverse)

tidal_t_table_print_reverse <- xtable(tidal_t_table_reverse,
                              caption = "P values of t-test comparing high tide and low tide mitochondrial gene expression")

print.xtable(tidal_t_table_print_reverse, include.rownames = FALSE)

##########
# graphs #
##########

## @knitr whole_time_course_reverse

# line graph of each contig average over biological replicates at each timepoint
raw_timecourse_reverse <- BaseLineGraph(arrays_average_reverse)
raw_timecourse_reverse + scale_color_brewer("gene", palette = "Set3")

## @knitr tidal_contig_graph_reverse

# line graph for each contig: tidal
tidal_reverse <- BaseLineGraph(tidal_average_reverse)
tidal_reverse + scale_color_brewer("gene", palette = "Set3")


## @knitr tidal_filter_graph_reverse

# filtered graph: tidal (plus error bars)
tidal_filtered_reverse <- BaseLineGraph(tidal_average_filtered_reverse)
tidal_filtered_reverse +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)


## @knitr circadian_contig_graph_reverse

# line graph for each contig: circadian
circa_reverse <- BaseLineGraph(circadian_average_reverse)
circa_reverse + scale_color_brewer("gene", palette = "Set3")

## @knitr circa_filter_graph_reverse

# filtered graph: circadian
circa_filtered_reverse<- BaseLineGraph(circadian_average_filtered_reverse)
circa_filtered_reverse +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

## @knitr
