#############
# filtering #
#############

#############################
# summarising and filtering #
#############################

# collapse levels to give tidal and circdian data
tidal <- arrays_long
levels(tidal$time) <- c("HW", "LW", "HW", "LW")
tidal_average <- ddply(tidal, c("Probeset.ID", "gene", "time"), 
                       summarise, 
                       mean_expression = mean(expression), 
                       sd = sd(expression), 
                       n = length(expression), 
                       se = sd/sqrt(n))

circadian <- arrays_long
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



##########
# graphs #
##########


# Which genes are present (barplots of counts for each gene)
raw <- BaseBarGraph(perfect_hits)
raw + ggtitle("Summary of mitochondrial representation on NimbleGen arrays") + 
    ylab("Count")


# line graph of each contig average over biological replicates at each timepoint
#pdf("mitoexpression%03d.pdf", width = 8.3, height = 11.7, onefile = FALSE)

mito_colour_palette <-  c("chartreuse", "darkmagenta","turquoise", "indianred", "purple", "deepskyblue","firebrick", "gold", "gray86", "moccasin")
    
raw_timecourse <- BaseLineGraph(arrays_average) 
raw_timecourse + scale_color_manual("gene", values = mito_colour_palette) +
    ggtitle("Mitochondrial gene expression over whole time course")

# line graph for each contig: tidal
tidal <- BaseLineGraph(tidal_average)    
tidal + scale_color_manual("gene", values = mito_colour_palette) +
    ggtitle("Mitochondrial gene expression: tidal time")


# filtered graph: tidal (plus error bars)
tidal_filtered <- BaseLineGraph(tidal_average_filtered)                          
tidal_filtered + 
    ggtitle("Mitochondrial gene expression: tidal filtered for >0.5 expression change") +
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

# line graph for each contig: circadian
circa <- BaseLineGraph(circadian_average)
circa + scale_color_manual("gene", values = mito_colour_palette) + 
    ggtitle("Mitochondrial gene expression: circadian time") 
    

# filtered graph: circadian
circa_filtered <- BaseLineGraph(circadian_average_filtered)
circa_filtered + 
    ggtitle("Mitochondrial gene expression: circadian filtered for >0.5 expression change") + 
    geom_errorbar(aes(ymin=mean_expression-se, ymax=mean_expression+se), width=.1)

