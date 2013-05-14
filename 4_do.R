##########
# graphs #
##########


# Which genes are present (barplots of counts for each gene)
raw <- BaseBarGraph(perfect_hits)
raw + ggtitle("Summary of mitochondrial representation on NimbleGen arrays") + 
    ylab("Count")





# line graph of each contig average over biological replicates at each timepoint
#pdf("mitoexpression%03d.pdf", width = 8.3, height = 11.7, onefile = FALSE)

mito_colour_palette <-  c("lightsteelblue", "chartreuse", "chocolate1", "darkmagenta",
                     "navy", "turquoise", "indianred", "olivedrab1", "deepskyblue",
                     "darkorchid", "firebrick", "gold", "gray86", "moccasin")
    
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

#dev.off()