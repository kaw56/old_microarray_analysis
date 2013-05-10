# graph functions for microarray analysis
library(ggplot2)

BaseBarGraph <- function(data.set) {
    ggplot(data.set, 
           aes(x=gene, fill=gene)) + 
        geom_bar() + 
        theme_bw() +
        guides(fill=FALSE) + 
        theme(panel.grid.major.x=element_blank(), 
              panel.grid.minor.x=element_blank(), 
              axis.title.x=element_blank(), 
              plot.title=element_text(face = "bold"))
}

BaseLineGraph <- function(data.set) {
    ggplot(data.set, aes(x=time, 
                              y=mean_expression, 
                              color = gene, 
                              group=Probeset.ID))  + 
        geom_line() + 
        theme_bw() + 
        theme(axis.title.x=element_blank(), 
              plot.title=element_text(face = "bold"))
}