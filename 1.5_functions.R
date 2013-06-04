# reshaping functions

library(reshape2)

## @knitr func1
#rename
renaming <- function(df) {
    names(df) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")
    return(df)
}

# make the data long format
make_long <- function(df) {
    melt(df, id.vars=c("Probeset.ID", "gene"), 
         variable.name="replicate", 
         value.name="expression")
}

