# reshaping functions

#separate out time points
renaming <- function(dataframe) {
    names(dataframe) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")    
}
    