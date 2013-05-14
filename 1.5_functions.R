# reshaping functions

#rename
renaming <- function(df) {
    names(df) <- c("Probeset.ID","rep1", "rep2", "rep3", "gene")
    return(df)
}