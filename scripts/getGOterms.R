library(GOSim)
library(GO.db)

# Purpose: To find all GO terms containing a given string and their offspring
# Input: Search string, all lowercase
# Output: A vector (unique) of the GO IDs of the identified terms and all their offspring
# By: Tovah Markowitz, PhD
# Date: 6/7/24

# Example usage: searchTerm <- "immun"
# goTerms <- getGOterms(searchTerm)

getGOterms <- function(searchTerm) {
# get all GO terms and their descriptions
goterms <- unlist(Term(GOTERM))
# find the indices of the terms whose descriptions include the term(s) of interest
Idx <- grep(searchTerm, goterms)

# get the GO terms of relevance
myGOterms <- names(goterms)[Idx]

# Use GOSim to find every offspring of all GO terms, default is human
GOoffspring <- getOffsprings()

# Find the offspring
# Note: not all terms in GO.bp will be found in GOSim
# Terminal GO terms will have offspring listed as NA
Idx2 <- which(names(GOoffspring) %in% myGOterms)
myOffspring <- GOoffspring[Idx2]

# Get a unique list of terms and offspring
myGOterms2 <- na.omit(unique(c(myGOterms, unlist(myOffspring))))
return(myGOterms2)
}
