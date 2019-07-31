## This code is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## 
## For a copy of the GNU General Public License, 
## see <http://www.gnu.org/licenses/>.

## (c) 2018 Sarah Humberg, Felix Schönbrodt, Maria Bloechl


# Function to identify redundant models
frm <- function(aic.results){
  
  # define list with all sets of models which are nested in each other, respectively
  nestinglist <- list(c("full", "somrr", "somfr", "omfr", "s7", "null"),
                      c("full", "somrr", "omrr", "omfr", "s7", "null"),
                      c("full", "s4", "pb", "null"), 
                      c("full", "s5", "s1", "null"), 
                      c("full", "s5", "pb", "null"), 
                      c("full", "s6", "pb", "null"), 
                      c("full", "s6", "s2", "null"), 
                      c("full", "s3", "s1", "null"), 
                      c("full", "s3", "s2", "null"))
  
  # create empty vector to be filled with names of models we might want to remove
  toremove <- c()
  
  # go through all rows of the table, starting with the second row
  for (i in 2:nrow(aic.results)){
    # compare this row i to each of the previous rows k
    for (k in 1:(i-1)){
      # Are the log-likelihoods of the two models essentially the same?
      similar_LL <- abs(aic.results[i,"LL"] - aic.results[k,"LL"]) < 1
      
      # Is model k nested in model i? --> test whether there is a vector in the 
      # nesting-list which contains both of the two models i and k
      modelnames <- c(as.character(aic.results[i,"Modnames"]), 
                      as.character(aic.results[k,"Modnames"]))
      matches <- c()
      
      for (j in 1:length(nestinglist)){
        matches = c(matches,sum(modelnames %in% nestinglist[[j]]))
      }
      
      nested <- any(matches==2)          
      
      # If both conditions are true (similar LL and nested), add the name of 
      # model i to the list of models we might want to remove
      if(similar_LL & nested){
        toremove <- c(toremove, as.character(aic.results[i,"Modnames"]))
      }
    }
  }
  
  toremove <- as.vector(unique(toremove))
  
  # Print warning if models were detected that we might want to remove
  if(length(toremove) > 0){
    warning(paste0("There were nested models with log-likelihood differences < 1. 
                   You might want to remove the following models from the set: ", 
                   "c(", paste0(toremove, collapse=", "), ")"))
  }
  
  res <- list(toremove = toremove)
  
}
