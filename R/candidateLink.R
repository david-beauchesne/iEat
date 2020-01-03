candLink <- function(similar, candidates) {
        # Candidate links
        for(l in names(similar)) { # extracting source candidates
            if(l %in% candidates) { # if candidate is already in candidate list, add source' with wt to its weight
              candidates[which(candidates %in% l), 'weight'] <- as.numeric(candidates[which(candidates %in% l), 'weight']) + similar[l]
            } else {
                  candidates <- rbind(candidates, c(l,similar[l])) # if candidate is not in the list, add it source' with wt to its weight
            }
        }
    return(candidates)
}
