KNN <- function(taxa, matSim, K, minSim) {
    # K nearest neighbour (KNN) majority vote selection to identify most similar taxa
    similar <- matSim[taxa, ] %>%
                .[!names(.) %in% taxa] %>% # removing i from most similar targetSim
                .[order(., decreasing = TRUE)] %>%
                {
                    if(.[K+1] == .[K]) # if K + 1 == K, randomly sample a most similar taxa and pick K most similar taxa
                        c(.[which(. > .[K])], .[sample(which(. == .[K]))])[1:K]
                        else .[1:K]
                } %>%
                .[!. == 0 & . < minSim] # remove all similarities == 0 and similarities below minSim
    return(similar)
}
