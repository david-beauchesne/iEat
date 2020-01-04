#' @name iEat
#'
#' @title Instance-based machine learning method to predict biotic interactions
#'
#' @description This method was published in Beauchesne et al. (2016) Life & Environment 66(3-4):333-342. The steps of the algorithm are visually presented in Figure 1 of the paper.
#'
#' @param S0 Matrix, catalogue of empirical data used to infer predictions
#' @param S1 Vector of taxa forming networking for which topology is predicted
#' @param S2 Vector of taxa in S1 for which we wish to predict resources (if unspecified, S2 == S1 and the whole network is predicted)
#' @param sourceSim Matrix (numeric), source similarity matrix between S1 taxa and the union of S0 and S1 taxa, structure sourceSim[unique(S0,S1), S1]
#' @param targetSim Matrix (numeric), source similarity matrix between S1 taxa and the union of S0 and S1 taxa, structure sourceSim[unique(S0,S1), S1] (if unspecified, targetSim == sourceSim and no similarity distinction between resources and consumers)
#' @param K Integer, how many neighbours for K nearest neighbour evaluation
#' @param minSim Integer, minimum similarity value accepted to consider taxa as similar (implemented to avoid unrealistic interactions)
#' @param minWt Integer, minimum weight for candidate source to become a predicted source
#' @param predict String, specifies whether the predictions are made from the "full algorithm", the "catalogue" or the "predictive" contribution. If unspecified, predict == 'full algorithm'. See Beauchesne et al. (2016) for more details. If predict == 'catalogue', the methodology corresponds to the approach presented by Gray et al. (2015).
#'
#' @return
#' A dataframe with source taxa for which target predictions are made, target infered from catalogue data (empirical) and target infered from KNN algorithm
#'
#' @author
#' David Beauchesne
#'
#' @importFrom magrittr %>%
#'
#' @example
#' # Simulated data for testing
#' ncat <- 100
#' npred <- 10
#' S0 <- paste0('Taxon_',1:ncat) %>%
#'         data.frame(taxon = .,
#'             target = replicate(n = length(.), expr = paste(sample(., round(runif(1,1,8))), collapse = ' | ')),
#'             source = replicate(n = length(.), expr = paste(sample(., round(runif(1,1,8))), collapse = ' | ')),
#'             row.names = .,
#'             stringsAsFactors = FALSE)
#' S1 <- as.character(sample(S0[,'taxon'], npred))
#' S2 <- S1
#' sourceSim <- targetSim <- matrix(nrow = nrow(S0), ncol = length(S1), data = runif(nrow(S0) * length(S1), min = 0, max = 1), dimnames = list(S0[,'taxon'],S1))
#' # Predictions on simulated data
#' iEat_bin(S0, S1, S2, sourceSim)
#'
#' @rdname iEat_bin
#'
#' @export

# /TODO: Tanimoto: NAs or ""?

iEat <- function(S0, S1, S2 = S1, sourceSim, targetSim = sourceSim, K = 5, minSim = 0.3, minWt = 0.5, predict = 'full algorithm') {

    #Checkups for data structure
    if (sum(!S0[, 'taxon'] %in% rownames(sourceSim)) > 0 | sum(!S0[, 'taxon'] %in% rownames(targetSim)) > 0)
        stop('Taxa in S0 have to be included as rows in the similarity matrices.')
    if (sum(!S1 %in% rownames(sourceSim)) > 0 | sum(!S1 %in% rownames(targetSim)) > 0)
        stop('Taxa in S1 have to be included as rows in the similarity matrices.')
    if (sum(!S2 %in% rownames(sourceSim)) > 0 | sum(!S2 %in% rownames(targetSim)) > 0)
        stop('Taxa in S2 have to be included as columns in the similarity matrices.')
    if (sum(!S2 %in% S1) > 0)
        stop('Taxa in S2 have to be included in S1')
    if (!is.numeric(sourceSim) | !is.numeric(targetSim))
        stop('Similarity matrices have to be numerical')


    # Empty matrix created to store algorithm predictions
    predictions <- data.frame(source = S1,
                              target_catalogue = character(length(S1)),
                              target_predictive = character(length(S1)),
                              row.names = S1,
                              stringsAsFactors = FALSE)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # STEPS S1-S3
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Catalogue contribution of the algorithm
    # if taxa are known to interact in the catalogue, they are assumed to interact in the inferred network
    if (predict == 'full algorithm' | predict == 'catalogue') {
        # Evaluate for all species in S2
        for(i in S2) {
            # STEP S1
            # Identify resources of S2[i] in S0
            targetS2 <- unlist(strsplit(S0[i, 'target'], " \\|\\ "))

            # STEPS S2-S3
            # Add interactions between species S2[i] and resources found in S0 that are also in S1
            # Add link with S1 taxa that are considered as linked in S0
            if (length(targetS2)) {
                predictions[i, 'target_catalogue'] <- paste(targetS2[targetS2 %in% S1], collapse = ' | ')
            }
        }
    }

    # Predictive contribution of the algorithm, KNN algorithm to infer interactions
    if (predict == 'full algorithm' | predict == 'predictive') {
        for(i in S2) {

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # STEPS S4-S7
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # Empty matrix for resource candidate list for S2[i] consumer
            candidates <- matrix(nrow = 0, ncol = 2, dimnames = list(c(), c('target','weight')))

            # Identify resources of S2[i] that are in S0, but not in S1
            targetS2 <- unlist(strsplit(S0[i, 'target'], " \\|\\ ")) %>%
                        .[!. %in% S1]

            # Extract K most similar resources to targetS2 in S1
            if (length(targetS2)) {
                for(j in targetS2) {
                    candidates <- KNN(taxa = j, matSim = targetSim[, S1], K = K, minSim = minSim) %>%
                                  candLink(similar = ., candidates = candidates)
                }
            }

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # STEP S8
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # Identify K most similar source (consumers) to i in S0
            simSource <- KNN(taxa = i, matSim = sourceSim, K = K, minSim = minSim)

            if(length(simSource)) {
                for(j in names(simSource)) {
                    # List of resources for source j
                    target <- unlist(strsplit(S0[j, 'target'], " \\|\\ "))

                    if(length(target)) {
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                        # STEPS S9-S12
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                        # If candidate target are in S1, add to candidate list with weight = 1
                        for(k in target[target %in% S1]) {
                            similar <- 1
                            names(similar) <- k
                            candidates <- candLink(similar = similar, candidates = candidates)
                        }

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                        # STEPS S13-S16
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                        # If candidate target are not in S1, identify K most similar resources in S1
                        for(k in target[!target %in% S1]) {
                            candidates <- KNN(taxa = k, matSim = targetSim[,S1], K = K, minSim = minSim) %>%
                                          candLink(similar = ., candidates = candidates)
                        }#k
                    }#if
                }#j
            }#if

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # STEP S17
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            predictions[i, 'target_predictive'] <- candidates %>%
                                                   .[.[, 'weight'] >= minWt, 'target'] %>%
                                                   paste(., collapse = ' | ')
        } #i
    }#if
    return(predictions)
}#iEat
