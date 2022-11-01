#' @name iEat_to_foodWeb
#'
#' @title transform iEat prediction result to matrix format
#'
#' @param  Matrix, catalogue of empirical data used to infer predictions
#' 
#' @return
#' A dataframe with source taxa for which target predictions are made, target infered from catalogue data (empirical) and target infered from KNN algorithm
#'
#' @author
#' David Beauchesne
#'
#' @importFrom magrittr %>%
#'
#' @export

iEat_to_foodWeb <- function(iEatResult) {
    library(stringr)
    foodWeb <- matrix(nrow = nrow(iEatResult), ncol = nrow(iEatResult), data = 0, dimnames = list(rownames(iEatResult), rownames(iEatResult)))

    for(i in 1:nrow(iEatResult)) {
        resources <- unique(c(unlist(str_split(iEatResult[i, 'target_catalogue'], ' \\| ')), unlist(str_split(iEatResult[i, 'target_predictive'], ' \\| '))))
        foodWeb[which(rownames(foodWeb) %in% resources), i] <- 1
    }

    return(foodWeb)
}
