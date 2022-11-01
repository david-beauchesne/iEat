#' @name iEat_to_foodWeb
#'
#' @title transform iEat prediction result to matrix format
#'
#' @param iEatResult output from iEat function
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

  for (i in 1:nrow(iEatResult)) {
    resources <- unique(c(unlist(str_split(iEatResult[i, "target_catalogue"], " \\| ")), unlist(str_split(iEatResult[i, "target_predictive"], " \\| "))))
    foodWeb[which(rownames(foodWeb) %in% resources), i] <- 1
  }

  return(foodWeb)
}
