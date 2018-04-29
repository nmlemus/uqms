#' Compute the lambda values of the GLD using the method of the Lmoments
#'
#' This function allows you to show image with scale.
#' @param data The raw data to fit.
#' @keywords lmoments
#' @examples
#' gldfitLmoments(data)
#' @export

gldfitLmoments <- function(data){
  # compute the lmoments
  lmr <- lmoms(data)
  params_gld <- pargld(lmr) # fit the GLD
  lambdas = params_gld$para;

  result <- list("lambdas" = lambdas)
  return(result)
}

#' Compute the lambda values of the GLD using the method of Moments
#'
#' This function allows you to show image with scale.
#' @param data The raw data to fit.
#' @keywords lmoments
#' @examples
#' gldfitMoments(data)
#' @export

gldfitMoments <- function(data){
  lambdas <- fun.RMFMKL.mm(data)
  result <- list("lambdas" = lambdas)
  return(result)
}

gldfitMaxlikelihood <- function(data){
  lambdas <- fun.RMFMKL.ml(data)
  result <- list("lambdas" = lambdas)
  return(result)
}
