#' Plot the GLD based in the Lambda values
#'
#' This function allows you to show image with scale.
#' @param L Lambda values of the GLD.
#' @keywords lmoments
#' @examples
#' L = c(0, 2, 0.25, 1.5)
#' gldPlot(L)
#'
#' @export
#'

gldPlot <- function(L){
  # Create a Y vector between 0 and 1.
  Y <- c(0.0001, 0.001*(1:9), 0.01*(1:99), .99+0.001*(1:9), .9999)
  # Compute x = Q(Y)
  x <- L[1] + (Y^L[3] - (1-Y)^L[4])/L[2]
  # Compute f(x) = PDF(x)
  y <- L[2]/(L[3]*Y^(L[3]-1) + L[4]*(1-Y)^(L[4]-1))
  plot(x, y, type="l")
}

#' Return (x, y) values of the GLD based in the Lambda values
#'
#' This function allows you to return the (x,y) values of the GLD, based in the Lambda values.
#' @param L Lambda values of the GLD.
#' @examples
#' L = c(0, 2, 0.25, 1.5)
#' gldToPlot(L)
#'
#' @export
#'

gldToPlot <- function(L){
  # Create a Y vector between 0 and 1.
  Y <- c(0.0001, 0.001*(1:9), 0.01*(1:99), .99+0.001*(1:9), .9999)
  # Compute x = Q(Y)
  x <- L[1] + (Y^L[3] - (1-Y)^L[4])/L[2]
  # Compute f(x) = PDF(x)
  y <- L[2]/(L[3]*Y^(L[3]-1) + L[4]*(1-Y)^(L[4]-1))

  result <- list("x" = x, "y" = y)
  return (result)
}

#' Compare if two GLDs belongs to the same distribution
#'
#' This function return D distance from the ks-test that compare if both GLDs belongs to the same distribution.
#' @param L1 Lambda values of the first GLD.
#' @param L2 Lambda values of the second GLD.
#' @examples
#' L1 = c(0, 2, 0.25, 1.5)
#' L2 = c(0, 2, 0.3, 1.75)
#' D <- gldComparison(L1, L2)
#'
#' @export
#'

gldComparison <- function(L1, L2, alternative = c("random", "quantiles")) {
  if (alternative == "random"){
    d1 <- rgl(2000, L1, param = "fmkl")
    d2 <- rgl(2000, L2, param = "fmkl")
    D12 <- fun.diag.ks.g(L1, d2, param = "fmkl")
    D21 <- fun.diag.ks.g(L2, d1, param = "fmkl")
  } else if (alternative == "quantiles"){
    d1 <- gldToPlot(L1)$x
    d2 <- gldToPlot(L2)$x
    D12 <- fun.diag.ks.g(L1, d2, param = "rs")
    D21 <- fun.diag.ks.g(L2, d1, param = "rs")
  }


  gof = 0
  if (D12 > 850 || D21 > 850){
    gof = 1
  }
  result <- list("D12" = D12, "D21" = D21, "gof" = gof)
  return (result)
}
