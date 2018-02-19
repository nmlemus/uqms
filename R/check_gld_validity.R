#' Check if one GLD is valid
#'
#' This function check if one GLD is valid.
#' @param lambdas lambda values of a single GLD function.
#' @keywords lambdas
#' @examples
#' checkGLDValid1(lambdas)
#' @export
#'
checkGLDValid1 <- function(lambdas){
  isVal = gl.check.lambda.alt(lambdas[1],
                              lambdas[2],
                              lambdas[3],
                              lambdas[4],
                              param = "fkml")
  return(isVal)
}

#' Check if many GLDs are valid
#'
#' This function check if many GLDs are valid.
#' @param lambdas lambda values of multiple GLD functions.
#' @keywords lambdas
#' @examples
#' checkGLDValidN(lambdas)
#' @export
#'
checkGLDValidN <- function(lambdas){
  n = length(lambdas[,1])
  isValid = array(0, dim = c(n, 1))

  for(i in 1:n){
    isVal = gl.check.lambda.alt(lambdas[i,1],
                                lambdas[i,2],
                                lambdas[i,3],
                                lambdas[i,4],
                                param = "fkml")
    if(isVal){
      isValid[i] = 1;
    }
  }

  result <- list("isValid" = isValid)
  return(result)
}
