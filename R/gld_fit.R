#' Create a multidimentional array with simulations form a set of .csv.
#'
#' This function read all the csv that are stored in a directory and create an NxMxS array, where N and M are sparial dimensions
#' and S is the number of simulations on each spatial point.
#' @author Noel Moreno Lemus
#' @param directory the directory where the .csv files are stored.
#' @param pattern a pattern to read the .csv (e.g. pattern="^[h]").
#' @param dimension the dimensions to be used, NxMxS.
#' @examples
#' directory = "~/PhD/thesis_phd/python_codes/datasets/20"
#' pattern="^[h]"
#' dimension = c(60, 60, 500)
#' mArray <- createDataset(directory, pattern, dimension)
#'
#' @export
#'

createDataset <- function(directory, pattern, dimension){
  # Set the working directory where the .csv are stored.
  setwd(directory)
  # Patter to read part of the .csv (e.g. pattern="^[h]")
  temp = list.files(pattern=pattern)

  # Read all the files
  for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

  # Create and array of 60 x 60 x 500 to store all the simulations at timestep = 20
  mArray = array(0, dim = c(dimension[1], dimension[2], dimension[3]));

  # Read all the .csv, convert into matrix and remove the first column.
  for (i in 1:length(temp)){
    print(i)
    h1 = as.matrix(read.csv(temp[i]))
    mArray[,,i] = h1
  }
  return(mArray)
}

#' Fit the GLD to a multidimensional array in parallel.
#'
#' This function compute the GLD that best fit each distribution on each spatial location.
#' @author Noel Moreno Lemus
#' @param nodes number of nodes to be used in the computation.
#' @param dimension the dimensions to be used, NxMxS.
#' @param mArray a multidimentional array with dimensions @param dimension
#' @examples
#' mArray <- createDataset(directory, pattern, dimension)
#' result <- gldFitParallel(4, dimension = c(60, 60), mArray = mArray)
#'
#' @export
#'

gldFitParallel <- function(nodes, dimension, mArray){
  library(doSNOW)
  # Initialize a cluster
  cl <- makeCluster(nodes)
  #registerDoParallel(cl)
  registerDoSNOW(cl)
  # Initialize de variables
  lambdas = array(0, dim = c(dimension[1], dimension[2], 4))
  statistics = array(0, dim = c(dimension[1], dimension[2], 4))

  myarray = mArray

  stime <- system.time({
    theRes <- foreach(i = 1:dimension[1],
                      .combine = rbind,
                      .packages = c("GLDEX")) %dopar% {
                        lambdas = NULL
                        statistics = NULL
                        for (j in 1:dimension[2]) {
                          mydata = myarray[i, j, ]
                          stackedLambdaMM <- fun.RMFMKL.mm(mydata[!is.na(mydata)])
                          lambdas = rbind(lambdas, stackedLambdaMM)
                          # Moments of fitted distribution:
                          st = fun.theo.mv.gld(
                            stackedLambdaMM[1],
                            stackedLambdaMM[2],
                            stackedLambdaMM[3],
                            stackedLambdaMM[4],
                            param = "fmkl",
                            normalise = "Y"
                          )
                          statistics = rbind(statistics, st)
                          # Moments of the original data:
                          fun.moments.r(mydata, normalise = "Y")
                        }
                        list(lambdas, statistics)
                      }
  })
  stime
  stopCluster(cl)
  result = list("theRes" = theRes, "stime" = stime)
}

extractLambdas <- function(theRes, dimension){
  library("stringr")

  # Initialize de variables
  lambdas = array(0, dim = c(dimension[1], dimension[2], 4))
  statistics = array(0, dim = c(dimension[1], dimension[2], 4))

  for(i in 1:dimension[1]){
    a1 = str_c("theRes[",i,",1]$result.",i)
    a2 = str_c("theRes[",i,",2]$result.",i)
    temp1 = eval(parse(text = a1))
    temp2 = eval(parse(text = a2))
    for(j in 1:dimension[2]){
      print(c(i, "->", j))
      lambdas[i,j,] = temp1[j,]
      statistics[i,j,] = temp2[j,]
    }
  }
  result = list("lambdas" = lambdas, "statistics" = statistics)
  return(result)
}

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
