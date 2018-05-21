#' Create a synthetic dataset to test the clustering algorithm using the GLD
#'
#' This function create a synthetic dataset using the approach proposed in B. Jiang, J. Pei, Y. Tao and X. Lin, "Clustering Uncertain Data Based on Probability Distribution Similarity," in IEEE Transactions on Knowledge and Data Engineering, vol. 25, no. 4, pp. 751-763, April 2013. doi: 10.1109/TKDE.2011.221
#' @author Noel Moreno Lemus
#' @param k number of clusters to test.
#' @param d dimension (e.g. if we have temperatute and preasure d = 2)
#' @param n the number of objects of the data
#' @param s the sample size
#' @examples
#' dataset <- createSyntheticData(3, 2, 1000, 500)
#'
#' @references
#' B. Jiang, J. Pei, Y. Tao and X. Lin, "Clustering Uncertain Data Based on Probability Distribution Similarity," in IEEE Transactions on Knowledge and Data Engineering, vol. 25, no. 4, pp. 751-763, April 2013.
#' doi: 10.1109/TKDE.2011.221
#' @export
#'

createSyntheticData <- function (k = 11, d = 1, n = 1000, s = 1000){
  no_gaussian = (k - 1)/2
  count = 1
  dataset = array(0, dim = c(n, 1, s))

  # Gaussian Distribution
  for (i in 1:no_gaussian) {
    for (j in  1:90) {
      print(count)
      dataset[count, 1, ] = rnorm(s, 0, 0.05*i)
      count = count + 1
    }
  }

  # Exponential Distribution
  for (i in 1:no_gaussian) {
    for (j in  1:90) {
      print(count)
      dataset[count, 1, ] = rexp(s, i)
      count = count + 1
    }
  }

  # Gamma Distribution
 # for (i in 1:no_gaussian) {
  #  for (j in  1:90) {
   #   print(count)
    #  dataset[count, 1, ] = rgamma(s, i)
     # count = count + 1
  #  }
#  }

  print(c("n-count", n-count))
  # Uniform Distribution
  for (i in 1:(n-count+1)) {
      print(count)
      dataset[count, 1, ] = runif(s)
      count = count + 1
  }
  dataset
}

gldFitTest <- function(dataset) {
  a = dim(dataset)
  library(doSNOW)
  # Initialize a cluster
  cl <- makeCluster(4)
  # registerDoParallel(cl)
  registerDoSNOW(cl)
  # Initialize de variables
  lambdas = array(0, dim = c(a[1], 1, 4))
  statistics = array(0, dim = c(a[1], 1, 4))

  stime <- system.time({
    theRes <- foreach(
      i = 1:a[1],
      .combine = rbind,
      .packages = c("GLDEX")
    ) %dopar% {
      lambdas = NULL
      statistics = NULL
      #for (j in 1:501) {
      mydata = dataset[i, 1,]
      stackedLambdaMM <- fun.RMFMKL.mm(mydata)
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
      #fun.moments.r(mydata, normalise = "Y")
      list(lambdas, statistics)
    }

    #}
  })
  print(c("Time -> ", stime))
  stopCluster(cl)
  theRes
}

gldGoFTest <- function(dataset, lambdas){
  library(doSNOW)
  # Initialize a cluster
  cl <- makeCluster(4)
  # registerDoParallel(cl)
  registerDoSNOW(cl)
  # Initialize de variables
  ksTestMatrix = array(0, dim = c(1000, 1, 1))
  pvaluesMatrix = array(0, dim = c(1000, 1, 1))

  stime <- system.time({
    gofRes <- foreach(i = 1:1000,
                      .combine = rbind,
                      .packages = c("GLDEX")) %dopar% {
                        ksTestMatrix = NULL
                        pvaluesMatrix = NULL

                        mydata = dataset[i, 1, ]
                        # Komogorov-Smirnoff (KS) resample test
                        ksTest = fun.diag.ks.g(result = lambdas[i, 1, ], data = dataset[i, 1, ], no.test = 1000,
                                               param = "fmkl")
                        ksTestMatrix = rbind(ksTestMatrix, ksTest)

                        # Kolmogorov-Smirnov test
                        ksTest1 = ks.gof(mydata, "pgl", lambdas[i, 1, ], param="fmkl")
                        pvaluesMatrix = rbind(pvaluesMatrix, ksTest1$p.value)

                        list(ksTestMatrix, pvaluesMatrix)
                      }

  })
  print(c("Time -> ", stime))
  stopCluster(cl)
  gofRes
}

extractGoF <- function(gofRes){
  D = array(0, dim = c(1000))
  p_values = array(0, dim = c(1000))
  for (i in 1:1000) {
    a1 = str_c("gofRes[",i,",1]$result.",i)
    a2 = str_c("gofRes[",i,",2]$result.",i)
    temp1 = eval(parse(text = a1))
    temp2 = eval(parse(text = a2))
    D[i] = temp1[1,]
    p_values[i] = temp2[1,]
  }

  result = list("D" = D, "p_values" = p_values)
  return(result)
}

extractLambdas <- function(theRes, dimDataSet){
  lambdas = array(0, dim = c(dimDataSet, 1, 4))
  statistics = array(0, dim = c(dimDataSet, 1, 4))

  for (i in 1:dimDataSet) {
    a1 = str_c("theRes[",i,",1]$result.",i)
    a2 = str_c("theRes[",i,",2]$result.",i)
    temp1 = eval(parse(text = a1))
    temp2 = eval(parse(text = a2))
    lambdas[i,1,] = temp1[1,]
  }
  lambdas
}

extractL1234 <- function(lambdas){

}

gldClustering1D <- function(lambdas, no_clusters, l234 = TRUE){
  lambda1 = lambdas[,,1]
  lambda2 = lambdas[,,2]
  lambda3 = lambdas[,,3]
  lambda4 = lambdas[,,4]

  dimension = length(lambda1)

  dimTotal = dimension


  if (l234){
    x = array(0, dim = c(dimTotal, 3))
  } else {
    x = array(0, dim = c(dimTotal, 2))
  }


  count = 0;
  for(i in 1:dimension){

    count = count + 1;
    #print(count)
    #x[count,1] = lambda1[i, j];
    if (l234){
      x[count,1] = lambda2[i];
      x[count,2] = lambda3[i];
      x[count,3] = lambda4[i];
    } else {
      x[count,1] = lambda3[i];
      x[count,2] = lambda4[i];
    }


  }

  cl <- kmeans(x, no_clusters)

  #cl <- dbscan(x, eps = 0.5)

  results = list("clusters" = cl, "x" = x)
  return(results)
}

plot.clusters.l3l4 <- function(cluster_number, color) {

  dim_clusters = length(clusters)

  for(i in 1:dim_clusters){
    #for(j in 1:dim_clusters[2]){
    if (clusters[i]==cluster_number){
      points(lambda3[i], lambda4[i], col = color, pch = 16)
    }
    #}
  }
}

plot.gld.by.cluster <- function(cluster_number, n){
  lista = list()

  # y = (clusters==cluster_number);
  count = 1;
  i = 1;

  dim_clusters = length(clusters)

  for(i in 1:dim_clusters){
    #for(j in 1:dim_clusters[2]){
      if (clusters[i]==cluster_number && count < n){
        lista[[length(lista)+1]] = rgl(2000, 0, 2, lambda3[i], lambda4[i])
        count = count + 1;
      }
    #}
  }


  #while(count < n) {
  #  if(y[i]){
  #    print(i)
  #    lista[[length(lista)+1]] = rgl(2000, lambda1[i,], lambda2[i,], lambda3[i,], lambda4[i,]);
  #    count = count + 1;
  #  }
  #  i = i + 1;
  # }

  plot.multi.dens(lista)
  #plot.multi.gld(lista)

}

fit.multi.dist <- function(data, dists){
  n = length(dists)

  stime <- system.time({
    for (i in 1:n) {
      print(dists[i])
      fit <- fitdist(data, dists[i])
    }
  })
  stime
}


