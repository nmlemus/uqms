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
#' @author Noel Moreno Lemus
#' @param L1 Lambda values of the first GLD.
#' @param L2 Lambda values of the second GLD.
#' @examples
#' L1 = c(0, 2, 0.25, 1.5)
#' L2 = c(0, 2, 0.3, 1.75)
#' D <- gldComparison(L1, L2)
#'
#' @export
#'

gldComparison <- function(L1, L2, param = "fmkl", no.test = 1000, len = floor(0.9*no.test), alpha = 0.05) {
    sample1 <- rgl(len * no.test, L1[1], L1[2], L1[3], L1[4], param)
    sample2 <- rgl(len * no.test, L2[1], L2[2], L2[3], L2[4], param)

    # test and sample.fitted to use the same names used in fun.diag.ks.g
    test <- split(sample1, 1:no.test)
    sample.fitted <- split(sample2, 1:no.test)

    result.o <- sum(sapply(1:no.test, function(i, test, sample.fitted)
    ks.gof(test[[i]], sample.fitted[[i]])$p.value, test, sample.fitted) > alpha)
    return(result.o)
}

distGLDComparison <- function(x, centers) {
  if (ncol(x) != ncol(centers))
    stop(sQuote("x")," and ", sQuote("centers")," must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))

  alpha = 0.05
  no.test = 1000
  len = floor(0.9*no.test)
  sample1 <- rgl(len * no.test, 0, 2, x[1], x[2], "fmkl")
  test <- split(sample1, 1:no.test)

  for (k in 1:nrow(centers)) {
    sample2 <- rgl(len * no.test, 0, 2, centers[k, 1], centers[k, 2], "fmkl")
    sample.fitted <- split(sample2, 1:no.test)
    result.o <- sum(sapply(1:no.test, function(i, test, sample.fitted)
    ks.gof(test[[i]], sample.fitted[[i]])$p.value, test, sample.fitted) > alpha)

    z[, k] <- (1000 - result.o)
  }
  z
}

#' Evaluate if the centroid of a cluster is a good representative of the other members of this cluster
#'
#' This function return a list of D distances from the ks-test that compare the centroid of a cluster wiht n members of this cluster.
#' @author Noel Moreno Lemus
#' @param cluster_number The cluster ID we are interesting in.
#' @param n Number of elements of the cluster to analize.
#' @param centroid The lambda values of the centroid of the cluster.
#' @examples
#'
#' Ds <- gldClusterComparison(3, 60, c(0, 2, 1.5, 1.3))
#'
#' @export
#'

gldClusterComparison <- function(cluster_number, n, centroid) {
  lista = list()
  count = 1;
  i = 1;
  dim_clusters = dim(clusters)
  for(i in 1:dim_clusters[1]){
    for(j in 1:dim_clusters[2]){
      if (clusters[i,j]==cluster_number && count <= n){
        D = gldComparison(centroid, c(0, lambda2[i,j], lambda3[i,j], lambda4[i,j]))
        lista[[length(lista)+1]] = D;
        count = count + 1;
      }
    }
  }
  as.numeric(lista)
}


#' Clustering of the GLDs in function of its l2, l3 and l4 values
#'
#' TODO.
#' @author Noel Moreno Lemus
#' @param lambdas A matrix of n x m x 4 of all the lambda values.
#' @param no_clusters Number of clusters.
#'
#' @examples
#'
#' gldClustering(lambdas)
#'
#' @export
#'
gldClustering <- function(lambdas, no_clusters){
  lambda1 = lambdas[,,1]
  lambda2 = lambdas[,,2]
  lambda3 = lambdas[,,3]
  lambda4 = lambdas[,,4]

  dimension = dim(lambda1)

  dimTotal = dimension[1]*dimension[2]

  lambdas34 = array(0, dim = c(dimTotal, 6))
  x = array(0, dim = c(dimTotal, 3))

  count = 0;
  for(i in 1:dimension[1]){
    for(j in 1:dimension[2]){
      count = count + 1;
      print(count)
      lambdas34[count,1] = i;
      lambdas34[count,2] = j;
      lambdas34[count,3] = lambda1[i, j];
      lambdas34[count,4] = lambda2[i, j];
      lambdas34[count,5] = lambda3[i, j];
      lambdas34[count,6] = lambda4[i, j];
      #x[count,1] = lambda1[i, j];
      x[count,1] = lambda2[i, j];
      x[count,2] = lambda3[i, j];
      x[count,3] = lambda4[i, j];
    }
  }

  cl <- kmeans(x, no_clusters)

  #cl <- dbscan(x, eps = 0.5)
  clusters = array(0, dim = c(dimension[1], dimension[2]))
  for(k in 1:dimTotal){
    clusters[lambdas34[k,1], lambdas34[k,2]] = cl$cluster[k]
  }

  image_display(clusters)

  results = list("clusters" = clusters, "x" = x)
  return(results)
}


#' Plot the clusters in l3-l4 space.
#'
#' TODO.
#' @author Noel Moreno Lemus
#' @param clusters An n x m matrix wiht the clusters by positions.
#' @param x Array used to create the clusters in function @method gldClustering.
#'
#' @examples
#'
#' gldClustersL3L4(clusters, x)
#'
#' @export
#'
gldClustersL3L4 <- function (clusters, x) {
  library(latex2exp)
  plot(x[,2], x[,3], type = "n", xlab = TeX('$\\lambda_{3}'), ylab = TeX('$\\lambda_{4}'), main = TeX('Clusters in $\\lambda_{3}-\\lambda_{4}'))  # setting up coord. system

  no_clusters = max(clusters)

  #legend_list = list()

  for (i in 1:no_clusters) {
    plot.clusters.l3l4(i, i)
    #legend_list[[length(legend_list)+1]] = paste("Cluster ", i)
  }
  #legend("right", 95, legend_list, col=c(1:no_clusters), lty=1, cex=0.8)
}
