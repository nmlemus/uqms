plot.gld.by.cluster <- function(cluster_number, n){
  lista = list()

 # y = (clusters==cluster_number);
  count = 1;
  i = 1;

  dim_clusters = dim(clusters)
  
  for(i in 1:dim_clusters[1]){
    for(j in 1:dim_clusters[2]){
      if (clusters[i,j]==cluster_number && count < n){
        lista[[length(lista)+1]] = rgl(2000, 0, lambda2[i,j], lambda3[i,j], lambda4[i,j])
        count = count + 1;
      }
    }
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

plot.gld.by.cluster2 <- function(cluster_number, n){
  lista = list()

  # y = (clusters==cluster_number);
  count = 1;
  i = 1;

  dim_clusters = dim(clusters)
  
  for(i in 1:dim_clusters[1]){
    for(j in 1:dim_clusters[2]){
      if (clusters[i,j]==cluster_number && count < n){
        if (count == 1){
          resultado = gldToPlot(c(lambda1[i,j], lambda2[i,j], lambda3[i,j], lambda4[i,j]));
          plot(resultado$x, resultado$y, type = 'l')
        } else {
          resultado = gldToPlot(c(lambda1[i,j], lambda2[i,j], lambda3[i,j], lambda4[i,j]));
          lines(resultado$x, resultado$y, col = count)
        }
        count = count + 1;
      }
    }
  }


  #while(count < n) {
  #  if(y[i]){
  #    print(i)
  #    lista[[length(lista)+1]] = rgl(2000, lambda1[i,], lambda2[i,], lambda3[i,], lambda4[i,]);
  #    count = count + 1;
  #  }
  #  i = i + 1;
  # }

  # plot.multi.dens(lista)

}

plot.clusters.l3l4 <- function(cluster_number, color) {
  
  dim_clusters = dim(clusters)
  
  for(i in 1:dim_clusters[1]){
    for(j in 1:dim_clusters[2]){
      if (clusters[i,j]==cluster_number){
        points(lambda3[i,j], lambda4[i,j], col = color, pch = 16)
      }
    }
  }
}
