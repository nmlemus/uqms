#' Color image dislay with scale
#'
#' This function allows you to show image with scale.
#' @param date The image matrix you like to display.
#' @keywords image
#' @export
#' @examples
#' image_display()
#' @export

image_display <- function(data){
  library('gplots')
  library('fields')
  library('RColorBrewer')

  # Original 11 colors
  cols = rev(brewer.pal(11,'Spectral'))

  # Repeat the 11th color an extra 5 times
  cols = c(cols,rep(cols[11],5))

  rf <- colorRampPalette(cols)   # make colors
  r <- rf(64)

  z =  data
  x <- (1:nrow(z))
  y <- (1:ncol(z))

  #first reverse, then transpose, it's the same as rotate 90 degrees
  rotate_clockwise         <- function(x) { t(     apply(x, 2, rev))}
  foo = rotate_clockwise(z)

  image.plot(y, x, foo, col = r)
}
