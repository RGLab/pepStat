library(grid)

textPlot <- function(text) {
  grid.newpage()
  grid.text(text)
}

parseTextField <- function(x) {
  if (identical(x, "")) return(NULL)
  return( scan( textConnection(x), what=character(), sep="," ) )
}
