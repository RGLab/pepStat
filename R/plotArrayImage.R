plotArrayImage <- function(peptideSet, array.index = NULL,
  low = "white", high = "steelblue",
  ask = dev.interactive() & 1 < length(array.index))
{
  if (is.null(array.index))
    stop("Must specify indicies of arrays to be plotted")
  
  layout <- preproc(peptideSet)$printer
  if (is.null(layout))
    stop("Cannot find printer layout in preproc info of peptideSet")
  
  gr <- layout$ngrid.r
  gc <- layout$ngrid.c
  sr <- layout$nspot.r
  sc <- layout$nspot.c
  
  plot.coord <- getPlotCoords(peptideSet)
  plot.grid <- getPlotGrid(peptideSet)
  y <- exprs(peptideSet)[, array.index, drop = FALSE]
  n.plots <- ncol(y)
  
  if (nrow(plot.coord) != nrow(y))
    stop("Slide layout does not match number of features in peptideSet")
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  for(i in 1:n.plots) {
    dev.hold()
    d.tmp <- data.frame(Intensity = y[,i], plot.coord)
    p <- ggplot() +
      geom_tile(data = d.tmp, aes(x = x, y = y, fill = Intensity), 
                colour = "white") +
      scale_fill_gradient(low = low, high = high) + 
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank()) +
      geom_segment(data = plot.grid, 
                   aes(y = y, x = x, xend = xend, yend = yend)) +
      ggtitle(paste("Sample Name:", sampleNames(peptideSet)[array.index[i]]))
    print(p)
    dev.flush()
  }
}

plotArrayResiduals <- function(peptideSet, array.index = NULL, smooth = FALSE,
  low = "blue", high = "red",
  ask = dev.interactive() & 1 < length(array.index))
{
  if (is.null(array.index))
    stop("Must specify indicies of arrays to be plotted")
  
  layout <- preproc(peptideSet)$printer
  if (is.null(layout))
    stop("Cannot find printer layout in preproc info of peptideSet")
  
  gr <- layout$ngrid.r
  gc <- layout$ngrid.c
  sr <- layout$nspot.r
  sc <- layout$nspot.c
  
  plot.coord <- getPlotCoords(peptideSet)
  plot.grid <- getPlotGrid(peptideSet)
  
  res <- makeResiduals(peptideSet, array.index, plot.coord)
  prefix = "Residuals"

  if (smooth){
    res <- smoothResiduals(res, layout)
    prefix <- paste("Smoothed", prefix)
  }
  
  ind <- which(names(res) %in% c("x", "y", "id"))
  res.mat <- as.matrix(res[,-ind])
  n.plots <- ncol(res.mat)
  plot.coord <- res[,c("x", "y")]
  
  if (nrow(plot.coord) != nrow(res.mat))
    stop("Slide layout does not match number of features in peptideSet")
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  for(i in 1:n.plots) {
    dev.hold()
    d.tmp <- data.frame(Intensity = res.mat[,i], plot.coord)
    p <- ggplot() +
      geom_tile(data = d.tmp, aes(x = x, y = y, fill = Intensity), 
                colour = "white") +
      scale_fill_gradient(low = low, high = high) + 
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank()) +
      geom_segment(data = plot.grid, 
                   aes(y = y, x = x, xend = xend, yend = yend)) +
      ggtitle(paste(prefix, "for Sample Name", 
                    sampleNames(peptideSet)[array.index[i]]))
    print(p)
    dev.flush()
  }
}

smoothResiduals <- function(resid.df, layout){
  gr <- layout$ngrid.r
  gc <- layout$ngrid.c
  sr <- layout$nspot.r
  sc <- layout$nspot.c
  nrow = gr * sr
  ncol = gc * sc
  
  ind <- which(names(resid.df) %in% c("x", "y", "id"))
  res.mat <- as.matrix(resid.df[, -ind])
  plot.coord <- resid.df[,c("x", "y")]
  sm <- apply(res.mat, 2, function(y){
    z <- smooth.2d(y, x = plot.coord, nrow = gr * sr, ncol = gc * sc)
    as.vector(z$z)
  })
  data.frame(sm, x = rep(1:nrow, ncol), y = rep(1:ncol, each = nrow))
}

makeResiduals <- function(peptideSet, array.index, plot.coord){
  tmp <- data.frame(int = exprs(peptideSet)[, array.index],
                    id = factor(featureID(peptideSet)), plot.coord)
  
  out <- ddply(.data = tmp, .variables = "id", .fun = function(v){
    ind <- which(names(v) %in% c("x", "y", "id"))
    tmp.int <- as.matrix(v[, -ind])
    m <- colMeans(tmp.int)
    res <- tmp.int - rep(m, each = nrow(v))
    data.frame(res = res, x = v$x, y = v$y)
  })
  
  out
}

getPlotCoords <- function(peptideSet)
{
  layout <- preproc(peptideSet)$printer
  if (is.null(layout))
    stop("Cannot find printer layout in preproc info of peptideSet")

  gr <- layout$ngrid.r
  gc <- layout$ngrid.c
  sr <- layout$nspot.r
  sc <- layout$nspot.c
  
  block <- rep(1:(gr * gc), each = sr * sc)
  column <- rep(1:sc, sr * gc * gr)
  row <- rep(rep(1:sr, each = sc), gc * gr)  

  y <- column + ((block - 1) %% gc) * sc
  x <- row + floor((block - 1) / gc) * sr
  grid.coord <- data.frame(x, y)  
  grid.coord
}

getPlotGrid <- function(peptideSet)
{
  layout <- preproc(peptideSet)$printer
  if (is.null(layout))
    stop("Cannot find printer layout in preproc info of peptideSet")
  
  gr <- layout$ngrid.r
  gc <- layout$ngrid.c
  sr <- layout$nspot.r
  sc <- layout$nspot.c
    
  x <- (0:gr) * sc + .5
  y <- (0:gc) * sr + .5
  A <- cbind(x, x, y[1], y[length(y)])
  B <- cbind(x[1], x[length(x)], y, y)
  D <- data.frame(rbind(A,B))
  names(D) <- c("x", "xend", "y", "yend")
  return(D)
}