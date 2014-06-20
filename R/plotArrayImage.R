#' Plot microarray images
#'
#' Plot a color image of the intensities on a slide. These plots are helpful to
#' diagnose spatial abnormalities.
#'
#'
#' @param peptideSet A \code{peptideSet} object. The object must contain all the
#' original probes. See details below.
#' @param array.index A vector subsetting \code{exprs(peptideSet)}, indicating
#' which slides to plot
#' @param smooth A \code{logical}, a 2D spatial smoother is applied to residuals,
#' the fitted residuals are plotted.
#' @param low A \code{character} string. The color of the lowest slide intensity.
#' passed to \code{scale_fill_gradient2}.
#' the fitted residuals are plotted.
#' @param high A \code{character} string. The color of the highest slide intensity.
#' passed to \code{scale_fill_gradient2}.
#' @param ask A \code{logical}. If TRUE, the user is asked before each plot. See
#' \code{par(ask=.)}.
#'
#' @details
#' The most coherent results are achieved when the \code{peptideSet} object is
#' read with \code{makePeptideSet} with empty.control.list = NULL and rm.control.list
#' = NULL
#'
#' @author Gregory Imholte
#'
#' @aliases plotArrayImage
#' @aliases plotArrayResiduals
#'
#' @importFrom ggplot2 ggplot ggtitle theme geom_tile aes_string element_blank
#' geom_segment scale_fill_gradient2
#' @importFrom fields smooth.2d
#' @importFrom plyr ddply
#' @export
#' @rdname plotArray
#' @example examples/pipeline.R
plotArrayImage <- function(peptideSet, array.index = NULL,
  low = "white", high = "steelblue",
  ask = dev.interactive(orNone = TRUE) & 1 < length(array.index)){
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
      geom_tile(data = d.tmp, aes_string(x = "x", y = "y", fill = "Intensity"),
                colour = "white") +
      scale_fill_gradient2(low = low, high = high) +
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank()) +
      geom_segment(data = plot.grid,
                   aes_string(y = "y", x = "x", xend = "xend", yend = "yend")) +
      ggtitle(paste("Sample Name:", sampleNames(peptideSet)[array.index[i]]))
    print(p)
    dev.flush()
  }
}

#' @rdname plotArray
#' @export
plotArrayResiduals <- function(peptideSet, array.index = NULL, smooth = FALSE,
  low = "blue", high = "red",
  ask = dev.interactive(orNone = TRUE) & 1 < length(array.index))
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
      geom_tile(data = d.tmp, aes_string(x = "x", y = "y", fill = "Intensity"),
                colour = "white") +
      scale_fill_gradient2(low = low, high = high) +
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank()) +
      geom_segment(data = plot.grid,
                   aes_string(y = "y", x = "x", xend = "xend", yend = "yend")) +
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

#   block <- as.numeric(peptideSet@featureRange@values[[1]]@listData$block)
#   column <- as.numeric(peptideSet@featureRange@values[[1]]@listData$column)
#   row <- as.numeric(peptideSet@featureRange@values[[1]]@listData$row)
  block <- as.numeric(values(peptideSet)$block)
  column <- as.numeric(values(peptideSet)$column)
  row <- as.numeric(values(peptideSet)$row)

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

  x <- (0:gr) * sr + .5
  y <- (0:gc) * sc + .5
  A <- cbind(x, x, y[1], y[length(y)])
  B <- cbind(x[1], x[length(x)], y, y)
  D <- data.frame(rbind(A,B))
  names(D) <- c("x", "xend", "y", "yend")
  return(D)
}
