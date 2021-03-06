% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sits_smooth_aux_functions.R
\name{.sits_map_layer_cluster}
\alias{.sits_map_layer_cluster}
\title{Parallel processing of classified images}
\usage{
.sits_map_layer_cluster(
  cube,
  cube_out,
  overlapping_y_size = 0,
  func,
  func_args = NULL,
  multicores = 1,
  memsize = 1,
  ...
)
}
\arguments{
\item{cube}{Probability data cube}

\item{cube_out}{Output probability data cube}

\item{func}{a function that receives RasterBrick and
returns any Raster*.}

\item{multicores}{Number of process to run the Bayesian smoothing in
snow subprocess.}

\item{memsize}{Maximum overall memory (in GB) to run the Bayesian
smoothing.}

\item{...}{optional arguments to merge final raster
(see \link[raster]{writeRaster} function)}

\item{overlapping_rows}{number of overlapping rows of each chunk.}

\item{args}{additional arguments to pass to \code{fun} function.}
}
\value{
RasterBrick object
}
\description{
Process chunks of raster bricks individually in parallel.
}
\keyword{internal}
