% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sits_cube.R
\name{sits_cube.raster_cube}
\alias{sits_cube.raster_cube}
\title{Backward compatibility}
\usage{
\method{sits_cube}{raster_cube}(type = "RASTER", name = NULL, ...)
}
\arguments{
\item{type}{Type of cube}

\item{name}{Name of the input data}

\item{...}{Other parameters to be passed for specific types}
}
\value{
A message
}
\description{
Using "RASTER" as a type for a data cube is deprecated in favour
             of "STACK" type.
}
