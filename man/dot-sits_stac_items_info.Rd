% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sits_stac.R
\name{.sits_stac_items_info}
\alias{.sits_stac_items_info}
\title{Format assets}
\usage{
.sits_stac_items_info(items, bands)
}
\arguments{
\item{items}{a \code{STACItemCollection} object returned by rstac package.}

\item{bands}{a \code{character} with the bands names to be filtered.}
}
\value{
a \code{tibble} with date, band and path information, arranged
 by the date.
}
\description{
Format assets
}
\keyword{internal}
