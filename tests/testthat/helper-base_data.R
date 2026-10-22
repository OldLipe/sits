# Helpers to build base data (e.g., digital elevation models) fixtures.
#
# Base data is normally produced by sits_add_base_cube() + sits_get_data(),
# which require a DEM collection from a remote provider. These helpers build
# equivalent fixtures from the local MOD13Q1 example files, so the base data
# support can be tested without network access.

#' Add synthetic base data to a set of training samples
#'
#' @param samples Training samples (class "sits").
#' @param sd      Standard deviation of the elevation of each label.
#' @return Samples of class "sits_base" with an ELEVATION base band.
test_samples_base <- function(samples, sd = 80) {
    labels <- sits_labels(samples)
    # one mean elevation per label, so the base encoder has a signal
    elevation <- seq(500, 1200, length.out = length(labels))
    names(elevation) <- labels
    samples[["base_data"]] <- purrr::map2(
        samples[["time_series"]], samples[["label"]],
        function(ts, label) {
            tibble::tibble(
                Index = as.Date(ts[["Index"]][[1L]]),
                ELEVATION = stats::rnorm(1L, elevation[[label]], sd)
            )
        }
    )
    class(samples) <- unique(c("sits_base", "sits", class(samples)))
    samples
}

#' Build a local data cube with base data
#'
#' Writes a synthetic elevation raster on the grid of the MOD13Q1 example
#' files and attaches it to the example cube as base data.
#'
#' @param output_dir Directory to store the elevation file.
#' @return A data cube of class "base_raster_cube".
test_base_cube <- function(output_dir) {
    cube <- sits_cube(
        source = "BDC",
        collection = "MOD13Q1-6.1",
        data_dir = system.file("extdata/raster/mod13q1", package = "sits"),
        progress = FALSE
    )
    # elevation as a smooth spatial gradient over the cube grid
    template <- .raster_open_rast(.tile_path(cube))
    elevation <- .raster_rast(template, nlayers = 1L)
    nrows <- .raster_nrows(template)
    ncols <- .raster_ncols(template)
    values <- outer(seq_len(nrows), seq_len(ncols), function(row, col) {
        500 + 2 * row + 1.5 * col
    })
    elevation <- .raster_set_values(elevation, as.integer(round(values)))
    .raster_write_rast(
        rast = elevation,
        file = file.path(
            output_dir, "TERRA_MODIS_012010_ELEVATION_2013-09-14.tif"
        ),
        data_type = "INT2S",
        overwrite = TRUE,
        missing_value = -32768L
    )
    dem_cube <- sits_cube(
        source = "MPC",
        collection = "COP-DEM-GLO-30",
        data_dir = output_dir,
        bands = "ELEVATION",
        progress = FALSE
    )
    sits_add_base_cube(cube, dem_cube)
}
