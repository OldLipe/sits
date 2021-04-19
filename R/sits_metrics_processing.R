#' @title Preprocess a set of values retrieved from a raster brick
#' @name  .sits_raster_metrics_data_preprocess
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param  cube             data cube being processed
#' @param  band_cube        band to be processed
#' @param  metric           metric to be process
#' @param  extent           extent to be read
#' @param  filter           smoothing filter to be applied.
#' @param  stats            normalization parameters.
#' @param  impute_fn        imputing function to be applied to replace NA
#' @param  multicores       number of cores to process the time series.
#' @param  .verbose         prints information about processing times
#' @return Matrix with pre-processed values.
.sits_raster_metrics_data_preprocess <- function(cube,
                                                 band_cube,
                                                 metric,
                                                 extent,
                                                 filter = NULL,
                                                 stats = NULL,
                                                 impute_fn,
                                                 multicores,
                                                 .verbose = FALSE) {

    # get the file information for the cube
    file_info <- cube$file_info[[1]]

    # define the input raster files for band
    bnd_files <- dplyr::filter(file_info, band == band_cube)$path

    # are there bands associated to the files?
    assertthat::assert_that(length(bnd_files) > 0,
                            msg = paste0(".sits_raster_data_preprocess:
                                         no files for band ", band_cube)
    )

    # read the values
    values <- .sits_raster_api_read_extent(bnd_files, extent)

    # get the missing values, minimum values and scale factors
    missing_value <- .sits_config_missing_values(cube$sensor, band_cube)
    minimum_value <- .sits_config_minimum_values(cube$sensor, band_cube)
    maximum_value <- .sits_config_maximum_values(cube$sensor, band_cube)

    # correct NA, minimum, maximum, and missing values
    values[values < minimum_value] <- NA
    values[values > maximum_value] <- NA
    values[values == missing_value] <- NA

    # does the cube have a cloud band?
    cld_band <- .sits_config_cloud_band(cube)

    if (cld_band %in% sits_bands(cube)) {
        cld_files <- dplyr::filter(file_info, band == cld_band)$path
        clouds <- .sits_raster_api_read_extent(cld_files, extent)
    }
    else {
        clouds <- NULL
    }

    # change the points under clouds to NA
    if (!purrr::is_null(clouds)) {
        cld_index <- .sits_config_cloud_values(cube)
        values[clouds %in% cld_index] <- NA
    }

    # remove NA pixels
    if (any(is.na(values))) {
        if (.verbose) task_start_time <- lubridate::now()

        values <- .sits_raster_data_na_remove(
            values = values,
            impute_fn = impute_fn,
            multicores = multicores
        )

        if (.verbose) {
            .sits_processing_task_time(
                "Impute NA",
                task_start_time
            )
        }
    }

    # scale the data set
    scale_factor <- .sits_config_scale_factors(cube$sensor, band_cube)
    values <- scale_factor * values

    # generates metrics bands
    values <- .sits_raster_metrics(values, metric, multicores)

    # filter the data
    if (!(purrr::is_null(filter))) {
        values <- .sits_raster_data_filter(values, filter, multicores)
    }
    # normalize the data
    # if (!purrr::is_null(stats)) {
    #     values <- .sits_normalize_matrix(values, stats, band_cube, multicores)
    # }

    values_dt <- data.table::as.data.table(values)
    return(values_dt)
}
#' @title Filter the time series values in the case of a matrix
#' @name .sits_raster_metrics
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description This function filters a matrix.
#'
#' @param  values         matrix of values.
#' @param  metric         metric to be process
#' @param  multicores     Number of cores.
#' @return Filtered matrix.
.sits_raster_metrics <- function(values, metric, multicores) {
    # TODO: criar uma lista de função
    if ("min" %in% metric)
        metric_function <- min_ts
    if ("max" %in% metric)
        metric_function <- max_ts
    if ("sum" %in% metric)
        metric_function <- sum_ts
    if ("mean" %in% metric)
        metric_function <- mean_ts
    if ("std" %in% metric)
        metric_function <- std_ts
    if ("skew" %in% metric)
        metric_function <- skew_ts
    if ("kurt" %in% metric)
        metric_function <- kurt_ts
    if ("amplitude" %in% metric)
        metric_function <- amplitude_ts
    if ("fslope" %in% metric)
        metric_function <- fslope_ts
    if ("abs_sum" %in% metric)
        metric_function <- abs_sum_ts
    if ("amd" %in% metric)
        metric_function <- amd_ts
    if ("mse" %in% metric)
        metric_function <- mse_ts
    if ("fqr" %in% metric)
        metric_function <- fqr_ts
    if ("sqr" %in% metric)
        metric_function <- sqr_ts
    if ("tqr" %in% metric)
        metric_function <- tqr_ts
    if ("iqr" %in% metric)
        metric_function <- iqr_ts

    #auxiliary function to scale a block of data
    metrics <- function(chunk , f) {

        values <- f(chunk)
        names(values) <- metric

        return(values)
    }

    if (multicores > 1) {
        chunks <- .sits_raster_data_split(values, multicores)
        rows <- parallel::mclapply(chunks,
                                   metrics,
                                   metric_function,
                                   mc.cores = multicores
        )
        values <- do.call(rbind, rows)
    }
    else {
        values <- metrics(values, metric_function)
    }

    return(values)
}
