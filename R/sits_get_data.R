#' @title Get time series from data cubes and cloud services
#' @name sits_get_data
#' @author Gilberto Camara
#'
#' @description Retrieve a set of time series from a data cube or from
#' a time series service. Data cubes and puts it in a "sits tibble".
#' Sits tibbles are the main structures of sits package.
#' They contain both the satellite image time series and their metadata.
#'
#' @note
#' There are four ways of specifying data to be retrieved using the
#' "samples" parameter:
#' \itemize{
#' \item{CSV file: }{Provide a CSV file with columns
#' "longitude", "latitude", "start_date", "end_date" and "label" for
#' each sample}
#' \item{SHP file: }{Provide a shapefile in POINT or POLYGON geometry
#' containing the location of the samples and an attribute to be
#' used as label. Also, provide start and end date for the time series.}
#' \item{sits object: }{A sits tibble.}
#' \item{sf object: }{An "sf" object with POINT or POLYGON geometry.}
#' \item{data.frame: }{A data.frame with with mandatory columns
#' "longitude", "latitude".}
#' }
#
#' @param cube            Data cube from where data is to be retrieved.
#' @param samples         Samples location (sits, sf, or data.frame).
#' @param ...             Specific parameters for specific cases.
#' @param start_date      Start of the interval for the time series
#'                        in "YYYY-MM-DD" format (optional).
#' @param end_date        End of the interval for the time series in
#'                        "YYYY-MM-DD" format (optional).
#' @param label           Label to be assigned to the time series (optional).
#' @param bands           Bands to be retrieved (optional).
#' @param crs             A coordinate reference system of samples.
#'                        The provided crs could be a character
#'                        (e.g, "EPSG:4326" or "WGS84" or a proj4string), or a
#'                        a numeric with the EPSG code (e.g. 4326).
#'                        This parameter only works for 'csv' or data.frame'
#'                        samples. Default is 4326.
#' @param impute_fn       Imputation function for NA values.
#' @param label_attr      Attribute in the shapefile or sf object to be used
#'                        as a polygon label.
#' @param n_sam_pol       Number of samples per polygon to be read
#'                        (for POLYGON or MULTIPOLYGON shapefile).
#' @param pol_avg         Summarize samples for each polygon?
#' @param pol_id          ID attribute for polygons.
#' @param output_dir      Directory where the time series will be saved as rds.
#'                        Default is the current path.
#' @param multicores      Number of threads to process the time series.
#' @param progress        A logical value indicating if a progress bar
#'                        should be shown. Default is \code{FALSE}.
#'
#' @return A tibble with the metadata and data for each time series
#' <longitude, latitude, start_date, end_date, label, cube, time_series>.
#'
#' @note
#' Please refer to the sits documentation available in
#' <https://e-sensing.github.io/sitsbook/> for detailed examples.
#'
#' @examples
#' if (sits_run_examples()) {
#'     # reading a lat/long from a local cube
#'     # create a cube from local files
#'     data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
#'     raster_cube <- sits_cube(
#'         source = "BDC",
#'         collection = "MOD13Q1-6",
#'         data_dir = data_dir,
#'         delim = "_",
#'         parse_info = c("X1", "tile", "band", "date")
#'     )
#'     samples <- tibble::tibble(longitude = -55.66738, latitude = -11.76990)
#'     point_ndvi <- sits_get_data(raster_cube, samples)
#'     #
#'     # reading samples from a cube based on a  CSV file
#'     csv_file <- system.file("extdata/samples/samples_sinop_crop.csv",
#'         package = "sits"
#'     )
#'     points <- sits_get_data(cube = raster_cube, samples = csv_file)
#'
#'     # reading a shapefile from BDC (Brazil Data Cube)
#'     # needs a BDC access key that can be obtained
#'     # for free by registering in the BDC website
#'     if (nchar(Sys.getenv("BDC_ACCESS_KEY")) > 0) {
#'         # create a data cube from the BDC
#'         bdc_cube <- sits_cube(
#'             source = "BDC",
#'             collection = "CB4_64_16D_STK-1",
#'             bands = c("NDVI", "EVI"),
#'             tiles = c("022024", "022025"),
#'             start_date = "2018-09-01",
#'             end_date = "2018-10-28"
#'         )
#'         # define a shapefile to be read from the cube
#'         shp_file <- system.file("extdata/shapefiles/bdc-test/samples.shp",
#'             package = "sits"
#'         )
#'         # get samples from the BDC based on the shapefile
#'         time_series_bdc <- sits_get_data(
#'             cube = bdc_cube,
#'             samples = shp_file)
#'     }
#' }
#'
#' @export
sits_get_data <- function(cube,
                          samples,
                          ...,
                          start_date = NULL,
                          end_date = NULL,
                          label = "NoClass",
                          bands = NULL,
                          crs = NULL,
                          impute_fn = sits_impute_linear(),
                          label_attr = NULL,
                          n_sam_pol = 30,
                          pol_avg = FALSE,
                          pol_id = NULL,
                          multicores = 2,
                          output_dir = getwd(),
                          progress = FALSE) {

    # Pre-conditions
    .check_is_raster_cube(cube)
    .check_is_regular(cube)
    # Get default cube bands if not provided
    bands <- .default(bands, .cube_bands(cube, add_cloud = TRUE))
    .check_bands_in_cube(bands = bands, cube = cube)
    # Get cube timeline
    tl <- .cube_timeline(cube)
    # Get start_date and end_date as default
    start_date <- .default(start_date, as.Date(tl[[1]]))
    end_date <- .default(end_date, as.Date(tl[length(tl)]))
    .check_dates_parameter(c(start_date, end_date))
    .check_lgl_parameter(pol_avg)
    .check_num_parameter(n_sam_pol)
    .check_multicores(multicores)
    .check_output_dir(output_dir)
    .check_progress(progress)
    .check_that(
        !(pol_avg && !.has(pol_id)),
        msg = "invalid 'pol_id' parameter."
    )
    crs <- .default(crs, 4326)
    .check_crs(crs)

    # Transform samples to tibble format
    samples <- .samples_format_tbl(
        samples    = samples,
        start_date = start_date,
        end_date   = end_date,
        label      = label,
        bands      = bands,
        crs        = crs,
        label_attr = label_attr,
        n_sam_pol  = n_sam_pol,
        pol_avg    = pol_avg,
        pol_id     = pol_id
    )
    # Check samples format
    .check_samples(samples)

    # Get time series
    data <- .samples_get_data(
        cube       = cube,
        samples    = samples,
        bands      = bands,
        impute_fn  = impute_fn,
        multicores = multicores,
        crs        = crs,
        output_dir = output_dir,
        progress   = progress
    )
    if (pol_avg && "polygon_id" %in% colnames(data)) {
        data <- .sits_avg_polygon(data = data)
    }
}

#' @rdname sits_get_data
#'
#' @export
sits_get_data.default <- function(cube, samples, ...) {
    stop("Invalid samples parameter for sits_get_data")
}


#' @title Dispatch function to get time series from data cubes and cloud
#' services
#' @name .samples_get_data
#' @author Gilberto Camara
#' @keywords internal
#' @noRd
#' @param cube            Data cube from where data is to be retrieved.
#' @param samples         Samples to be retrieved.
#' @param crs             A coordinate reference system of samples.
#'                        The provided crs could be a character
#'                        (e.g, "EPSG:4326" or "WGS84" or a proj4string), or a
#'                        a numeric with the EPSG code (e.g. 4326).
#'                        This parameter only works for 'csv' or data.frame'
#'                        samples. Default is 4326.
#' @param bands           Bands to be retrieved (optional).
#' @param impute_fn       Imputation function for NA values.
#' @param multicores      Number of threads to process the time series.
#' @param progress        A logical value indicating if a progress bar
#'                        should be shown. Default is \code{FALSE}.
#'
#' @return                A tibble with a set of time series retrieved
#'                        from a data cube.
.samples_get_data <- function(cube,
                              samples,
                              bands,
                              crs,
                              impute_fn,
                              multicores,
                              output_dir,
                              progress) {

    # Dispatch
    UseMethod(".samples_get_data", cube)
}

#' @name .samples_get_data
#' @keywords internal
#' @noRd
#' @export
.samples_get_data.raster_cube <- function(cube,
                                          samples,
                                          bands,
                                          crs,
                                          impute_fn,
                                          multicores,
                                          output_dir,
                                          progress) {
    # sits forces samples to WGS84, so we reproject otherwise
    if (crs != 4326) {
        samples <- .proj_transform_samples(samples, crs = crs)
    }

    # Is the cloud band available?
    cld_band    <- .source_cloud()
    if (cld_band %in% bands) {
        bands <- bands[bands != cld_band]
    } else {
        cld_band <- NULL
    }

    # get cubes timeline
    tl <- sits_timeline(cube)

    tiles_bands <- tidyr::expand_grid(tile = .cube_tiles(cube),
                                      band = bands) %>%
        purrr::pmap(function(tile, band) {
            return(list(tile, band))
        })

    # prepare parallelization
    .sits_parallel_start(workers = multicores, log = FALSE)
    on.exit(.sits_parallel_stop(), add = TRUE)

    samples_tiles_bands <- .sits_parallel_map(tiles_bands, function(tile_band) {
        tile_id <- tile_band[[1]]
        band <- tile_band[[2]]

        tile <- sits_select(cube, bands = c(band, cld_band), tiles = tile_id)

        hash_bundle <- digest::digest(list(tile, samples), algo = "md5")

        filename <- .file_path(
            "samples", hash_bundle,
            ext = ".rds",
            output_dir = output_dir
        )

        if (file.exists(filename)) {
            tryCatch(
                {
                    # ensuring that the file is not corrupted
                    timeseries <- readRDS(filename)

                    return(timeseries)
                },
                error = function(e) {
                    unlink(filename)
                    gc()
                }
            )
        }

        # get XY
        xy_tb <- .proj_from_latlong(
            longitude = samples[["longitude"]],
            latitude  = samples[["latitude"]],
            crs       = .cube_crs(tile)
        )
        # join lat-long with XY values in a single tibble
        samples <- dplyr::bind_cols(samples, xy_tb)
        # filter the points inside the data cube space-time extent
        samples <- dplyr::filter(
            samples,
            .data[["X"]] > tile$xmin & .data[["X"]] < tile$xmax &
                .data[["Y"]] > tile$ymin & .data[["Y"]] < tile$ymax &
                .data[["start_date"]] <= as.Date(tl[length(tl)]) &
                .data[["end_date"]] >= as.Date(tl[1])
        )
        # are there points to be retrieved from the cube?
        if (nrow(samples) == 0) {
            return(NULL)
        }
        # create a matrix to extract the values
        xy <- matrix(
            c(samples[["X"]], samples[["Y"]]),
            nrow = nrow(samples),
            ncol = 2
        )
        colnames(xy) <- c("X", "Y")


        # build the sits tibble for the storing the points
        samples_tbl <- slider::slide_dfr(samples, function(point) {

            # get the valid timeline
            dates <- .timeline_during(
                timeline   = tl,
                start_date = as.Date(point[["start_date"]]),
                end_date   = as.Date(point[["end_date"]])
            )
            sample <- tibble::tibble(
                longitude  = point[["longitude"]],
                latitude   = point[["latitude"]],
                start_date = dates[[1]],
                end_date   = dates[[length(dates)]],
                label      = point[["label"]],
                cube       = tile[["collection"]],
                polygon_id = point[["polygon_id"]]
            )
            # store them in the sample tibble
            sample$time_series <- list(tibble::tibble(Index = dates))
            # return valid row of time series
            return(sample)
        })

        # extract time series
        ts <- .raster_data_get_ts(
            tile = tile,
            points = samples_tbl,
            bands = band,
            xy = xy,
            cld_band = cld_band,
            impute_fn = impute_fn,
            output_dir = output_dir
        )

        ts[["tile"]] <- tile_id
        ts[["#..id"]] <- seq_len(nrow(ts))

        saveRDS(ts, filename)

        return(ts)
    }, progress = progress)

    ts_tbl <- samples_tiles_bands %>%
        dplyr::bind_rows() %>%
        tidyr::unnest("time_series") %>%
        dplyr::group_by(
            .data[["longitude"]], .data[["latitude"]],
            .data[["start_date"]], .data[["end_date"]],
            .data[["label"]], .data[["cube"]],
            .data[["Index"]], .data[["tile"]], .data[["#..id"]]
        )

    if ("polygon_id" %in% colnames(ts_tbl)) {
        ts_tbl <- dplyr::group_by(ts_tbl, .data[["polygon_id"]], .add = TRUE)
    }

    ts_tbl <- ts_tbl %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(bands), stats::na.omit)) %>%
        dplyr::arrange(.data[["Index"]]) %>%
        dplyr::ungroup() %>%
        tidyr::nest(time_series = !!c("Index", bands)) %>%
        dplyr::select(-c("tile", "#..id"))


    # get the first point that intersect more than one tile
    # eg sentinel 2 mgrs grid
    ts_tbl <- ts_tbl %>%
        dplyr::group_by(
            .data[["longitude"]], .data[["latitude"]],
            .data[["start_date"]], .data[["end_date"]],
            .data[["label"]], .data[["cube"]]) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()

    # recreate hash values
    hash_bundle <- purrr::map_chr(tiles_bands, function(tile_band) {
        tile_id <- tile_band[[1]]
        band <- tile_band[[2]]
        tile <- sits_select(cube, bands = c(band, cld_band), tiles = tile_id)
        digest::digest(list(tile, samples), algo = "md5")
    })

    # recreate file names to delete them
    # samples will be recycled for each hash_bundle
    temp_timeseries <- .file_path(
        "samples", hash_bundle,
        ext = "rds",
        output_dir = output_dir
    )

    # delete temporary rds
    unlink(temp_timeseries)
    gc()

    # check if data has been retrieved
    .sits_get_data_check(nrow(samples), nrow(ts_tbl))

    if (!inherits(ts_tbl, "sits")) {
        class(ts_tbl) <- c("sits", class(ts_tbl))
    }

    return(ts_tbl)
}

.samples_data_read <- function(tile, chunk, points, xy, cld_band, impute_fn,
                               output_dir) {
    block <- .block(chunk)
    # For cubes that have a time limit to expire (MPC cubes only)
    tile <- .cube_token_generator(tile)
    # Read and preprocess values of cloud
    # Get cloud values (NULL if not exists)
    cloud_mask <- .tile_cloud_read_block(tile = tile, block = block)
    # Read and preprocess values of each band

    bands <- .cube_bands(tile, add_cloud = FALSE)
    values <- purrr::map_dfc(bands, function(band) {

        # get the scale factors, max, min and missing values
        band_params   <- .tile_band_conf(tile, band)
        missing_value <- .miss_value(band_params)
        minimum_value <- .min_value(band_params)
        maximum_value <- .max_value(band_params)
        scale_factor  <- .scale(band_params)
        offset_value  <- .offset(band_params)

        # Get band values (stops if band not found)
        values <- .tile_read_block(tile = tile, band = band, block = block)

        # Remove cloud masked pixels
        if (.has(cloud_mask)) {
            values[cloud_mask] <- NA
        }
        # Remove NA pixels
        if (.has(impute_fn)) {
            values <- impute_fn(values)
        }

        # Create a new raster
        bbox <- .bbox(chunk)
        r_obj <- .raster_new_rast(
            nrows = block[["nrows"]], ncols = block[["ncols"]],
            xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
            ymin = bbox[["ymin"]], ymax = bbox[["ymax"]],
            nlayers = length(.tile_timeline(tile)), crs = bbox[["crs"]]
        )
        # Copy values
        r_obj <- .raster_set_values(
            r_obj = r_obj,
            values = values
        )

        # Return values
        as.data.frame(values)
    })
    # Compose final values
    values <- as.matrix(values)
    # Return values
    values
}

#' @name .samples_get_data
#' @keywords internal
#' @noRd
#' @export
.samples_get_data.class_cube <- function(cube,
                                         samples,
                                         bands,
                                         crs,
                                         impute_fn,
                                         multicores,
                                         output_dir,
                                         progress) {

    # Filter only tiles that intersects with samples
    cube <- .cube_filter_spatial(
        cube = cube,
        roi = .point_as_sf(point = .point(x = samples, crs = crs))
    )

    # pre-condition - check bands
    if (is.null(bands)) {
        bands <- .cube_bands(cube)
    }

    .check_cube_bands(cube, bands = bands)

    # get cubes timeline
    tl <- sits_timeline(cube)

    tiles_bands <- tidyr::expand_grid(tile = .cube_tiles(cube),
                                      band = bands) %>%
        purrr::pmap(function(tile, band) {
            return(list(tile, band))
        })

    # prepare parallelization
    .sits_parallel_start(workers = multicores, log = FALSE)
    on.exit(.sits_parallel_stop(), add = TRUE)

    samples_tiles_bands <- .sits_parallel_map(tiles_bands, function(tile_band) {
        tile_id <- tile_band[[1]]
        band <- tile_band[[2]]

        tile <- sits_select(cube, bands = band, tiles = tile_id)

        hash_bundle <- digest::digest(list(tile, samples), algo = "md5")

        filename <- .file_path(
            "samples", hash_bundle,
            ext = ".rds",
            output_dir = output_dir
        )

        if (file.exists(filename)) {
            tryCatch(
                {
                    # ensuring that the file is not corrupted
                    timeseries <- readRDS(filename)

                    return(timeseries)
                },
                error = function(e) {
                    unlink(filename)
                    gc()
                }
            )
        }

        # get XY
        xy_tb <- .proj_from_latlong(
            longitude = samples[["longitude"]],
            latitude  = samples[["latitude"]],
            crs       = .cube_crs(tile)
        )
        # join lat-long with XY values in a single tibble
        samples <- dplyr::bind_cols(samples, xy_tb)
        # filter the points inside the data cube space-time extent
        samples <- dplyr::filter(
            samples,
            .data[["X"]] > tile$xmin & .data[["X"]] < tile$xmax &
                .data[["Y"]] > tile$ymin & .data[["Y"]] < tile$ymax &
                .data[["start_date"]] <= as.Date(tl[length(tl)]) &
                .data[["end_date"]] >= as.Date(tl[1])
        )
        # are there points to be retrieved from the cube?
        if (nrow(samples) == 0) {
            return(NULL)
        }
        # create a matrix to extract the values
        xy <- matrix(
            c(samples[["X"]], samples[["Y"]]),
            nrow = nrow(samples),
            ncol = 2
        )
        colnames(xy) <- c("X", "Y")
        # build the sits tibble for the storing the points
        samples_tbl <- slider::slide_dfr(samples, function(point) {

            # get the valid timeline
            dates <- .timeline_during(
                timeline   = tl,
                start_date = as.Date(point[["start_date"]]),
                end_date   = as.Date(point[["end_date"]])
            )
            sample <- tibble::tibble(
                longitude  = point[["longitude"]],
                latitude   = point[["latitude"]],
                start_date = dates[[1]],
                end_date   = dates[[length(dates)]],
                label      = point[["label"]],
                cube       = tile[["collection"]],
                polygon_id = point[["polygon_id"]]
            )
            # store them in the sample tibble
            sample$predicted <- list(tibble::tibble(
                from = dates[[1]], to = dates[[2]])
            )
            # return valid row of time series
            return(sample)
        })
        ts <- .raster_class_get_ts(
            tile = tile,
            points = samples_tbl,
            band = "class",
            xy = xy,
            output_dir = output_dir
        )

        ts[["tile"]] <- tile_id
        ts[["#..id"]] <- seq_len(nrow(ts))

        saveRDS(ts, filename)

        return(ts)
    }, progress = progress)

    ts_tbl <- samples_tiles_bands %>%
        dplyr::bind_rows() %>%
        tidyr::unnest("predicted") %>%
        dplyr::group_by(
            .data[["longitude"]], .data[["latitude"]],
            .data[["start_date"]], .data[["end_date"]],
            .data[["label"]], .data[["cube"]],
            .data[["from"]], .data[["to"]], .data[["tile"]],
            .data[["#..id"]]
        )

    if ("polygon_id" %in% colnames(ts_tbl)) {
        ts_tbl <- dplyr::group_by(ts_tbl, .data[["polygon_id"]], .add = TRUE)
    }

    ts_tbl <- ts_tbl %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(bands), stats::na.omit)) %>%
        dplyr::arrange(.data[["from"]]) %>%
        dplyr::ungroup() %>%
        tidyr::nest(predicted = !!c("from", "to", bands)) %>%
        dplyr::select(-c("tile", "#..id"))

    # get the first point that intersect more than one tile
    # eg sentinel 2 mgrs grid
    ts_tbl <- ts_tbl %>%
        dplyr::group_by(.data[["longitude"]], .data[["latitude"]],
                        .data[["start_date"]], .data[["end_date"]],
                        .data[["label"]], .data[["cube"]]) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()

    # recreate hash values
    hash_bundle <- purrr::map_chr(tiles_bands, function(tile_band) {
        tile_id <- tile_band[[1]]
        band <- tile_band[[2]]
        tile <- sits_select(cube, bands = band, tiles = tile_id)
        digest::digest(list(tile, samples), algo = "md5")
    })

    # recreate file names to delete them
    # samples will be recycled for each hash_bundle
    temp_timeseries <- .file_path(
        "samples", hash_bundle,
        ext = "rds",
        output_dir = output_dir
    )

    # delete temporary rds
    unlink(temp_timeseries)
    gc()

    # check if data has been retrieved
    .sits_get_data_check(nrow(samples), nrow(ts_tbl))

    class(ts_tbl) <- unique(c("predicted", "sits", class(ts_tbl)))

    return(ts_tbl)
}

#' @title Check if all points have been retrieved
#' @name .sits_get_data_check
#' @keywords internal
#' @noRd
#' @param n_rows_input     Number of rows in input.
#' @param n_rows_output    Number of rows in output.
#'
#' @return No return value, called for side effects.
#'
.sits_get_data_check <- function(n_rows_input, n_rows_output) {

    # Have all input rows being read?
    if (n_rows_output == 0) {
        message("No points have been retrieved")
        return(invisible(FALSE))
    }

    if (n_rows_output < n_rows_input) {
        message("Some points could not be retrieved")
    } else {
        message("All points have been retrieved")
    }
}

#' @title Extracts the time series average by polygon.
#' @name .sits_avg_polygon
#' @keywords internal
#' @noRd
#' @description This function extracts the average of the automatically
#' generated points for each polygon in a shapefile.
#'
#' @param data A sits tibble with points time series.
#'
#' @return A sits tibble with the average of all points by each polygon.
.sits_avg_polygon <- function(data) {
    bands <- sits_bands(data)
    columns_to_avg <- c(bands, "latitude", "longitude")

    data_avg <- data %>%
        tidyr::unnest(cols = "time_series") %>%
        dplyr::group_by(
            .data[["Index"]],
            .data[["start_date"]],
            .data[["end_date"]],
            .data[["label"]],
            .data[["cube"]],
            .data[["polygon_id"]]
        ) %>%
        dplyr::summarise(dplyr::across(!!columns_to_avg, mean, na.rm = TRUE),
                         .groups = "drop"
        ) %>%
        tidyr::nest("time_series" = c("Index", dplyr::all_of(bands))) %>%
        dplyr::select(!!colnames(data))

    class(data_avg) <- class(data)

    return(data_avg)
}

# ---- Samples API ----

.is_sf <- function(data) {
    inherits(data, "sf")
}

.are_sits_samples <- function(samples) {
    inherits(samples, "sits")
}

.are_samples_polygons <- function(samples, pol_avg) {
    pol_avg && "polygon_id" %in% colnames(samples)
}

are_samples_empty <- function(samples) {
    nrow(samples) == 0
}

.samples_type <- function(samples) {
    if (.is_sf(samples)) {
        return("sf")
    }
    if (is.data.frame(samples)) {
        return("data.frame")
    }
    if (is.character(samples)) {
        .check_file(samples, extensions = c("csv", "shp"))
        return(tolower(.file_ext(samples)))
    }
    if (.are_sits_samples(samples)) {
        return("sits")
    }
    else {
        stop("invalid 'samples' parameter.")
    }
}

.samples_switch <- function(samples, ...) {
    switch(.samples_type(samples), ...)
}

.samples_transform <- function(samples, crs, as_crs) {
    pts_coords <- .point_coords(
        .point(x = samples, crs = crs, as_crs = as_crs)
    )
    samples[, c("longitude", "latitude")] <- pts_coords[, c("X", "Y")]
    samples
}

.samples_format_tbl <- function(samples,
                                start_date,
                                end_date,
                                label,
                                bands,
                                crs,
                                label_attr,
                                n_sam_pol,
                                pol_avg,
                                pol_id) {
    .samples_switch(
        samples,
        "data.frame" = .samples_from_df(
            samples = samples,
            label = label,
            start_date = start_date,
            end_date = end_date
        ),
        "sits" = samples,
        "csv" = .samples_from_csv(samples),
        "shp" = .samples_from_shp(
            shp_file    = samples,
            label       = label,
            shp_attr    = label_attr,
            start_date  = start_date,
            end_date    = end_date,
            n_sam_pol   = n_sam_pol,
            pol_avg     = pol_avg,
            pol_id      = pol_id
        ),
        "sf" = .samples_from_sf(
            sf_object     = samples,
            label         = label,
            label_attr    = label_attr,
            start_date    = start_date,
            end_date      = end_date,
            n_sam_pol     = n_sam_pol,
            pol_id        = pol_id
        )
    )
}

`.set_samples_label<-` <- function(samples, value) {
    samples[["label"]] <- value
    samples
}

`.set_samples_start_date<-` <- function(samples, value) {
    samples[["start_date"]] <- value
    samples
}

`.set_samples_end_date<-` <- function(samples, value) {
    samples[["end_date"]] <- value
    samples
}

.samples_from_df <- function(samples, label, start_date, end_date) {
    # check if samples contains all the required columns
    .check_chr_contains(
        x = colnames(samples),
        contains = c("latitude", "longitude"),
        discriminator = "all_of",
        msg = "missing lat/long information in data frame"
    )
    # Fill missing columns
    if (!("label" %in% colnames(samples))) {
        .set_samples_label(samples) <- label
    }
    if (!("start_date" %in% colnames(samples))) {
        .set_samples_start_date(samples) <- start_date
    }
    if (!("end_date" %in% colnames(samples))) {
        .set_samples_end_date(samples) <- end_date
    }
    class(samples) <- c("sits", class(samples))
    samples
}
