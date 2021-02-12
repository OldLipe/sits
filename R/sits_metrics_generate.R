#' @title Clean data cube to improve quality
#' @name  sits_metrics_generate
#' @description Interpolate data over time to fill cloud pixels.
#'
#' @param cube       input data cube
#' @param data_dir   data directory where output data is written
#' @param name       name of the output data cube
#' @param n_metrics  number of metrics to be generate
#' @param impute_fn  imputing function to be applied to replace NA
#' @param memsize    size of memory
#' @param multicores number of cores
#'
#' @return           new data cube with interpolated cloud data
#'
#' @examples
#' # define a data cube of CBERS-4 AWFI data
#' data_dir <- system.file("extdata/raster/cbers", package = "sits")
#'
#' cbers_022024 <- sits_cube(
#'     type = "STACK",
#'     name = "cbers_022024",
#'     satellite = "CBERS-4",
#'     sensor = "AWFI",
#'     resolution = "64m",
#'     data_dir = data_dir,
#'     parse_info = c("X1", "X2", "band", "date")
#' )
#'
#' cbers_022024_no_clds <- sits_cloud_remove(
#'     cube = cbers_022024,
#'     data_dir = tempdir(),
#'     name = "cbers_022024_no_cld"
#' )
#' @export
sits_metrics_generate <- function(cube,
                                  data_dir,
                                  name,
                                  n_metrics = 8,
                                  impute_fn = sits_impute_linear(),
                                  memsize = 8,
                                  multicores = 2) {

    # precondition - is there cloud information available?
    cloud_band <- .sits_config_cloud_band(cube)
    assertthat::assert_that(cloud_band %in% sits_bands(cube),
                            msg = "cloud information band not available in cube"
    )
    # precondition - does the data directory exists?
    assertthat::assert_that(assertthat::is.dir(data_dir),
                            msg = "invalid data directory"
    )

    # total number of instances
    n_instances <- length(sits_timeline(cube))

    # estimate the blocks to be read
    blocks <- .metrics_blocks_estimate(
        cube = cube,
        n_metrics = n_metrics,
        memsize = memsize
    )

    # interpolate the cloud bricks
    files <- .sits_create_metrics_bands(
        cube = cube,
        data_dir = data_dir,
        blocks = blocks,
        impute_fn = impute_fn,
        multicores = multicores
    )

    # find out what are the bands of the cube
    bands <- sits_bands(cube)
    bands <- bands[bands != cloud_band]

    # create the output cube
    cube_new <- .sits_raster_brick_cube(
        satellite = cube$satellite,
        sensor = cube$sensor,
        name = name,
        timeline = sits_timeline(cube),
        bands = bands,
        files = files
    )

    class(cube_new) <- c("brick_cube", "raster_class", class(cube_new))

    return(cube_new)
}
#' @title Create the output of the cloud estimation procedure
#' @name .sits_create_metrics_bands
#' @keywords internal
#'
#' @param cube        input data cube
#' @param data_dir    directory where data is to be stored
#' @param blocks      block information
#' @param impute_fn   imputation function to remove NA
#' @param multicores  number of cores to use
#'
#' @return            a tibble with date, band and path information.
.sits_create_metrics_bands <- function(cube,
                                       data_dir,
                                       blocks,
                                       impute_fn,
                                       multicores) {

    # get initial time for classification
    start_time <- lubridate::now()
    message(sprintf("Starting creating metrics bands at %s", start_time))

    # define the bands
    cloud_band <- .sits_config_cloud_band(cube)
    bands <- sits_bands(cube)
    # ensure that the cloud band is available
    assertthat::assert_that(cloud_band %in% sits_bands(cube),
                            msg = ".sits_clouds_interpolate: no cloud band"
    )
    # define the bands that are not associated to clouds
    bands_no_cloud <- bands[bands != cloud_band]

    # get the file information from the cube
    file_info <- cube$file_info[[1]]

    # define the parameters for the output files
    params = .sits_raster_api_params_cube(cube)
    # process the bands
    file_list <- purrr::map(bands_no_cloud, function(bnd) {
        message(paste0("Creating metrics from band: ", bnd))
        start_task_time <- lubridate::now()
        # find out the information about the band
        info_band <- dplyr::filter(file_info, band == bnd)
        # what is the number of layers?
        num_layers <- nrow(info_band)

        # what are the start and end date?
        start_date <- info_band[1, ]$date
        end_date <- info_band[num_layers, ]$date

        # define the output filename
        filename <- paste0(
            data_dir, "/",
            cube$satellite, "_",
            cube$sensor, "_",
            start_date, "_", end_date, "_",
            bnd, "MAX_METRIC", ".tif"
        )

        # create a raster object
        r_obj <- suppressWarnings(
            terra::rast(
                nrows = params$nrows,
                ncols = params$ncols,
                nlyrs = 8,
                xmin = params$xmin,
                xmax = params$xmax,
                ymin = params$ymin,
                ymax = params$ymax,
                crs = params$crs
            )
        )

        assertthat::assert_that(terra::nrow(r_obj) == params$nrows,
                                msg = ".sits_raster_api_write: unable to create raster object"
        )
        # open the file with writeStart
        suppressWarnings(terra::writeStart(
            r_obj,
            filename = filename,
            overwrite = TRUE,
            wopt = list(
                gdal = c("COMPRESS = LZW"),
                filetype = "GTiff",
                datatype = "INT2U"
            )
        ))

        # read the blocks
        bs <- purrr::map(c(1:blocks$n), function(b) {

            # measure performance
            start_block_time <- lubridate::now()
            # define the extent
            extent <- c(
                blocks$row[b], blocks$nrows[b],
                blocks$col, blocks$ncols
            )
            names(extent) <- (c("row", "nrows", "col", "ncols"))

            # preprocess the input data
            # produce a matrix
            values_block <- .sits_raster_metrics_data_preprocess(
                cube = cube,
                band_cube = bnd,
                extent = extent,
                impute_fn = impute_fn,
                multicores = multicores
            )
            #browser()
            # rescale the data
            mult_factor <- 1/as.numeric(cube$scale_factors[[1]][bnd])
            values_block <- mult_factor * values_block
            # write a block of values

            terra::writeValues(r_obj,
                               as.matrix(values_block),
                               start = blocks$row[b],
                               nrows = blocks$nrows[b])

            task <- paste0("process block ", b, " of band ", bnd)
            .sits_processing_task_time(task, start_block_time)

            return(b)
        })
        # close the file
        terra::writeStop(r_obj)

        task <- paste0("Removed clouds from band ", bnd)
        .sits_processing_task_time(task, start_task_time)
        return(filename)
    })

    # report on time used for processing
    task <- paste0("Removed clouds from all bands")
    .sits_processing_task_time(task, start_time)

    # return the file info
    files <- unlist(file_list)
    return(files)
}

#' @title ...
#' @description ...
#'
#' @param cube ...
#' @param n_metrics ...
#' @param memsize ...
#'
#' @return a list ...
.metrics_blocks_estimate <- function(cube, n_metrics, memsize) {

    #n_bands <- length(unique(cube$file_info[[1]]$band))
    n_bands <- 4
    n_instances <- length(sits_timeline(cube))

    # number of bytes per pixel
    nbytes <- 8
    # estimated memory bloat
    bloat <- as.numeric(sits:::.sits_config_memory_bloat())

    # number of rows and cols
    nrows <- as.numeric(cube$nrows)
    ncols <- as.numeric(cube$ncols)

    # the full output band

    # single instance depends on the number of bands
    single_data_size <- nrows * ncols * nbytes * n_bands * 1

    # estimated full size of the data
    full_size <- as.numeric(n_instances) * single_data_size

    # estimated size of memory required
    mem_required <- (full_size + as.numeric(sits:::.sits_mem_used())) * bloat

    # number of passes to read the full data sets
    nblocks <- ceiling(mem_required / (memsize * 1e+09))

    # Cloud processing uses the whole image
    sub_image <- sits:::.sits_raster_sub_image_default(cube)

    # calculate the blocks
    blocks <- sits:::.sits_raster_block_list(
        nblocks = nblocks,
        sub_image = sub_image
    )

    return(blocks)
}
