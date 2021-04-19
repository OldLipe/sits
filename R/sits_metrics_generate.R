#' @title Generate bands metrics from cube images
#' @name  sits_metrics_generate
#' @description Generate a set of metrics by band
#'
#' @param cube          input data cube
#' @param output_dir    data directory where output data is written
#' @param name          name of the output data cube
#' @param list_metrics  number of metrics to be generate
#' @param impute_fn     imputing function to be applied to replace NA
#' @param memsize       size of memory
#' @param multicores    number of cores
#'
#' @return              new data cube of metrics
#'
#' @examples
#' # define a data cube of CBERS-4 AWFI data
#' data_dir <- system.file("extdata/raster/cbers", package = "sits")
#'
#' cbers_022024 <- sits_cube(
#'     source = "LOCAL",
#'     name = "cbers_022024",
#'     satellite = "CBERS-4",
#'     sensor = "AWFI",
#'     resolution = "64m",
#'     data_dir = data_dir,
#'     parse_info = c("X1", "X2", "band", "date")
#' )
#'
#' cbers_022024_no_clds <- sits_metrics_generate(
#'     cube = cbers_022024,
#'     output_dir = tempdir(),
#'     name = "cbers_022024_no_cld",
#'     list_metrics = list("B13" = c("max", "min")),
#'    impute_fn = sits_impute_linear(),
#'    memsize = 8,
#'    multicores = 2
#' )
#' @export
sits_metrics_generate <- function(cube,
                                  output_dir,
                                  name,
                                  list_metrics = list("EVI" = c("min"),
                                                      "NDVI" = "min"),
                                  impute_fn = sits_impute_linear(),
                                  memsize = 8,
                                  multicores = 2) {

    # precondition - is there cloud information available?
    cloud_band <- .sits_config_cloud_band(cube)
    assertthat::assert_that(cloud_band %in% sits_bands(cube),
                            msg = "cloud information band not available in cube"
    )
    # precondition - does the data directory exists?
    assertthat::assert_that(assertthat::is.dir(output_dir),
                            msg = "invalid data directory"
    )

    cube_metrics <- purrr::map(seq_len(nrow(cube)), function(row) {
        # total number of instances
        n_instances <- length(sits_timeline(cube[row,]))

        # estimate the blocks to be read
        blocks <- .metrics_blocks_estimate(
            cube = cube[row,],
            #n_metrics = sum(lengths(list_metrics)),
            nbands = length(names(list_metrics)),
            memsize = memsize
        )

        # interpolate the cloud bricks
        files <- .sits_create_metrics_bands(
            cube = cube[row,],
            output_dir = output_dir,
            blocks = blocks,
            list_metrics = list_metrics,
            impute_fn = impute_fn,
            multicores = multicores
        )

        # find out what are the bands of the cube
        bands <- sits_bands(cube[row,])
        bands <- bands[bands != cloud_band]

        # transform bands name
        bands <- .sits_config_bands_convert(satellite = cube[row,]$satellite,
                                            sensor = cube[row,]$sensor,
                                            bands_files = bands)

        bands <- as.vector(.transform_bands_to_metrics(bands, list_metrics))

        # create the output cube
        cube_new <- .sits_metrics_local_cube(
            satellite = cube[row,]$satellite,
            sensor = cube[row,]$sensor,
            name = name,
            timeline = sits_timeline(cube[row,]),
            bands = bands,
            files = files,
            cube = cube[row,])

        cube_new
    })

    cube_metrics <- dplyr::bind_rows(cube_metrics)

    class(cube_metrics) <- c("raster_cube", class(cube_metrics))

    return(cube_metrics)
}

#' @title Compare sits bands with metrics bands
#' @name .transform_bands_to_metrics
#' @keywords internal
#'
#' @param bands a \code{character} with bands name supported by sits
#' @param list_metrics a named \code{list} with bands name and metrics to be
#' generated
#'
#' @return a \code{character} vector with  merged bands name and metrics.
.transform_bands_to_metrics <- function(bands, list_metrics) {

    if (!all(names(list_metrics) %in% bands))
        stop("The provided bands doesnt match with supported bands in SITS.")

    bands <- sapply(seq_len(length(list_metrics)), function(i) {
        paste0(names(list_metrics[i]), list_metrics[[i]])
    }) %>% unlist()

    bands
}

#' @title Create the images metrics
#' @name .sits_create_metrics_bands
#' @keywords internal
#'
#' @param cube         input data cube
#' @param output_dir   directory where data is to be stored
#' @param blocks       block information
#' @param list_metrics a named \code{list} with bands name and metrics to be
#' generated
#' @param impute_fn    imputation function to remove NA
#' @param multicores   number of cores to use
#'
#' @return             a tibble with date, band and path information.
.sits_create_metrics_bands <- function(cube,
                                       output_dir,
                                       blocks,
                                       list_metrics,
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

    # get the subset of bands selected
    bands_metrics <- bands_no_cloud[bands_no_cloud %in% names(list_metrics)]

    # get the file information from the cube
    file_info <- cube$file_info[[1]]

    # define the parameters for the output files
    params = .sits_raster_api_params_cube(cube)

    # process the bands
    file_list <- purrr::map(bands_metrics, function(bnd) {
        message(paste0("Creating metrics from band: ", bnd))
        files <- purrr::map_chr(list_metrics[[bnd]], function(metric) {
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
                output_dir, "/",
                cube$satellite, "_",
                cube$sensor, "_",
                cube$tile, "_",
                start_date, "_", end_date, "_",
                bnd, toupper(metric), ".tif"
            )

            # create a raster object
            r_obj <- suppressWarnings(
                terra::rast(
                    nrows = params$nrows,
                    ncols = params$ncols,
                    nlyrs = 1,
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
                    datatype = "INT2S"
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
                    metric = metric,
                    extent = extent,
                    impute_fn = impute_fn,
                    multicores = multicores
                )
                # rescale the data
                scl_factor <- .sits_config_scale_factors(cube$sensor, bnd)[[1]]
                mult_factor <- 1/as.numeric(scl_factor)
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

            task <- paste0("Create metrics from band: ", bnd)
            .sits_processing_task_time(task, start_task_time)
            return(filename)
        })
    })
    # report on time used for processing
    task <- paste0("All metrics were created.")
    .sits_processing_task_time(task, start_time)

    # return the file info
    files <- unlist(file_list)
    return(files)
}

#' @title Estimate the size of blocks
#' @name .metrics_blocks_estimate
#' @description This function estimates the number of blocks using the number of
#' bands and the memory size
#'
#' @param cube    A data cube
#' @param nbands  Number of bands
#' @param memsize Memory size
#'
#' @return a list with information about the beginning and end of each block
.metrics_blocks_estimate <- function(cube, nbands, memsize) {

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
    single_data_size <- nrows * ncols * nbytes * nbands * 1

    # estimated full size of the data
    full_size <- as.numeric(n_instances) * single_data_size

    # estimated size of memory required
    mem_required <- (full_size + as.numeric(sits:::.sits_mem_used())) * bloat

    # number of passes to read the full data sets
    nblocks <- ceiling(mem_required / (memsize * 1e+09))

    # Cloud processing uses the whole image
    sub_image <- sits:::.sits_raster_sub_image_default(cube)

    # calculate the blocks
    blocks <- .sits_raster_metrics_block_list(
        nblocks = nblocks,
        sub_image = sub_image
    )

    return(blocks)
}

#' @title Calculate a list of blocks to be read from disk to memory
#' @name .sits_raster_metrics_block_list
#' @keywords internal
#' @param  nblocks   number of blocks to read from each image
#' @param  sub_image area of interest in the image
#' @return a list with n (number of blocks), row (vector of starting rows),
#' nrow (vector with number of rows for each block)
#' and size (vector with size of each block)
.sits_raster_metrics_block_list <- function(nblocks, sub_image) {
    # number of rows per block
    block_rows <- ceiling(sub_image["nrows"] / nblocks)

    first_row <- unname(sub_image["first_row"])
    last_row <- first_row + unname(sub_image["nrows"]) - 1

    # initial row of each block
    row_vec <- seq.int(
        from = first_row,
        to = last_row,
        by = block_rows
    )

    # number of rows in each block
    n_rows <- length(row_vec)
    assertthat::assert_that(n_rows > 0, msg = "empty row vector")
    nrows_vec <- rep.int(block_rows, n_rows)

    # check that total number of rows is the same as the sum of all blocks
    # correct the last block for overflow
    if (sum(nrows_vec) != sub_image["nrows"]) {
        nrows_vec[length(nrows_vec)] <-
            sub_image["nrows"] - sum(nrows_vec[1:(length(nrows_vec) - 1)])
    }

    # find out the size of the block in pixels
    size_vec <- nrows_vec * sub_image["ncols"]

    # elements of the block list
    # n          number of blocks
    # row        starting row in each block (vector)
    # nrows      number of rows in each block (vector)
    # col        first col
    # ncols      number of cols in each block

    blocks <- list(
        n = length(row_vec),
        row = row_vec,
        nrows = nrows_vec,
        col = sub_image["first_col"],
        ncols = sub_image["ncols"],
        size = size_vec
    )

    message(
        "Using ", blocks$n, " blocks of size ",
        blocks$nrows[1], " x ", blocks$ncols
    )

    return(blocks)
}
