#' @title Determine the cube params to write in the metadata
#' @name .sits_raster_api_params_cube
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description    Based on the R object associated to a raster object,
#'                 determine its params
#' @param file     A valid raster image
#' @return A tibble with the cube params
.sits_raster_api_params_cube <- function(cube) {

    params <- tibble::tibble(
        nrows = cube$nrows,
        ncols = cube$ncols,
        xmin  = cube$xmin,
        xmax  = cube$xmax,
        ymin  = cube$ymin,
        ymax  = cube$ymax,
        xres  = cube$xres,
        yres  = cube$yres,
        crs   = cube$crs
    )
    return(params)
}

#' @title Read a part of a raster file and return a matrix
#' @name .sits_raster_api_read_extent
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param  r_files        Files associated to the raster object
#' @param  extent         Image extent to be read.
#' @return                Data.table of values
.sits_raster_api_read_extent <- function(r_files, extent = NULL) {

    # create terra objects
    t_obj <- terra::rast(r_files)
    # start reading
    terra::readStart(t_obj)
    if (terra::nlyr(t_obj) == 1) {
        values <- matrix(
            as.matrix(
                terra::readValues(x      = t_obj,
                                  row    = extent["row"],
                                  nrows  = extent["nrows"],
                                  col    = extent["col"],
                                  ncols  = extent["ncols"])
            ), nrow = extent["nrows"], byrow = TRUE
        )
    } else {
        values <- terra::readValues(x      = t_obj,
                                    row    = extent["row"],
                                    nrows  = extent["nrows"],
                                    col    = extent["col"],
                                    ncols  = extent["ncols"],
                                    mat = TRUE)
    }
    terra::readStop(t_obj)

    return(values)
}

#' @title Remove cloud pixels and NA values by imputation
#' @name  .sits_raster_data_na_remove
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param  values           matrix of values retrieved from a raster object
#' @param  impute_fn        imputing function to be applied to replace NA
#' @param  multicores       number of cores to process the time series.
#' @return Data.table with pre-processed values.
.sits_raster_data_na_remove <- function(values,
                                        impute_fn,
                                        multicores) {
    cld_remove_block <- function(block) {
        # interpolate NA
        block <- impute_fn(block)
    }
    # use multicores to speed up filtering
    if (multicores > 1) {
        chunks <- .sits_raster_data_split(values, multicores)
        rows <- parallel::mclapply(chunks,
                                   cld_remove_block,
                                   mc.cores = multicores
        )
        values <- do.call(rbind, rows)
    }
    else {
        values <- impute_fn(values)
    }

    return(values)
}

#' @title Filter the time series values in the case of a matrix
#' @name .sits_raster_data_filter
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description This function filters a matrix.
#'
#' @param  values         matrix of values.
#' @param  filter         Filter function to apply to matrix.
#' @param  multicores     Number of cores.
#' @return Filtered matrix.
.sits_raster_data_filter <- function(values, filter, multicores) {

    # auxiliary function to scale a block of data
    filter_matrix_block <- function(chunk) {
        filtered_block <- filter(chunk)
        return(filtered_block)
    }
    # use multicores to speed up filtering
    if (multicores > 1) {
        chunks <- .sits_raster_data_split(values, multicores)
        rows <- parallel::mclapply(chunks,
                                   filter_matrix_block,
                                   mc.cores = multicores
        )
        values <- do.call(rbind, rows)
    }
    else {
        values <- filter(values)
    }

    return(values)
}

#' @title Split a data.table or a matrix for multicore processing
#' @name .sits_raster_data_split
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description This function splits a data.table into a
#'              list of chunks for multicore processing.
#'
#' @param data             Data (data.table or matrix).
#' @param ncores           Number of cores for processing.
#' @return                 List of pairs of positions (first row, last row)
#'                         to be assigned to each core.
#'
.sits_raster_data_split <- function(data, ncores) {
    # number of rows in the data
    nrows <- nrow(data)
    # find the number of rows per core
    step <- ceiling(nrows / ncores)

    # create a vector with the initial rows per block
    blocks <- seq(from = 1, to = nrows, by = step)

    # fill the list with the initial and final row per block
    block_lst <- purrr::map2(blocks, 1:ncores, function(blk, i) {
        start <- blk
        end <- start + step - 1
        if (i == ncores) {
            end <- nrows
        }
        return(data[start:end, ])
    })
    return(block_lst)
}

#' @title Create a raster brick data cube
#' @name .sits_metrics_local_cube
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description  Builds a BRICK cube
#'
#' @param  satellite             Name of satellite
#' @param  sensor                Name of sensor
#' @param  name                  Name of the data cube.
#' @param  timeline              Vector of dates with the timeline of the bands.
#' @param  bands                 Vector of bands contained in the Raster Brick
#'                               set (in the same order as the files).
#' @param  files                 Vector with the file paths of the raster files.
#' @param cube                   A sits cube
#' @return A tibble with metadata information about a raster data set.
#'
.sits_metrics_local_cube <- function(satellite,
                                     sensor,
                                     name,
                                     timeline,
                                     bands,
                                     files,
                                     cube = NULL) {
    # obtain the parameters
    params <- .sits_raster_api_params_file(files[1])
    assertthat::assert_that(nrow(params) > 0,
                            msg = paste(".sits_raster_brick_cube: error in",
                            "retrieving raster params")
    )

    # bands in SITS are uppercase
    bands <- toupper(bands)
    # get scale factors

    times_brick <- rep(timeline[1], time = length(files))

    # get tile name
    tile <- NA
    if (!is.null(cube))
        tile <- cube$tile

    # get the file information
    file_info <- .sits_raster_api_file_info(bands, times_brick, files)

    # creating a local cube
    cube <- .sits_cube_create(
        name = name,
        source = "LOCAL",
        satellite = satellite,
        sensor = sensor,
        tile = tile,
        bands = bands,
        nrows = params$nrows,
        ncols = params$ncols,
        xmin = params$xmin,
        xmax = params$xmax,
        ymin = params$ymin,
        ymax = params$ymax,
        xres = params$xres,
        yres = params$yres,
        crs = params$crs,
        file_info = file_info
    )
    return(cube)
}

#' @title Create a tibble with file information to include in the cube
#' @name  .sits_raster_api_file_info
#' @keywords internal
#'
#' @param  bands    List of spectral bands
#' @param  timeline Cube timeline
#' @param  files    List of files associated to the
.sits_raster_api_file_info <- function(bands, timeline, files) {
    # create a tibble to store the file info
    # iterate through the list of bands and files
    assertthat::assert_that(length(bands) == length(timeline) &
                                length(files) == length(timeline),
                            msg = paste(".sits_raster_api_file_info: unmatched",
                            "bands, files and timeline"))

    # create the file info
    file_info_lst <- purrr::pmap(
        list(bands, timeline, files),
        function(b, t, f) {
            fil <- tibble::tibble(
                band = b,
                date = lubridate::as_date(t),
                path = f
            )
            return(fil)
        })
    # join the list into a tibble
    file_info <- dplyr::bind_rows(file_info_lst)

    return(file_info)
}
