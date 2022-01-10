.bg_raster_preprocess <- function(cube,
                                  file,
                                  file_cld,
                                  block,
                                  missing_value,
                                  minimum_value,
                                  maximum_value,
                                  scale_factor,
                                  offset_value, ...,
                                  impute_fn = sits_impute_linear(),
                                  filter_fn = NULL,
                                  normalize_fn = NULL) {

    # pre-conditions
    .check_num(missing_value, len_min = 1, len_max = 1,
               msg = "invalid missing_value parameter")

    .check_num(minimum_value, len_min = 1, len_max = 1,
               msg = "invalid minimum_value parameter")

    .check_num(maximum_value, len_min = 1, len_max = 1,
               msg = "invalid maximum_value parameter")

    .check_num(scale_factor, len_min = 1, len_max = 1,
               msg = "invalid scale_factor parameter")

    .check_num(offset_value, len_min = 1, len_max = 1,
               msg = "invalid offset_value parameter")

    # does the cube have a cloud band?
    cld_band <- .source_cloud()
    if (!is.null(file_cld)) {

        cld_index <- .source_cloud_interp_values(
            source = .cube_source(cube = cube),
            collection = .cube_collection(cube = cube)
        )

        clouds <- .raster_read_stack(files = file_cld, block = block)

        # get information about cloud bitmask
        if (.source_cloud_bit_mask(
            source = .cube_source(cube = cube),
            collection = .cube_collection(cube = cube))) {

            clouds <- as.matrix(clouds)
            cld_rows <- nrow(clouds)
            clouds <- matrix(bitwAnd(clouds, sum(2 ^ cld_index)),
                             nrow = cld_rows) > 0
        } else {

            clouds <- clouds %in% cld_index
        }

    } else {
        clouds <- NULL
    }

    # read the values
    values <- .raster_read_stack(file, block = block)

    # correct NA, minimum, maximum, and missing values
    values[values == missing_value] <- NA
    values[values < minimum_value] <- NA
    values[values > maximum_value] <- NA

    # change the points under clouds to NA
    if (!purrr::is_null(clouds)) {
        values[clouds] <- NA
    }

    # impute NA pixels
    if (!is.null(impute_fn) && any(is.na(values))) {

        .check_that(inherits(impute_fn, "function"))

        values <- impute_fn(values)
    }

    return(values)
}

.bg_get_block_size <- function(cube, roi, sub_image, memsize, multicores) {

    sub_image <-  .sits_raster_sub_image(cube = cube, roi = roi)

    blocks_estimate <- function(tile,
                                sub_image,
                                memsize,
                                multicores) {

        # total number of instances in the time
        ninstances <- length(sits_timeline(tile))

        # get the number of bands
        nbands <- length(.cube_bands(tile, add_cloud = TRUE))

        # number of instances per classification interval
        #ninterval <- length(.cube_timeline(tile))

        # number of bytes per pixel
        nbytes <- 8

        # estimated processing bloat
        proc_bloat <- as.numeric(.config_processing_bloat())
        if (proc_bloat == 0) proc_bloat <- multicores

        # number of rows and cols
        nrows <- sub_image[["nrows"]]
        ncols <- sub_image[["ncols"]]

        # single instance size
        single_data_size <- nrows * ncols * nbytes

        # total size including all bands
        nbands_data_size <- single_data_size * nbands
        nbands_data_size_no_cloud <- single_data_size * (nbands - 1)

        # estimated size of the data for classification
        input_data_size <- as.numeric(ninstances) * nbands_data_size
        #output_data_size <- as.numeric(ninterval) * nbands_data_size_no_cloud
        #class_data_size <- (input_data_size + output_data_size) * proc_bloat
        class_data_size <- (input_data_size) * proc_bloat

        # number of passes to read the full data sets
        nblocks <- ceiling(class_data_size * 1e-09 / memsize * multicores)

        return(nblocks)
    }

    # get the number of blocks
    nblocks <- blocks_estimate(
        tile = cube,
        sub_image = sub_image,
        memsize = memsize,
        multicores = multicores
    )

    blocks <- .sits_raster_block_list(
        nblocks = nblocks,
        sub_image = sub_image
    )

    return(blocks)
}

.preprocess_images <- function(cube,
                               roi,
                               impute_fn = sits_impute_linear(),
                               filter_fn = NULL,
                               normalize_fn = NULL,
                               multicores,
                               memsize,
                               output_dir,
                               bands = NULL) {

    stopifnot("CLOUD" %in% sits_bands(cube))

    # 0 - criar um subimage
    sub_image <- .sits_raster_sub_image(cube = cube, roi = roi)

    # 1 - definir a quantidade de chunks
    blocks <- .bg_get_block_size(
        cube = cube,
        roi = roi,
        sub_image = sub_image,
        memsize = memsize,
        multicores = multicores
    )

    stopifnot(all(bands %in% sits_bands(cube)))

    # 2 - um loop sobre as bandas
    if (is.null(bands))
        bands <- .cube_bands(cube, add_cloud = FALSE)

    for (band_cube in bands) {

        # get the missing values, minimum values and scale factors
        missing_value <- .cube_band_missing_value(cube = cube, band = band_cube)
        minimum_value <- .cube_band_minimum_value(cube = cube, band = band_cube)
        maximum_value <- .cube_band_maximum_value(cube = cube, band = band_cube)
        scale_factor <- .cube_band_scale_factor(cube, band = band_cube)
        offset_value <- .cube_band_offset_value(cube = cube, band = band_cube)

        # selecionar as bandas do cubo
        file_info_band <- .cube_file_info(cube)

        band_files <- dplyr::filter(file_info_band, band == band_cube)[["path"]]
        cld_files <- dplyr::filter(file_info_band, band == "CLOUD")[["path"]]

        # 3 - um loop sobre os blocos
        .sits_parallel_start(workers = multicores, log = FALSE)
        on.exit(.sits_parallel_stop(), add = TRUE)

        # read the blocks and compute the probabilities
        filenames <- .sits_parallel_map(blocks, function(b) {

            # read the data
            distances <- .bg_raster_preprocess(
                cube       = cube,
                file       = band_files,
                file_cld   = cld_files,
                block      = b,
                missing_value = missing_value,
                minimum_value = minimum_value,
                maximum_value = maximum_value,
                scale_factor  = scale_factor,
                offset_value = offset_value,
                impute_fn = impute_fn,
                filter_fn = filter_fn,
                normalize_fn = normalize_fn
            )

            for (i in seq_along(sits_timeline(cube))) {

                # define the file name of the raster file to be written
                filename_block <- paste0(
                    output_dir,
                    "block_", band_cube, "_", sits_timeline(cube)[[i]], "_",
                    b[["first_row"]], "_", b[["nrows"]], ".tif"
                )

                # compute block spatial parameters
                params <- .cube_params_block(cube, b)

                # create a new raster
                r_obj <- .raster_new_rast(
                    nrows   = params$nrows,
                    ncols   = params$ncols,
                    xmin    = params$xmin,
                    xmax    = params$xmax,
                    ymin    = params$ymin,
                    ymax    = params$ymax,
                    nlayers = 1,
                    crs     = params$crs
                )

                # copy values
                r_obj <- .raster_set_values(r_obj  = r_obj,
                                            values = distances[, i])

                # write the probabilities to a raster file
                .raster_write_rast(
                    r_obj        = r_obj,
                    file         = filename_block,
                    format       = "GTiff",
                    data_type    = "INT2S",
                    gdal_options = .config_gtiff_default_options(),
                    overwrite    = TRUE
                )
            }
            return(NULL)
        }, progress = TRUE)

        for (i in seq_along(sits_timeline(cube))) {

            files <- list.files(path = output_dir,
                                pattern = paste0("_", band_cube,
                                                 "_", sits_timeline(cube)[[i]]),
                                full.names = TRUE)

            pattern_name <- paste(cube$satellite[[1]],
                                  cube$sensor[[1]],
                                  band_cube, sits_timeline(cube)[[i]], ".tif",
                                  sep = "_")

            merge_files(filenames = files,
                        pattern_name = pattern_name,
                        output_dir = output_dir)
        }
        gc()
    }
}


merge_files <- function(filenames, pattern_name, output_dir) {

    # join predictions
    .raster_merge(
        in_files = filenames,
        out_file = paste0(output_dir, pattern_name),
        format = "GTiff",
        gdal_datatype = "Int16",
        gdal_options = .config_gtiff_default_options(),
        overwrite = TRUE
    )
}
