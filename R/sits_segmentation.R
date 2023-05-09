#' @title Segment an image using supercells
#'
#' @name sits_supercells
#'
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#'
#' @description
#' Apply a segmentation on a data cube based on the "supercells" package. This is an
#' adaptation and extention to remote sensing data of the SLIC superpixels
#' algorithm proposed by Achanta et al. (2012). See references for more details.
#'
#' @param cube          Regular data cube
#' @param tiles         Tiles to be segmented
#' @param bands         Bands to include in the segmentation
#' @param date          Date to select the image to be segmented
#' @param step          Distance (in number of cells) between initial supercells' centers.
#' @param compactness   A compactness value. Larger values cause clusters to
#'                      be more compact/even (square).
#' @param iter          Number of iterations to create the output.
#' @param minarea       Specifies the minimal size of a supercell (in cells).
#' @param chunks        Should the input (x) be split into chunks before deriving
#'                      supercells? Either FALSE,
#'                      TRUE (default - only large input objects are split),
#'                      or a numeric value (representing the side length of the chunk
#'                      in the number of cells).
#' @param future        Should the future package be used for parallelization
#'                      of the calculations?
#' @param multicores    Number of cores for parallel processing
#'
#' @references
#'                      Achanta, Radhakrishna, Appu Shaji, Kevin Smith, Aurelien Lucchi,
#'                      Pascal Fua, and Sabine Süsstrunk. 2012. “SLIC Superpixels Compared
#'                      to State-of-the-Art Superpixel Methods.” IEEE Transactions on
#'                      Pattern Analysis and Machine Intelligence 34 (11): 2274–82.
#'
#'                      Nowosad, Jakub, and Tomasz F. Stepinski. 2022. “Extended SLIC
#'                      Superpixels Algorithm for Applications to Non-Imagery Geospatial
#'                      Rasters.” International Journal of Applied Earth Observation
#'                      and Geoinformation 112 (August): 102935.
#'
#' @examples
#' # example code
#' if (sits_run_examples()) {
#'     # Example of classification of a data cube
#'     # create a data cube from local files
#'     data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
#'     cube <- sits_cube(
#'         source = "BDC",
#'         collection = "MOD13Q1-6",
#'         data_dir = data_dir
#'     )
#'     # segment the image
#'     segments <- sits_supercells(
#'         cube = cube,
#'         tile = "012010",
#'         bands = "NDVI",
#'         date = sits_timeline(cube)[1],
#'         step = 10
#'     )
#' }
#' @export
sits_supercells <- function(
        cube,
        tiles = cube[1, "tile"],
        bands,
        date,
        step = 50,
        compactness = 1,
        iter = 10,
        minarea = 30,
        chunks = TRUE,
        future = FALSE,
        multicores = 4
){
    # check package availability
    .check_require_packages(c("supercells", "future"))
    # check input parameters
    # cube is regular
    .check_is_regular(cube)
    # tile belongs to the cube
    .check_chr_within(tiles, .cube_tiles(cube),
                      msg = "tiles not available in the cube")
    # bands are OK
    .check_chr_within(bands, .cube_bands(cube),
                      msg = "bands not available in the cube")
    # date is ok
    .check_that(as.Date(date) %in% .cube_timeline(cube)[[1]],
                      msg = "date not available in the cube")
    # step is OK
    .check_int_parameter(step, min = 1, max = 500)
    # compactness is OK
    .check_int_parameter(compactness, min = 1, max = 50)
    # iter is OK
    .check_int_parameter(iter, min = 10, max = 100)
    # minarea is OK
    .check_int_parameter(minarea, min = 10, max = 100)
    # chunks
    .check_lgl_parameter(chunks)
    # future
    .check_lgl_parameter(future)
    # multicores
    .check_int_parameter(multicores, min = 1, max = 1000)

    # obtain the image files to perform the segmentation
    # get the tile
    tile_rows <- .cube_filter_tiles(cube, tiles)

    cells_tile <- slider::slide(tile_rows, function(row){
        # filter tile by band and date
        row <- row %>%
            .tile_filter_bands(bands) %>%
            .tile_filter_dates(date)
        # get the paths of required image files
        files <- purrr::map_chr(bands, function(band){
            file <- .tile_path(row, band, date)
            return(file)
        })
        # obtain the SpatRaster (terra) object
        rast <- terra::rast(files)
        # prepare the future plan
        if (future)
            future::plan(strategy = "future::multisession",
                         workers = multicores)

        # segment the terra object
        cells_sf <- supercells::supercells(
            x = rast,
            compactness = compactness,
            step = step,
            iter = iter,
            chunks = chunks,
            future = future
        )
        class(cells_sf) <- c("supercells", class(cells_sf))
        return(cells_sf)
    })
    # returns a named list
    names(cells_tile) <- tiles
    class(cells_tile) <- c("segments", class(cells_tile))
    return(cells_tile)
}

st_segmentation <- function(samples = NULL, index = NULL, gridx = 12,gridy=12, ...) {

    # function that returns `randomForest::randomForest` model
    train_fun <- function(samples) {

        train_data <- .sits_distances(samples)
        if (!is.null(index)) {
            train_data <- train_data[index]
        }
        mygrid = kohonen::somgrid(gridx, gridy, "rectangular")

        sits:::.check_require_packages("kohonen")

        trainingdata <- list(band = as.matrix(train_data[,3:ncol(train_data)]))

        model <- kohonen::som(trainingdata$band, grid = mygrid,...)

        # construct model predict enclosure function and returns
        predict_fun <- function(values) {

            sits:::.check_require_packages("kohonen")


            prediction <- predict(model,
                                  ,newdata = values)
            browser()
            cluster <- as.factor(prediction$unit.classif)
            dt <- data.table::data.table(cluster)
            dt[, id := .I]
            df <- data.table::dcast.data.table(dt, id ~ cluster, fun.aggregate = length)
            df$id <- NULL
            matrix(nrow = nrow(df), ncol = ncol(df))
            df <- as.matrix(df)
            df <- matrix(as.double(df),nrow = nrow(df), ncol = ncol(df))
            colnames(df) <- levels(cluster)
            return(df)

        }
        # Set model class
        predict_fun <- sits:::.set_class(
            predict_fun, "st_segmentation", "sits_model", class(predict_fun)
        )
        return(predict_fun)
    }
    result <- sits:::.sits_factory_function(samples, train_fun)
    return(result)
}

.sits_distances <- function(data) {
    # get bands order
    bands <- names(data$time_series[[1]][-1])

    # create a tibble with the time series transposed from columns to rows
    # and create sample_id and reference columns as the first two
    # columns for training
    distances_tbl <- data %>%
        dplyr::mutate(
            sample_id = seq_len(nrow(data))
        ) %>%
        tidyr::unnest("time_series") %>%
        dplyr::select("sample_id", "label", dplyr::all_of(bands)) %>%
        dplyr::group_by(.data[["sample_id"]]) %>%
        dplyr::mutate(temp_index = seq_len(dplyr::n())) %>%
        dplyr::ungroup()

    if (length(bands) > 1) {
        distances_tbl <- tidyr::pivot_wider(distances_tbl,
                                            names_from = "temp_index",
                                            values_from = !!bands,
                                            names_sep = ""
        )
    } else {
        distances_tbl <- tidyr::pivot_wider(distances_tbl,
                                            names_from = "temp_index",
                                            values_from = !!bands,
                                            names_prefix = bands,
                                            names_sep = ""
        )
    }

    distances_tbl <- data.table::data.table(distances_tbl)
    # postcondition
    sits:::.check_na(distances_tbl)

    return(distances_tbl)
}
