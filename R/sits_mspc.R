#' @title Get information from items
#' @name .sits_mspc_items
#' @keywords internal
#'
#' @param url        a \code{character} representing a URL for the MSPC catalog.
#' @param collection a \code{character} with the collection to be searched.
#' @param tiles      a \code{character} with the names of the tiles.
#' @param roi        selects images (tiles) that intersect according to the
#'  region of interest provided. Expressed either as an \code{sfc} or \code{sf}
#'  object from sf package, a \code{character} with GeoJSON following the rules
#'  from RFC 7946, or a \code{vector} with bounding box named XY values in
#'  WGS 84 ("xmin", "xmax", "ymin", "ymax").
#' @param start_date a \code{character} corresponds to the initial date
#'                   when the cube will be created.
#' @param end_date   a \code{character} corresponds to the final date when the
#'                   cube will be created.
#' @param bands      a \code{character} vector with the bands name.
#' @param ...        other parameters to be passed for specific types.
#'
#' @return           a \code{STACItemCollection} object representing the search
#'                   by rstac.
.sits_mspc_items <- function(url, collection, tiles, roi,
                             start_date, end_date, bands, ...){

    # obtain the datetime parameter for STAC like parameter
    datetime <- .sits_stac_datetime(start_date, end_date)

    # obtain the bbox and intersects parameters
    if (!is.null(roi)) {
        roi <- .sits_stac_roi(roi)
    } else {
        roi[c("bbox", "intersects")] <- list(NULL, NULL)
    }

    # get the limit items to be returned in each page
    limit_items <- .sits_config_rstac_limit()

    # creating a rstac object
    rstac_query <- rstac::stac(url) %>%
        rstac::stac_search(collection = collection,
                           bbox       = roi$bbox,
                           intersects = roi$intersects,
                           datetime   = datetime,
                           limit      = limit_items)

    # making the request and filtering images by interval
    items_info <- rstac_query %>%
        rstac::post_request()

    # if specified, a filter per tile
    if (!is.null(tiles)) {
        if (collection == "landsat-8-c2-l2") {
            sep_tile <- .sits_usgs_format_tiles(tiles)
            items_info <- .sits_mspc_search_tiles_landsat(items_info, sep_tile)
        } else
            items_info <- .sits_mspc_search_tiles_sentinel(items_info, sep_tile)
    }

    # checks if the collection returned any items
    assertthat::assert_that(
        !(rstac::items_length(items_info) == 0),
        msg = ".sits_mspc_items: the provided search returned zero items."
    )

    # progress bar status
    pgr_fetch  <- FALSE

    # if more than 1000 items are found the progress bar is displayed
    if (rstac::items_matched(items_info) > 1000)
        pgr_fetch <- TRUE

    # fetching all the metadata and updating to upper case instruments
    items_info <- rstac::items_fetch(items_info,
                                     progress = pgr_fetch)

    # getting sensor name
    platform <- toupper(items_info$features[[1]]$properties$platform)
    sensor <- .sits_config_sensors(platform)

    # getting bands name
    items_info <- .sits_stac_bands(items = items_info,
                                   bands = bands,
                                   source = "MSPC",
                                   sensor = sensor)

    # store tile info in items object
    items_info$features <- purrr::map(items_info$features, function(features) {
        if (collection == "landsat-8-c2-l2")
            features$properties$tile <- paste0(
                features$properties[["landsat:wrs_path"]],
                features$properties[["landsat:wrs_row"]])
        else
            features$properties$tile <- paste0(
                features$properties[["s2:mgrs_tile"]])
        features
    })

    return(items_info)

}

#' @title Search items tiles in landsat collection
#' @name .sits_mspc_search_tiles_landsat
#' @keywords internal
#'
#' @param items a \code{STACItemCollection} object returned by rstac package.
#'  grouped.
#' @param tiles a \code{character} with the names of the tiles.
#'
#' @return      a \code{STACItemCollection} object representing the search
#'              by rstac.
.sits_mspc_search_tiles_landsat <- function(items, tiles) {

    #' TODO: improve this function

    # checks if the supplied tiles are in the searched items
    index_features <- purrr::map_lgl(items$features, function(feature) {
        wrs_path <- feature[["properties"]][["landsat:wrs_path"]]
        wrs_row <- feature[["properties"]][["landsat:wrs_row"]]

        if (tiles$wrs_path %in% wrs_path && tiles$wrs_row %in% wrs_row)
            return(TRUE)
        return(FALSE)
    })

    # selects the tiles found in the search
    items$features <- items$features[index_features]

    # checks if the search return zero items
    assertthat::assert_that(
        rstac::items_length(items) != 0,
        msg = paste(".sits_mspc_search_tiles: the supplied tile(s) were",
                    "not found.")
    )

    return(items)
}

#' @title Search items tiles in sentinel collection
#' @name .sits_mspc_search_tiles_sentinel
#' @keywords internal
#'
#' @param items a \code{STACItemCollection} object returned by rstac package.
#'  grouped.
#' @param tiles a \code{character} with the names of the tiles.
#'
#' @return      a \code{STACItemCollection} object representing the search
#'              by rstac.
.sits_mspc_search_tiles_sentinel <- function(items, tiles) {

    #' TODO: improve this function

    # checks if the supplied tiles are in the searched items
    index_features <- purrr::map_lgl(items$features, function(feature) {
        mgrs_tile <- feature[["properties"]][["s2:mgrs_tile"]]

        if (tiles %in% mgrs_tile)
            return(TRUE)
        return(FALSE)
    })

    # selects the tiles found in the search
    items$features <- items$features[index_features]

    # checks if the search return zero items
    assertthat::assert_that(
        rstac::items_length(items) != 0,
        msg = paste(".sits_mspc_search_tiles: the supplied tile(s) were",
                    "not found.")
    )

    return(items)
}

#' @title Get the STAC information corresponding to a tile.
#' @name .sits_mspc_tile_cube
#' @keywords internal
#'
#' @param name        Name of output data cube.
#' @param items       \code{STACItemCollection} object returned by rstac.
#' @param collection  Image collection in USGS
#' @param file_info   file information with date/time.
#'
#' @return          a \code{tibble} with metadata information about a
#'                  raster data set.
.sits_mspc_tile_cube <- function(name, items, collection, file_info){
    # store items properties attributes
    item_prop <- items$features[[1]]$properties

    # get stac bands
    bands <- items[["bands"]]

    params <- .sits_raster_api_params_file(file_info$path[1])

    # format stac crs
    item_prop[["proj:epsg"]] <- .sits_format_crs(item_prop[["proj:epsg"]])

    # obtain bbox extent
    bbox <- .sits_stac_get_bbox(items, item_prop[["proj:epsg"]])

    # get resolution
    # TODO: this hardcode will be removed
    if (collection == "landsat-8-c2-l2") {
        item_prop[["eo:gsd"]] <- 30
        sensor <- "OLI"
        satellite <- "LANDSAT-8"
    }
    else {
        # TODO: just test, this will be removed
        sensor <- "MSI"
        satellite <- "SENTINEL-2"
        item_prop[["eo:gsd"]] <- 10
    }

    res <- list(xres = item_prop[["eo:gsd"]], yres = item_prop[["eo:gsd"]])

    # add resolution to file_info
    file_info <- dplyr::mutate(file_info,
                               res = as.numeric(res[["xres"]]),
                               .before = path)

    # create a tibble to store the metadata
    tile <- .sits_cube_create(
        name       = name,
        source     = "MSPC",
        collection = collection,
        satellite  = satellite,
        sensor     = sensor,
        tile       = item_prop[["tile"]],
        bands      = bands,
        labels     = labels,
        nrows      = params$nrows,
        ncols      = params$ncols,
        xmin       = params$xmin,
        xmax       = params$xmax,
        ymin       = params$ymin,
        ymax       = params$ymax,
        xres       = res[["xres"]],
        yres       = res[["yres"]],
        crs        = item_prop[["proj:epsg"]],
        file_info  = file_info)

    tile <- .sits_config_bands_stac_write(tile)

    return(tile)

}
