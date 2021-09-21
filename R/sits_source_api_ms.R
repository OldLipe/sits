#' @keywords internal
#' @export
.source_access_test.ms_cube <- function(source, collection, ..., bands,
                                          dry_run = TRUE) {
    # require package
    if (!requireNamespace("rstac", quietly = TRUE)) {
        stop(paste("Please install package rstac from CRAN:",
                   "install.packages('rstac')"), call. = FALSE
        )
    }

    items_query <- .stac_items_query(source = source,
                                     collection = collection,
                                     limit = 1, ...)

    # assert that service is online
    tryCatch({
        items <- rstac::post_request(items_query)
    }, error = function(e) {
        stop(paste(".source_access_test.stac_cube: service is unreachable\n",
                   e$message), call. = FALSE)
    })

    # assign items
    items <- rstac::items_sign(items = items,
                                    sign_fn = rstac::sign_planetary_computer())

    items <- .source_items_bands_select(source = source,
                                        collection = collection,
                                        items = items,
                                        bands = bands[[1]], ...)

    href <- .source_item_get_hrefs(source = source,
                                   item = items$feature[[1]], ...,
                                   collection = collection)

    # assert that token and/or href is valid
    if (dry_run)
        tryCatch({
            .raster_open_rast(href)
        }, error = function(e) {
            stop(paste(".source_access_test.stac_cube: cannot open url\n",
                       href, "\n", e$message), call. = FALSE)
        })


    return(invisible(NULL))
}

#' @keywords internal
#' @export
.source_item_get_date.ms_cube <- function(source,
                                          item, ...,
                                          collection = NULL) {
    item[[c("properties", "datetime")]]
}

#' @keywords internal
#' @export
.source_item_get_hrefs.ms_cube <- function(source,
                                           item, ...,
                                           collection = NULL) {

    href <- unname(purrr::map_chr(item[["assets"]], `[[`, "href"))

    # add gdal vsi in href urls
    return(.stac_add_gdal_vsi(href))
}

#' @keywords internal
#' @export
.source_item_get_bands.ms_cube <- function(source,
                                           item, ...,
                                           collection = NULL) {
    names(item[["assets"]])
}

#' @keywords internal
#' @export
.source_item_get_resolutions.ms_cube <- function(source,
                                                 item, ...,
                                                 collection = NULL) {

    # TODO: change this value
    return(10)
}

#' @keywords internal
#' @export
.source_items_new.ms_cube <- function(source,
                                      collection, ...,
                                      stac_query,
                                      tiles = NULL) {

    # set caller to show in errors
    .check_set_caller(".source_items_new.ms_cube")

    # if specified, a filter per tile is added to the query
    if (!is.null(tiles))
        stac_query <- rstac::ext_query(q = stac_query,
                                       "s2:mgrs_tile" == tiles)

    # making the request
    items_info <- rstac::post_request(q = stac_query, ...)

    # assign items
    items_info <- rstac::items_sign(items = items_info,
                                    sign_fn = rstac::sign_planetary_computer())

    # check if matched items
    .check_that(
        x = rstac::items_matched(items_info) > 0,
        msg = "no items matched the query criteria."
    )

    # fetching all the metadata
    items_info <- rstac::items_fetch(items = items_info)

    return(items_info)
}

#' @keywords internal
#' @export
.source_items_tiles_group.ms_cube <- function(source,
                                              items, ...,
                                              collection = NULL) {

    rstac::items_group(items, field = c("properties", "s2:mgrs_tile"))
}

#' @keywords internal
#' @export
.source_items_tile_get_crs.ms_cube <- function(source,
                                               tile_items, ...,
                                               collection = NULL) {

    # making request to collection endpoint to get crs info
    # format collection crs
    crs <- .sits_proj_format_crs(
        tile_items[["features"]][[1]][[c("properties", "proj:epsg")]]
    )

    return(crs)
}

#' @keywords internal
#' @export
.source_items_tile_get_name.ms_cube <- function(source,
                                                tile_items, ...,
                                                collection = NULL) {

    tile_items[["features"]][[1]][[c("properties", "s2:mgrs_tile")]]
}

#' @keywords internal
#' @export
.source_items_tile_get_bbox.ms_cube <- function(source,
                                                tile_items, ...,
                                                collection = NULL) {
    # get collection crs
    crs <- .source_items_tile_get_crs(source = source,
                                      tile_items = tile_items,
                                      collection = collection)

    # get bbox by geometry attributei in tile_items
    bbox <- .stac_get_bbox(tile_items, crs)

    return(bbox)
}

#' @keywords internal
#' @export
.source_items_tile_get_size.ms_cube <- function(source,
                                                tile_items, ...,
                                                collection = NULL) {

    href <- .source_item_get_hrefs(source = source,
                                   item = tile_items[["features"]][[1]], ...,
                                   collection = collection)

    # read the first image and obtain the size parameters
    params <- .raster_params_file(href)

    size <- c(nrows = params[["nrows"]], ncols = params[["ncols"]])

    return(size)
}
