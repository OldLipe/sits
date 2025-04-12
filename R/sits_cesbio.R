sits_cesbio <- function(samples = NULL,
                        start_date = NULL,
                        end_date = NULL,
                        xa = 3,
                        xb = 5,
                        shadow_value = 0.40,
                        neigh_value = 0.50) {
    # Training function
    train_fun <- function(samples) {
        # Create a stats tibble
        # TODO: verify  if  data has only one polarization

        detect_change_fun <- function(values, ...) {
            dots <- list(...)
            # Extract tile
            tile <- dots[["tile"]]
            bbox <- dots[["bbox"]]
            block <- dots[["block"]]

            # Get the number of dates in the timeline
            tile_tl <- .tile_timeline(tile)
            n_times <- length(tile_tl)

            # Get the start and end time of the detection period
            start_detection <- xb
            end_detection <- n_times - xa
            if (.has(start_date) && .has(end_date)) {
                filt_idxs <- which(tile_tl >= start_date & tile_tl <= end_date)
                start_detection <- min(filt_idxs) - 1
                end_detection <- max(filt_idxs)
            }
            # TODO: transform to db
            # Calculate the radar change ratio
            values <- C_cesbio_calc_rcr(
                values = values,
                xa     = xa,
                xb     = xb
            )

            # Filter the shadow pixels
            shadow_values <- C_cesbio_detect_shadow(
                rcr = values,
                shadow_value = shadow_value
            )

            # Filter the neigh pixels
            neigh_values <- C_cesbio_detect_neigh(
                rcr = values,
                neigh_value = neigh_value
            )

            # Polygonize values
            browser()
            shadow_values <- .dc_as_polygon(
                values = shadow_values,
                block = block,
                bbox = bbox
            )
            # Get date that corresponds to the index value
            #shadow_values <- tile_tl[shadow_values]
            # Polygonize values
            neigh_values <- .dc_as_polygon(
                values = neigh_values,
                block = block,
                bbox = bbox
            )
            # Get date that corresponds to the index value
            #neigh_values <- tile_tl[.as_chr(neigh_values)]

        }
        # Set model class
        predict_fun <- .set_class(
            detect_change_fun, "cesbio_model", "sits_model",
            class(detect_change_fun)
        )
        return(predict_fun)
    }
    # If samples is informed, train a model and return a predict function
    # Otherwise give back a train function to train model further
    result <- .factory_function(samples, train_fun)
    return(result)

}

#' @title Detect change as a polygon
#' @name .dc_as_polygon
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @keywords internal
#' @noRd
#' @param values     Matrix of values for a raster (time series)
#' @param block      Data block that is being processed
#' @param bbox       Bounding box of the block
#' @return           Vector object with polygons
.dc_as_polygon <- function(values, block, bbox) {
    # Create a template raster
    template_raster <- .raster_new_rast(
        nrows = block[["nrows"]], ncols = block[["ncols"]],
        xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
        ymin = bbox[["ymin"]], ymax = bbox[["ymax"]],
        nlayers = 1, crs = bbox[["crs"]]
    )
    # Set values and NA value in template raster
    values <- .raster_set_values(template_raster, values)
    values <- .raster_set_na(values, 0)
    names(values) <- "date"
    # Extract polygons raster and convert to sf object
    values <- .raster_as_polygon(values)
    if (nrow(values) == 0) {
        return(values)
    }
    # Get only polygons segments
    values <- suppressWarnings(sf::st_collection_extract(values, "POLYGON"))
    # Return the segment object
    return(values)
}

.raster_as_polygon <- function(rast, ...) {
    tem_rast <- tempfile(fileext = ".tif")
    .raster_write_rast(
        r_obj = rast,
        file = tem_rast,
        overwrite = FALSE,
        data_type = "INT1U"
    )

    temp_vect <- tempfile(fileext = ".gpkg")
    gdalraster::polygonize(
        raster_file = tem_rast,
        out_dsn = temp_vect,
        out_layer = "pol",
        out_fmt = "gpkg"
    )
    sf::st_read(temp_vect)
}

