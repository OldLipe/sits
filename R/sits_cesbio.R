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
            browser()
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

            # Get date that corresponds to the index value
            shadow_values <- tile_tl[shadow_values]
            # Polygonize values
            shadow_values <- .detect_change_as_polygon(
                values = shadow_values,
                block = block,
                bbox = bbox
            )

            # Get date that corresponds to the index value
            neigh_values <- tile_tl[.as_chr(neigh_values)]
            # Polygonize values
            neigh_values <- .detect_change_as_polygon(
                values = neigh_values,
                block = block,
                bbox = bbox
            )

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
