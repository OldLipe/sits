sits_bocd <- function(samples = NULL,
                      start_date = NULL,
                      end_date = NULL) {
    # Training function
    train_fun <- function(samples) {
        # Create a stats tibble
        detect_change_fun <- function(values, ...) {
            dots <- list(...)
            # Extract tile
            tile <- dots[["tile"]]
            bbox <- dots[["bbox"]]
            block <- dots[["block"]]

            # Get the number of dates in the timeline
            tile_tl <- .tile_timeline(tile)
            n_times <- length(tile_tl)

            change_values <- matrix(nrow = nrow(values), ncol = 1)
            for (i in seq_len(nrow(values))) {
                cp_pts <- ocp::onlineCPD(t(values[i, ]), multivariate = TRUE)
                v <- NULL
                if (.has(cp_pts$changepoint_list$maxCPs[[1]])) {
                    v <- min(cp_pts$changepoint_list$maxCPs[[1]][-1])
                }
                change_values[i,] <- v
            }

            # Polygonize values
            values <- .dc_as_polygon(
                values = v,
                block = block,
                bbox = bbox
            )
            # Get date that corresponds to the index value
            values[["date"]] <- tile_tl[unlist(values[["date"]])]
            values
        }
        # Set model class
        predict_fun <- .set_class(
            detect_change_fun, "bocd_model", "sits_model",
            class(detect_change_fun)
        )
        return(predict_fun)
    }
    # If samples is informed, train a model and return a predict function
    # Otherwise give back a train function to train model further
    result <- .factory_function(samples, train_fun)
    return(result)
}
