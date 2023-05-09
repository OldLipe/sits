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
