#' @title ...
#' @name sits_metrics_cube
#'
#' @description a lazy-eval metrics creation cube
#'
#' @param cube    input data cube.
#' @param metrics metrics list. melhorar a documentacao
#'
#' @return the same input data cube
#' @export
sits_metrics_cube <- function(cube, metrics) {

    .check_lst(metrics,
               min_len = 1,
               msg = "invalid 'metrics' parameters.")

    m_bands <- names(metrics)

    .check_chr_within(x = m_bands,
                      within = sits_bands(cube),
                      msg = "invalid bands in metrics list.")

    # change o sits_config on-the-fly
    .set_source_metrics(cube, metrics)

    # add attributes in cube metrics
    attr(cube, "metrics") <- metrics

    return(cube)
}

#' @title ...
#' @name .set_source_metrics
#'
#' @keywords internal
#'
#' @description change config on the fly
#'
#' @param cube    input data cube.
#' @param metrics metrics list. melhorar a documentacao
#'
#' @return an invisible null
.set_source_metrics <- function(cube, metrics) {

    # aqui eu pego a coleção
    col <- sits:::.config_get(key = c("sources", sits:::.cube_source(cube),
                                      "collections", sits:::.cube_collection(cube)))

    bands <- col$bands[names(metrics)]

    rep_bands <- rep(bands, lengths(metrics))

    metrics_vec <- unlist(
        sapply(names(metrics), function(x) {paste(x, metrics[[x]], sep = "_")})
    )

    rep_bands <- lapply(seq_along(rep_bands), function(i) {
        rep_bands[[i]]$band_name <-  metrics_vec[[i]]
        # TODO: ver como tratar os valores altos das metricas
        # polares
        rep_bands[[i]]
    })

    names(rep_bands) <- unname(metrics_vec)

    col$bands <- c(rep_bands, col$bands)
    new_source <-  sits:::.config_get(key = "sources")
    new_source[[c(sits:::.cube_source(cube),
                  "collections",
                  sits:::.cube_collection(cube))]] <- col
    sits:::.config_set_options(sources = new_source)

    return(invisible(NULL))
}

#' @title ...
#' @name sits_get_metrics
#'
#' @description materialize metrics in sits_tibble and sits_cube
#'
#' @param data documentar
#' @param metrics documentar
#' @param ... documentar
#'
#' @return documentar
#' @export
sits_get_metrics <- function(data, metrics, ...) {


    # set caller to show in errors
    .check_set_caller("sits_generate_metrics")

    # Dispatch
    UseMethod("sits_get_metrics", data)
}

#' @rdname sits_get_metrics
#' @export
sits_get_metrics.sits <- function(data, metrics, ...) {

    sits_ts <- sits::sits_select(data, names(metrics))

    metrics_list <- purrr::map(sits_bands(sits_ts), function(band) {
        # select the chosen bands for the time series
        sits_timeseries <- sits::sits_select(sits_ts, band)

        ts <- purrr::map(sits_timeseries$time_series, function(ts) {
            as.matrix(t(unlist(ts[, band])))
        })

        ts <- do.call(rbind, ts)

        .extract_metric(ts, metrics[[band]])
    })

    names(metrics_list) <- sits::sits_bands(sits_ts)
    add_prefix_column <- function(time_series) {

        time_series <- lapply(names(time_series), function(name) {
            tibble::as_tibble(time_series[[name]]) %>%
                dplyr::rename_with(~ toupper(paste0(name, "_", .x))) %>%
                as.matrix()
        })

        metrics_tibble <- do.call(cbind, lapply(time_series, as.data.frame)) %>%
            tibble::as_tibble()

        metrics_tibble$Index <- max(sits::sits_timeline(data))

        return(metrics_tibble)
    }

    metrics_tibble <- add_prefix_column(metrics_list)

    sits_ts$time_series <- purrr::map(1:nrow(metrics_tibble), function(x) {
        no_index <- colnames(metrics_tibble[x, ])
        no_index <- no_index[!no_index %in% "Index"]
        metrics_tibble[x, c("Index", no_index)]
    })

    return(sits_ts)
}

#' @title ...
#' @name sits_get_metrics
#' @keywords internal
#'
#' @param ts documentar
#' @param metrics documentar
#'
#' @return documentar
.extract_metric <- function(ts, metrics) {

    # TODO: create an eval character
    metrics_list <- list(
        if ("max" %in% metrics) sitsfeats::max_ts(ts),
        if ("min" %in% metrics) sitsfeats::min_ts(ts),
        if ("mean" %in% metrics) sitsfeats::mean_ts(ts),
        if ("median" %in% metrics) sitsfeats::median_ts(ts),
        if ("sum" %in% metrics) sitsfeats::sum_ts(ts),
        if ("std" %in% metrics) sitsfeats::std_ts(ts),
        if ("skew" %in% metrics) sitsfeats::skew_ts(ts),
        if ("kurt" %in% metrics) sitsfeats::kurt_ts(ts),
        if ("amplitude" %in% metrics) sitsfeats::amplitude_ts(ts),
        if ("fslope" %in% metrics) sitsfeats::fslope_ts(ts),
        if ("abs_sum" %in% metrics) sitsfeats::abs_sum_ts(ts),
        if ("amd" %in% metrics) sitsfeats::amd_ts(ts),
        if ("mse" %in% metrics) sitsfeats::mse_ts(ts),
        if ("fqr" %in% metrics) sitsfeats::fqr_ts(ts),
        if ("sqr" %in% metrics) sitsfeats::sqr_ts(ts),
        if ("tqr" %in% metrics) sitsfeats::tqr_ts(ts),
        if ("iqr" %in% metrics) sitsfeats::iqr_ts(ts),
        if ("area_q1" %in% metrics) sitsfeats::area_q1(ts),
        if ("area_q2" %in% metrics) sitsfeats::area_q2(ts),
        if ("area_q3" %in% metrics) sitsfeats::area_q3(ts),
        if ("area_q4" %in% metrics) sitsfeats::area_q4(ts),
        if ("polar_balance" %in% metrics) sitsfeats::polar_balance(ts),
        if ("angle" %in% metrics) sitsfeats::angle(ts),
        if ("area_ts" %in% metrics) sitsfeats::area_ts(ts),
        if ("ecc_metric" %in% metrics) sitsfeats::ecc_metric(ts),
        if ("gyration_radius" %in% metrics) sitsfeats::gyration_radius(ts),
        if ("csi" %in% metrics) sitsfeats::csi(ts)
    )

    names(metrics_list) <- c("max",
                             "min",
                             "mean",
                             "median",
                             "sum",
                             "std",
                             "skew",
                             "kurt",
                             "amplitude",
                             "fslope",
                             "abs_sum",
                             "amd",
                             "mse",
                             "fqr",
                             "sqr",
                             "tqr",
                             "iqr",
                             "area_q1",
                             "area_q2",
                             "area_q3",
                             "area_q4",
                             "polar_balance",
                             "angle",
                             "area_ts",
                             "ecc_metric",
                             "gyration_radius",
                             "csi")

    metrics_list <- Filter(Negate(is.null), metrics_list)

    metrics_list
}
