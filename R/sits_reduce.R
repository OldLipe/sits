sits_reduce <- function(data, ...) {
    .check_na(data)
    .check_null(data)
    UseMethod("sits_reduce", data)
}

sits_reduce.sits <- function(data, ...) {
    data <- .check_samples(data)
    .check_set_caller("sits_apply.sits")

    .reduce(data, col = "time_series", fn = dplyr::mutate, ...)
}
