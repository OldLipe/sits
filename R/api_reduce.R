.reduce <- function(data, col, fn, ...) {
    # Pre-condition
    .check_chr_within(
        col,
        within = names(data),
        msg = "invalid column name"
    )
    # Select data do unpack
    x <- data[col]
    # Prepare to unpack
    x[["#.."]] <- seq_len(nrow(data))
    # unpack
    x <- tidyr::unnest(x, cols = dplyr::all_of(col))
    x <- dplyr::group_by(x, .data[["#.."]], .drop = FALSE)
    # apply user function
    x <- fn(x, ...)
    # pack
    x <- dplyr::ungroup(x)
    x <- tidyr::nest(x, `..unnest_col` = -"#..")
    # remove garbage
    x[["#.."]] <- NULL
    names(x) <- col
    # prepare result
    data[[col]] <- x[[col]]
    return(data)
}
