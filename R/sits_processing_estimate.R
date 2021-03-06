#' @title Estimate the processing time
#' @name .sits_est_class_time
#' @keywords internal
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#'
#' @description This function normalizes image values.
#'
#' @param  start_time     initial processing time.
#' @param  n_blocks       number of blocks
#' @param  block          current block.
#' @return Expected time to complete the function.
.sits_est_class_time <- function(start_time,
                                 n_blocks,
                                 block) {
    # compute current time
    current_time <- lubridate::now()

    # compute elapsed time and estimates remaining time
    elaps_time <- lubridate::time_length(current_time - start_time,
        unit = "minute"
    )
    if (block < n_blocks) {
        message(sprintf(
            "Elapsed time %s minute(s).
             Estimated total process time %s minute(s)...",
             round(as.numeric(elaps_time), 1),
             round(as.numeric((n_blocks / block) * elaps_time), 1)
        ))
    }
    else {
        message(paste0(
            "Classification finished at ", current_time,
            ". Total elapsed time: ",
            round(as.numeric(elaps_time), 1),
            "minute(s)."
        ))
    }
    return(invisible(TRUE))
}

#' @title Estimate the processing time of a task
#' @name .sits_processing_task_time
#' @keywords internal
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#'
#' @description This function normalizes image values.
#'
#' @param  task           Description of the task
#' @param  start_time     Initial processing time.

#' @return message about processing time for a task
.sits_processing_task_time <- function(task, start_time) {
    # compute current time
    current_time <- lubridate::now()

    # compute elapsed time and estimates remaining time
    elapsed_time <- lubridate::time_length(current_time - start_time,
        unit = "minute"
    )
    message(paste0(
        "Elapsed time for task ", task, " : ",
        round(elapsed_time, digits = 2), " minutes"
    ))
    return(invisible(TRUE))
}
