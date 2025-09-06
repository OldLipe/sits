.elects_classification_head <- torch::nn_module(
    classname = ".elects_classification_head",
    initialize = function(hidden_dims, nclasses) {
        self$projection <- torch::nn_sequential(
            torch::nn_linear(hidden_dims, nclasses, bias = TRUE),
            torch::nn_log_softmax(3)
        )
    },
    forward = function(x) {
        self$projection(x)
    }
)

.elects_decision_head <- torch::nn_module(
    classname = ".elects_decision_head",
    initialize = function(hidden_dims) {
        self$projection = torch::nn_sequential(
            torch::nn_linear(hidden_dims, 1, bias = TRUE),
            torch::nn_sigmoid()
        )

        # initialize bias to predict late in first epochs
        torch::nn_init_normal_(self$projection[[1]]$bias, mean = 2e-1, std = 1e-1)
    },
    forward = function(x) {
        self$projection(x)$squeeze(3)
    }
)

.elects_early_detection <- function(res) {
    logprobabilities <- res[[1]]
    deltas <- res[[2]]
    batchsize <- logprobabilities$shape[1]
    sequencelength <- logprobabilities$shape[2]
    nclasses <- logprobabilities$shape[3]


    sample_stop_decision <- function(delta) {
        dist <- torch::torch_stack(c(1 - delta, delta), dim = 2)
        torch::distr_categorical(dist)$sample()$bool()
    }

    stop_list <- list()

    for (i in seq_len(sequencelength)) {
        if (i < sequencelength - 1) {
            stop_now = sample_stop_decision(deltas[, t])
            stop_list[[i]] <- stop_now
        } else {
            # make sure to stop last
            last_stop <- torch::ones(stop_now$shape)$bool()
            stop_list[[i]] <- last_stop
        }
    }

    # stack over the time dimension (multiple stops possible)
    stopped <- torch::torch_stack(stop_list, dim = 2)$bool()

    # is only true if stopped for the first time
    first_stops <- (stopped$cumsum(2) == 1) && stopped

    # time of stopping
    t_stop <-  first_stops$long()$argmax(2)

    # all predictions
    predictions <- logprobabilities$argmax(-1)

    # predictions at time of stopping
    predictions_at_t_stop <- torch::torch_masked_select(predictions, first_stops)

    return(list(logprobabilities, deltas, predictions_at_t_stop, t_stop))
}

.elects_calculate_probability_making_decision <- function(deltas) {
    batchsize <- deltas$shape[1]
    sequencelength <- deltas$shape[2]

    pts <- list()

    initial_budget <- torch::torch_ones(batchsize, device = "cuda")

    budget <- list(initial_budget)
    for (i in seq_len(sequencelength)) {
        pt <- deltas[, i] * budget[[length(budget)]]
        budget[[i]] <- budget[[length(budget)]] - pt
        pts[[i]] <- pt
    }

    # last time
    pt <- budget[[length(budget)]]
    pts[[length(pts)]] <- pt

    return(torch::torch_stack(pts, dim = -1))
}

.elects_probability_correct_class <- function(logprobabilities, targets) {
    batchsize <- logprobabilities$shape[1]
    seqquencelength <- logprobabilities$shape[2]
    nclasses <- logprobabilities$shape[3]
    eye <- torch::torch_eye(nclasses, dtype = "int", device = "cuda")

    targets_one_hot <- eye[targets]

    # implement the y*\hat{y} part of the loss function
    y_haty <- torch::torch_masked_select(
        logprobabilities, targets_one_hot$bool()$view(c(batchsize, 1, nclasses))
    )
    return(y_haty$view(c(batchsize, seqquencelength))$exp())
}

.acc_metric <- luz::luz_metric(
    abbrev = "Acc",
    initialize = function() {
        self$correct <- 0
        self$total <- 0
    },
    update = function(pred, target) {
        #pred <- torch::torch_argmax(preds, dim = 2)
        self$correct <- self$correct + (pred == target)$
            to(dtype = torch::torch_float())$
            sum()$
            item()
        self$total <- self$total + pred$numel()
    },
    compute = function() {
        self$correct/self$total
    }
)

.early_reward_loss <- torch::nn_module(
    classname = "early_reward_loss",
    initialize = function(alpha = 0.5, epsilon = 5, weight = NULL) {
        self$negative_log_likelihood <- torch::nn_nll_loss(
            weight = weight, reduction = "none"
        )
        self$alpha <- alpha
        self$epsilon <- epsilon
    },
    forward = function(log_class_probabilities,
                       probability_stopping,
                       y_true,
                       return_stats = FALSE) {
        N <- log_class_probabilities$shape[1]
        T <- log_class_probabilities$shape[2]
        C <- log_class_probabilities$shape[3]

        # equation 3
        Pt <- .elects_calculate_probability_making_decision(
            probability_stopping
        )

        # equation 7 additive smoothing
        Pt <- Pt + self$epsilon / T

        # equation 6, right term
        t <- torch::torch_ones(N, T, device = "cuda") * torch::torch_arange(
            T, dtype = "float", device = "cuda")

        earliness_reward <- Pt * .elects_probability_correct_class(
            log_class_probabilities, y_true
        ) * (1 - t / T)
        earliness_reward <- earliness_reward$sum(2)$mean(1)

        y_true <- y_true$to(device = "cpu")$repeat_interleave(repeats = as.integer(T))$to(device = "cuda")
        # equation 6 left term
        cross_entropy <- self$negative_log_likelihood(
            log_class_probabilities$view(c(N*T, C)), y_true$view(c(N*T)))$view(c(N, T))

        classification_loss <- (cross_entropy * Pt)$sum(2)$mean(1)

        # equation 6
        loss <- self$alpha * classification_loss - (1 - self$alpha) * earliness_reward

        # if (return_stats) {
        #     stats_lst <- list(
        #         classification_loss = classification_loss$cpu()$detach(),
        #         earliness_reward = earliness_reward$cpu()$detach(),
        #         probability_making_decision = Pt$cpu()$detach()
        #     )
        #     return(list(loss,stats_lst))
        # }

        return(loss)
    }
)

.elects_predict <- function(model, data, ..., callbacks = list(),
                            accelerator = NULL, verbose = NULL,
                            dataloader_options = NULL) {
    luz:::enable_mps_fallback()
    ctx <- luz:::predict_context$new(
        model = model$model,
        newdata = data,
        callbacks = callbacks,
        accelerator = accelerator,
        verbose = verbose,
        dataloader_options = dataloader_options,
        callbacks_default = luz:::default_predict_callbacks
    )

    pars <- rlang::list2(...)
    if (is.null(pars$stack))
        stack <- TRUE
    else
        stack <- pars$stack

    predict_fn <- if (is.null(ctx$model$predict)) ctx$model else ctx$model$predict
    on.exit({
        e <- rlang::current_env()
        rm(list = rlang::env_names(e), envir = e)
    }, add = TRUE)
    torch::with_no_grad({
        ctx$call_callbacks("on_predict_begin")
        luz:::with_handlers(
            !!! ctx$handlers,
            .expr = {
                coro::loop(for(batch in ctx$data) {
                    ctx$batch <- batch
                    ctx$call_callbacks("on_predict_batch_begin")
                    res <- do.call(predict_fn, list(ctx$input))
                    res <- torch::torch_stack(res[c(3, 4)], dim = 2)
                    ctx$pred[[length(ctx$pred) + 1]] <- res
                    #ctx$pred[[length(ctx$pred) + 1]] <- torchres[[3]]
                    #idx[[length(idx) + 1]] <- res[[4]]
                    ctx$call_callbacks("on_predict_batch_end")
                })
            }
        )
        ctx$call_callbacks("on_predict_end")
    })

    if (stack) {
        ctx$pred <- torch::torch_cat(ctx$pred)
    }

    ctx$pred
}
