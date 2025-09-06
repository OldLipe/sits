#' @title Train a model using LSTM for early detection
#' @name sits_elects
#'
#' @description Implementation of ELECTS
#'
#' @note
#' DOCUMENT
#'
#' @references
#' ADD ref
#'
#' @param samples            Time series with the training samples
#'                           (tibble of class "sits").
#' @param samples_validation Time series with the validation samples
#'                           (tibble of class "sits").
#'                           If \code{samples_validation} parameter is provided,
#'                           \code{validation_split} is ignored.
#' @param epochs             Number of iterations to train the model
#'                           (integer, min = 1, max = 20000).
#' @param batch_size         Number of samples per gradient update
#'                           (integer, min = 16L, max = 2048L)
#' @param validation_split   Fraction of training data
#'                           to be used as validation data.
#' @param optimizer          Optimizer function to be used.
#' @param opt_hparams        Hyperparameters for optimizer:
#'                           \code{lr} : Learning rate of the optimizer
#'                           \code{eps}: Term added to the denominator
#'                                to improve numerical stability.
#'                           \code{weight_decay}:       L2 regularization rate.
#' @param lr_decay_epochs    Number of epochs to reduce learning rate.
#' @param lr_decay_rate      Decay factor for reducing learning rate.
#' @param patience           Number of epochs without improvements until
#'                           training stops.
#' @param min_delta	         Minimum improvement in loss function
#'                           to reset the patience counter.
#' @param verbose            Verbosity mode (TRUE/FALSE). Default is FALSE.
#'
#' @return A fitted model to be used for classification of data cubes.
#'
#'
#' @examples
#' if (sits_run_examples()) {
#'     # create a lightTAE model
#'     torch_model <- sits_train(samples_modis_ndvi, sits_elects())
#'     # plot the model
#'     plot(torch_model)
#'     # create a data cube from local files
#'     data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
#'     cube <- sits_cube(
#'         source = "BDC",
#'         collection = "MOD13Q1-6.1",
#'         data_dir = data_dir
#'     )
#'     # classify a data cube
#'     probs_cube <- sits_classify(
#'         data = cube, ml_model = torch_model, output_dir = tempdir()
#'     )
#'     # plot the probability cube
#'     plot(probs_cube)
#'     # smooth the probability cube using Bayesian statistics
#'     bayes_cube <- sits_smooth(probs_cube, output_dir = tempdir())
#'     # plot the smoothed cube
#'     plot(bayes_cube)
#'     # label the probability cube
#'     label_cube <- sits_label_classification(
#'         bayes_cube,
#'         output_dir = tempdir()
#'     )
#'     # plot the labelled cube
#'     plot(label_cube)
#' }
#' @export
sits_elects <- function(samples = NULL,
                        samples_validation = NULL,
                        hidden_dims = 64L,
                        num_rnn_layers = 2L,
                        dropout = 0.2,
                        epochs = 150L,
                        batch_size = 128L,
                        validation_split = 0.2,
                        optimizer = torch::optim_adamw,
                        opt_hparams = list(
                            lr = 0.0005,
                            eps = 1e-08,
                            weight_decay = 7e-04
                        ),
                        lr_decay_epochs = 50L,
                        lr_decay_rate = 1.0,
                        patience = 20L,
                        min_delta = 0.01,
                        verbose = FALSE) {
    # set caller for error msg
    .check_set_caller("sits_elects")
    # Verifies if 'torch' and 'luz' packages is installed
    .check_require_packages(c("torch", "luz"))
    # documentation mode? verbose is FALSE
    verbose <- .message_verbose(verbose)
    # Function that trains a torch model based on samples
    train_fun <- function(samples) {
        # does not support working with DEM or other base data
        if (inherits(samples, "sits_base")) {
            stop(.conf("messages", "sits_train_base_data"), call. = FALSE)
        }
        # Avoid add a global variable for 'self'
        self <- NULL
        # Check validation_split parameter if samples_validation is not passed
        if (is.null(samples_validation)) {
            .check_num_parameter(validation_split,
                                 exclusive_min = 0.0, max = 0.5
            )
        }
        # Pre-conditions
        .check_pre_sits_lighttae(
            samples = samples, epochs = epochs,
            batch_size = batch_size,
            lr_decay_epochs = lr_decay_epochs,
            lr_decay_rate = lr_decay_rate,
            patience = patience, min_delta = min_delta,
            verbose = verbose
        )

        # Check opt_hparams
        # Get parameters list and remove the 'param' parameter
        optim_params_function <- formals(optimizer)[-1L]
        .check_opt_hparams(opt_hparams, optim_params_function)
        optim_params_function <- utils::modifyList(
            x = optim_params_function,
            val = opt_hparams
        )
        # Samples labels
        labels <- .samples_labels(samples)
        # Samples bands
        bands <- .samples_bands(samples)
        # Samples timeline
        timeline <- .samples_timeline(samples)
        # Create numeric labels vector
        code_labels <- seq_along(labels)
        names(code_labels) <- labels
        # Number of labels, bands, and number of samples (used below)
        n_labels <- length(labels)
        n_bands <- length(bands)
        n_times <- .samples_ntimes(samples)
        # Data normalization
        ml_stats <- .samples_stats(samples)
        # Organize train and the test data
        train_test_data <- .torch_train_test_samples(
            samples = samples,
            samples_validation = samples_validation,
            ml_stats = ml_stats,
            labels = labels,
            code_labels = code_labels,
            timeline = timeline,
            bands = bands,
            validation_split = validation_split
        )
        # Obtain the train and the test data
        train_samples <- train_test_data[["train_samples"]]
        test_samples <- train_test_data[["test_samples"]]
        n_samples_train <- nrow(train_samples)
        n_samples_test <- nrow(test_samples)

        # Organize data for model training
        train_x <- array(
            data = as.matrix(.pred_features(train_samples)),
            dim = c(n_samples_train, n_times, n_bands)
        )
        train_y <- unname(code_labels[.pred_references(train_samples)])
        # Create the test data
        test_x <- array(
            data = as.matrix(.pred_features(test_samples)),
            dim = c(n_samples_test, n_times, n_bands)
        )
        test_y <- unname(code_labels[.pred_references(test_samples)])
        # Set torch seed
        torch::torch_manual_seed(sample.int(10000L, 1L))
        # Define the ELECTS architecture
        model_elects <- torch::nn_module(
            classname = "model_elects",
            initialize = function(input_dim,
                                  hidden_dims,
                                  num_rnn_layers,
                                  nclasses,
                                  dropout) {
                # normalization over D-dimension. T-dimension is untouched
                # project to hidden_dims length
                self$intransforms <- torch::nn_sequential(
                    torch::nn_layer_norm(input_dim),
                    torch::nn_linear(input_dim, hidden_dims)
                )

                self$backbone <- torch::nn_lstm(
                    input_size = hidden_dims,
                    hidden_size = hidden_dims,
                    num_layers = num_rnn_layers,
                    bias = FALSE,
                    batch_first = TRUE,
                    dropout = dropout,
                    bidirectional = FALSE
                )

                # Heads
                self$classification_head <- .elects_classification_head(
                    hidden_dims, nclasses
                )
                self$stopping_decision_head <- .elects_decision_head(
                    hidden_dims
                )
            },
            forward = function(x) {
                x <- self$intransforms(x)
                res <- self$backbone(x)
                outputs <- res[[1]]
                last_state_list <- res[[2]]

                log_class_probabilities <- self$classification_head(outputs)
                probabilitiy_stopping <- self$stopping_decision_head(outputs)

                return(list(log_class_probabilities, probabilitiy_stopping))
            },
            step = function() {
                ctx$loss <- list()
                for (opt_name in names(ctx$optimizers)) {
                    return_stats <- FALSE
                    if (ctx$training) {
                        model_pred <- ctx$model(ctx$input)
                        ctx$call_callbacks("on_train_batch_after_pred")

                        preds <- self$predict(x = NULL, probs = model_pred)
                        ctx$pred <- preds[[3]]
                        opt <- ctx$optimizers[[opt_name]]

                        loss_fn <- .early_reward_loss()
                        loss <- loss_fn(model_pred[[1]], model_pred[[2]], ctx$target, return_stats)

                        ctx$call_callbacks("on_train_batch_after_loss")

                        ctx$call_callbacks("on_train_batch_before_step")
                        opt$zero_grad()
                        loss$backward()
                        opt$step()
                        ctx$call_callbacks("on_train_batch_after_step")
                        ctx$loss[[opt_name]] <- loss$detach()
                    } else {
                        return_stats <- TRUE
                        model_pred <- ctx$model(ctx$input)
                        preds <- ctx$model$predict(x = NULL, probs = model_pred)
                        ctx$pred <- preds[[3]]
                        ctx$call_callbacks("on_valid_batch_after_pred")

                        loss_fn <- .early_reward_loss()
                        loss <- loss_fn(model_pred[[1]], model_pred[[2]], ctx$target, return_stats)
                        ctx$loss[[opt_name]] <- loss
                        ctx$call_callbacks("on_valid_batch_after_loss")
                    }
                }
            },
            predict = function(x, probs = NULL) {
                torch::local_no_grad()
                res <- probs
                if (is.null(res)) {
                    res <- self$forward(x)
                }

                logprobabilities <- res[[1]]
                deltas <- res[[2]]
                batchsize <- logprobabilities$shape[1]
                sequencelength <- logprobabilities$shape[2]
                nclasses <- logprobabilities$shape[3]


                sample_stop_decision <- function(delta) {
                    dist <- torch::torch_stack(c(1 - delta, delta), dim = 2)
                    res <- torch::distr_categorical(dist)$sample()
                    res == 2
                }

                stop_list <- list()

                for (i in seq_len(sequencelength)) {
                    if (i < sequencelength - 1) {
                        stop_now <- sample_stop_decision(deltas[, i])
                        stop_list[[i]] <- stop_now
                    } else {
                        # make sure to stop last
                        last_stop <- torch::torch_ones(stop_now$shape)$bool()$to(device = "cuda")
                        stop_list[[i]] <- last_stop
                    }
                }
                # stack over the time dimension (multiple stops possible)
                stopped <- torch::torch_stack(stop_list, dim = 2)$bool()

                # is only true if stopped for the first time
                first_stops <- (stopped$cumsum(2) == 1) & stopped

                # time of stopping
                t_stop <-  first_stops$to("long")$argmax(2)

                # all predictions
                predictions <- logprobabilities$argmax(-1)

                # predictions at time of stopping
                predictions_at_t_stop <- torch::torch_masked_select(predictions, first_stops)

                return(list(logprobabilities, deltas, predictions_at_t_stop, t_stop))
            }
        )
        # torch 12.0 with luz not working with Apple MPS
        cpu_train <- .torch_cpu_train()
        # Train the model using luz
        torch_model <-
            luz::setup(
                module = model_elects,
                loss = .early_reward_loss,
                metrics = list(.acc_metric()),
                optimizer = optimizer
            ) |>
            luz::set_hparams(
                input_dim = n_bands,
                nclasses = n_labels,
                hidden_dims = hidden_dims,
                num_rnn_layers = num_rnn_layers,
                dropout = dropout
            ) |>
            luz::set_opt_hparams(
                !!!optim_params_function
            ) |>
            luz::fit(
                data = list(train_x, train_y),
                epochs = epochs,
                valid_data = list(test_x, test_y),
                callbacks = list(
                    luz::luz_callback_early_stopping(
                        monitor = "valid_loss",
                        mode = "min",
                        patience = patience,
                        min_delta = min_delta
                    ),
                    luz::luz_callback_lr_scheduler(
                        torch::lr_step,
                        step_size = lr_decay_epochs,
                        gamma = lr_decay_rate
                    )
                ),
                accelerator = luz::accelerator(cpu = cpu_train),
                dataloader_options = list(batch_size = batch_size),
                verbose = verbose
            )
        # Serialize model
        serialized_model <- .torch_serialize_model(torch_model[["model"]])

        # Retrieve attention mask
        # Get the encoder
        # encoder <- torch_model$model$temporal_encoder
        # Retrieve the attention mask from the encoder
        # attn_mask <- encoder$attention_heads$attention$attention_mask

        # Function that predicts labels of input values
        predict_fun <- function(values) {
            # Verifies if torch package is installed
            .check_require_packages("torch")
            # Set torch threads to 1
            # Note: function does not work on MacOS
            suppressWarnings(torch::torch_set_num_threads(1L))
            # Unserialize model
            torch_model[["model"]] <- .torch_unserialize_model(serialized_model)
            # Transform input into a 3D tensor
            # Reshape the 2D matrix into a 3D array
            n_samples <- nrow(values)
            n_times <- .samples_ntimes(samples)
            n_bands <- length(bands)
            # Performs data normalization
            values <- .pred_normalize(pred = values, stats = ml_stats)
            values <- array(
                data = as.matrix(values), dim = c(n_samples, n_times, n_bands)
            )
            # CPU or GPU classification?
            if (.torch_gpu_classification()) {
                # Get batch size
                batch_size <- sits_env[["batch_size"]]
                # transform the input array to a dataset
                values <- .torch_as_dataset(values)
                # Transform data set to dataloader to use the batch size
                values <- torch::dataloader(values, batch_size = batch_size)
                # GPU classification
                values <- .try(
                    .elects_predict(torch_model, values),
                    .msg_error = .conf("messages", ".check_gpu_memory_size")
                )
            } else {
                #  CPU classification
                values <- .elects_predict(torch_model, values)
            }
            # Convert from tensor to array
            values <- torch::as_array(values)
            # Update the columns names to labels
            #colnames(values) <- labels
            values
        }
        # Set model class
        predict_fun <- .set_class(
            predict_fun, "elects_model", "torch_model", "sits_model", class(predict_fun)
        )
        predict_fun
    }
    # If samples is informed, train a model and return a predict function
    # Otherwise give back a train function to train model further
    .factory_function(samples, train_fun)
}

