#' @title Train a model with two temporal encoders
#' @name sits_multitemporal
#'
#' @author Felipe Souza, \email{lipecaso@@gmail.com}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Trains a model with two temporal encoders, each one reading a
#' subset of the bands of the training samples. The encoded features of both
#' encoders are concatenated and passed to a multi-layer perceptron, which
#' produces the values used by the loss function.
#'
#' Any temporal deep learning model of \code{sits} can be used as an encoder:
#' \code{\link[sits]{sits_tempcnn}}, \code{\link[sits]{sits_lighttae}},
#' \code{\link[sits]{sits_tae}}, \code{\link[sits]{sits_resnet}} and
#' \code{\link[sits]{sits_lstm_fcn}}. The two encoders can be of the same type
#' (e.g., two temporal CNN with different hyperparameters) or of different
#' types (e.g., a temporal CNN and a lightweight temporal attention encoder).
#' Feature based models, such as \code{\link[sits]{sits_mlp}}, cannot be used
#' as encoders, since they do not read the samples as time series.
#'
#' Bands are assigned to each encoder with \code{bands1} and \code{bands2}.
#' A band can be assigned to both encoders, and by default each encoder reads
#' all bands of the training samples.
#'
#' @note
#' \code{sits} provides a set of default values for all classification models.
#' These settings have been chosen based on testing by the authors.
#' Nevertheless, users can control all parameters for each model.
#' Novice users can rely on the default values,
#' while experienced ones can fine-tune deep learning models
#' using \code{\link[sits]{sits_tuning}}.
#'
#' Each encoder is built with its own hyperparameters, including its own
#' \code{seed}. The \code{seed} of \code{sits_multitemporal} controls the
#' fusion layers and the training process.
#'
#' The design of two encoders merged by a multi-layer perceptron follows the
#' layer-level fusion strategy described by Ofori-Ampofo et al. referenced
#' below. If you use this method, please cite the fusion paper and the papers
#' of the encoders selected.
#'
#' @references Stella Ofori-Ampofo, Charlotte Pelletier and Stefan Lang,
#' "Crop Type Mapping from Optical and Radar Time Series Using
#' Attention-Based Deep Learning",
#' Remote Sensing, 13, 4668, 2021. \doi{10.3390/rs13224668}.
#'
#' @param samples            Time series with the training samples.
#' @param samples_validation Time series with the validation samples. if the
#'                           \code{samples_validation} parameter is provided,
#'                           the \code{validation_split} parameter is ignored.
#' @param encoder1           Temporal model used as the first encoder.
#' @param bands1             Bands read by the first encoder
#'                           (default: all bands).
#' @param encoder2           Temporal model used as the second encoder.
#' @param bands2             Bands read by the second encoder
#'                           (default: all bands).
#' @param encoder_dim        Dimension of the features produced by each
#'                           encoder. A single value is used by both encoders;
#'                           a vector of two values sets each one.
#' @param fusion_layers      Number of nodes of the fusion perceptron.
#' @param fusion_dropout_rates Dropout rates (0,1) of the fusion perceptron.
#' @param epochs             Number of iterations to train the model.
#' @param batch_size         Number of samples per gradient update.
#' @param validation_split   Fraction of training data to be used for
#'                           validation.
#' @param optimizer          Optimizer function to be used.
#' @param opt_hparams        Hyperparameters for optimizer:
#'                           lr : Learning rate of the optimizer
#'                           eps: Term added to the denominator
#'                                to improve numerical stability.
#'                           weight_decay:       L2 regularization
#' @param lr_decay_epochs    Number of epochs to reduce learning rate.
#' @param lr_decay_rate      Decay factor for reducing learning rate.
#' @param patience           Number of epochs without improvements until
#'                           training stops.
#' @param min_delta	         Minimum improvement in loss function
#'                           to reset the patience counter.
#' @param seed               Seed for random values.
#' @param verbose            Verbosity mode (TRUE/FALSE). Default is FALSE.
#'
#' @return A fitted model to be used for classification.
#'
#' @examples
#' if (sits_run_examples()) {
#'     # a temporal CNN and a lightweight temporal attention encoder,
#'     # each one reading one band of the training samples
#'     torch_model <- sits_train(
#'         samples_l8_rondonia_2bands,
#'         sits_multitemporal(
#'             encoder1 = sits_tempcnn(),
#'             bands1 = "NDVI",
#'             encoder2 = sits_lighttae(),
#'             bands2 = "EVI",
#'             epochs = 20,
#'             verbose = TRUE
#'         )
#'     )
#'     # plot the model
#'     plot(torch_model)
#'     # classify a point
#'     point_class <- sits_classify(
#'         data = samples_l8_rondonia_2bands[1:5, ],
#'         ml_model = torch_model
#'     )
#'     plot(point_class)
#' }
#' @family machine learning
#' @export
sits_multitemporal <- function(
  samples = NULL,
  samples_validation = NULL,
  encoder1 = sits_tempcnn(),
  bands1 = NULL,
  encoder2 = sits_lighttae(),
  bands2 = NULL,
  encoder_dim = 128L,
  fusion_layers = c(128L, 64L),
  fusion_dropout_rates = c(0.30, 0.30),
  epochs = 150L,
  batch_size = 64L,
  validation_split = 0.2,
  optimizer = torch::optim_adamw,
  opt_hparams = list(
    lr = 5.0e-04,
    eps = 1.0e-08,
    weight_decay = 1.0e-06
  ),
  lr_decay_epochs = 1L,
  lr_decay_rate = 0.95,
  patience = 20L,
  min_delta = 0.01,
  seed = NULL,
  verbose = FALSE
) {
  # set caller for error msg
  .check_set_caller("sits_multitemporal")
  # Verifies if 'torch' and 'luz' packages is installed
  .check_require_packages(c("torch", "luz"))
  # documentation mode? verbose is FALSE
  verbose <- .message_verbose(verbose)
  # Function that trains a torch model based on samples
  train_fun <- function(samples, embedding_dim = NULL) {
    # does not support working with DEM or other base data
    if (inherits(samples, "sits_base")) {
      stop(.conf("messages", "sits_train_base_data"), call. = FALSE)
    }
    # Verifies if torch package is installed
    .check_require_packages("torch")
    # Avoid add a global variable for 'self'
    self <- NULL
    # Check validation_split parameter if samples_validation is not passed
    if (is.null(samples_validation)) {
      .check_num_parameter(
        validation_split,
        exclusive_min = 0.0,
        max = 0.5
      )
    }
    # Preconditions
    .check_pre_sits_multitemporal(
      samples = samples,
      encoder1 = encoder1,
      bands1 = bands1,
      encoder2 = encoder2,
      bands2 = bands2,
      encoder_dim = encoder_dim,
      fusion_layers = fusion_layers,
      fusion_dropout_rates = fusion_dropout_rates,
      epochs = epochs,
      batch_size = batch_size,
      lr_decay_epochs = lr_decay_epochs,
      lr_decay_rate = lr_decay_rate,
      patience = patience,
      min_delta = min_delta,
      embedding_dim = embedding_dim,
      verbose = verbose
    )
    # Other pre-conditions:
    .check_int_parameter(seed, allow_null = TRUE)
    # Get optimizer hyperparameters
    optim_params_function <- .torch_optim_params(optimizer, opt_hparams)
    # Extract sample metadata and normalization statistics
    # Sample labels
    labels <- .samples_labels(samples)
    # Sample bands
    bands <- .samples_bands(samples)
    # Sample timeline
    timeline <- .samples_timeline(samples)
    # Create numeric labels vector
    code_labels <- seq_along(labels)
    names(code_labels) <- labels
    # Number of labels, bands, and number of samples (used below)
    n_labels <- length(labels)
    n_bands <- length(bands)
    n_times <- .samples_ntimes(samples)
    # Bands read by each encoder (by default, all bands)
    bands_enc <- list(
      .default(bands1, bands),
      .default(bands2, bands)
    )
    # Samples of each encoder, used to build it
    samples_enc <- purrr::map(bands_enc, function(bands_sel) {
      .samples_select_bands(samples, bands_sel)
    })
    # Position of the bands of each encoder in the samples bands, which
    # is the order of the third dimension of the input tensor
    index_enc <- purrr::map(samples_enc, function(samples_sel) {
      match(.samples_bands(samples_sel), bands)
    })
    # Dimension of the features produced by each encoder
    dim_enc <- rep_len(encoder_dim, 2L)
    # Post condition: do the encoders read the samples as time series?
    purrr::walk2(
      list(encoder1, encoder2),
      seq_len(2L),
      function(encoder, i) {
        .check_duo_temporal_encoder(
          encoder = encoder,
          samples = samples_enc[[i]],
          dim = dim_enc[[i]]
        )
      }
    )
    # Dimension of the fused features
    dim_fusion <- sum(dim_enc)
    # Data normalization
    ml_stats <- .samples_stats(samples)

    # Organize train and the test data
    train_samples <- .predictors(samples)
    feats <- .pred_features_normalize(
      pred = train_samples,
      stats = ml_stats
    )
    .pred_features(train_samples) <- feats
    # Post condition: is predictor data valid?
    .check_predictors(pred = train_samples, samples = samples)
    # Are there samples for validation?
    if (!is.null(samples_validation)) {
      .check_samples_validation(
        samples_validation = samples_validation,
        labels = labels,
        timeline = timeline,
        bands = bands
      )
      # Test samples are extracted from validation data
      test_samples <- .predictors(samples_validation)
      feats <- .pred_features_normalize(
        pred = test_samples,
        stats = ml_stats
      )
      .pred_features(test_samples) <- feats
    } else {
      # Split the data into training and validation data sets
      # Create partitions different splits of the input data
      test_samples <- .pred_sample(
        pred = train_samples,
        frac = validation_split
      )
      # Remove the lines used for validation
      sel <- !train_samples[["sample_id"]] %in%
        test_samples[["sample_id"]]
      train_samples <- train_samples[sel, ]
    }
    # Shuffle the data
    train_samples <- train_samples[
      sample(
        nrow(train_samples),
        nrow(train_samples)
      ),
    ]
    test_samples <- test_samples[
      sample(
        nrow(test_samples),
        nrow(test_samples)
      ),
    ]

    # Organize data for model training
    n_samples_train <- nrow(train_samples)
    n_samples_test <- nrow(test_samples)
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
    # luz drops the last incomplete batch when training, so a batch
    # larger than the training set yields no gradient updates at all
    .check_that(
      x = nrow(train_x) >= batch_size,
      msg = .conf("messages", "sits_train_batch_size")
    )
    # Create a torch seed (we define a new variable to allow users
    # to access this seed number from the model environment)
    torch_seed <- .torch_set_seed(seed)

    # Define the architecture of the two encoders and their fusion
    duo_model <- torch::nn_module(
      classname = "model_duo_temporal",
      initialize = function(
        n_labels,
        index_enc,
        dim_enc,
        dim_fusion,
        fusion_layers,
        fusion_dropout_rates
      ) {
        # Index of the bands read by each encoder
        self$index1 <- as.integer(index_enc[[1L]])
        self$index2 <- as.integer(index_enc[[2L]])
        # Build each encoder from its own samples, so that it is
        # created with the number of bands it reads
        self$encoder1 <- encoder1(
          samples_enc[[1L]],
          embedding_dim = dim_enc[[1L]]
        )
        self$encoder2 <- encoder2(
          samples_enc[[2L]],
          embedding_dim = dim_enc[[2L]]
        )
        # Building the encoders sets the torch seed to the seed of
        # each one, so restore ours to keep the fusion layers and
        # the training process reproducible
        torch::torch_manual_seed(torch_seed)
        # Fusion of the encoded features
        self$fusion <- .torch_multi_linear_batch_norm_relu_dropout(
          input_dim = dim_fusion,
          hidden_dims = fusion_layers,
          dropout_rates = fusion_dropout_rates
        )
        # Dimension of the last fusion layer
        dim_last <- fusion_layers[[length(fusion_layers)]]
        if (!.has(embedding_dim)) {
          # reduce to linear tensor with n_labels
          # softmax is done externally
          self$nn_linear <- torch::nn_linear(dim_last, n_labels)
        } else {
          self$nn_linear <- torch::nn_linear(
            dim_last,
            embedding_dim
          )
        }
      },
      forward = function(x) {
        # input is 3D n_samples x n_times x n_bands
        # each encoder reads the bands assigned to it. drop = FALSE
        # keeps the band dimension, which a single band would lose,
        # and the encoders read contiguous tensors
        x1 <- x[,, self$index1, drop = FALSE]
        x2 <- x[,, self$index2, drop = FALSE]
        # encode each set of bands
        e1 <- self$encoder1(x1$contiguous())
        e2 <- self$encoder2(x2$contiguous())
        # fuse the encoded features and classify
        torch::torch_cat(list(e1, e2), dim = 2L) |>
          self$fusion() |>
          self$nn_linear()
      }
    )
    # return encoder model
    if (.has(embedding_dim)) {
      return(duo_model(
        n_labels = n_labels,
        index_enc = index_enc,
        dim_enc = dim_enc,
        dim_fusion = dim_fusion,
        fusion_layers = fusion_layers,
        fusion_dropout_rates = fusion_dropout_rates
      ))
    }
    # train with CPU or GPU? some encoders are not supported by all
    # GPU backends, so the device is chosen by running each one on it
    cpu_train <- .torch_cpu_train_encoders(
      encoders = list(encoder1, encoder2),
      samples = samples_enc,
      dims = dim_enc
    )
    # Train the model using luz
    torch_model <-
      luz::setup(
        module = duo_model,
        loss = torch::nn_cross_entropy_loss(),
        metrics = list(luz::luz_metric_accuracy()),
        optimizer = optimizer
      ) |>
      luz::set_opt_hparams(
        !!!optim_params_function
      ) |>
      luz::set_hparams(
        n_labels = n_labels,
        index_enc = index_enc,
        dim_enc = dim_enc,
        dim_fusion = dim_fusion,
        fusion_layers = fusion_layers,
        fusion_dropout_rates = fusion_dropout_rates
      ) |>
      luz::fit(
        data = list(train_x, train_y),
        epochs = epochs,
        valid_data = list(test_x, test_y),
        callbacks = list(
          luz::luz_callback_early_stopping(
            monitor = "valid_loss",
            patience = patience,
            min_delta = min_delta,
            mode = "min"
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
    # remove data used for training
    force(rm(
      train_samples,
      test_samples,
      train_y,
      train_x,
      test_y,
      test_x
    ))
    gc()
    # Serialize model
    serialized_model <- force(.torch_serialize_model(torch_model$model))
    # Function that predicts labels of input values
    predict_fun <- function(values) {
      # Verifies if torch package is installed
      .check_require_packages("torch")

      # Unserialize model
      torch_model$model <- .torch_unserialize_model(
        model = torch_model$model,
        raw = serialized_model
      )
      # Use dataloader?
      if (.torch_use_dataloader(values)) {
        # Predict!
        values <- .torch_predict_chunks(
          torch_model = torch_model,
          dataset = values[["dataset"]],
          callback = values[["callback"]]
        )
        # Prepare results
        values <- unlist(values)
      } else {
        # Transform input into a 3D tensor
        n_samples <- nrow(values)
        # Performs data normalization
        values <- .pred_features_normalize(values, stats = ml_stats)
        # Represent matrix values as array
        C_as_array_inplace(values, c(n_samples, n_times, n_bands))
        # CPU classification
        values <- stats::predict(
          object = torch_model,
          newdata = values,
          accelerator = luz::accelerator(cpu = TRUE)
        )
        # Convert from tensor to array
        values <- torch::as_array(values)
        # Update the columns names to labels
        colnames(values) <- labels
      }
      # Return!
      values
    }
    # Set model class
    predict_fun <- .set_class(
      predict_fun,
      "torch_model",
      "sits_model",
      class(predict_fun)
    )
    predict_fun
  }
  # If samples is informed, train a model and return a predict function
  # Otherwise give back a train function to train model further
  .factory_function(samples, train_fun)
}
