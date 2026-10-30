#' @title Train a L-TAE model with base data
#' @name sits_lighttae_base
#'
#' @author Felipe Souza, \email{lipecaso@@gmail.com}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Charlotte Pelletier, \email{charlotte.pelletier@@univ-ubs.fr}
#'
#' @description Implementation of a Light Temporal Attention Encoder (L-TAE)
#' combined with base data (e.g., digital elevation models). The model has
#' two encoders: the first one encodes the satellite image time series with
#' a pixel-set encoder followed by a lightweight temporal self-attention
#' encoder, while the second one encodes the base data using a multi-layer
#' perceptron. The encoded features are then concatenated and passed to the
#' decoder, which produces the values used by the loss function.
#'
#' For the L-TAE model without base data, please see
#' \code{\link[sits]{sits_lighttae}}.
#'
#' Base data is attached to a data cube with
#' \code{\link[sits]{sits_add_base_cube}} and extracted with
#' \code{\link[sits]{sits_get_data}}. Since base data is stable in time
#' (one value per band per location), it is not suitable for the temporal
#' encoder; a dedicated encoder is used instead.
#'
#' @note
#' \code{sits} provides a set of default values for all classification models.
#' These settings have been chosen based on testing by the authors.
#' Nevertheless, users can control all parameters for each model.
#' Novice users can rely on the default values,
#' while experienced ones can fine-tune deep learning models
#' using \code{\link[sits]{sits_tuning}}.
#'
#' The time series encoder is based on the papers by Vivien Garnot
#' referenced below. The two-encoder design with feature-level fusion
#' follows the layer-level fusion strategy described by Ofori-Ampofo et al.
#' If you use this method, please cite the original TAE, L-TAE and fusion
#' papers.
#'
#' @references
#' Vivien Garnot, Loic Landrieu, Sebastien Giordano, and Nesrine Chehata,
#' "Satellite Image Time Series Classification with Pixel-Set Encoders
#' and Temporal Self-Attention",
#' 2020 Conference on Computer Vision and Pattern Recognition.
#' pages 12322-12331.
#' \doi{10.1109/CVPR42600.2020.01234}
#'
#' Vivien Garnot, Loic Landrieu,
#' "Lightweight Temporal Self-Attention  for Classifying
#' Satellite Images Time Series",
#' arXiv preprint arXiv:2007.00586, 2020.
#'
#' Stella Ofori-Ampofo, Charlotte Pelletier and Stefan Lang,
#' "Crop Type Mapping from Optical and Radar Time Series Using
#' Attention-Based Deep Learning",
#' Remote Sensing, 13, 4668, 2021. \doi{10.3390/rs13224668}.
#'
#' @param samples            Time series with the training samples
#'                           (tibble of class "sits_base").
#' @param samples_validation Time series with the validation samples
#'                           (tibble of class "sits_base").
#'                           If \code{samples_validation} parameter is provided,
#'                           \code{validation_split} is ignored.
#' @param base_layers        Number of nodes of the base data encoder.
#' @param base_dropout_rates Dropout rates (0,1) of the base data encoder.
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
#' @param seed               Seed for random values.
#' @param verbose            Verbosity mode (TRUE/FALSE). Default is FALSE.
#'
#' @return A fitted model to be used for classification of data cubes.
#'
#' @examples
#' if (sits_run_examples()) {
#'     # load samples and adjust their dates to the data cube
#'     samples <- read.csv(
#'         system.file("extdata/samples/samples_sinop_crop.csv",
#'             package = "sits"
#'         )
#'     )
#'     samples[["start_date"]] <- "2019-06-01"
#'     samples[["end_date"]] <- "2019-08-30"
#'     # region of interest covering the samples
#'     roi <- c(
#'         xmin = min(samples[["longitude"]]),
#'         xmax = max(samples[["longitude"]]),
#'         ymin = min(samples[["latitude"]]),
#'         ymax = max(samples[["latitude"]])
#'     )
#'     output_dir <- paste0(tempdir(), "/base_cube")
#'     if (!dir.exists(output_dir)) {
#'         dir.create(output_dir)
#'     }
#'     # create and regularize an optical data cube
#'     s2_cube <- sits_cube(
#'         source = "MPC",
#'         collection = "SENTINEL-2-L2A",
#'         bands = c("B05", "CLOUD"),
#'         start_date = "2019-06-01",
#'         end_date = "2019-08-30",
#'         roi = roi,
#'         crs = 4326
#'     )
#'     s2_reg <- sits_regularize(
#'         cube = s2_cube,
#'         period = "P16D",
#'         res = 232,
#'         roi = roi,
#'         multicores = 2,
#'         output_dir = output_dir
#'     )
#'     # create and regularize a digital elevation model cube
#'     dem_cube <- sits_cube(
#'         source = "MPC",
#'         collection = "COP-DEM-GLO-30",
#'         roi = roi,
#'         crs = 4326
#'     )
#'     dem_reg <- sits_regularize(
#'         cube = dem_cube,
#'         res = 232,
#'         roi = roi,
#'         multicores = 2,
#'         output_dir = output_dir
#'     )
#'     # add the elevation model as base data
#'     base_cube <- sits_add_base_cube(s2_reg, dem_reg)
#'     # extract training samples with base data
#'     samples_base <- sits_get_data(
#'         cube = base_cube,
#'         samples = samples,
#'         multicores = 2
#'     )
#'     # create a L-TAE model with a base data encoder
#'     # this example has few samples, so batch_size is reduced to ensure
#'     # that gradient updates are performed (luz drops incomplete batches)
#'     torch_model <- sits_train(
#'         samples_base,
#'         sits_lighttae_base(epochs = 20, batch_size = 8, verbose = TRUE)
#'     )
#'     # plot the model
#'     plot(torch_model)
#'     # classify the data cube
#'     probs_cube <- sits_classify(
#'         data = base_cube,
#'         ml_model = torch_model,
#'         output_dir = output_dir
#'     )
#'     # plot the probability cube
#'     plot(probs_cube)
#'     # label the probability cube
#'     label_cube <- sits_label_classification(
#'         probs_cube,
#'         output_dir = output_dir
#'     )
#'     # plot the labelled cube
#'     plot(label_cube)
#' }
#' @family machine learning
#' @export
sits_lighttae_base <- function(
  samples = NULL,
  samples_validation = NULL,
  base_layers = c(64L, 32L),
  base_dropout_rates = c(0.20, 0.20),
  epochs = 150L,
  batch_size = 64L,
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
  seed = NULL,
  verbose = FALSE
) {
  # set caller for error msg
  .check_set_caller("sits_lighttae_base")
  # Verifies if 'torch' and 'luz' packages is installed
  .check_require_packages(c("torch", "luz"))
  # documentation mode? verbose is FALSE
  verbose <- .message_verbose(verbose)
  # Function that trains a torch model based on samples
  train_fun <- function(samples, embedding_dim = NULL) {
    # Verifies if torch package is installed
    .check_require_packages("torch")
    # Avoid add a global variable for 'self'
    self <- NULL
    # Check validation_split parameter if samples_validation is not passed
    if (is.null(samples_validation)) {
      .check_num_parameter(validation_split, exclusive_min = 0.0, max = 0.5)
    }
    # Pre-conditions
    .check_pre_sits_lighttae_base(
      samples = samples,
      base_layers = base_layers,
      base_dropout_rates = base_dropout_rates,
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
    # Samples labels
    labels <- .samples_labels(samples)
    # Samples bands (base bands are encoded separately)
    bands <- .samples_bands(samples, include_base = FALSE)
    # Samples base bands
    base_bands <- .samples_base_bands(samples)
    # Samples timeline
    timeline <- .samples_timeline(samples)
    # Create numeric labels vector
    code_labels <- seq_along(labels)
    names(code_labels) <- labels
    # Number of labels, bands, and number of samples (used below)
    n_labels <- length(labels)
    n_bands <- length(bands)
    n_times <- .samples_ntimes(samples)
    # Number of observations of the base data (one for DEM-like data)
    n_base_times <- nrow(samples[["base_data"]][[1L]])
    # Size of each block of the flat predictors layout
    n_ts_feats <- n_bands * n_times
    n_base_feats <- length(base_bands) * n_base_times
    # Data normalization (covers time series and base bands)
    ml_stats <- .samples_stats(samples)

    # obtain training samples
    train_samples <- .predictors(samples)
    # [n, n_times*n_bands + n_base_feats]
    feats <- .pred_features_normalize(train_samples, stats = ml_stats)
    .pred_features(train_samples) <- feats
    # Post condition: is predictor data valid?
    .check_predictors(pred = train_samples, samples = samples)
    # Are there samples for validation?
    if (.has(samples_validation)) {
      # Validation samples must also contain base data
      .check_that(
        x = inherits(samples_validation, "sits_base"),
        msg = .conf("messages", "sits_train_base_samples")
      )
      .check_samples_validation(
        samples_validation = samples_validation,
        labels = labels,
        timeline = timeline,
        bands = .samples_bands(samples)
      )
      # Test samples are extracted from validation data
      test_samples <- .predictors(samples_validation)
      feats <- .pred_features_normalize(test_samples, stats = ml_stats)
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
    # The model receives the flat predictors layout and splits the
    # time series from the base features internally
    train_x <- as.matrix(.pred_features(train_samples))
    train_y <- unname(code_labels[.pred_references(train_samples)])
    # Create the test data
    test_x <- as.matrix(.pred_features(test_samples))
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

    # Define the L-TAE architecture with a base data encoder
    light_tae_base_model <- torch::nn_module(
      classname = "model_ltae_base",
      initialize = function(
        n_bands,
        n_times,
        n_labels,
        n_ts_feats,
        n_base_feats,
        timeline,
        base_layers,
        base_dropout_rates,
        layers_spatial_encoder = c(32L, 64L, 128L),
        n_heads = 16L,
        n_neurons = c(256L, 128L),
        dropout_rate = 0.2,
        dim_input_decoder = 128L,
        dim_layers_decoder = c(64L, 32L)
      ) {
        self$n_bands <- n_bands
        self$n_times <- n_times
        self$n_ts_feats <- n_ts_feats
        self$n_base_feats <- n_base_feats
        # define an spatial encoder
        self$spatial_encoder <-
          .torch_pixel_spatial_encoder(
            n_bands = n_bands,
            layers_spatial_encoder = layers_spatial_encoder
          )
        # number of input channels == last layer of mlp2
        in_channels <-
          layers_spatial_encoder[length(layers_spatial_encoder)]
        # define a temporal encoder
        self$temporal_encoder <-
          .torch_light_temporal_attention_encoder(
            timeline = timeline,
            in_channels = in_channels,
            n_heads = n_heads,
            n_neurons = n_neurons,
            dropout_rate = dropout_rate
          )
        # base data encoder - a multi-layer perceptron
        self$base_encoder <-
          .torch_multi_linear_batch_norm_relu_dropout(
            input_dim = n_base_feats,
            hidden_dims = base_layers,
            dropout_rates = base_dropout_rates
          )
        # the decoder takes the fused features, so its input is the
        # temporal encoding concatenated with the base encoding
        dim_fusion <- dim_input_decoder +
          base_layers[[length(base_layers)]]
        # add a final layer to the decoder
        # with a dimension equal to the number of layers
        if (!.has(embedding_dim)) {
          dim_layers_decoder[length(dim_layers_decoder) + 1L] <-
            n_labels
          self$decoder <- .torch_multi_linear_batch_norm_relu(
            dim_fusion,
            dim_layers_decoder
          )
          self$embedding_layer <- NULL
        } else {
          # Intermediate layers with BatchNorm + ReLU
          self$decoder <- .torch_multi_linear_batch_norm_relu(
            dim_fusion,
            dim_layers_decoder
          )
          # Final plain linear layer for embeddings
          # (no BatchNorm or ReLU — preserves full
          # representational capacity)
          last_hidden <-
            dim_layers_decoder[length(dim_layers_decoder)]
          self$embedding_layer <- torch::nn_linear(
            last_hidden,
            embedding_dim
          )
        }
      },
      forward = function(input) {
        # input is 2D n_samples x (n_ts_feats + n_base_feats)
        # time series features come first, base features last
        # (see .predictors.sits_base)
        # narrow returns views, so no data is copied here
        x_ts <- input$narrow(
          dim = 2L,
          start = 1L,
          length = self$n_ts_feats
        )
        x_base <- input$narrow(
          dim = 2L,
          start = self$n_ts_feats + 1L,
          length = self$n_base_feats
        )
        # predictors are band-major and time-minor; the pixel-set
        # encoder expects (batch, n_times, n_bands) and calls view(),
        # so the result must be made contiguous
        x_ts <- x_ts$reshape(
          c(x_ts$shape[[1L]], self$n_bands, self$n_times)
        )$permute(c(1L, 3L, 2L))$contiguous()
        # encode the time series
        e_ts <- x_ts |>
          self$spatial_encoder() |>
          self$temporal_encoder()
        # encode the base data
        e_base <- x_base$reshape(
          c(x_base$shape[[1L]], self$n_base_feats)
        ) |>
          self$base_encoder()
        # fuse the encoded features and decode
        out <- torch::torch_cat(list(e_ts, e_base), dim = 2L) |>
          self$decoder()
        if (!is.null(self$embedding_layer)) {
          out <- self$embedding_layer(out)
        }
        out
        # softmax is done externally
        # by .ml_normalize.torch_model function
      }
    )

    # return encoder model
    if (.has(embedding_dim)) {
      return(light_tae_base_model(
        n_bands = n_bands,
        n_times = n_times,
        n_labels = n_labels,
        n_ts_feats = n_ts_feats,
        n_base_feats = n_base_feats,
        timeline = timeline,
        base_layers = base_layers,
        base_dropout_rates = base_dropout_rates
      ))
    }
    # train with CPU or GPU?
    cpu_train <- .torch_cpu_train()
    # Train the model using luz
    torch_model <-
      luz::setup(
        module = light_tae_base_model,
        loss = torch::nn_cross_entropy_loss(),
        metrics = list(luz::luz_metric_accuracy()),
        optimizer = optimizer
      ) |>
      luz::set_hparams(
        n_bands = n_bands,
        n_times = n_times,
        n_labels = n_labels,
        n_ts_feats = n_ts_feats,
        n_base_feats = n_base_feats,
        timeline = timeline,
        base_layers = base_layers,
        base_dropout_rates = base_dropout_rates
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
    serialized_model <- .torch_serialize_model(torch_model$model)
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
        # Performs data normalization
        values <- .pred_features_normalize(values, stats = ml_stats)
        # CPU classification - the model splits the flat layout
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
