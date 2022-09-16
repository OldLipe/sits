test_that("Combine predictions", {
    # create a data cube from local files
    data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
    cube <- sits_cube(
        source = "BDC",
        collection = "MOD13Q1-6",
        data_dir = data_dir,
        delim = "_",
        parse_info = c("X1", "tile", "band", "date")
    )
    # select a set of samples
    samples_ndvi <- sits_select(samples_modis_4bands, bands = c("NDVI"))
    # create a random forest model
    rfor_model <- sits_train(samples_ndvi, sits_rfor())
    # classify a data cube using rfor model
    probs_rfor_cube <- sits_classify(data = cube, ml_model = rfor_model,
                                     output_dir = tempdir(),
                                     version = "rfor")
    # create an XGBoost model
    xgb_model <- sits_train(samples_ndvi, sits_xgboost())
    # classify a data cube using xgboost model
    probs_xgb_cube <- sits_classify(data = cube, ml_model = xgb_model,
                                    output_dir = tempdir(),
                                    version = "xgb")
    # create a list of predictions to be combined
    pred_cubes <- list(probs_rfor_cube, probs_xgb_cube)
    # combine predictions
    comb_probs_cube <- sits_combine_predictions(cubes = pred_cubes,
                                                output_dir = tempdir(),
                                                version = "comb_rfor_xgb")
    expect_equal(sits_labels(comb_probs_cube), sits_labels(probs_xgb_cube))
    expect_equal(sits_bbox(comb_probs_cube), sits_bbox(probs_xgb_cube))
    expect_equal(nrow(comb_probs_cube), nrow(probs_xgb_cube))
})