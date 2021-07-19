context("Config")

test_that("Internal", {
    cubes <- sits:::.sits_config_satveg_cubes()
    expect_true(length(cubes) > 1)

    bbox <- sits:::.sits_config_satveg_bbox(cubes[1])
    expect_true(bbox["xmin"] < bbox["xmax"])

    expect_equal(sits:::.sits_config_color("NoClass"), "#737373")

    expect_true(sits:::.sits_config_memory_bloat() > 1)

    data_dir <- system.file("extdata/raster/cbers", package = "sits")
    cbers_cube <- sits_cube(
        source = "LOCAL",
        name = "022024",
        satellite = "CBERS-4",
        sensor = "AWFI",
        data_dir = data_dir,
        delim = "_",
        parse_info = c("X1", "X2", "tile", "band", "date")
    )

    bands <- sits:::sits_bands(cbers_cube)

    expect_true(sits:::.sits_config_minimum_values(cbers_cube, bands)[1] > -100000)
    expect_true(sits:::.sits_config_maximum_values(cbers_cube, bands)[1] < 100000)
})

test_that("Show", {
    con <- file(paste0(tempdir(), "config.txt"))
    writeLines(capture.output(sits_config_show()), con)
    close(con)

    lin <- readLines(paste0(tempdir(), "config.txt"))

    expect_equal(lin[1], "default:")
    expect_true(grepl("bloat", lin[5]))
    expect_true(grepl("rstac", lin[8]))
})
