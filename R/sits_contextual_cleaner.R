sits_contextual_cleaner <- function(cube,
                                    window_size = 3L,
                                    target_class = 1L,
                                    mode_class = 2L,
                                    memsize = 4L,
                                    multicores = 2L,
                                    output_dir,
                                    version = "v1",
                                    progress = TRUE) {
    # Overlapping pixels
    overlap <- ceiling(window_size / 2L) - 1L
    # Get block size
    block <- .raster_file_blocksize(.raster_open_rast(.tile_path(cube)))
    # Check minimum memory needed to process one block
    job_block_memsize <- .jobs_block_memsize(
        block_size = .block_size(block = block, overlap = overlap),
        npaths = 1L, nbytes = 8L,
        proc_bloat = .conf("processing_bloat")
    )

    # Get input band
    band <- .cube_bands(cube)
    # Update multicores parameter
    multicores <- .jobs_max_multicores(
        job_block_memsize = job_block_memsize,
        memsize = memsize,
        multicores = multicores
    )
    # Update block parameter
    block <- .jobs_optimal_block(
        job_block_memsize = job_block_memsize,
        block = block,
        image_size = .tile_size(.tile(cube)),
        memsize = memsize,
        multicores = multicores
    )
    # Prepare parallelization
    .parallel_start(workers = multicores)
    on.exit(.parallel_stop(), add = TRUE)

    # Process each tile sequentially
    clean_cube <- .cube_foreach_tile(cube, function(tile) {
        # Process the data
        .contextual_cleaner_tile(
            tile = tile,
            block = block,
            band = band,
            window_size = window_size,
            target_class = target_class,
            mode_class = mode_class,
            overlap = overlap,
            output_dir = output_dir,
            version = version,
            progress = progress
        )
    })
    # Update cube class and return
    .set_class(clean_cube, "class_cube", class(clean_cube))
}

.contextual_cleaner_tile <- function(tile,
                                     block,
                                     band,
                                     window_size,
                                     target_class,
                                     mode_class,
                                     overlap,
                                     output_dir,
                                     version,
                                     progress) {
    # Output file
    out_file <- .file_derived_name(
        tile = tile, band = band, version = version, output_dir = output_dir
    )
    # Resume tile
    if (.raster_is_valid(out_file, output_dir = output_dir)) {
        # recovery message
        .check_recovery()
        # Create tile based on template
        tile <- .tile_derived_from_file(
            file = out_file, band = band,
            base_tile = tile, derived_class = .tile_derived_class(tile),
            labels = .tile_labels(tile),
            update_bbox = FALSE
        )
        return(tile)
    }
    # Create chunks as jobs
    chunks <- .tile_chunks_create(tile = tile, overlap = overlap, block = block)
    # Process jobs sequentially
    block_files <- .jobs_map_parallel_chr(chunks, function(chunk) {
        # Get job block
        block <- .block(chunk)
        # Block file name for each fraction
        block_files <- .file_block_name(
            pattern = .file_pattern(out_file),
            block = block,
            output_dir = output_dir
        )
        # Resume processing in case of failure
        if (.raster_is_valid(block_files)) {
            return(block_files)
        }
        # Read bands data
        values <- .clean_data_read(
            tile = tile, block = block, band = band
        )
        # Apply kernel modal
        values <- context_cleaner(
            x = as.matrix(values),
            ncols = block[["ncols"]],
            nrows = block[["nrows"]],
            band = 0L,
            window_size = window_size,
            target_class = target_class,
            mode_class = mode_class
        )
        # Prepare fractions to be saved
        band_conf <- .tile_band_conf(tile = tile, band = band)
        # Job crop block
        crop_block <- .block(.chunks_no_overlap(chunk))
        # Prepare and save results as raster
        .raster_write_block(
            files = block_files, block = block, bbox = .bbox(chunk),
            values = values, data_type = .data_type(band_conf),
            missing_value = .miss_value(band_conf),
            crop_block = crop_block
        )
        # Free memory
        gc()
        # Returned block files for each fraction
        block_files
    }, progress = progress)
    # Merge blocks into a new class_cube tile
    .tile_derived_merge_blocks(
        file = out_file,
        band = band,
        labels = .tile_labels(tile),
        base_tile = tile,
        derived_class = .tile_derived_class(tile),
        block_files = block_files,
        multicores = 1L,
        update_bbox = FALSE
    )
}
