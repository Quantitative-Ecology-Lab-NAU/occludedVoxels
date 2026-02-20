# Script for quantifying occlusion for TLS scans
# adapted from Andrew's script check_occlusion5
# Written 01/15/2026 by Nate Beine
# Updated 06/10/2024 - added detailed comments, cleaned up code and made it hella faster,
# added more diagnostic outputs, and improved visualization

# ========================================
# LIDAR OCCLUSION ANALYSIS - SINGLE AND BATCH PROCESSING
# ========================================

# BASE DIRECTORY ----
# Set the base directory for all file paths
BASE_DIR <- "d:/OneDrive - Northern Arizona University/Students/NBeine/NateOcclusionCode"

# ========================================
# USER PARAMETERS - CONFIGURE HERE
# ========================================

# Analysis parameters
VOXEL_SIZE <- 0.01          # Fine voxel size for vegetation detail (0.01 = 1cm)
OCCLUSION_SAMPLE_RES <- 0.03 # Visualization grid resolution (coarse display only; statistics use true 1cm voxels)
NCORES <- parallelly::availableCores(logical=FALSE) - 1  # CPU cores to use (availableCores() - 1 uses all but one)
SCANNER_HEIGHT <- 1.37      # Scanner height above ground in meters
BUFFER_DISTANCE <- 0.1      # Buffer around fitted quadrat square in meters (XY only)
HEIGHT_LIMIT_DEFAULT <- 1.5 # Normalized height limit in meters (default quadrats)
HEIGHT_LIMIT_IN <- 2.1      # Normalized height limit in meters (_IN quadrats)
QUADRAT_SIZE <- 1.0         # Target quadrat side length in meters (Hough square fit)

# Output options
WRITE_STATS_FILE <- TRUE     # Write occlusion statistics to a text file (TRUE/FALSE)
WRITE_HTML_FILE <- FALSE     # Write 3D visualization to an HTML file (TRUE/FALSE)

# Display configuration
cat("\n========================================
     OCCLUSION ANALYSIS PARAMETERS
     ========================================
     Vegetation voxel size: ", VOXEL_SIZE, " m (analysis + statistics)
     Visualization grid resolution: ", OCCLUSION_SAMPLE_RES, " m (display only)
     CPU cores: ", NCORES, " of ", parallelly::availableCores(), " available
     Scanner height: ", SCANNER_HEIGHT, " m
     Quadrat buffer: ", BUFFER_DISTANCE, " m
     Quadrat size (Hough fit): ", QUADRAT_SIZE, " m x ", QUADRAT_SIZE, " m
     Height limits: ", HEIGHT_LIMIT_DEFAULT, " m (default), ", HEIGHT_LIMIT_IN, " m (IN quadrats)
     ========================================\n\n", sep="")

# SINGLE SCAN VERSION ----

# Setup ----
library(lidR)
library(plotly)
library(Rcpp)
library(parallelly)

Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")
sourceCpp(file.path(BASE_DIR, "check_occlusion6_enhanced.cpp"))

# File paths ----
file_tls <- file.path(BASE_DIR, "ABOR_GCNP_TP_02_25_registered.las")
file_mls <- file.path(BASE_DIR, "ABOR_GCNP_TP_02_25_Q1_BOTH.las")

# Function: get_mls_hull ----
# Reads the MLS point cloud, builds a concave hull to locate the quadrat,
# then uses a Hough-based dominant-orientation search to fit a square of
# exactly QUADRAT_SIZE x QUADRAT_SIZE metres aligned to the quadrat edges.
# The fitted square (and its buffer) replace the raw concave hull throughout.
get_mls_hull <- function(file_path, buffer_distance, quadrat_size = QUADRAT_SIZE) {
  cat("Reading MLS file and building concave hull for orientation detection...\n")
  mls <- readTLSLAS(file_path, select = "xyz")
  mls_crs <- st_crs(mls)

  # Concave hull of MLS points (used only for orientation + centroid)
  raw_hull <- st_concave_hull(mls)
  raw_hull <- sf::st_set_crs(raw_hull, mls_crs)
  hull_coords <- sf::st_coordinates(raw_hull)[, 1:2]
  n_pts <- nrow(hull_coords) - 1L  # last == first for closed ring

  # ---- Hough-style dominant orientation --------------------------------
  # Accumulate edge angles (weighted by edge length) in [0, 90) degree bins.
  # A square has 90-degree rotational symmetry so folding to [0,90) captures
  # both pairs of parallel sides in the same peak.
  bin_width  <- 0.5          # degrees per bin
  bins       <- seq(0, 90 - bin_width, by = bin_width)
  accumulator <- numeric(length(bins))

  for (i in seq_len(n_pts)) {
    dx <- hull_coords[i + 1L, 1] - hull_coords[i, 1]
    dy <- hull_coords[i + 1L, 2] - hull_coords[i, 2]
    edge_len <- sqrt(dx*dx + dy*dy)
    if (edge_len < 1e-10) next
    angle_deg <- atan2(dy, dx) * 180 / pi   # [-180, 180]
    angle_deg <- angle_deg %% 90            # fold to [0, 90)
    bin_idx   <- findInterval(angle_deg, bins, rightmost.closed = TRUE)
    bin_idx   <- max(1L, min(bin_idx, length(accumulator)))
    accumulator[bin_idx] <- accumulator[bin_idx] + edge_len
  }

  best_bin   <- which.max(accumulator)
  best_angle <- bins[best_bin] + bin_width / 2   # bin midpoint (degrees)
  cat("  Hough dominant orientation: ", round(best_angle, 2), " deg\n", sep="")

  # ---- Build 1m x 1m square centred on the raw hull centroid -----------
  centroid   <- sf::st_centroid(raw_hull)
  cx         <- sf::st_coordinates(centroid)[1, 1]
  cy         <- sf::st_coordinates(centroid)[1, 2]

  half  <- quadrat_size / 2
  # Local corners (unrotated)
  local <- matrix(c(-half, -half,
                     half, -half,
                     half,  half,
                    -half,  half,
                    -half, -half),   # close the ring
                  ncol = 2, byrow = TRUE)

  theta <- best_angle * pi / 180
  R     <- matrix(c(cos(theta), -sin(theta),
                    sin(theta),  cos(theta)), nrow = 2)
  world <- t(R %*% t(local))
  world[, 1] <- world[, 1] + cx
  world[, 2] <- world[, 2] + cy

  square_poly <- sf::st_sfc(sf::st_polygon(list(world)), crs = mls_crs)
  square_buf  <- sf::st_buffer(square_poly, buffer_distance)

  cat("  Fitted square: ", quadrat_size, "m x ", quadrat_size,
      "m  |  centroid (", round(cx,3), ", ", round(cy,3),
      ")  |  angle ", round(best_angle,2), " deg\n",
      "  Hull buffered by ", buffer_distance, " m\n", sep="")

  return(list(
    hull            = square_poly,
    hull_buffered   = square_buf,
    hough_angle     = best_angle,
    quadrat_centroid = c(cx, cy)
  ))
}

# Function: process_tls_file ----
process_tls_file <- function(file_path, file_label, mls_hull_data, occlusion_sample_res = OCCLUSION_SAMPLE_RES) {
  cat("\n========================================
     Processing: ", file_label, "
     ========================================
     Step 1: Loading data...", sep="")

  # Load and preprocess ----
  las <- readTLSLAS(file_path, filter = "-thin_with_voxel 0.001")
  cat("\n  Points after loading: ", nrow(las), "\n", sep="")

  cat("Step 2: Clipping to ray transect area...\n")
  scanner_location <- c(0, 0, SCANNER_HEIGHT)

  # Ensure hulls have same CRS as TLS data
  las_crs <- st_crs(las)
  hull_buffered_corrected <- sf::st_set_crs(mls_hull_data$hull_buffered, las_crs)

  # Get bounding box of buffered hull
  hull_bbox <- sf::st_bbox(hull_buffered_corrected)

  # Extend bounding box to include scanner position and ray transect
  # Create a box from scanner (0,0) to the quadrat area, with some margin
  x_min_clip <- min(0, hull_bbox$xmin) - 0.5
  x_max_clip <- max(0, hull_bbox$xmax) + 0.5
  y_min_clip <- min(0, hull_bbox$ymin) - 0.5
  y_max_clip <- max(0, hull_bbox$ymax) + 0.5

  # Clip to expanded bounding box (includes transect vegetation)
  las <- filter_poi(las, X >= x_min_clip & X <= x_max_clip &
                         Y >= y_min_clip & Y <= y_max_clip)
  cat("  Points after transect area clipping: ", nrow(las), "\n", sep="")

  cat("Step 3: Classifying ground...\n")
  las <- classify_ground(las, algorithm=pmf(ws=.013, th=0.02))
  cat("  Points after ground classification: ", nrow(las), "
       Step 4: Normalizing height...
         (Converts absolute elevation to height above ground)
         (Without normalization: Z would be elevation, making quadrat height vary with terrain)\n", sep="")
  las <- normalize_height(las, tin())

  cat("\n     Step 5: Clipping to normalized height < 1.5 m...\n")
  las <- filter_poi(las, Z < 1.5)

  cat("Step 6: Classifying noise...
       Step 7: Removing noise...\n")
  las <- classify_noise(las, ivf(0.1, 1))
  las <- remove_noise(las)

  cat("\n     DIAGNOSTIC - Points available for voxelization:
         Total points: ", nrow(las), "
         Z range: ", min(las$Z), " to ", max(las$Z), "\n", sep="")

  # Hierarchical Occlusion Analysis ----
  cat("\n     Step 8: Creating voxel structure...
         Fine voxel size (vegetation + analysis): ", VOXEL_SIZE, " m
         Visualization grid resolution: ", occlusion_sample_res, " m", sep="")

  # Create fine voxels for populated areas (exact vegetation measurements)
  myVoxels_populated <- voxel_metrics(las, ~length(Z), res = VOXEL_SIZE, all_voxels = FALSE)
  myVoxels_populated$V1[is.na(myVoxels_populated$V1)] <- 0
  cat("\n       Fine populated voxels: ", format(nrow(myVoxels_populated), big.mark = ","), "\n", sep="")

  # Create sampling grid from AOI extent
  aoi_bbox <- sf::st_bbox(mls_hull_data$hull_buffered)

  cat("\n     Step 9: Creating sampling grid from AOI extent...
         AOI extent:
           X: ", round(aoi_bbox$xmin, 2), " to ", round(aoi_bbox$xmax, 2), "
           Y: ", round(aoi_bbox$ymin, 2), " to ", round(aoi_bbox$ymax, 2), "\n", sep="")

  # Create coarse sampling grid for VISUALIZATION only (axis-aligned)
  cat("\n     Step 10: Creating visualization sampling grid...\n")
  x_sample <- seq(aoi_bbox$xmin, aoi_bbox$xmax, by = occlusion_sample_res)
  y_sample <- seq(aoi_bbox$ymin, aoi_bbox$ymax, by = occlusion_sample_res)
  z_sample <- seq(0, HEIGHT_LIMIT_DEFAULT - occlusion_sample_res/2, by = occlusion_sample_res)

  sample_grid <- expand.grid(X = x_sample, Y = y_sample, Z = z_sample)
  sample_grid$V1 <- 0

  # Filter sample grid to only keep points within buffered hull
  sample_points_sf <- sf::st_as_sf(sample_grid, coords = c("X", "Y"), crs = las_crs)
  hull_buffered_corrected <- sf::st_set_crs(mls_hull_data$hull_buffered, las_crs)
  in_buffered_hull <- sf::st_within(sample_points_sf, hull_buffered_corrected, sparse = FALSE)[,1]

  # Get original coordinates back
  sample_grid <- sample_grid[in_buffered_hull, ]
  cat("       Initial sample points (axis-aligned): ", format(nrow(sample_grid) + sum(!in_buffered_hull), big.mark = ","), "
         Sample points within buffer: ", format(nrow(sample_grid), big.mark = ","), "
         (Statistics will use only samples within unbuffered hull)\n", sep="")

  # Mark sample points that overlap with populated voxels
  cat("\n       Identifying occupied sample points...\n")
  sample_grid$key <- paste(
    round(sample_grid$X / occlusion_sample_res) * occlusion_sample_res,
    round(sample_grid$Y / occlusion_sample_res) * occlusion_sample_res,
    round(sample_grid$Z / occlusion_sample_res) * occlusion_sample_res,
    sep = "_"
  )

  # Match populated voxels to sampling grid (a populated voxel marks its sample cell as occupied)
  pop_coarse_key <- paste(
    round(myVoxels_populated$X / occlusion_sample_res) * occlusion_sample_res,
    round(myVoxels_populated$Y / occlusion_sample_res) * occlusion_sample_res,
    round(myVoxels_populated$Z / occlusion_sample_res) * occlusion_sample_res,
    sep = "_"
  )

  occupied_cells <- unique(pop_coarse_key)
  sample_grid$V1[sample_grid$key %in% occupied_cells] <- 1
  sample_grid$key <- NULL

  cat("         Occupied samples: ", sum(sample_grid$V1 > 0), "
         Empty samples to check: ", sum(sample_grid$V1 == 0), "\n", sep="")

  # Coarse occlusion for visualization ----
  # Uses check_occlusion_accurate with ALL populated voxels (including transect
  # vegetation between scanner and quadrat) so the visualization correctly
  # shows occlusion caused by vegetation outside the AOI.
  #
  # Penetration-depth model: rays extend past each target voxel to find the
  # furthest populated voxel.  Targets beyond the last return = occluded.
  # Targets before the last return = empty (scanner saw through).
  cat("\n     Step 11: Determining occlusion for visualization grid...
         Sample points to process: ", format(nrow(sample_grid), big.mark = ","), "
         Obstacle voxels (full transect): ", format(nrow(myVoxels_populated), big.mark = ","), "
         Using ", NCORES, " CPU cores\n", sep="")

  obstacles_matrix <- as.matrix(myVoxels_populated[, c("X", "Y", "Z")])
  voxels_matrix <- as.matrix(sample_grid[, c("X", "Y", "Z", "V1")])
  scanner_location <- as.numeric(c(0, 0, SCANNER_HEIGHT))
  occlusion_results <- check_occlusion_accurate(
    voxels_matrix, obstacles_matrix,
    scanner_location, VOXEL_SIZE, NCORES
  )
  cat("       Occlusion check complete!\n")

  sample_grid$occluded <- occlusion_results

  # ============================================================
  # ACCURATE 1CM VOXEL ANALYSIS
  # Builds a true 1cm grid and ray-traces every empty voxel
  # using ALL populated voxels as obstacles (DDA algorithm).
  # No scaling or multiplication - exact voxel counts.
  # ============================================================
  cat("\n     Step 12: Building accurate 1cm analysis grid...\n")

  # Create 1cm grid aligned to voxel_metrics coordinate system
  x_fine <- seq(floor(aoi_bbox$xmin / VOXEL_SIZE) * VOXEL_SIZE,
                ceiling(aoi_bbox$xmax / VOXEL_SIZE) * VOXEL_SIZE,
                by = VOXEL_SIZE)
  y_fine <- seq(floor(aoi_bbox$ymin / VOXEL_SIZE) * VOXEL_SIZE,
                ceiling(aoi_bbox$ymax / VOXEL_SIZE) * VOXEL_SIZE,
                by = VOXEL_SIZE)
  z_fine <- seq(0, HEIGHT_LIMIT_DEFAULT - VOXEL_SIZE / 2, by = VOXEL_SIZE)

  cat("       Grid dimensions: ", length(x_fine), " x ", length(y_fine), " x ", length(z_fine), "\n", sep="")

  # Efficiently filter XY positions to buffered hull (check unique XY pairs only)
  xy_unique <- expand.grid(X = x_fine, Y = y_fine)
  xy_sf <- sf::st_as_sf(xy_unique, coords = c("X", "Y"), crs = las_crs)
  in_hull_xy <- sf::st_within(xy_sf, hull_buffered_corrected, sparse = FALSE)[,1]
  xy_valid <- xy_unique[in_hull_xy, ]
  cat("       XY positions in buffered hull: ", nrow(xy_valid), " of ", nrow(xy_unique), "\n", sep="")

  # Build full 3D grid for valid XY positions
  analysis_grid_1cm <- data.frame(
    X = rep(xy_valid$X, each = length(z_fine)),
    Y = rep(xy_valid$Y, each = length(z_fine)),
    Z = rep(z_fine, times = nrow(xy_valid))
  )
  cat("       Total 1cm analysis voxels: ", format(nrow(analysis_grid_1cm), big.mark=","), "\n", sep="")

  # Mark populated voxels by exact 1cm index matching (integer-based, no float issues)
  cat("       Matching populated voxels to analysis grid...\n")
  pop_ix <- as.integer(floor(myVoxels_populated$X / VOXEL_SIZE))
  pop_iy <- as.integer(floor(myVoxels_populated$Y / VOXEL_SIZE))
  pop_iz <- as.integer(floor(myVoxels_populated$Z / VOXEL_SIZE))
  pop_keys <- sprintf("%d,%d,%d", pop_ix, pop_iy, pop_iz)
  pop_key_set <- unique(pop_keys)

  grid_ix <- as.integer(floor(analysis_grid_1cm$X / VOXEL_SIZE))
  grid_iy <- as.integer(floor(analysis_grid_1cm$Y / VOXEL_SIZE))
  grid_iz <- as.integer(floor(analysis_grid_1cm$Z / VOXEL_SIZE))
  grid_keys <- sprintf("%d,%d,%d", grid_ix, grid_iy, grid_iz)

  analysis_grid_1cm$V1 <- as.integer(grid_keys %in% pop_key_set)
  cat("       Populated 1cm cells: ", format(sum(analysis_grid_1cm$V1 > 0), big.mark=","), "\n", sep="")
  cat("       Empty 1cm cells to ray-trace: ", format(sum(analysis_grid_1cm$V1 == 0), big.mark=","), "\n", sep="")

  # Step 13: Accurate 1cm occlusion with DDA ray tracing ----
  # Process in CHUNKS to avoid C stack overflow when passing
  # multi-million-row matrices to Rcpp, with a progress bar.
  cat("\n     Step 13: Running accurate 1cm occlusion analysis...
       Single-return model with ALL ", format(nrow(myVoxels_populated), big.mark=","),
      " populated voxels as obstacles\n", sep="")
  cat("       (includes transect vegetation between scanner and quadrat)\n")

  # --- Chunked processing with progress bar ---
  CHUNK_SIZE <- 2000000L  # 2M voxels per chunk (flat grid makes this fast)
  n_total <- nrow(analysis_grid_1cm)
  n_chunks <- ceiling(n_total / CHUNK_SIZE)
  cat("       Total analysis voxels: ", format(n_total, big.mark=","),
      "  |  Processing in ", n_chunks, " chunk(s)\n", sep="")
  cat("       Obstacle voxels: ", format(nrow(obstacles_matrix), big.mark=","),
      "  |  Using ", NCORES, " cores\n", sep="")
  flush.console()

  accurate_results <- integer(n_total)

  # Create a text progress bar (like lidR uses) — after summary so it's not interrupted
  pb <- txtProgressBar(min = 0, max = n_total, style = 3, char = "=", width = 60)

  for (chunk_i in seq_len(n_chunks)) {
    idx_start <- (chunk_i - 1L) * CHUNK_SIZE + 1L
    idx_end   <- min(chunk_i * CHUNK_SIZE, n_total)
    chunk_rows <- idx_start:idx_end

    chunk_matrix <- as.matrix(analysis_grid_1cm[chunk_rows, c("X", "Y", "Z", "V1")])

    chunk_result <- check_occlusion_accurate(
      chunk_matrix, obstacles_matrix,
      as.numeric(scanner_location), VOXEL_SIZE, NCORES
    )

    accurate_results[chunk_rows] <- chunk_result

    # Update progress bar to the end of this chunk
    setTxtProgressBar(pb, idx_end)

    # Free chunk memory
    rm(chunk_matrix, chunk_result)
    gc(verbose = FALSE)
  }

  close(pb)

  analysis_grid_1cm$occluded <- accurate_results
  rm(accurate_results)
  gc(verbose = FALSE)
  cat("\n       1cm occlusion check complete!\n")

  # Step 14: Exact 1cm statistics on unbuffered quadrat ----
  cat("\n     Step 14: Computing exact 1cm statistics...\n")

  # Filter to unbuffered hull using efficient XY matching
  xy_unique_grid <- unique(analysis_grid_1cm[, c("X", "Y")])
  xy_sf_grid <- sf::st_as_sf(xy_unique_grid, coords = c("X", "Y"), crs = las_crs)
  hull_corrected_1cm <- sf::st_set_crs(mls_hull_data$hull, las_crs)
  in_unbuff_xy <- sf::st_within(xy_sf_grid, hull_corrected_1cm, sparse = FALSE)[,1]
  xy_in_quadrat <- xy_unique_grid[in_unbuff_xy, ]

  # Filter analysis grid to unbuffered hull via integer XY keys
  unbuff_xy_keys <- sprintf("%d,%d",
    as.integer(floor(xy_in_quadrat$X / VOXEL_SIZE)),
    as.integer(floor(xy_in_quadrat$Y / VOXEL_SIZE)))
  grid_xy_keys <- sprintf("%d,%d",
    as.integer(floor(analysis_grid_1cm$X / VOXEL_SIZE)),
    as.integer(floor(analysis_grid_1cm$Y / VOXEL_SIZE)))
  in_quadrat_mask <- grid_xy_keys %in% unique(unbuff_xy_keys)
  analysis_quadrat <- analysis_grid_1cm[in_quadrat_mask, ]

  # Exact 1cm counts - no estimation or multiplication!
  fine_voxel_volume <- VOXEL_SIZE^3
  fine_total_count <- nrow(analysis_quadrat)
  fine_occupied_count <- sum(analysis_quadrat$occluded == 1)
  fine_empty_count <- sum(analysis_quadrat$occluded == 0)
  fine_occluded_count <- sum(analysis_quadrat$occluded == 2)

  cat("       1cm voxels in quadrat: ", format(fine_total_count, big.mark=","), "\n", sep="")
  cat("       Populated:     ", format(fine_occupied_count, big.mark=","),
      " (", round(fine_occupied_count/fine_total_count*100, 2), "%)\n", sep="")
  cat("       Empty visible: ", format(fine_empty_count, big.mark=","),
      " (", round(fine_empty_count/fine_total_count*100, 2), "%)\n", sep="")
  cat("       Occluded:      ", format(fine_occluded_count, big.mark=","),
      " (", round(fine_occluded_count/fine_total_count*100, 2), "%)\n", sep="")
  cat("       Consistency:   ", format(fine_occupied_count + fine_empty_count + fine_occluded_count, big.mark=","),
      " = ", format(fine_total_count, big.mark=","), " total\n", sep="")

  # Volume calculations from exact 1cm counts
  total_volume <- fine_total_count * fine_voxel_volume
  occupied_volume <- fine_occupied_count * fine_voxel_volume
  empty_visible_volume <- fine_empty_count * fine_voxel_volume
  occluded_volume <- fine_occluded_count * fine_voxel_volume

  # Coarse sample statistics (for visualization reference)
  sample_points_for_stats <- sf::st_as_sf(sample_grid, coords = c("X", "Y"), crs = las_crs, remove = FALSE)
  hull_corrected_stats <- sf::st_set_crs(mls_hull_data$hull, las_crs)
  in_unbuffered_hull <- sf::st_within(sample_points_for_stats, hull_corrected_stats, sparse = FALSE)[,1]
  sample_grid_quadrat <- sample_grid[in_unbuffered_hull, ]
  sample_total <- nrow(sample_grid_quadrat)
  sample_occupied <- sum(sample_grid_quadrat$V1 > 0 & sample_grid_quadrat$occluded == 1)
  sample_empty_visible <- sum(sample_grid_quadrat$V1 == 0 & sample_grid_quadrat$occluded == 0)
  sample_occluded <- sum(sample_grid_quadrat$occluded == 2)

  # Separate voxels inside vs outside the actual quadrat (for visualization)
  las_crs <- st_crs(las)
  hull_corrected <- sf::st_set_crs(mls_hull_data$hull, las_crs)
  voxel_points <- sf::st_as_sf(myVoxels_populated, coords = c("X", "Y"), crs = las_crs)
  in_quadrat <- sf::st_within(voxel_points, hull_corrected, sparse = FALSE)[,1]

  myVoxels_in_quadrat <- myVoxels_populated[in_quadrat, ]
  myVoxels_outside_quadrat <- myVoxels_populated[!in_quadrat, ]

  # Separate outside voxels into ray transect vs other areas
  if (nrow(myVoxels_outside_quadrat) > 0) {
    cat("\n       Classifying transect vegetation...\n")
    # Calculate if voxel is in ray transect (cone from scanner to quadrat)
    hull_coords <- sf::st_coordinates(hull_corrected)[,1:2]
    quadrat_center_x <- mean(hull_coords[,1])
    quadrat_center_y <- mean(hull_coords[,2])

    # Get buffered hull extent for ray cone
    hull_buffered_corrected <- sf::st_set_crs(mls_hull_data$hull_buffered, las_crs)
    hull_buffered_coords <- sf::st_coordinates(hull_buffered_corrected)[,1:2]
    max_extent <- max(sqrt((hull_buffered_coords[,1] - scanner_location[1])^2 +
                           (hull_buffered_coords[,2] - scanner_location[2])^2))

    # Direction from scanner to quadrat center
    dir_x <- quadrat_center_x - scanner_location[1]
    dir_y <- quadrat_center_y - scanner_location[2]
    dist_to_quadrat <- sqrt(dir_x^2 + dir_y^2)

    # For each outside voxel, check if it's in the ray transect cone
    vox_dx <- myVoxels_outside_quadrat$X - scanner_location[1]
    vox_dy <- myVoxels_outside_quadrat$Y - scanner_location[2]
    vox_dist <- sqrt(vox_dx^2 + vox_dy^2)

    # Project onto ray direction
    along_ray <- (vox_dx * dir_x + vox_dy * dir_y) / dist_to_quadrat

    # Perpendicular distance from ray axis
    perp_dist <- abs(vox_dx * dir_y - vox_dy * dir_x) / dist_to_quadrat

    # Calculate maximum perpendicular distance allowed (expands from scanner to quadrat)
    quadrat_radius <- max(sqrt((hull_coords[,1] - quadrat_center_x)^2 +
                               (hull_coords[,2] - quadrat_center_y)^2))
    # Cone expands linearly from scanner (0 width) to quadrat (quadrat_radius + buffer)
    max_perp_at_distance <- pmax(0, (quadrat_radius + BUFFER_DISTANCE) * (along_ray / dist_to_quadrat))

    # In ray transect if: in front of scanner, within expanding cone
    # No upper distance limit — cone extends to the full scene extent beyond the quadrat
    in_ray_transect <- along_ray > 0 &
                       perp_dist <= max_perp_at_distance

    myVoxels_in_ray <- myVoxels_outside_quadrat[in_ray_transect, ]
    myVoxels_other <- myVoxels_outside_quadrat[!in_ray_transect, ]

    cat("         Voxels in ray transect (cone area): ", format(nrow(myVoxels_in_ray), big.mark=","), "
           Other TLS voxels: ", format(nrow(myVoxels_other), big.mark=","), "\n", sep="")
  } else {
    myVoxels_in_ray <- data.frame()
    myVoxels_other <- data.frame()
  }

  # Color coding for visualization ----
  # Voxels in quadrat - colored by density using viridis palette
  myVoxels_in_quadrat$density_category <- cut(myVoxels_in_quadrat$V1,
                                   breaks = c(-Inf, 0, 1, 5, 10, 20, 50, 100, Inf),
                                   labels = c("empty", "1 point", "2-5 points", "6-10 points",
                                             "11-20 points", "21-50 points", "50+ points", "50+ points"),
                                   ordered_result = TRUE)

  # Define viridis color mapping
  viridis_palette <- c(
    "1 point" = "#440154",      # dark purple
    "2-5 points" = "#3B528B",   # blue-purple
    "6-10 points" = "#21908C",  # teal
    "11-20 points" = "#5DC863", # green
    "21-50 points" = "#FDE725", # yellow
    "50+ points" = "#FDE725"    # bright yellow
  )

  myVoxels_in_quadrat$color <- viridis_palette[as.character(myVoxels_in_quadrat$density_category)]
  myVoxels_in_quadrat$occluded <- 1  # All populated voxels are visible (they have points)
  myVoxels_in_quadrat$location <- "quadrat"

  # Voxels in ray transect (outside quadrat) - light green
  if (nrow(myVoxels_in_ray) > 0) {
    myVoxels_in_ray$color <- "limegreen"  # Bright green for ray transect
    myVoxels_in_ray$occluded <- 1
    myVoxels_in_ray$location <- "ray_transect"
  }

  # Other TLS voxels (outside ray) - dark green
  if (nrow(myVoxels_other) > 0) {
    myVoxels_other$color <- "darkgreen"  # Dark green for other TLS areas
    myVoxels_other$occluded <- 1
    myVoxels_other$location <- "other_tls"
  }

  # Sample grid for occlusion visualization
  sample_grid$color <- "lightgrey"
  sample_grid$color[sample_grid$V1 > 0] <- "green"  # Occupied samples
  sample_grid$color[sample_grid$occluded == 2] <- "red"  # Occluded samples

  # Calculate volume proportions
  prop_occupied <- occupied_volume / total_volume
  prop_empty_visible <- empty_visible_volume / total_volume
  prop_occluded <- occluded_volume / total_volume

  return(list(
    fine_voxels = myVoxels_in_quadrat,
    transect_voxels = myVoxels_in_ray,
    other_tls_voxels = myVoxels_other,
    sample_grid = sample_grid,
    total_volume = total_volume,
    occupied_volume = occupied_volume,
    prop_occupied = prop_occupied,
    empty_visible_volume = empty_visible_volume,
    prop_empty_visible = prop_empty_visible,
    occluded_volume = occluded_volume,
    prop_occluded = prop_occluded,
    fine_occupied_count = fine_occupied_count,
    fine_empty_count = fine_empty_count,
    fine_occluded_count = fine_occluded_count,
    fine_total_count = fine_total_count,
    scanner_location = scanner_location,
    hull = mls_hull_data$hull,
    sample_total = sample_total,
    sample_occupied = sample_occupied,
    sample_empty_visible = sample_empty_visible,
    sample_occluded = sample_occluded,
    file_label = file_label,
    occlusion_sample_res = occlusion_sample_res
  ))
}

# SINGLE SCAN EXECUTION ----
cat("\n######### SINGLE SCAN - GETTING MLS HULL #########\n")
mls_hull_data <- get_mls_hull(file_mls, BUFFER_DISTANCE)

cat("\n######### SINGLE SCAN - PROCESSING TLS FILE #########\n")
tls_results <- process_tls_file(file_tls, "TLS", mls_hull_data, occlusion_sample_res = OCCLUSION_SAMPLE_RES)

# Visualization ----
cat("\n     Creating full transect visualization...
       (Showing entire ray transect from scanner to quadrat)\n")

# Get data
fine_voxels <- tls_results$fine_voxels  # Quadrat voxels
transect_voxels <- tls_results$transect_voxels  # Ray transect vegetation (light green)
other_tls_voxels <- tls_results$other_tls_voxels  # Other TLS vegetation (dark green)
sample_grid <- tls_results$sample_grid
scanner_location <- tls_results$scanner_location

# Split 3cm sample grid into the three visualization categories
empty_samples     <- sample_grid[sample_grid$V1 == 0 & sample_grid$occluded == 0, ]
populated_samples <- sample_grid[sample_grid$V1 > 0  & sample_grid$occluded != 2, ]
occluded_samples  <- sample_grid[sample_grid$occluded == 2, ]

cat("\n       Ray transect vegetation (1cm): ", format(nrow(transect_voxels), big.mark=","), "
         Other TLS vegetation (1cm): ", format(nrow(other_tls_voxels), big.mark=","), "
         Empty samples (3cm): ", format(nrow(empty_samples), big.mark=","), "
         Populated samples (3cm): ", format(nrow(populated_samples), big.mark=","), "
         Occluded samples (3cm): ", format(nrow(occluded_samples), big.mark=","), "\n", sep="")

plot_quadrat <- plot_ly()

# Add scanner position
scanner_df <- data.frame(
  X = scanner_location[1],
  Y = scanner_location[2],
  Z = scanner_location[3]
)
plot_quadrat <- plot_quadrat %>%
  add_markers(data = scanner_df, x = ~X, y = ~Y, z = ~Z,
              marker = list(size = 10, color = "maroon", symbol = "circle"),
              name = "Scanner Position",
              showlegend = TRUE)

# Add ray transect vegetation (light green)
if (nrow(transect_voxels) > 0) {
  plot_quadrat <- plot_quadrat %>%
    add_markers(data = transect_voxels, x = ~X, y = ~Y, z = ~Z,
                color = I("limegreen"),
                marker = list(size = 1.5, opacity = 0.5),
                name = "Ray Transect Vegetation",
                showlegend = TRUE)
}

# Add other TLS vegetation (dark green)
if (nrow(other_tls_voxels) > 0) {
  plot_quadrat <- plot_quadrat %>%
    add_markers(data = other_tls_voxels, x = ~X, y = ~Y, z = ~Z,
                color = I("darkgreen"),
                marker = list(size = 1.5, opacity = 0.3),
                name = "Other TLS Vegetation",
                showlegend = TRUE)
}

# Add empty samples (3cm — light grey squares)
if (nrow(empty_samples) > 0) {
  plot_quadrat <- plot_quadrat %>%
    add_markers(data = empty_samples, x = ~X, y = ~Y, z = ~Z,
                color = I("lightgrey"),
                marker = list(size = 4, opacity = 0.3, symbol = 'square'),
                name = "Empty visible (3cm)")
}

# Add populated samples (3cm — viridis colored squares)
if (nrow(populated_samples) > 0) {
  plot_quadrat <- plot_quadrat %>%
    add_markers(data = populated_samples, x = ~X, y = ~Y, z = ~Z,
                color = I("#440154"),
                marker = list(size = 4, opacity = 0.7, symbol = 'square'),
                name = "Populated (3cm)")
}

# Add occluded samples (3cm — red squares)
if (nrow(occluded_samples) > 0) {
  plot_quadrat <- plot_quadrat %>%
    add_markers(data = occluded_samples, x = ~X, y = ~Y, z = ~Z,
                color = I("red"),
                marker = list(size = 4, opacity = 0.7, symbol = 'square'),
                name = "Occluded (3cm)")
}

# Add quadrat boundary from hull

# Legend for visualization:
# - Maroon sphere: Scanner position at origin
# - Light green (limegreen): Ray transect vegetation (scanner → quadrat cone area)
# - Dark green: Other TLS vegetation (outside ray transect)
# - Colored points: Quadrat vegetation (blue=sparse, purple=dense) - inside unbuffered hull
# - Light grey: Empty visible space
# - Red: Occluded space (blocked by vegetation)
# - Black outline: Quadrat boundary (MLS concave hull)

# Get hull boundary coordinates
hull <- tls_results$hull
hull_coords <- sf::st_coordinates(hull)[,1:2]

z_min_box <- 0
z_max_box <- HEIGHT_LIMIT_DEFAULT

# Draw bottom boundary
plot_quadrat <- plot_quadrat %>%
  add_trace(x = hull_coords[,1],
            y = hull_coords[,2],
            z = rep(z_min_box, nrow(hull_coords)),
            type = "scatter3d",
            mode = "lines",
            line = list(color = "black", width = 4),
            showlegend = TRUE,
            name = "Quadrat Boundary",
            hoverinfo = "skip")

# Draw top boundary
plot_quadrat <- plot_quadrat %>%
  add_trace(x = hull_coords[,1],
            y = hull_coords[,2],
            z = rep(z_max_box, nrow(hull_coords)),
            type = "scatter3d",
            mode = "lines",
            line = list(color = "black", width = 4),
            showlegend = FALSE,
            hoverinfo = "skip")

plot_quadrat <- plot_quadrat %>%
  layout(scene = list(aspectmode = "data",
                      xaxis = list(title = 'X (m)'),
                      yaxis = list(title = 'Y (m)'),
                      zaxis = list(title = 'Z - Normalized Height (m)')),
         title = paste0("Ray Transect Occlusion Analysis: Scanner → Quadrat AOI"))

print(plot_quadrat)

# Optionally save the visualization as an HTML file based on the TLS input file name
if (WRITE_HTML_FILE) {
  html_file <- file.path(dirname(file_tls),
                         paste0(tools::file_path_sans_ext(basename(file_tls)),
                                "_occlusion_plot.html"))
  htmlwidgets::saveWidget(plot_quadrat, html_file, selfcontained = TRUE)
  cat("\n     Visualization saved to: ", html_file, "\n", sep="")
}

# SUMMARY STATISTICS ----
cat("\n========================================
     OCCLUSION STATISTICS FOR ", tls_results$file_label, "
     ========================================\n", sep="")

# Volume proportions (from exact 1cm analysis)
prop_occupied <- tls_results$occupied_volume / tls_results$total_volume
prop_empty_visible <- tls_results$empty_visible_volume / tls_results$total_volume
prop_occluded <- tls_results$occluded_volume / tls_results$total_volume

cat("\n     EXACT 1cm Voxel Analysis (DDA ray tracing):
         Total 1cm voxels in quadrat: ", format(tls_results$fine_total_count, big.mark = ","), "
         Populated voxels: ", format(tls_results$fine_occupied_count, big.mark = ","),
      " (", round(prop_occupied * 100, 2), "%)
         Empty visible voxels: ", format(tls_results$fine_empty_count, big.mark = ","),
      " (", round(prop_empty_visible * 100, 2), "%)
         Occluded voxels: ", format(tls_results$fine_occluded_count, big.mark = ","),
      " (", round(prop_occluded * 100, 2), "%)\n", sep="")

cat("\n     Volume Analysis (exact from 1cm voxels):
         Total quadrat volume: ", round(tls_results$total_volume, 4), " m³
         Occupied volume: ", round(tls_results$occupied_volume, 6), " m³ (", round(prop_occupied * 100, 2), "%)
         Empty visible volume: ", round(tls_results$empty_visible_volume, 4), " m³ (", round(prop_empty_visible * 100, 2), "%)
         Occluded volume: ", round(tls_results$occluded_volume, 4), " m³ (", round(prop_occluded * 100, 2), "%)
         Average points per populated voxel: ", round(mean(tls_results$fine_voxels$V1), 1), "\n", sep="")

cat("\n     Visualization Grid Reference (", tls_results$occlusion_sample_res*100, "cm coarse samples):
         Sample points: ", format(tls_results$sample_total, big.mark = ","), "
         Occupied samples: ", tls_results$sample_occupied, "
         Empty visible samples: ", tls_results$sample_empty_visible, "
         Occluded samples: ", tls_results$sample_occluded, "\n", sep="")

# Show color distribution
cat("\n     Voxel color distribution (by point density):\n")
color_counts <- table(tls_results$fine_voxels$color)
for (col in names(color_counts)) {
  cat("         ", col, ": ", color_counts[col], "\n", sep="")
}

# Optionally write statistics to a text file based on the TLS input file name
if (WRITE_STATS_FILE) {
  stats_file <- file.path(dirname(file_tls),
                          paste0(tools::file_path_sans_ext(basename(file_tls)),
                                 "_occlusion_stats.txt"))

  stats_lines <- c(
    "========================================",
    paste0("OCCLUSION STATISTICS FOR ", tls_results$file_label),
    paste0("Generated: ", Sys.time()),
    "========================================",
    "",
    "EXACT 1cm Voxel Analysis (DDA ray tracing):",
    paste0("  Total 1cm voxels in quadrat: ", format(tls_results$fine_total_count, big.mark = ",")),
    paste0("  Populated voxels: ", format(tls_results$fine_occupied_count, big.mark = ","),
           " (", round(prop_occupied * 100, 2), "%)"),
    paste0("  Empty visible voxels: ", format(tls_results$fine_empty_count, big.mark = ","),
           " (", round(prop_empty_visible * 100, 2), "%)"),
    paste0("  Occluded voxels: ", format(tls_results$fine_occluded_count, big.mark = ","),
           " (", round(prop_occluded * 100, 2), "%)"),
    "",
    "Volume Analysis (exact from 1cm voxels):",
    paste0("  Total quadrat volume: ", round(tls_results$total_volume, 4), " m³"),
    paste0("  Occupied volume: ", round(tls_results$occupied_volume, 6), " m³ (", round(prop_occupied * 100, 2), "%)"),
    paste0("  Empty visible volume: ", round(tls_results$empty_visible_volume, 4), " m³ (", round(prop_empty_visible * 100, 2), "%)"),
    paste0("  Occluded volume: ", round(tls_results$occluded_volume, 4), " m³ (", round(prop_occluded * 100, 2), "%)"),
    paste0("  Average points per populated voxel: ", round(mean(tls_results$fine_voxels$V1), 1)),
    "",
    paste0("Visualization Grid Reference (", tls_results$occlusion_sample_res*100, "cm coarse samples):"),
    paste0("  Sample points: ", format(tls_results$sample_total, big.mark = ",")),
    paste0("  Occupied samples: ", tls_results$sample_occupied),
    paste0("  Empty visible samples: ", tls_results$sample_empty_visible),
    paste0("  Occluded samples: ", tls_results$sample_occluded),
    "",
    "Voxel color distribution (by point density):"
  )

  color_counts_tbl <- table(tls_results$fine_voxels$color)
  for (col in names(color_counts_tbl)) {
    stats_lines <- c(stats_lines, paste0("  ", col, ": ", color_counts_tbl[col]))
  }

  stats_lines <- c(stats_lines, "========================================")

  writeLines(stats_lines, stats_file)
  cat("\n     Occlusion statistics written to: ", stats_file, "\n", sep="")
}

cat("\n     ========================================\n     SINGLE SCAN ANALYSIS COMPLETE\n     ========================================\n")
