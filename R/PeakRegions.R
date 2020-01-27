#' mz points to frequency regions
#'
#' Given M/Z point data in a data.frame, create IRanges based point "regions" of
#' width 1, using the `frequency_multiplier` argument to convert from the floating
#' point double to an `integer`.
#'
#' @param data a `data.frame` containing `mz`
#' @param frequency_fit_description the description for frequency ~ mz
#' @param mz_fit_description the description for mz ~ frequency
#' @param frequency_multiplier a value used to convert to integers.
#'
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @export
mz_points_to_frequency_regions <- function(mz_data, frequency_fit_description = c(0, -1/2, -1/3),
                                           mz_fit_description = c(0, -1, -2, -3), frequency_multiplier = 400){
  frequency_list = mz_scans_to_frequency(mz_data, frequency_fit_description, mz_fit_description)

  # this is a check to make sure we will be able to convert
  # to multiply and still get integers out
  range_freq = range(frequency_list$frequency$frequency)

  if (any(is.na(range_freq))) {
    stop("NA entries in conversion from M/Z to frequency! Stopping!")
  }

  init_conversion = suppressWarnings(as.integer(range_freq * frequency_multiplier))
  na_conversion = any(is.na(init_conversion))
  reduce_value = 1
  if (na_conversion) {
    reduce_value = reduce_value - 0.1
    while (na_conversion && (reduce_value >= 0.3)) {
      new_multiplier = (frequency_multiplier * reduce_value)
      new_conversion = suppressWarnings(as.integer(range_freq * new_multiplier))
      na_conversion = any(is.na(new_conversion))
    }

    if (!na_conversion && (reduce_value >= 0.3)) {
      message("Frequency point multiplier modified to avoid NA values ...")
      frequency_multiplier = new_multiplier
    } else {
      stop("No good point multiplier found, stopping.")
    }
  }

  frequency_data = frequency_list$frequency

  frequency_regions = frequency_points_to_frequency_regions(frequency_data, multiplier = frequency_multiplier)
  frequency_list$frequency = NULL
  frequency_list$frequency_multiplier = frequency_multiplier
  frequency_regions@metadata = frequency_list
  frequency_regions
}

#' frequency points to frequency_regions
#'
#' Given a set of frequency points in a data.frame, create IRanges based point "regions"
#' of width 1, useing the `multiplier` to convert from a floating point double
#' to an `integer`
#'
#' @param frequency_data a `data.frame`
#' @param frequency_variable which column is the `frequency` stored in
#' @param multiplier value used to convert to integers
#'
#' @export
frequency_points_to_frequency_regions = function(frequency_data, frequency_variable = "frequency", multiplier = 400){
  frequency_regions <- IRanges::IRanges(start = round(frequency_data[, frequency_variable] * multiplier), width = 1)
  if (is.null(frequency_data$point)) {
    frequency_data$point <- seq(1, nrow(frequency_data))
  }
  S4Vectors::mcols(frequency_regions) <- frequency_data
  frequency_regions@metadata = list(multiplier = multiplier)
  frequency_regions
}

#' create frequency regions
#'
#' Given a point-point spacing and a frequency range, create IRanges based regions
#' of specified width. Overlapping sliding regions can be creating by specifying
#' a `region_size` bigger than `delta`, adjacent tiled regions can be created
#' by specifying a `region_size` == `delta`.
#'
#' @param point_spacing how far away are subsequent points.
#' @param frequency_range the range of frequency to use
#' @param n_point how many points you want to cover
#' @param delta_point the *step* size between the beginning of each subsequent region
#' @param multiplier multiplier to convert from frequency to integer space
#'
#' @details For Fourier-transform mass spec, points are equally spaced in
#'   frequency space, which will lead to unequal spacing in M/Z space. Therefore,
#'   we create regions using the point-point differences in frequency space.
#'
#'   What will be returned is an `IRanges` object, where the widths are constantly
#'   increasing over M/Z space.
#'
#' @export
#' @return IRanges
create_frequency_regions <- function(point_spacing = 0.5, frequency_range = NULL,
                              n_point = 10, delta_point = 1,
                              multiplier = 500){
  if (is.null(frequency_range)) {
    stop("A valid range of frequencies must be supplied!")
  }

  region_size = n_point * point_spacing
  step_size = point_spacing * delta_point
  use_frequencies = round(frequency_range)

  frequency_start = round(seq(use_frequencies[1], use_frequencies[2] - region_size, by = step_size) * multiplier)
  frequency_end = round(seq(use_frequencies[1] + region_size, use_frequencies[2], by = step_size) * multiplier)

  regions <- IRanges::IRanges(start = frequency_start, end = frequency_end)

  #S4Vectors::mcols(regions) <- list(mz_start = mz_start, mz_end = mz_end)
  regions@metadata <- list(multiplier = multiplier,
                           point_spacing = point_spacing,
                           n_point = n_point,
                           delta_point = delta_point)
  regions

}

#' @export
PeakRegions <- R6::R6Class("PeakRegions",
  public = list(
    frequency_point_regions = NULL,
    frequency_fit_description = NULL,
    mz_fit_description = NULL,

    peak_regions = NULL,
    sliding_regions = NULL,
    tiled_regions = NULL,

    peak_region_list = NULL,

    frequency_multiplier = NULL,
    scan_peaks = NULL,
    peak_data = NULL,
    scan_level_arrays = NULL,
    is_normalized = FALSE,
    normalization_factors = NULL,

    n_scan = NULL,
    scans_per_peak = NULL,
    scan_perc = NULL,
    min_scan = NULL,
    max_subsets = NULL,
    scan_subsets = NULL,

    frequency_range = NULL,

    scan_correlation = NULL,
    keep_peaks = NULL,
    peak_index = NULL,

    scan_indices = NULL,
    instrument = NULL,

    set_min_scan = function(){
      if (!is.null(self$frequency_point_regions)) {
        self$n_scan <- length(unique(self$frequency_point_regions@elementMetadata$scan))
        self$min_scan <- min(c(round(self$scan_perc * self$n_scan), 25))
        if (self$min_scan <= 3) {
          self$min_scan <- 3
        }
      } else {
        warning("No frequency points to pull scan data from!")
      }

      invisible(self)
    },

    add_data = function(raw_ms){
      if (!is.null(raw_ms)) {
        if (inherits(raw_ms, "RawMS")) {
          raw_mz_data = raw_ms$extract_raw_data()
          self$instrument = raw_ms$get_instrument()
        } else if (inherits(raw_ms, "data.frame")) {
          raw_mz_data = raw_ms
        }
        self$frequency_point_regions = mz_points_to_frequency_regions(mz_data = raw_mz_data,
                                                                      frequency_fit_description = self$frequency_fit_description,
                                                                      mz_fit_description = self$mz_fit_description,
                                                                      frequency_multiplier = self$frequency_multiplier)
        self$frequency_multiplier = self$frequency_point_regions@metadata$frequency_multiplier
        self$frequency_range = range(S4Vectors::mcols(self$frequency_point_regions)$frequency)
        self$set_min_scan()

      }
      invisible(self)
    },

    initialize = function(raw_ms = NULL,
                          frequency_fit_description = c(0, -1/2, -1/3),
                          mz_fit_description = c(0, -1, -2, -3),
                          frequency_multiplier = 400,
                          scan_perc = 0.1, max_subsets = 100){
      #browser(expr = TRUE)
      self$frequency_multiplier <- frequency_multiplier
      self$frequency_fit_description = frequency_fit_description
      self$mz_fit_description = mz_fit_description
      self$scan_perc <- scan_perc
      self$max_subsets <- max_subsets

      self$add_data(raw_ms)

      invisible(self)
    }
  )
)


#' @export
PeakRegionFinder <- R6::R6Class("PeakRegionFinder",
  public = list(
    run_time = NULL,
    start_time = NULL,
    stop_time = NULL,
    peak_regions = NULL,

    sliding_region_size = NULL,
    sliding_region_delta = NULL,
    quantile_multiplier = NULL,

    tiled_region_size = NULL,
    tiled_region_delta = NULL,

    region_percentage = NULL,

    peak_method = NULL,
    min_points = NULL,
    sample_id = NULL,

    zero_normalization = NULL,

    progress = NULL,

    add_regions = function(){
      if (self$progress) {
        message("Adding sliding and tiled regions ...")
      }
      log_message("Addling sliding and tiled regions ...")
      sliding_regions <- function(self){
        create_frequency_regions(frequency_range = self$peak_regions$frequency_range, n_point = self$sliding_region_size,
                                           delta_point = self$sliding_region_delta,
                                           multiplier = self$peak_regions$frequency_multiplier)
      }
      tiled_regions <- function(self){
        create_frequency_regions(frequency_range = self$peak_regions$frequency_range, n_point = self$tiled_region_size,
                          delta_point = self$tiled_region_delta,
                          multiplier = self$peak_regions$frequency_multiplier)
      }
      run_regions <- list(sliding = sliding_regions,
                          tiled = tiled_regions)
      new_regions <- internal_map$map_function(names(run_regions), function(region){
        run_regions[[region]](self)
      })
      names(new_regions) <- names(run_regions)
      self$peak_regions$sliding_regions <- new_regions[["sliding"]]
      self$peak_regions$tiled_regions <- new_regions[["tiled"]]
      invisible(self)
    },



    reduce_sliding_regions = function(){
      if (self$progress) {
        message("Finding initial signal regions ...")
      }
      log_message("Finding initial signal regions ...")
      self$peak_regions$peak_regions <- find_signal_regions(self$peak_regions$sliding_regions, self$peak_regions$frequency_point_regions, self$quantile_multiplier)
      log_memory()
    },

    split_peak_regions = function(use_regions = NULL){
      if (self$progress) {
        message("Splitting signal regions by peaks ...")
      }
      log_message("Splitting signal regions by peaks ...")
      if (is.null(use_regions)) {
        use_regions <- seq_len(length(self$peak_regions$peak_regions))
      }
      peak_data <- split_regions(self$peak_regions$peak_regions[use_regions], self$peak_regions$frequency_point_regions, self$peak_regions$tiled_regions, self$peak_regions$min_scan, peak_method = self$peak_method, min_points = self$min_points)

      self$peak_regions$peak_region_list = peak_data

      self$peak_regions$peak_index <- seq_len(length(peak_data))
      #self$peak_regions$peak_regions <- subset_signal_regions(self$)
    },

    remove_double_peaks_in_scans = function(){
      peak_region_list = self$peak_regions$peak_region_list

      scan_peaks <-  internal_map$map_function(peak_region_list, function(in_list){
        tmp_peaks = in_list$peaks
        dup_scans <- tmp_peaks[, "scan"][duplicated(tmp_peaks[, "scan"])]
        tmp_peaks = tmp_peaks[!(tmp_peaks[, "scan"] %in% dup_scans), ]
        tmp_points = in_list$points
        point_scans = tmp_points@elementMetadata$scan
        keep_scan = point_scans %in% unique(tmp_peaks$scan)
        in_list$points = tmp_points[keep_scan]
        in_list
      })

      n_remain <- purrr::map_int(scan_peaks, ~ nrow(.x$peaks))
      keep_remain <- n_remain > 0
      scan_peaks <- scan_peaks[keep_remain]
      self$peak_regions$peak_region_list <- scan_peaks

      self$peak_regions$scan_correlation <- self$peak_regions$scan_correlation[keep_remain, ]
      self$peak_regions$peak_index <- self$peak_regions$peak_index[keep_remain]

      if (!is.null(self$peak_regions$scans_per_peak)) {
        self$peak_regions$scans_per_peak <- self$peak_regions$scans_per_peak[keep_remain]
      }
      invisible(self)
    },

    normalize_data = function(which_data = "both"){
      if (self$progress) {
        message("Normalizing scans ...")
      }
      log_message("Normalizing scans ...")

      if (!self$zero_normalization) {
        self$peak_regions <- two_pass_normalization(self$peak_regions, summary_function = median,
                                                    normalize_peaks = which_data)
      } else {
        self$peak_regions = zero_normalization(self$peak_regions, summary_function = median,
                                               normalize_peaks = which_data)
      }

    },

    find_peaks_in_regions = function(which_data = "raw"){
      if (self$progress) {
        message("Finding peaks in regions ...")
      }
      log_message("Finding peaks in regions ...")
      self$peak_regions <- characterize_peaks(self$peak_regions)
    },

    model_mzsd = function(){
      self$peak_regions$peak_data$Log10ObservedMZSDModel <- mz_sd_model(self$peak_regions$peak_data)
      invisible(self)
    },

    model_heightsd = function(){
      self$peak_regions$peak_data$Log10HeightSDModel <-
        int_sd_model(self$peak_regions$peak_data)
      invisible(self)
    },

    indicate_high_frequency_sd = function() {
      peak_region = self$peak_regions

      sd_cutoff = max(boxplot.stats(peak_region$peak_data$ObservedFrequencySD)$stats)

      high_peaks = peak_region$peak_data$ObservedFrequencySD > sd_cutoff

      self$peak_regions$peak_data$HighSD = high_peaks
      invisible(self)
    },

    add_data = function(raw_ms) {
      if (inherits(raw_ms, "RawMS")) {
        self$peak_regions$add_data(raw_ms)
      }
      invisible(self)
    },

    summarize_peaks = function(){
      list(TIC = sum(self$peak_regions$peak_data$Height),
           Sample = self$sample_id,
           Peaks = self$peak_regions$peak_data,
           ScanLevel = self$peak_regions$scan_level_arrays)
    },

    add_offset = function(){
      peak_data = self$peak_regions$peak_data
      frequency_offset = self$peak_regions$frequency_point_regions@metadata$difference_range$most_common
      mz_coefficients = self$peak_regions$frequency_point_regions@metadata$mz_coefficients
      mz_description = self$peak_regions$frequency_point_regions@metadata$mz_fit_description
      mz_1 = predict_exponentials(peak_data$ObservedFrequency, mz_coefficients, mz_description)
      mz_2 = predict_exponentials(peak_data$ObservedFrequency + frequency_offset, mz_coefficients, mz_description)
      peak_data$Offset = (mz_1 - mz_2) * self$offset_multiplier
      self$peak_regions$peak_data = peak_data
      invisible(self)
    },

    sort_ascending_mz = function(){
      self$peak_regions$peak_data = self$peak_regions$peak_data[order(self$peak_regions$peak_data$ObservedMZ), ]
      invisible(self)
    },

    characterize_peaks = function(){
      self$add_regions()
      self$reduce_sliding_regions()
      self$split_peak_regions()
      self$remove_double_peaks_in_scans()
      self$normalize_data()
      self$find_peaks_in_regions()
      self$indicate_high_frequency_sd()
      self$add_offset()
      self$sort_ascending_mz()
      if (nrow(self$peak_regions$peak_data) == 0) {
        stop("No peaks meeting criteria!")
      }
      #self$add_offsets()
      self$model_mzsd()
      #self$model_heightsd()
    },

    summarize = function(package_used = "package:FTMS.peakCharacterization"){
      self$stop_time <- Sys.time()
      self$run_time <- as.numeric(difftime(self$stop_time, self$start_time, units = "s"))
      # generate information about our objects
      pkg_description <- utils::packageDescription(substring(package_used, 9))

      if (!is.null(pkg_description$RemoteSha)) {
        pkg_sha <- pkg_description$RemoteSha
      } else {
        pkg_sha <- NA
      }

      p_finder <- as.list(self)
      p_finder[[".__enclos_env__"]] <- NULL
      p_finder$clone <- NULL

      p_regions <- as.list(self$peak_regions)
      p_finder$peak_regions <- NULL

      p_regions[[".__enclos_env__"]] <- NULL
      p_regions$clone <- NULL
      p_regions$frequency_point_regions <- NULL
      #p_regions$mz_model <- NULL
      p_regions$sliding_regions <- NULL
      p_regions$tiled_regions <- NULL
      p_regions$scan_peaks <- NULL
      p_regions$peak_data <- NULL
      p_regions$peak_regions <- NULL
      p_regions$scan_level_arrays <- NULL
      p_regions$keep_peaks <- NULL
      p_regions$normalization_factors <- NULL
      p_regions$peak_index <- NULL
      p_regions$scan_correlation <- NULL
      p_regions$peak_region_list <- NULL

      processing_info <- list(Package = package_used,
                              Version = pkg_description$Version,
                              Sha = pkg_sha,
                              RunTime = self$run_time,
                              PeakFinder = p_finder,
                              PeakRegions = p_regions
      )
      # and then about the peaks we had
      list(processing_metadata.json = processing_info,
           peak_list.json = self$summarize_peaks())

    },

    peak_meta = function(){
      mz_point_data <- as.data.frame(self$peak_regions$frequency_point_regions@elementMetadata)
      mz_point_data <- split(mz_point_data, mz_point_data$scan)
      mz_point_data <- mz_point_data[names(mz_point_data) %in% as.character(self$peak_regions$normalization_factors$scan)]

      tmp_ms_info = purrr::map_df(mz_point_data, function(in_scan){
        data.frame(scan = in_scan[1, "scan"],
                   tic = sum(in_scan[, "intensity"]),
                   raw_tic = sum(in_scan[, "RawIntensity"]))
      })

      list(run_time_info = list(
              run_time = as.numeric(difftime(self$stop_time, self$start_time, units = "s")),
              n_peak_regions = length(self$peak_regions$peak_region_list),
              n_scans = length(self$peak_regions$normalization_factors$scan)
              ),
           ms_info = tmp_ms_info,
           frequency_mz = self$peak_regions$frequency_point_regions@metadata,
           peaks = list(n_peaks = nrow(self$peak_regions$peak_data),
                        mz_range = range(self$peak_regions$peak_data$ObservedMZ, na.rm = TRUE)),
           other_info = list(
             n_scans = length(self$peak_regions$normalization_factors$scan),
             min_scan = self$peak_regions$min_scan,
             n_peaks = nrow(self$peak_regions$peak_data),
             intensity_range = range(self$peak_regions$peak_data$Height),
             dynamic_range = max(self$peak_regions$peak_data$Height) / min(self$peak_regions$peak_data$Height),
             median_intensity = median(self$peak_regions$peak_data$Height),
             instrument = self$peak_regions$instrument
           ))
    },

    initialize = function(raw_ms = NULL, sliding_region_size = 10, sliding_region_delta = 1, tiled_region_size = 1, tiled_region_delta = 1,
                          region_percentage = 0.99, offset_multiplier = 1, frequency_multiplier = 400,
                          quantile_multiplier = 1.5, peak_method = "lm_weighted", min_points = 4,
                          zero_normalization = FALSE, frequency_fit_description = c(0, -1/2, -1/3),
                          mz_fit_description = c(0, -1, -2, -3), progress = FALSE){
      if (inherits(raw_ms, "RawMS")) {
        self$peak_regions <- PeakRegions$new(raw_ms = raw_ms$extract_raw_data(), frequency_fit_description = frequency_fit_description,
                                             mz_fit_description = mz_fit_description, frequency_multiplier = frequency_multiplier)
      } else if (inherits(raw_ms, "PeakRegions")) {
        self$peak_regions <- raw_ms
      } else {
        self$peak_regions <- PeakRegions$new(raw_ms = NULL, frequency_fit_description = frequency_fit_description,
                                             mz_fit_description = mz_fit_description, frequency_multiplier = frequency_multiplier)
      }

      self$sliding_region_size <- sliding_region_size
      self$sliding_region_delta <- sliding_region_delta
      self$tiled_region_size <- tiled_region_size
      self$tiled_region_delta <- tiled_region_delta
      self$region_percentage <- region_percentage
      self$quantile_multiplier = quantile_multiplier

      self$peak_method = peak_method
      self$min_points = min_points
      self$zero_normalization = zero_normalization
      self$progress = progress
      self$offset_multiplier = offset_multiplier


      invisible(self)
    }

  ),
  lock_objects = FALSE,
  lock_class = FALSE
)

count_overlaps <- function(regions, point_regions){
  zero_intensities <- point_regions@elementMetadata$intensity == 0

  nonzero_counts <- IRanges::countOverlaps(regions, point_regions[!zero_intensities])
  zero_counts <- IRanges::countOverlaps(regions, point_regions[zero_intensities])
  count_ratio <- log10(nonzero_counts + 1) - log10(zero_counts + 1)

  count_data = data.frame(nonzero_counts = nonzero_counts,
                          zero_counts = zero_counts,
                          count_ratio = count_ratio)
  S4Vectors::mcols(regions) = count_data
  regions
}

#' find signal
#'
#' Given some regions and point_regions, find the regions that actually should
#' contain **real** data. See *details* for an explanation of what is considered
#' **real**.
#'
#' @param regions the regions we want to query
#' @param point_regions the individual points
#' @param multiplier how much above base quantiles to use (default = 1.5)
#'
#' @export
#' @return IRanges
find_signal_regions <- function(regions, point_regions, multiplier = 1.5){
  regions = count_overlaps(regions, point_regions)
  nz_counts = regions@elementMetadata$nonzero_counts

  chunk_indices = seq(1, length(nz_counts), by = 10000)

  chunk_perc = purrr::map_dbl(chunk_indices, function(in_index){
    use_counts = nz_counts[seq(in_index, min(in_index + 9999, length(nz_counts)))]
    if (max(use_counts) > 0) {
      return(stats::quantile(use_counts, 0.99))
    } else {
      return(0)
    }
  })

  cutoff_value = ceiling(median(chunk_perc) * multiplier)

  regions = regions[nz_counts > cutoff_value]
  IRanges::reduce(regions)
}

create_na_peak <- function(peak_method = "lm_weighted"){
  data.frame(ObservedCenter = as.numeric(NA),
             Height = as.numeric(NA),
             Area = as.numeric(NA),
             SSR = as.numeric(NA),
             type = peak_method,
             stringsAsFactors = FALSE)
}

get_reduced_peaks <- function(in_range, peak_method = "lm_weighted", min_points = 4,
                              which = c("mz", "frequency")){
  range_point_data <- in_range@elementMetadata

  possible_peaks <- pracma::findpeaks(range_point_data$log_int, nups = round(min_points / 2))

  if (!is.null(possible_peaks)) {
    n_peak <- nrow(possible_peaks)
    peaks <- purrr::map_df(seq(1, n_peak), function(in_peak){
      #print(in_peak)
      "!DEBUG Peak `in_peak`"
      peak_loc <- seq(possible_peaks[in_peak, 3], possible_peaks[in_peak, 4])
      peak_data <- range_point_data[peak_loc, ]
      weights <- peak_data$intensity / max(peak_data$intensity)
      out_peak = purrr::map_dfc(which, function(in_which){
        tmp_peak = get_fitted_peak_info(peak_data, use_loc = in_which, w = weights)
        names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
        tmp_peak
      })
      #out_peak <- get_fitted_peak_info(peak_data, w = weights)
      #out_peak <- get_peak_info(range_data[peak_loc, ], peak_method = peak_method, min_points = min_points)
      out_peak$points <- I(list(IRanges::start(in_range)[peak_loc]))
      out_peak$scan <- range_point_data[peak_loc[1], "scan"]
      out_peak
    })
  } else {
    peaks <- purrr::map_dfc(which, function(in_which){
      tmp_peak = create_na_peak()
      names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
      tmp_peak
    })
    peaks$points <- NA
    peaks$scan <- range_point_data$scan[1]
  }
  peaks
}

#' split signal region
#'
#' Given a region that should contain signal, and the point data within it,
#' find the peaks, and return the region, and the set of points that make
#' up each point from each scan.
#'
#' @param region_list a list with points and tiles IRanges objects
#' @param peak_method the method for getting the peaks
#' @param min_points how many points are needed for a peak
#'
#' @export
#' @return list
split_region_by_peaks <- function(region_list, peak_method = "lm_weighted", min_points = 4){
  log_memory()
  frequency_point_regions = region_list$points
  tiled_regions = region_list$tiles
  frequency_point_regions@elementMetadata$log_int <- log(frequency_point_regions@elementMetadata$intensity + 1e-8)

  scan_runs = rle(frequency_point_regions@elementMetadata$scan)
  #if (min(scan_runs$lengths) < min_points + 2) {
  frequency_point_splitscan <- split(frequency_point_regions, frequency_point_regions@elementMetadata$scan)
  reduced_peaks <- purrr::map_df(names(frequency_point_splitscan), function(in_scan){
    get_reduced_peaks(frequency_point_splitscan[[in_scan]], peak_method = peak_method, min_points = min_points)
  })
  # } else {
  #   reduced_peaks = get_reduced_peaks(frequency_point_regions, peak_method = peak_method,
  #                                     min_points = min_points)
  # }

  reduced_peaks <- reduced_peaks[!is.na(reduced_peaks$ObservedCenter.frequency), ]

  if (nrow(reduced_peaks) > 0) {
    #reduced_peaks = convert_found_peaks(as.data.frame(S4Vectors::mcols(frequency_point_regions)), reduced_peaks)
    reduced_points <- frequency_points_to_frequency_regions(reduced_peaks, "ObservedCenter.frequency", frequency_point_regions@metadata$frequency_multiplier)

    secondary_regions = split_reduced_points(reduced_points, tiled_regions, n_zero = 1)

    n_region = length(secondary_regions$region)

    out_regions = purrr::map(seq_len(n_region), function(in_region){
      f_points = IRanges::IRanges(start = unlist(secondary_regions$peaks[[in_region]]$points),
                                  width = 1)
      f_data = IRanges::subsetByOverlaps(frequency_point_regions, f_points)
      t_data = IRanges::subsetByOverlaps(tiled_regions, secondary_regions$region[[in_region]])

      list(points = f_data, tiles = t_data, region = secondary_regions$region[[in_region]],
           peaks = rename_peak_data(secondary_regions$peaks[[in_region]]))
    })

  } else {
    out_regions = vector("list", length = 1)
    out_regions[[1]] = list(points = NULL, tiles = NULL,
                            region = NULL, peaks = NULL)
  }
  out_regions
}


convert_peaks_with_consensus_model = function(reduced_points){
  point_data = as.data.frame(S4Vectors::mcols(reduced_points))
  all_models = point_data[, c("intercept", "slope")]
  use_model = data.frame(intercept = median(all_models$intercept, na.rm = TRUE),
                         slope = median(all_models$slope, na.rm = TRUE))
  new_frequency = purrr::map_df(.x = point_data$ObservedCenter.mz,
                                .f = ~ mz_frequency_interpolation(.x, model = use_model))
  new_frequency$frequency = new_frequency$predicted_frequency
  new_frequency$predicted_frequency = NULL
  point_data[, c("frequency", "intercept", "slope")] = new_frequency[, c("frequency", "intercept", "slope")]

  new_points = mz_points_to_frequency_regions(point_data, reduced_points@metadata$frequency_multiplier)
  new_points
}

split_reduced_points = function(reduced_points, tiled_regions, n_zero = 2){
  overlap_counts <- data.frame(peak_count = IRanges::countOverlaps(tiled_regions, reduced_points))
  rle_counts = rle(overlap_counts$peak_count)
  rle_df = data.frame(lengths = rle_counts$lengths, values = rle_counts$values,
                      row = seq(1, length(rle_counts$lengths)))
  rle_index = vector("list", nrow(rle_df))
  start_index = 1
  for (irow in seq(1, nrow(rle_df))) {
    rle_index[[irow]] = seq(from = start_index, length.out = rle_df[irow, "lengths"])
    start_index = max(rle_index[[irow]]) + 1
  }

  mask_values = dplyr::filter(rle_df, values == 0, lengths >= n_zero)
  exclude_indices = unlist(rle_index[mask_values$row])

  sub_tiles <- IRanges::reduce(tiled_regions[-exclude_indices])
  S4Vectors::mcols(sub_tiles) <- list(region = seq(1, length(sub_tiles)))

  reduced_points <- IRanges::mergeByOverlaps(reduced_points, sub_tiles)

  sub_point <- as.list(S4Vectors::split(reduced_points, reduced_points$region))
  names(sub_point) <- NULL
  sub_region <- IRanges::IRangesList()

  for (iregion in seq_len(length(sub_point))) {
    all_points <- unique(unlist(sub_point[[iregion]]$points))
    tmp_range <- IRanges::IRanges(start = min(all_points), end = max(all_points))
    S4Vectors::mcols(tmp_range) <- I(list(as.numeric(sub_point[[iregion]]$scan)))
    sub_region[[iregion]] <- tmp_range
  }

  return(list(region = sub_region, peaks = sub_point))
}



split_regions <- function(signal_regions, frequency_point_regions, tiled_regions, min_scan, peak_method = "lm_weighted", min_points = 4) {
  # I'm really curious if we could reduce the memory footprint by splitting this into a list first, and then doing the
  # multiprocessing on it. If so, that would be a good thing to do.
  signal_list = as.list(split(signal_regions, seq(1, length(signal_regions))))
  min_scan2 = floor(min_scan / 2)
  point_regions_list = purrr::map(signal_list, function(in_region){
    points = IRanges::subsetByOverlaps(frequency_point_regions, in_region)

    nonzero = as.data.frame(points@elementMetadata) %>%
      dplyr::filter(intensity > 0) %>%
      dplyr::pull(scan) %>% unique(.) %>% length(.)

    if (nonzero >= min_scan2) {
      tiles = IRanges::subsetByOverlaps(tiled_regions, in_region)
      return(list(points = points, tiles = tiles, region = in_region))
    } else {
      return(NULL)
    }
  })

  not_null = purrr::map_lgl(point_regions_list, ~ !is.null(.x))
  point_regions_list = point_regions_list[not_null]
  point_regions_list = point_regions_list[sample(length(point_regions_list))]

  # split_data = purrr::imap(point_regions_list, function(.x, .y){
  #   message(.y)
  #   split_region_by_peaks(.x, peak_method = peak_method, min_points = min_points)
  # })
  split_data = purrr::imap(point_regions_list, function(in_list, in_id){
    message(in_id)
    split_region_by_peaks(in_list, peak_method = peak_method, min_points = min_points)
  })
  # split_data = internal_map$map_function(point_regions_list, split_region_by_peaks,
  #                                         peak_method = peak_method, min_points = min_points)

  null_regions = purrr::map_lgl(split_data, ~ is.null(.x[[1]]$points))
  split_data = split_data[!null_regions]

  out_data = unlist(split_data, recursive = FALSE, use.names = FALSE)
  return(out_data)
}

two_pass_normalization <- function(peak_regions, intensity_measure = c("RawHeight", "Height"), summary_function = median, normalize_peaks = "both"){
  scan_peaks <- purrr::map(peak_regions$peak_region_list, "peaks")

  normalization_factors <- single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function)

  normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  keep_after_round1 = which(purrr::map_lgl(normed_peaks, ~!is.null(.x)))

  normed_scan_cor <- purrr::map_dbl(normed_peaks[keep_after_round1], intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] <- 0
  low_cor <- abs(normed_scan_cor) <= 0.5

  normalization_factors <- single_pass_normalization(scan_peaks[keep_after_round1], intensity_measure = intensity_measure, summary_function = summary_function, use_peaks = low_cor)

  normed_peaks2 <- internal_map$map_function(scan_peaks[keep_after_round1], normalize_scan_peaks, normalization_factors)
  logical_round2 = purrr::map_lgl(normed_peaks2, ~!is.null(.x))
  keep_after_round2 = keep_after_round1[logical_round2]

  normed_list_regions = internal_map$map_function(peak_regions$peak_region_list[keep_after_round2], function(in_region){
    in_region$points = normalize_raw_points(in_region$points, normalization_factors)
    in_region$peaks = normalize_scan_peaks(in_region$peaks, normalization_factors)
    in_region
  })
  normed_raw <- normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

  peak_regions$frequency_point_regions <- normed_raw
  peak_regions$is_normalized <- "both"
  peak_regions$normalization_factors <- normalization_factors
  peak_regions$peak_index = peak_regions$peak_index[keep_after_round2]
  peak_regions$peak_region_list = normed_list_regions

  normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor[logical_round2])
  peak_regions$scan_correlation <- normed_scan_cor
  log_memory()
  peak_regions

}

calculate_number_of_scans <- function(in_peak){
  length(unique(in_peak$scan))
}

calculate_number_of_scans_normalized <- function(in_peak, normalized_scans){
  sum(unique(in_peak$scan) %in% normalized_scans)
}

intensity_scan_correlation <- function(scan_peak){
  if (nrow(scan_peak) < 3) {
    return(NA)
  }
  cor(scan_peak$Height, scan_peak$scan, method = "spearman", use = "complete.obs")
}

zero_normalization = function(peak_regions, intensity_measure = c("RawHeight", "Height"), summary_function = median, normalize_peaks = "both"){
  scan_peaks <- purrr::map(peak_regions$peak_region_list, "peaks")

  all_scans = unique(unlist(purrr::map(scan_peaks, ~ .x$scan)))
  normalization_factors <- data.frame(scan = sort(all_scans), normalization = 0)

  normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  normed_list_regions = internal_map$map_function(peak_regions$peak_region_list, function(in_region){
    in_region$points = normalize_raw_points(in_region$points, normalization_factors)
    in_region$peaks = normalize_scan_peaks(in_region$peaks, normalization_factors)
    in_region
  })


  normed_scan_cor <- purrr::map_dbl(normed_peaks, intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] <- 0
  low_cor <- abs(normed_scan_cor) <= 0.5

  normed_raw <- normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

  peak_regions$peak_region_list = normed_list_regions
  peak_regions$frequency_point_regions <- normed_raw
  peak_regions$is_normalized <- "both"
  peak_regions$normalization_factors <- normalization_factors

  normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                HighCor = !low_cor)
  n_scans <- purrr::map_int(scan_peaks, calculate_number_of_scans)
  peak_regions$scan_correlation <- normed_scan_cor
  peak_regions
}

single_pass_normalization <- function(scan_peaks, intensity_measure = c("RawHeight", "Height"), summary_function = median, use_peaks = NULL, min_ratio = 0.7){
  n_scan_per_peak <- purrr::map_int(scan_peaks, function(x){
    if (sum(duplicated(x$scan)) == 0) {
      return(length(x$scan))
    } else {
      return(0L)
    }
  })

  use_measure <- NULL
  for (imeasure in intensity_measure){
    if (imeasure %in% names(scan_peaks[[1]])) {
      use_measure <- imeasure
      break()
    }
  }

  if (is.null(use_measure)) {
    stop("The desired intensity measure for normalization is not present!")
  }

  scan_cutoff <- quantile(n_scan_per_peak, 0.95)

  if (is.null(use_peaks)) {
    use_peaks <- rep(TRUE, length(scan_peaks))
  }

  normalize_peaks <- which((n_scan_per_peak >= scan_cutoff) & use_peaks)

  all_scans <- data.frame(scan = unique(unlist(purrr::map(scan_peaks, function(x){unique(x$scan)}))))

  peak_intensity <- purrr::map_dfc(scan_peaks[normalize_peaks], function(x){
    tmp_data <- dplyr::left_join(all_scans, as.data.frame(x[, c("scan", use_measure)]), by = "scan")
    log(tmp_data[, use_measure])
  })

  all_na = purrr::map_lgl(seq_len(nrow(peak_intensity)), ~ sum(is.na(peak_intensity[.x, ])) == ncol(peak_intensity))

  all_scans = all_scans[!all_na, , drop = FALSE]
  peak_intensity = peak_intensity[!all_na, , drop = FALSE]

  intensity_ratio <- purrr::map_dfr(seq_len(nrow(peak_intensity)), function(x){
    #message(x)
    peak_intensity[x, ] / max(peak_intensity[x, ], na.rm = TRUE)
  })
  peak_intensity[intensity_ratio < min_ratio] <- NA

  intensity_scans <- purrr::map_int(seq_len(ncol(peak_intensity)), function(x){
    sum(!is.na(peak_intensity[, x]))
  })

  peak_intensity <- peak_intensity[, intensity_scans >= scan_cutoff, drop = FALSE]

  notna_scans <- rowSums(!is.na(as.matrix(peak_intensity)))
  keep_scans <- notna_scans >= 25

  if (sum(keep_scans) == 0) {
    stop("No scans left to use in normalization!")
  }

  all_scans <- all_scans[keep_scans, , drop = FALSE]
  peak_intensity <- peak_intensity[keep_scans, ]

  scan_distances <- purrr::map_dbl(seq(1, nrow(peak_intensity)), function(in_scan){
    scan_peaks <- peak_intensity[in_scan, ,drop = FALSE]
    scan_peaks_matrix <- matrix(unlist(scan_peaks), nrow = nrow(peak_intensity) - 1, ncol = ncol(scan_peaks), byrow = TRUE)

    other_matrix <- as.matrix(peak_intensity[-in_scan, , drop = FALSE])
    scan_other_diff <- scan_peaks_matrix - other_matrix
    scan_distances <- purrr::map_dbl(seq(1, nrow(scan_other_diff)), function(x){
      sqrt(sum(scan_other_diff[x, ]^2, na.rm = TRUE))
    })
    sum(scan_distances)
  })

  normalize_scan <- which.min(scan_distances)

  scan_norm_matrix <- matrix(unlist(peak_intensity[normalize_scan, , drop = FALSE]),
                             nrow = nrow(peak_intensity), ncol = ncol(peak_intensity), byrow = TRUE)

  diff_matrix <- as.matrix(peak_intensity) - scan_norm_matrix

  normalization_factors <- purrr::map_dbl(seq_len(nrow(diff_matrix)), function(in_row){
    summary_function(diff_matrix[in_row, ], na.rm = TRUE)
  })
  data.frame(scan = all_scans$scan, normalization = normalization_factors)
}

normalize_raw_points <- function(raw_points, normalization_factors){

  raw_points@elementMetadata$index <- seq_len(length(raw_points))

  if (is.null(raw_points@elementMetadata$RawIntensity)) {
    raw_points@elementMetadata$RawIntensity <- raw_points@elementMetadata$intensity
  }

  point_height_scan <- as.data.frame(raw_points@elementMetadata[, c("index", "RawIntensity", "scan")])

  point_height_scan <- dplyr::right_join(point_height_scan, normalization_factors, by = "scan")
  point_height_scan <- point_height_scan[!is.na(point_height_scan$RawIntensity), ]
  point_height_scan <- point_height_scan[order(point_height_scan$index), ]
  raw_points <- raw_points[point_height_scan$index]

  stopifnot(raw_points@elementMetadata$index == point_height_scan$index)

  raw_points@elementMetadata$intensity <- exp(log(point_height_scan$RawIntensity) -
                                                point_height_scan$normalization)

  raw_points
}

normalize_scan_peaks <- function(in_peak, normalization_factors){
  #split_normalization <- split(normalization_factors$normalization, normalization_factors$scan)

  #normed_peaks <- internal_map$map_function(scan_peaks, function(in_peak){
  if (is.null(in_peak$RawHeight)) {
    in_peak$RawHeight <- in_peak$Height
  }
  in_peak$index <- seq(1, nrow(in_peak))
  height_scan <- as.data.frame(in_peak[, c("index", "RawHeight", "scan")])
  height_scan <- dplyr::right_join(height_scan, normalization_factors, by = "scan")
  height_scan <- height_scan[!is.na(height_scan$RawHeight), ]
  if (nrow(height_scan) > 0) {
    height_scan <- height_scan[order(height_scan$index), ]
    in_peak <- in_peak[height_scan$index, ]
    stopifnot(in_peak@elementMetadata$index == height_scan$index)
    in_peak$Height <- exp(log(height_scan$RawHeight) - height_scan$normalization)
  } else {
    in_peak = NULL
  }
  in_peak
}

apply_normalization_peak_regions <- function(peak_regions, normalization_factors, which_data = "both") {
  to_normalize <- switch(which_data,
                            both = c("frequency_point_regions", "scan_peaks"),
                            raw = "frequency_point_regions",
                            scan_peaks = "scan_peaks")

  split_normalization <- split(normalization_factors$normalization, normalization_factors$scan)

  if ("frequency_point_regions" %in% to_normalize) {
    split_points <- as.list(split(peak_regions$frequency_point_regions, peak_regions$frequency_point_regions@elementMetadata$scan))
    split_points <- split_points[names(split_points) %in% names(split_normalization)]

    normed_points <- purrr::map2(split_points, split_normalization, function(in_points, in_norm){
      in_points@elementMetadata$intensity <- exp(log(in_points@elementMetadata$intensity) - in_norm)
      in_points
    })

    peak_regions$frequency_point_regions <- unlist(IRanges::IRangesList(normed_points))
  }

  if ("scan_peaks" %in% to_normalize) {
    scan_peaks <- peak_regions$scan_peaks

    normed_peaks <- internal_map$map_function(scan_peaks, function(in_peak){
      peak_index <- seq_len(nrow(in_peak))
      split_peak_by_scan <- split(peak_index, in_peak$scan)

      keep_scans <- names(split_peak_by_scan)[names(split_peak_by_scan) %in% names(split_normalization)]
      in_peak <- in_peak[in_peak$scan %in% as.numeric(keep_scans), ]

      for (iscan in keep_scans) {
        in_peak[in_peak$scan %in% as.numeric(iscan), "Height"] <- exp(log(in_peak[in_peak$scan %in% as.numeric(iscan), "Height"]) - split_normalization[[iscan]])
      }
      in_peak
    })
    peak_regions$scan_peaks <- normed_peaks
  }
  peak_regions$is_normalized <- to_normalize
  peak_regions$normalization_factors <- normalization_factors

  return(peak_regions)
}

get_merged_peak_info <- function(point_data, peak_method = "lm_weighted", min_points = 4){

  point_data <- point_data[point_data$intensity > 0, ]
  point_data$log_int <- log(point_data$intensity + 1e-8)
  weights <- point_data$intensity / max(point_data$intensity)
  mz_peak_info <- get_fitted_peak_info(point_data, use_loc = "mz", w = weights)
  freq_peak_info <- get_fitted_peak_info(point_data, use_loc = "frequency", w = weights)
  mz_peak_info$ObservedMZ <- mz_peak_info$ObservedCenter
  mz_peak_info$ObservedCenter <- NULL
  mz_peak_info$ObservedMZMean <- mean(point_data[, "mz"])
  mz_peak_info$ObservedMZMedian <- median(point_data[, "mz"])
  mz_peak_info$ObservedFrequency <- freq_peak_info$ObservedCenter
  mz_peak_info$ObservedFrequencyMean <- mean(point_data[, "frequency"])
  mz_peak_info$ObservedFrequencyMedian <- median(point_data[, "frequency"])

  mz_peak_info

}


#' characterize peaks from points and picked peaks
#'
#' @param peak_region the PeakRegion object to work on
#'
#' @return list
#' @export
characterize_peaks <- function(peak_region){

  stopifnot(peak_region$is_normalized == "both")

  #peak_ranges <- peak_region$peak_regions
  #picked_peaks <- peak_region$scan_peaks
  #frequency_point_regions <- peak_region$frequency_point_regions
  use_scans <- peak_region$normalization_factors$scan
  peak_region$n_scan <- n_scan <- length(use_scans)
  peak_region$set_min_scan()

  peak_region_list = peak_region$peak_region_list
  peak_nscan = purrr::map_int(peak_region_list, ~ length(base::intersect(use_scans, .x$peaks$scan)))
  keep_peaks = which(peak_nscan >= floor(peak_region$min_scan / 2))
  use_region_list = peak_region_list[keep_peaks]
  log_memory()

  # peak_data <- purrr::imap(use_region_list, function(.x, .y){
  #   message(.y)
  #   characterize_mz_points(.x, peak_scans = use_scans)
  # })
  peak_data <- internal_map$map_function(
    use_region_list, characterize_mz_points, peak_scans = use_scans)
  log_memory()

  individual_peak_heights <- log10(purrr::map_dbl(peak_data, function(x){x$peak_info$Height}))
  scan_peak_heights <- purrr::map(peak_data, function(x){x$scan_data$LogHeight})

  individual_peak_nscan <- purrr::map_int(peak_data, function(x){x$peak_info$NScan})

  # update the minimum number of scans required, and then filter
  # all of the relevant bits
  #### This needs to be updated to account for the fact that we filtered above too.
  passing_data = which(individual_peak_nscan >= peak_region$min_scan)
  peak_region$keep_peaks <- keep_peaks[passing_data]

  peak_data <- peak_data[passing_data]
  individual_peak_heights <- individual_peak_heights[passing_data]
  scan_peak_heights <- scan_peak_heights[passing_data]
  individual_peak_nscan <- individual_peak_nscan[passing_data]
  peak_index <- peak_region$peak_index[passing_data]

  corrected_data <- correct_peak_sd_height(individual_peak_heights,
                                           scan_peak_heights,
                                           individual_peak_nscan,
                                           n_scan)

  template_scan_data <- matrix(NA, nrow = 1, ncol = n_scan)
  colnames(template_scan_data) <- sort(use_scans)

  correction_ratios <- corrected_data$OriginalHeight - corrected_data$CorrectedHeight

  corrected_peak_info <- purrr::map(seq_len(length(peak_data)),
                                      function(in_peak){
      original_data <- peak_data[[in_peak]]
      scan_data <- original_data$scan_data
      corrected_scan_data <- scan_data
      corrected_scan_data$LogHeight <- scan_data$LogHeight - correction_ratios[in_peak]
      peak_info <- original_data$peak_info
      peak_info$CorrectedLog10Height <- corrected_data[in_peak, "CorrectedHeight"]
      peak_info$CorrectedLog10HeightSD <- corrected_data[in_peak, "CorrectedSD"]
      peak_info$PeakID <- peak_index[in_peak]
      rownames(peak_info) <- NULL

      frequency_scans <- mz_scans <- original_scan_heights <- corrected_scan_heights <- template_scan_data
      original_scan_heights[1, as.character(scan_data$Scan)] <- scan_data$LogHeight
      rownames(original_scan_heights) <- peak_index[in_peak]
      corrected_scan_heights[1, as.character(corrected_scan_data$Scan)] <- corrected_scan_data$LogHeight
      rownames(corrected_scan_heights) <- peak_index[in_peak]

      mz_scans[1, as.character(scan_data$Scan)] <- scan_data$ObservedMZ
      rownames(mz_scans) <- peak_index[in_peak]

      frequency_scans[1, as.character(scan_data$Scan)] <- scan_data$ObservedFrequency

      list(peak = peak_info, original_scan = original_scan_heights,
           corrected_scan = corrected_scan_heights,
           mz_scan = mz_scans,
           frequency_scan = frequency_scans)
                                      })

  peak_info <- purrr::map_df(corrected_peak_info, "peak")
  #peak_info <- add_offset(peak_info, peak_region$mz_model)
  peak_info$ScanCorrelation <- peak_region$scan_correlation[peak_region$keep_peaks, "ScanCorrelation"]
  peak_info$HighScan = peak_info$NScan >= quantile(peak_info$NScan, 0.9)

  original_height <- do.call(rbind, purrr::map(corrected_peak_info, "original_scan"))
  corrected_height <- do.call(rbind, purrr::map(corrected_peak_info, "corrected_scan"))
  mz_scan <- do.call(rbind, purrr::map(corrected_peak_info, "mz_scan"))
  frequency_scan <- do.call(rbind, purrr::map(corrected_peak_info, "frequency_scan"))

  peak_region$peak_data <- peak_info
  peak_region$scan_level_arrays <- list(Log10Height = original_height,
                                        CorrectedLog10Height = corrected_height,
                                        ObservedMZ = mz_scan,
                                        ObservedFrequency = frequency_scan,
                                        Scan = colnames(original_height),
                                        PeakID = rownames(original_height))
  log_memory()

  invisible(peak_region)
}

characterize_mz_points <- function(in_region, peak_scans = NULL){
  log_memory()
  in_points = in_region$points
  scan_peaks = in_region$peaks
  if (is.null(peak_scans)) {
    peak_scans <- unique(in_points@elementMetadata$scan)
  }

  # first trim to the scans actually available from the scan peaks
  peak_scans <- base::intersect(peak_scans, scan_peaks$scan)

  if ((nrow(scan_peaks) == 0) || (length(peak_scans) == 0)) {
    peak_info <- data.frame(ObservedMZ = NA,
                            Height = NA,
                            Area = NA,
                            SSR = NA,
                            type = "NA",
                            ObservedMZMean = NA,
                            ObservedMZMedian = NA,
                            ObservedMZSD = NA,
                            Log10ObservedMZSD = NA,
                            Log10Height = NA,
                            HeightSD = NA,
                            Log10HeightSD = NA,
                            Start = NA,
                            Stop = NA,
                            NScan = 0L,
                            NPoint = NA)
    scan_heights <- data.frame(Scan = NA,
                               LogHeight = NA,
                               ObservedMZ = NA)
  } else {
    in_points <- in_points[in_points@elementMetadata$scan %in% peak_scans]

    peak_info <- get_merged_peak_info(as.data.frame(S4Vectors::mcols(in_points)))
    #peak_info$ObservedMZ = peak_info$ObservedCenter
    #peak_info$ObservedCenter = NULL

    peak_info$ObservedMZSD <- sd(scan_peaks$ObservedMZ)
    peak_info$ObservedFrequencySD <- sd(scan_peaks$ObservedFrequency)
    peak_info$Log10ObservedMZSD <- sd(log10(scan_peaks$ObservedMZ))
    peak_info$Log10Height <- log10(peak_info$Height)
    peak_info$HeightSD <- sd(scan_peaks$Height)
    peak_info$Log10HeightSD <- sd(log10(scan_peaks$Height))

    point_data <- as.data.frame(S4Vectors::mcols(in_points))
    peak_start <- min(point_data[point_data$intensity > 0, "mz"])
    peak_stop <- max(point_data[point_data$intensity > 0, "mz"])

    peak_info$Start <- peak_start
    peak_info$Stop <- peak_stop
    peak_info$NScan <- length(peak_scans)

    point_by_scan <- split(point_data, point_data$scan)
    peak_info$NPoint <- median(purrr::map_int(point_by_scan, nrow))

    scan_heights <- data.frame(Scan = scan_peaks$scan, LogHeight = log10(scan_peaks$Height), ObservedMZ = scan_peaks$ObservedMZ,
                               ObservedFrequency = scan_peaks$ObservedFrequency)
  }



  list(peak_info = peak_info, scan_data = scan_heights)
}


characterize_peaks_in_regions <- function(frequency_point_regions, peak_regions, n_scans, max_subsets = 100, peak_index = NULL){
  peak_region_scans <- peak_regions@elementMetadata$X



  peak_data <- dplyr::bind_rows(peak_data)
  peak_data$NScan <- n_scans

  if (is.null(peak_index)) {
    peak_data$PeakID <- seq_len(nrow(peak_data))
  } else {
    peak_data$PeakID <- peak_index
  }
  peak_data
}

characterize_picked_peaks <- function(scan_peaks, n_scans, peak_index = NULL){

  peak_data <- purrr::map_df(scan_peaks, function(in_peak){
    n_point <- purrr::map_int(in_peak$points, length)
    data.frame(ObservedMZ = mean(in_peak$ObservedMZ),
               Height = mean(in_peak$Height),
               Log10Height = mean(log10(in_peak$Height)),
               Area = mean(in_peak$Area),
               Log10Area = mean(log10(in_peak$Area)),
               ObservedMZSD = sd(in_peak$ObservedMZ),
               Log10ObservedMZSD = sd(log10(in_peak$ObservedMZ)),
               HeightSD = sd(in_peak$Height),
               Log10HeightSD = sd(log10(in_peak$Height)),
               AreaSD = sd(in_peak$Area),
               Log10AreaSD = sd(log10(in_peak$Area)),
               Start = min(in_peak$ObservedMZ),
               Stop = max(in_peak$ObservedMZ),
               NPoint = mean(n_point))
  })
  peak_data$NSubset <- n_scans
  peak_data$NScan <- n_scans
  peak_data$PeakID <- seq_len(nrow(peak_data))
  peak_data
}

add_offset <- function(peak_data, mz_model){
  offsets <- purrr::map_df(peak_data$PeakID, function(in_id){
    #print(in_id)
    peak_mz <- peak_data$ObservedMZ[peak_data$PeakID %in% in_id]
    offset_value <- mz_model$y[which.min(abs(peak_mz - mz_model$x))]
    data.frame(PeakID = in_id, Offset = offset_value)
  })
  dplyr::left_join(peak_data, offsets, by = "PeakID")
}

mz_sd_model <- function(in_data){
  min_scan <- max(in_data$NScan)
  fit_data <- in_data[(in_data$NScan >= min_scan) & (in_data$Height > 1), ]
  if ("ScanCorrelated" %in% names(fit_data)) {
    fit_data <- fit_data[!fit_data$ScanCorrelated, ]
  }

  lm_fit <- lm(Log10ObservedMZSD ~ ObservedMZ + CorrectedLog10Height + CorrectedLog10HeightSD + NPoint, data = fit_data)
  sd_pred <- abs(predict(lm_fit, newdata = in_data))
  sd_pred
}

int_sd_model <- function(in_data){
  min_scan <- quantile(in_data$NScan, 0.75)
  in_data$logHeight <- log10(in_data$Height)
  fit_data <- in_data[in_data$NScan >= min_scan, ]
  if ("Ignore" %in% names(fit_data)) {
    fit_data <- fit_data[!fit_data$Ignore, ]
  }
  fit_data <- fit_data[!is.na(fit_data$Log10HeightSD), ]
  fit_data <- fit_data[fit_data$Height > 10, ]
  fit_data$width <- fit_data$Stop - fit_data$Start
  fit_data$weight_2 <- 1 - (fit_data$width / max(fit_data$width))
  sd_fit <- lm(I(1/Log10HeightSD) ~ I(1/logHeight), data = fit_data, weights = fit_data$weight_2)
  sd_pred <- 1/predict(sd_fit, newdata = in_data)
  sd_pred
}


model_sds <- function(values, sds, loess_span = 0.75){
  sd_frame <- data.frame(x = values, y = sds)
  loess_fit <- stats::loess(y ~ x, data = sd_frame, span = loess_span, control = loess.control(surface = "direct"))
  loess_pred <- predict(loess_fit, values)
  loess_pred
}


#' generate bootstrap samples
#'
#' Generates `n_bootstrap` samples, each with `n_sample` entries, from the provided
#' `indices`. See **Details** for more information.
#'
#' @param n_indices the indices to take the bootstrap sample from
#' @param n_bootstrap how many bootstrap samples to generate? (default is 1000)
#' @param n_sample how many items should be in each sample? See **Details**.
#' @param min_indices the minimum number of indices to be willing to work with
#'
#' @details A *bootstrap* sample is a sample where the entries have been sampled
#'  **[with replacement](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))**.
#'  In general, a *bootstrap* sample has the same number of entries as the original
#'  sample. However, there may be cases where *upsampling* to specific number
#'  may be useful. Therefore, the parameter `n_sample` is provided to enable
#'  *upsampling* the original indices to a larger number. If `n_sample` is NULL,
#'  then no *upsampling* will be done.
#'
#' @return list of bootstrap sampled indices
#' @export
bootstrap_samples <- function(n_indices, n_bootstrap = 100, n_sample = NULL, min_indices = 4){
  if (n_indices < min_indices) {
    return(NULL)
  }

  if (is.null(n_sample)) {
    n_sample <- n_indices
  }
  purrr::map(seq_len(n_bootstrap), function(x){
    sample(n_indices, n_sample, replace = TRUE)
  })
}

rename_peak_data = function(data_frame){
  name_convert = c("ObservedCenter.mz" = "ObservedMZ",
                   "Height.mz" = "Height",
                   "Area.mz" = "Area",
                   "ObservedCenter.frequency" = "ObservedFrequency"
  )
  data_names = names(data_frame)
  for (iname in names(name_convert)) {
    name_match = data_names %in% iname
    if (sum(name_match) == 1) {
      data_names[which(name_match)] = name_convert[iname]
    }
  }
  names(data_frame) = data_names
  data_frame
}
