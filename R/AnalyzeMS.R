#' analyze ftms mass-spec data
#'
#' This class allows you to analyze mass spec data, and controls the execution
#' of reading in the mass
#'
#' @export
AnalyzeMS <- R6::R6Class("AnalyzeMS",
  public = list(
   load_file = function(){
     self$zip_ms <- ZipMS$new(self$in_file, self$metadata_file, self$out_file, temp_loc = self$temp_loc)
   },
   found_peaks = NULL,

   find_peaks = function(...){
     if ("R6" %in% class(self$peak_finder)) {
       self$peak_finder$raw_data <- self$zip_ms$raw_ms
       self$peak_finder$out_file <- self$zip_ms$out_file
       self$peak_finder$run_correspondence()
       self$found_peaks <- self$peak_finder$export_data()
       self$peak_finder$raw_data <- NULL
       self$peak_finder$multi_scan <- NULL
       peak_finder <- self$peak_finder
       save(peak_finder, file = file.path(self$zip_ms$temp_directory, "peak_finder.rds"))
       scans_to_json(peak_finder, file_output = self$zip_ms$temp_directory)
     } else if ("function" %in% class(self$peak_finder)) {
       self$found_peaks <- self$peak_finder(self$zip_ms$raw_ms, ...)
     }
     self$zip_ms$add_peak_list(self$found_peaks)
   },

   write_zip = function(){
     if (!is.null(self$out_file)) {
       self$zip_ms$write_zip(out_file = self$out_file)
     } else {
       self$zip_ms$write_zip()
     }
   },

   set_peak_finder = function(in_function){
     self$peak_finder <- in_function
   },

   peak_finder = NULL,

   zip_ms = NULL,
   in_file = NULL,
   metadata_file = NULL,
   out_file = NULL,
   temp_loc = NULL,

   run_all = function(){
     self$load_file()
     self$find_peaks()
     self$write_zip()
     self$zip_ms$cleanup()
   },

   initialize = function(in_file, metadata_file = NULL, out_file = NULL, peak_finder = NULL, temp_loc = NULL){
     self$in_file <- in_file

     if (!is.null(metadata_file)) {
       self$metadata_file <- metadata_file
     }

     if (!is.null(out_file)) {
       self$out_file <- out_file
     }


     if (!is.null(peak_finder)) {
       self$peak_finder <- peak_finder
     }

     if (!is.null(temp_loc)) {
       self$temp_loc <- temp_loc
     }
   }
  )
)


#' R6 peak_finder
#'
#' an R6 based peak_finder class that controls execution in a more fine grained
#' manner, and allows the possibility of easily checking each step and saving
#' intermediate outputs when and if desired.
#'
#' @export
PeakFinder <- R6::R6Class("PeakFinder",
  public = list(
    run_time = NULL,
    # Options for multi_scans ----
    raw_data = NULL,
    peak_method = NULL,
    min_points = 4,
    n_peak = Inf,
    flat_cut = 0.98,


    # Options for peak correspondence ----
    mz_range = c(-Inf, Inf),
    digital_resolution_multiplier = 1,
    initial_rmsd_multiplier = 3,

    # functions needed by others
    sd_fit_function = NULL,
    sd_predict_function = NULL,
    offset_predict_function = NULL,
    offset_fit_function = NULL,
    offset_correction_function = NULL,
    noise_function = NULL,

    # control options
    max_iteration = 30,
    scan_fraction = 0.1,
    collapse_peaks = FALSE,
    remove_low_ic_scans = TRUE,
    notify_progress = FALSE,
    max_failures = 5,
    keep_all_master_peak_lists = FALSE,
    keep_intermediates = FALSE,

    # Now the various bits and pieces ----
    raw_filter = NULL,
    apply_raw_filter = function(){
      if (!is.null(self$raw_filter)) {
        self$raw_data <- self$raw_filter(self$raw_data)
      }
    },
    out_file = NULL,
    report_function = NULL,
    create_report = function(){
      if (!is.null(self$report_function)) {
        self$report_function(self)
      }
    },

    vocal = FALSE,


    multi_scan = NULL,
    create_multi_scan = function(){
      if (self$vocal) {
        message("Creating Multi-Scans ....")
      }
      self$multi_scan <- SIRM.FTMS.peakCharacterization::MultiScans$new(self$raw_data, peak_method = self$peak_method,
                                                                        min_points = self$min_points,
                                                                        n_peak = self$n_peak, flat_cut = self$flat_cut,
                                                                        sd_fit_function = self$sd_fit_function,
                                                                        sd_predict_function = self$sd_predict_function)
      invisible(self)
    },
    multi_scan_peaklist = NULL,
    create_multi_scan_peaklist = function(){
      if (self$vocal) {
        message("Creating Multi-Scan Peaklist ....")
      }
      self$multi_scan_peaklist <- SIRM.FTMS.peakCharacterization::MultiScansPeakList$new(self$multi_scan, noise_function = self$noise_function,
                                                                                         sd_predict_function = self$sd_predict_function)
    },

    scan_start = NULL,
    scan_dr_filter = NULL,
    filter_dr_models = function(){
      if (self$vocal) {
        message("Filtering based on digital resolution models ....")
      }
      self$scan_start <- self$multi_scan_peaklist$scan_indices
      self$multi_scan_peaklist$remove_bad_resolution_scans()
      self$scan_dr_filter <- self$multi_scan_peaklist$scan_indices
    },

    correspondent_peaks = NULL,
    create_correspondent_peaks = function(...){
      if (self$vocal) {
        message("Peak Correspondence ....")
      }

      use_peaklist <- self$multi_scan_peaklist$clone(deep = TRUE)

      self$correspondent_peaks <-
        SIRM.FTMS.peakCharacterization::FindCorrespondenceScans$new(use_peaklist,
                                                                    peak_calc_type = self$peak_type,
                                                                    max_iteration = self$max_iteration,
                                                                    digital_resolution_multiplier = self$digital_resolution_multiplier,
                                                                    rmsd_multiplier = self$initial_rmsd_multiplier,
                                                                    max_failures = self$max_failures,
                                                                    mz_range = self$mz_range,
                                                                    notify_progress = self$vocal,
                                                                    noise_function = self$noise_function,
                                                                    sd_fit_function = self$sd_fit_function,
                                                                    sd_predict_function = self$sd_predict_function,
                                                                    offset_fit_function = self$offset_fit_function,
                                                                    offset_predict_function = self$offset_predict_function,
                                                                    offset_correction_function = self$offset_correction_function,
                                                                    collapse_peaks = self$collapse_peaks,
                                                                    keep_all_master_peak_lists = self$keep_all_master_peak_lists,
                                                                    keep_intermediates = self$keep_intermediates
                                                                    )
      self$correspondent_peaks$iterative_correspondence()
    },

    scan_information_content = NULL,

    collapse_correspondent_peaks = function(){
      if (self$vocal) {
        message("Collapsing correspondent peaks ....")
      }
      self$correspondent_peaks$master_peak_list <- collapse_correspondent_peaks(self$correspondent_peaks$master_peak_list)
    },

    scan_normalized = NULL,
    normalize_scans_by_correspondent_peaks = function(intensity_measure = "Height", summary_function = mean){
      if (self$vocal) {
        message("Normalizing scans ....")
      }
      mpl <- self$correspondent_peaks$master_peak_list

      intensity_options <- c(Height = "scan_height", Area = "scan_area", NormalizedArea = "scan_normalizedarea")
      intensity_internal <- intensity_options[intensity_measure]
      n_corresponding <- mpl$count_notna()

      normalizing_peaks <- n_corresponding >= quantile(n_corresponding, 0.95)

      n_scan <- ncol(mpl$scan_mz)

      peak_intensities <- log(mpl[[intensity_internal]][normalizing_peaks, ])

      # need to check that there are at least 25 correspondent peaks in each of
      # the scans, if not we will drop the scan and not bother to normalize it
      n_peak_scan <- vapply(seq(1, ncol(peak_intensities)), function(in_scan){
        sum(!is.na(peak_intensities[, in_scan]))
      }, numeric(1))

      keep_scans <- n_peak_scan >= 25

      if (sum(keep_scans) != n_scan){
        self$multi_scan_peaklist <- self$multi_scan_peaklist$reorder(keep_scans)
        mpl <- mpl$reorder(keep_scans)
        peak_intensities <- peak_intensities[, keep_scans]
        n_scan <- sum(keep_scans)
      }

      self$scan_normalized <- self$multi_scan_peaklist$scan_indices

      # this is getting the difference of a scan to all the other scans.
      # It uses the fact that the peaks are the rows, and the scans are the columns.
      # Take out the scan you want, make it a matrix that will be same size as the other.
      # Remove that scan from the full matrix.
      # Take their differences, this is the difference of that scan to all other scans
      # for all of the peaks.
      # Average the difference in each scan across the peaks
      # Sum them after to find the total difference
      scan_diffs <- vapply(seq(1, n_scan), function(in_scan){
        scan_peaks <- peak_intensities[, in_scan, drop = FALSE]
        scan_peaks_matrix <- matrix(scan_peaks, nrow = nrow(scan_peaks), ncol = n_scan - 1)

        other_matrix <- peak_intensities[, -in_scan, drop = FALSE]
        scan_other_diff <- scan_peaks_matrix - other_matrix
        peak_means <- colMeans(scan_other_diff, na.rm = TRUE)
        sum(peak_means)
      }, numeric(1))

      normalize_scan <- which.min(abs(scan_diffs))

      scan_norm_matrix <- matrix(peak_intensities[, normalize_scan, drop = FALSE],
                                 nrow = nrow(peak_intensities), ncol = n_scan)

      diff_matrix <- peak_intensities - scan_norm_matrix
      normalization_factors <- apply(diff_matrix, 2, summary_function, na.rm = TRUE)
      normalization_matrix <- matrix(normalization_factors, nrow = nrow(mpl[[intensity_internal]]),
                                     ncol = ncol(mpl[[intensity_internal]]), byrow = TRUE)

      # we normalize both the height and the area with the same factors to see what
      # makes the biggest difference
      mpl$scan_height <- exp(log(mpl$scan_height) - normalization_matrix)
      mpl$scan_area <- exp(log(mpl$scan_area) - normalization_matrix)
      mpl$scan_normalizedarea <- exp(log(mpl$scan_normalizedarea) - normalization_matrix)

      mpl$is_normalized <- TRUE
      mpl$normalization_factors <- normalization_factors
      mpl$normalized_by <- intensity_measure
      self$correspondent_peaks$master_peak_list <- mpl
      invisible(self)
    },

    intermediates = FALSE,
    save_intermediates = function(filename = NULL){
      if (self$intermediates) {
        if (is.null(filename)) {
          filename <- paste0(self$raw_data$raw_metadata$run$id, ".RData")
        }
        peakfinder <- self
        peakfinder$multi_scan <- NULL
        save(peakfinder, file = filename, compress = FALSE)
      }

    },

    peak_data = NULL,
    get_mz = function(scan_mz, sd_model){
      mean_mz <- mean(scan_mz)
      median_mz <- median(scan_mz)
      sd_mz <- sd(scan_mz)

      model_sd = self$sd_predict_function(sd_model, mean_mz)
      values = scan_mz

      list(Mean = mean_mz,
           Median = median_mz,
           SD = sd_mz,
           ModelSD = model_sd,
           Values = values)
    },
    get_height_area = function(scan_height_area){
      mean_h <- mean(scan_height_area)
      median_h <- median(scan_height_area)
      sd_h <- sd(scan_height_area)
      rsd_h <- sd_h / mean_h
      values <- scan_height_area

      list(Mean = mean_h,
           Median = median_h,
           SD = sd_h,
           RSD = rsd_h,
           Values = values)
    },
    create_peak_data = function(){
      if (self$vocal) {
        message("Generating the final peak data ....")
      }

      master_peaks <- self$correspondent_peaks$master_peak_list
      if (is.null(master_peaks$sd_model)) {
        master_peaks$calculate_sd_model()
      }
      sd_model <- master_peaks$sd_model
      n_peak <- length(master_peaks$master)

      peak_data <- lapply(seq(1, n_peak), function(in_peak){
        tmp_index <- !is.na(master_peaks$scan_mz[in_peak, ])

        mz <- self$get_mz(master_peaks$scan_mz[in_peak, tmp_index], sd_model)
        height <- self$get_height_area(master_peaks$scan_height[in_peak, tmp_index])
        area <- self$get_height_area(master_peaks$scan_area[in_peak, tmp_index])
        norm_area <- self$get_height_area(master_peaks$scan_normalizedarea[in_peak, tmp_index])

        list(N = sum(tmp_index),
             Scans = master_peaks$scan[tmp_index],
             ObservedMZ = mz,
             Height = height,
             Area = area,
             NormalizedArea = norm_area)
      })
      self$peak_data <- peak_data
    },
    processing_info = NULL,
    create_processing_info = function(){
      function_pkg <- "package:SIRM.FTMS.peakCharacterization"
      pkg_description <- utils::packageDescription(substring(function_pkg, 9))

      if (!is.null(pkg_description$RemoteSha)) {
        pkg_sha <- pkg_description$RemoteSha
      } else {
        pkg_sha <- ""
      }

      document_peakfinder <- as.list(PeakFinder$new(self$peak_method, self$noise_function,
                                                    self$raw_filter, self$create_report))
      document_peakfinder[[".__enclos_env__"]] <- NULL
      document_peakfinder$clone <- NULL

      dr_model <- self$correspondent_peaks$sd_models[[1]]

      if (class(dr_model) == "loess") {
        dr_model <- self$remove_loess_class(dr_model)
      }

      sd_model <- self$correspondent_peaks$sd_models[[length(self$correspondent_peaks$sd_models)]]

      if (class(sd_model) == "loess") {
        sd_model <- self$remove_loess_class(sd_model)
      }

      self$processing_info <- list(Package = function_pkg,
                                   Version = pkg_description$Version,
                                   Sha = pkg_sha,
                                   FunctionCalled = document_peakfinder,
                                   Parameters = list(Method = self$peak_method,
                                                     Scans = self$raw_data$scan_range),
                                   Models = list(DigitalResolutionModel = dr_model,
                                                 SDModel = sd_model)
      )
    },

    remove_loess_class = function(sd_model){
      attr(sd_model, "class") <- NULL
      sd_model$call <- NULL
      sd_model$terms <- NULL
      sd_model
    },

    run_correspondence = function(){
      # only do this stuff if needed, otherwise, start from the already saved
      # data.
      tictoc::tic()
      if (is.null(self$multi_scan) & (!is.null(self$raw_data))) {
        self$apply_raw_filter()
        self$create_multi_scan()
        self$create_multi_scan_peaklist()
      } else if (is.null(self$multi_scan_peaklist)) {
        stop("Need a MultiScanPeakList to work with!", call. = TRUE)
      }
      self$filter_dr_models()
      self$create_correspondent_peaks(median_corrected = FALSE)
      self$collapse_correspondent_peaks()
      self$normalize_scans_by_correspondent_peaks()
      self$save_intermediates()
      self$create_report()
      self$create_peak_data()
      self$create_processing_info()
      tmp_time <- tictoc::toc()
      self$run_time <- tmp_time$toc - tmp_time$tic
      invisible(self)
    },

    export_data = function(){
      SIRM.FTMS.peakCharacterization::PeakPickingAnalysis$new(self$peak_data, self$processing_info)
    },

    initialize = function(peak_method = "lm_weighted", noise_function = noise_sorted_peaklist, raw_filter = NULL,
                          report_function = NULL, intermediates = FALSE, sd_fit_function = NULL,
                          sd_predict_function = NULL, offset_fit_function = NULL, offset_predict_function = NULL,
                          offset_correction_function = NULL){
      if (!is.null(peak_method)) {
        self$peak_method <- peak_method
      }

      if (!is.null(noise_function)) {
        self$noise_function = noise_function
      }

      if (!is.null(raw_filter)) {
        self$raw_filter <- raw_filter
      }

      if (!is.null(report_function)) {
        self$report_function <- report_function
      }

      if (!is.null(sd_fit_function)) {
        self$sd_fit_function <- sd_fit_function
      } else {
        self$sd_fit_function <- default_sd_fit_function
      }

      if (!is.null(sd_predict_function)) {
        self$sd_predict_function <- sd_predict_function
      } else {
        self$sd_predict_function <- default_sd_predict_function
      }

      if (!is.null(offset_fit_function)) {
        self$offset_fit_function <- offset_fit_function
      } else {
        self$offset_fit_function <- default_offset_fit_function
      }

      if (!is.null(offset_predict_function)) {
        self$offset_predict_function <- offset_predict_function
      } else {
        self$offset_predict_function <- default_offset_predict_function
      }

      if (!is.null(offset_correction_function)) {
        self$offset_correction_function <- offset_correction_function
      } else {
        self$offset_correction_function <- default_correct_offset_function
      }

      self$intermediates <- intermediates

      invisible(self)
    }
  )
)

default_sd_fit_function <- function(x, y){
  loess_frame <- data.frame(x = x, y = y)
  loess_fit <- stats::loess(y ~ x, data = loess_frame, control = stats::loess.control(surface = "direct"))
  loess_fit
}

default_sd_predict_function <- function(model, x){
  stats:::predict.loess(model, newdata = x)
}


#' peak finding and reporting
#'
#' Given a RawMS object, actually does peak finding at the scan level and then
#' peak correspondence.
#'
#' @return list
#' @export
peak_finder <- function(raw_data, method = "lm_weighted", noise_function = noise_sorted_peaklist){
  # example data:
  # load("zip_ms_example.RData")
  # raw_data <- zip_ms$raw_ms
  multi_scan <- SIRM.FTMS.peakCharacterization::MultiScans$new(raw_data, peak_method = method)
  multi_scan_peak_list <- SIRM.FTMS.peakCharacterization::MultiScansPeakList$new(multi_scan, noise_function = noise_function)

  correspondent_peaks <- SIRM.FTMS.peakCharacterization::FindCorrespondenceScans$new(multi_scan_peak_list, multiplier = 3)
  correspondent_peaks$master_peak_list <- SIRM.FTMS.peakCharacterization::normalize_scans(correspondent_peaks$master_peak_list)

  get_mz <- function(scan_mz, sd_model){
    mean_mz <- mean(scan_mz)
    median_mz <- median(scan_mz)
    sd_mz <- sd(scan_mz)

    model_sd = exponential_predict(sd_model, mean_mz)[1]
    values = scan_mz

    list(Mean = mean_mz,
         Median = median_mz,
         SD = sd_mz,
         ModelSD = model_sd,
         Values = values)
  }

  get_height_area <- function(scan_height_area){
    mean_h <- mean(scan_height_area)
    median_h <- median(scan_height_area)
    sd_h <- sd(scan_height_area)
    rsd_h <- sd_h / mean_h
    values <- scan_height_area

    list(Mean = mean_h,
         Median = median_h,
         SD = sd_h,
         RSD = rsd_h,
         Values = values)
  }



  create_peak_data <- function(correspondent_peaks){
    sd_model <- correspondent_peaks$sd_models[[1]] # grab the digital resolution model
    master_peaks <- correspondent_peaks$master_peak_list
    n_peak <- length(master_peaks$master)

    peak_data <- lapply(seq(1, n_peak), function(in_peak){
      tmp_index <- !is.na(master_peaks$scan_mz[in_peak, ])

      mz <- get_mz(master_peaks$scan_mz[in_peak, tmp_index], sd_model)
      height <- get_height_area(master_peaks$scan_height[in_peak, tmp_index])
      area <- get_height_area(master_peaks$scan_area[in_peak, tmp_index])
      norm_area <- get_height_area(master_peaks$scan_area[in_peak, tmp_index] / mz$ModelSD)

      list(N = sum(tmp_index),
           ObservedMZ = mz,
           Height = height,
           Area = area,
           NormalizedArea = norm_area)
    })
    peak_data
  }

  peak_data <- create_peak_data(correspondent_peaks)

  function_call <- "peak_finder"
  function_pkg <- find(function_call)
  pkg_description <- utils::packageDescription(substring(function_pkg, 9))

  if (!is.null(pkg_description$RemoteSha)) {
    pkg_sha <- pkg_description$RemoteSha
  } else {
    pkg_sha <- ""
  }

  processing_meta <- list(Package = function_pkg,
                          Version = pkg_description$Version,
                          Sha = pkg_sha,
                          FunctionCalled = peak_finder,
                          Parameters = list(Method = method,
                                            Scans = raw_data$scan_range),
                          Models = list(DigitalResolutionModel = correspondent_peaks$sd_models[[1]],
                                        SDModel = correspondent_peaks$sd_models[[length(correspondent_peaks$sd_models)]])
  )

  PeakPickingAnalysis$new(peak_data, processing_meta)
}


#' scans to json
#'
#' write the scan level peaks to a JSON format
#'
#' @param peak_finder the peak_finder object that has everything
#' @param exclude_noise should noise peaks be excluded?
#' @param file_output where should the json go?
#'
#' @export
#' @return character
scans_to_json <- function(peak_finder, exclude_noise = TRUE, file_output = NULL){
  sd_model <- peak_finder$correspondent_peaks$master_peak_list$sd_model
  sd_pred_function <- peak_finder$correspondent_peaks$master_peak_list$sd_predict_function

  peak_lists <- peak_finder$multi_scan_peaklist$get_scan_peak_lists()
  peak_frames <- purrr::map_df(peak_lists, function(in_list){
    out_frame <- in_list$peak_list
    out_frame$scan <- in_list$scan
    if (exclude_noise) {
      out_frame <- out_frame[out_frame$not_noise, ]
      out_frame$noise_rank <- 0
    } else {
      noise_frame <- out_frame[!out_frame$not_noise, ]
      nonoise_frame <- out_frame[out_frame$not_noise, ]

      noise_frame <- noise_frame[order(noise_frame$Height), ]
      noise_frame$noise_rank <- seq(nrow(noise_frame), 1, -1)

      nonoise_frame$noise_rank <- 0
      out_frame <- rbind(nonoise_frame, noise_frame)
    }
    out_frame
  })
  Peaks <- purrr::map(seq(1, nrow(peak_frames)), function(in_row){
    tmp_data <- peak_frames[in_row, ]
    list(N = 1,
         Scans = tmp_data$scan,
         ObservedMZ = list(Mean = tmp_data$ObservedMZ,
                           Median = tmp_data$ObservedMZ,
                           SD = sd_pred_function(sd_model, tmp_data$ObservedMZ),
                           ModelSD = sd_pred_function(sd_model, tmp_data$ObservedMZ),
                           Values = tmp_data$ObservedMZ),
         Height = list(Mean = tmp_data$Height,
                       Median = tmp_data$Height,
                       SD = NA,
                       RSD = NA,
                       Values = tmp_data$Height),
         Area = list(Mean = tmp_data$Area,
                     Median = tmp_data$Area,
                     SD = NA,
                     RSD = NA,
                     Values = tmp_data$Area),
         NormalizedArea = list(Mean = tmp_data$NormalizedArea,
                               Median = tmp_data$NormalizedArea,
                               SD = NA,
                               RSD = NA,
                               Values = tmp_data$NormalizedArea),
         NoiseInfo = list(IsNoise = !tmp_data$not_noise,
                          NoiseRank = tmp_data$noise_rank)
    )
  })

  if (!is.null(file_output)) {
    scans_json <- peak_list_2_json(list(Peaks = Peaks))
    cat(scans_json, file = file.path(file_output, "scans_peaklist.json"))
  }
  invisible(scans_json)
}


#' add scans to file
#'
#' adds the jsonified scan peaks to an already existing zip file
#'
#' @param zip_file the zip file to work with
#' @param out_file if desired, a new file name to generate
#'
#' @export
#' @return NULL
add_scans_to_file <- function(zip_file, out_file = zip_file, exclude_noise = TRUE){
  zip_data <- zip_ms(zip_file, out_file = out_file, load_raw = TRUE, load_peak_list = FALSE)

  tmp_env <- new.env()
  if (file.exists(file.path(zip_data$temp_directory, "peak_finder.rds"))) {
    load(file.path(zip_data$temp_directory, "peak_finder.rds"), envir = tmp_env)
    scans_to_json(tmp_env$peak_finder, exclude_noise = exclude_noise, file_output = zip_data$temp_directory)
    zip_data$write_zip()
    zip_data$cleanup()

  } else {
    warning("peak_finder.rds not present!")
  }
}
