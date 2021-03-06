#' import json
#'
#' import json from a file correctly given some things where things get written
#' differently
#'
#' @param json_file the json file to read
#' @return list
import_json <- function(json_file){
  json_data <- jsonlite::fromJSON(json_file, simplifyVector = FALSE)
  if (length(json_data) == 1) {
    out_list <- json_data[[1]]
  } else {
    out_list <- json_data
  }
  out_list
}


#' get raw metadata
#'
#' When raw files are copied, we also generated metadata about their original locations
#' and new locations, and some other useful info. We would like to capture it, and
#' keep it along with the metadata from the mzml file. So, given a list of mzml
#' files, and a location for the raw files, this function creates metadata json
#' files for the mzml files.
#'
#' @param mzml_files the paths to the mzml files
#' @param raw_file_loc the directory holding raw files and json metadata files
#'
#' @importFrom purrr map_lgl
#' @importFrom R.utils isAbsolutePath getAbsolutePath
#'
#' @export
raw_metadata_mzml <- function(mzml_files, raw_file_loc, recursive = TRUE){
  # mzml_files <- dir("/home/rmflight/data/test_json_meta/mzml_data", full.names = TRUE)
  # raw_file_loc <- "/home/rmflight/data/test_json_meta"
  # recursive <- TRUE
  #
  if (!isAbsolutePath(mzml_files[1])) {
    mzml_files <- getAbsolutePath(mzml_files)
    names(mzml_files) <- NULL
  }


  raw_json_files <- dir(raw_file_loc, pattern = "json$", full.names = TRUE, recursive = recursive)
  if (!isAbsolutePath(raw_json_files[1])) {
    raw_json_files <- getAbsolutePath(raw_json_files)
    names(raw_json_files) <- NULL
  }

  raw_json_data <- data.frame(json_file = raw_json_files, id = basename_no_file_ext(raw_json_files), stringsAsFactors = FALSE)
  mzml_data <- data.frame(mzml_file = mzml_files, id = basename_no_file_ext(mzml_files), stringsAsFactors = FALSE)
  mzml_data <- mzml_data[file.exists(mzml_files), ]


  json_mzml_match <- dplyr::inner_join(raw_json_data, mzml_data, by = "id")
  json_mzml_match$mzml_meta <- FALSE

  did_write_mzml_meta <- purrr::map_lgl(seq(1, nrow(json_mzml_match)), function(in_row){
    raw_meta <- import_json(json_mzml_match[in_row, "json_file"])
    #print(json_mzml_match[in_row, "mzml_file"])
    mzml_meta <- try(get_mzml_metadata(json_mzml_match[in_row, "mzml_file"]))
    if (!inherits(mzml_meta, "try-error")) {
      source_file_data <- mzml_meta$fileDescription$sourceFileList$sourceFile
      tmp_model <- as.character(mzml_meta$referenceableParamGroupList$referenceableParamGroup[[1]]$name)
      tmp_serial <- as.character(mzml_meta$referenceableParamGroupList$referenceableParamGroup[[2]]$value)
      tmp_sha1 <- as.character(mzml_meta$fileDescription$sourceFileList$sourceFile[[3]]$value)
      mzml_meta$run$instrument <- list(model = tmp_model,
                                       serial = tmp_serial)
      if (tmp_sha1 == raw_meta$sha1) {
        mzml_file <- json_mzml_match[in_row, "mzml_file"]
        mzml_meta_file <- list(file = basename(mzml_file),
                               saved_path = mzml_file,
                               sha1 = digest::digest(mzml_file, algo = "sha1", file = TRUE))
        mzml_meta$file <- list(raw = raw_meta,
                               mzml = mzml_meta_file)

        outfile <- paste0(tools::file_path_sans_ext(json_mzml_match[in_row, "mzml_file"]), ".json")
        cat(jsonlite::toJSON(mzml_meta, pretty = TRUE, auto_unbox = TRUE), file = outfile)
        did_write <- TRUE
      } else {
        warning("SHA-1 of Files does not match! Not writing JSON metadata!")
        did_write <- FALSE
      }


    } else {
      did_write <- FALSE
    }
    did_write
  })
  json_mzml_match$mzml_meta <- did_write_mzml_meta
  json_mzml_match
}

basename_no_file_ext <- function(in_files){
  file_no_ext <- tools::file_path_sans_ext(in_files)
  basename(file_no_ext)
}

zip_ms_from_zip <- function(in_file){
  ZipMS$new(in_file)
  ZipMS
}

zip_ms_from_mzml <- function(in_file, out_dir){
  message("creating zip file from mzML and populating raw metadata")
  zip_file <- mzml_to_zip(in_file, out_dir)
  if (!is.null(zip_file)) {
    zip_ms_from_zip(zip_file)
  }
}

#' Represents the zip mass spec file
#'
#' This reference class represents the zip mass spec file. It does this by
#' providing objects for the zip file, the metadata, as well as various bits
#' underneath such as the raw data and peak lists, and their
#' associated metadata. Although it is possible to work with the ZipMS object directly, it
#' is heavily recommended to use the AnalyzeMS object
#' for carrying out the various steps of an analysis, including peak finding.
#'
#' @section `ZipMS` Methods:
#'
#'  use `?method-name` to see more details about each individual method
#'
#'  \describe{
#'   \item{`zip_ms`}{make a new `ZipMS`}
#'   \item{`show_temp_dir`}{show where files are stored}
#'   \item{`save`}{save the file}
#'   \item{`cleanup`}{unlink the `temp_directory`}
#'   \item{`add_peak_list`}{add a `Peaks` to the data}
#'  }
#'
#'
#' @section `ZipMS` Data Members:
#'  \describe{
#'    \item{`zip_file`}{the zip file that was read in}
#'    \item{`metadata`}{the actual metadata for the file}
#'    \item{`metadata_file`}{the metadata file}
#'    \item{`raw_ms`}{a `RawMS` holding the raw data}
#'    \item{`peaks`}{a `Peaks` holding the peak analysis}
#'    \item{`id`}{the sample id}
#'    \item{`out_file`}{the file where data will be saved}
#'  }
#'
#' @seealso AnalyzeMS
#' @return ZipMS object
#'
"ZipMS"


#' make a new ZipMS
#'
#' @param in_file the file to use (either .zip or .mzML)
#' @param mzml_meta_file metadata file (.json)
#' @param out_file the file to save to at the end
#' @param load_raw logical to load the raw data
#' @param load_peak_list to load the peak list if it exists
#'
#' @export
#' @return ZipMS
zip_ms <- function(in_file, mzml_meta_file = NULL, out_file = NULL, load_raw = TRUE,
                   load_peak_list = TRUE){
  ZipMS$new(in_file, mzml_meta_file = mzml_meta_file, out_file = out_file, load_raw = load_raw, load_peak_list = load_peak_list)
}

#' ZipMS - save
#'
#' @name save
#' @param out_file the file to save to
#'
#' @details `out_file`, if it is `NULL`, will be taken from when the
#'  object was generated, and by default will be set to the same as the `in_file`.
#'  If not `NULL`, then it is checked that the `id` is part of the
#'  `out_file`, and if not, the `id` is added to the actual file name.
#'
#' @examples
#' \dontrun{
#'  new_ms <- zip_ms("in_file")
#'  new_ms$save()
#'  new_ms$save("out_file")
#' }
NULL

#' ZipMS - show_temp_dir
#'
#' shows where the temp directory `ZipMS` is using is
#'
#' @name show_temp_dir
#' @usage ZipMS$show_temp_dir()
#'
NULL

#' ZipMS - cleanup
#'
#' cleans up after things are done
#'
#' @name cleanup
#' @usage ZipMS$cleanup()
#'
NULL

#' ZipMS - add_peak_list
#'
#' adds a peak list to the `ZipMS`
#'
#' @name add_peak_list
#' @usage ZipMS$add_peak_list(peak_list_data)
#' @param peak_list_data a [Peaks()] object
#'
NULL

#' @export
ZipMS <- R6::R6Class("ZipMS",
  public = list(
    zip_file = NULL,
    zip_metadata = NULL,
    metadata = NULL,
    metadata_file = NULL,
    raw_ms = NULL,
    peaks = NULL,
    peak_finder = NULL,
    json_summary = NULL,
    id = NULL,
    out_file = NULL,
    temp_directory = NULL,

    load_raw = function(){
      self$raw_ms <- RawMS$new(file.path(self$temp_directory, self$metadata$raw$raw_data),
                file.path(self$temp_directory, self$metadata$raw$metadata))

    },

    load_peak_finder = function(){
      if (file.exists(file.path(self$temp_directory, "peak_finder.rds"))) {
        peak_finder <- try(readRDS(file.path(self$temp_directory, "peak_finder.rds")))

        if (inherits(peak_finder, "try-error")) {
          peak_finder <- try({
            tmp_env <- new.env()
            load(file.path(self$temp_directory, "peak_finder.rds"), envir = tmp_env)
            tmp_env$peak_finder
          })
        }
        if (inherits(peak_finder, "PeakRegionFinder")) {
          self$peak_finder <- peak_finder
          rm(peak_finder)
          message("Peak Finder Binary File Loaded!")
        } else {
          stop("peak_finder.rds is not valid!")
        }

      }
    },

    save_json = function(){
      lists_2_json(self$json_summary, temp_dir = self$temp_directory)
    },

    save_peak_finder = function(){
      peak_finder <- self$peak_finder
      saveRDS(peak_finder, file.path(self$temp_directory, "peak_finder.rds"))
      invisible(self)
    },

    load_peak_list = function(){
      if (file.exists(file.path(self$temp_directory, "peak_finder.rds"))) {
        tmp_env <- new.env()
        load(file.path(self$temp_directory, "peak_finder.rds"), envir = tmp_env)
        peak_data <- tmp_env$peak_finder$correspondent_peaks$master_peak_list$clone(deep = TRUE)
        rm(tmp_env)
      } else {
        warning("No peak_finder.rds found, not returning peaks!")
        peak_data <- NULL
      }
      peak_data
    },

    compare_raw_corresponded_densities = function(mz_range = c(150, 1600), window = 1, delta = 0.1){
      if (!is.null(self$raw_ms)) {
        raw_peak_mz <- raw_peaks(self$raw_ms)
        raw_peak_density <- calculate_density(raw_peak_mz, use_range = mz_range, window = window, delta = delta)
        raw_peak_density$type <- "raw"
      } else {
        warning("No raw data to get peaks from!")
        raw_peak_density <- data.frame(window = NA, density = NA, type = "raw", stringsAsFactors = FALSE)
      }
      if (!is.null(self$peaks)) {
        correspondent_peak_mz <- self$peaks$master
        correspondent_peak_density <- calculate_density(correspondent_peak_mz, use_range = mz_range, window = window, delta = delta)
        correspondent_peak_density$type <- "correspondent"
      } else {
        warning("No correspondent peaks to get peaks from!")
        correspondent_peak_density <- data.frame(window = NA, density = NA, type = "correspondent", stringsAsFactors = FALSE)
      }
      peak_densities <- rbind(raw_peak_density, correspondent_peak_density)
      peak_densities$type <- forcats::fct_relevel(peak_densities$type, "raw", "correspondent")

      peak_densities
    },

    initialize = function(in_file, mzml_meta_file = NULL, out_file = NULL, load_raw = TRUE,
                          load_peak_list = TRUE,
                          temp_loc = NULL){
      private$do_load_raw <- load_raw
      private$do_load_peak_list <- load_peak_list

      if (is.null(temp_loc)) {
        temp_loc <- tempfile(pattern = "zipms_tmp")
      } else {
        temp_loc <- tempfile(pattern = "zipms_tmp", tmpdir = temp_loc)
      }

      dir.create(temp_loc)
      self$temp_directory <- temp_loc

      in_file <- path.expand(in_file)
      is_zip <- regexpr("*.zip", in_file)
      if (is_zip != -1) {
        in_zip <- in_file
        self$zip_file <- in_zip
        unzip(in_zip, exdir = self$temp_directory)

      } else {
        file.copy(in_file, file.path(self$temp_directory, basename(in_file)))
        if (!is.null(mzml_meta_file)) {
          file.copy(mzml_meta_file, file.path(self$temp_directory, basename(mzml_meta_file)))
        }
        initialize_zip_metadata(self$temp_directory)
        self$zip_file <- in_file
      }

      get_zip_raw_metdata(self)

      check_zip_file(self$temp_directory)

      self$metadata_file <- "metadata.json"
      self$metadata <- load_metadata(self$temp_directory, self$metadata_file)
      self$id <- self$metadata$id

      if (load_raw && (!is.null(self$metadata$raw$raw_data))) {
        self$raw_ms <- self$load_raw()
      }

      if (load_peak_list && (!is.null(self$metadata$peakpicking_analysis$output))) {
        self$peaks <- self$load_peak_list()
      }

      private$calc_md5_hashes()

      self$out_file <- private$generate_filename(out_file)

      invisible(self)
    },

    show_temp_dir = function(){
      print(self$temp_directory)
    },

    write_zip = function(out_file = NULL){
      if (is.null(out_file)) {
        out_file <- self$out_file
      } else {
        out_file <- private$generate_filename(out_file)
        self$out_file <- out_file
      }
      zip(out_file, list.files(self$temp_directory, full.names = TRUE), flags = "-jq")
      write_zip_file_metadata(self)
    },

    cleanup = function(){
      unlink(self$temp_directory, recursive = TRUE, force = TRUE)
      #file.remove(self$temp_directory)
    },

    finalize = function(){
      unlink(self$temp_directory, recursive = TRUE)
    },

    add_peak_list = function(peak_list_data){
      json_peak_meta <- jsonlite::toJSON(peak_list_data$peakpicking_parameters,
                                         pretty = TRUE, auto_unbox = TRUE)
      cat(json_peak_meta, file = file.path(self$temp_directory,
                                      "peakpicking_parameters.json"))
      self$metadata$peakpicking_analysis <- list(parameters =
                                                   "peakpicking_parameters.json",
                                                 output = "raw_peaklist.json")

      json_meta <- jsonlite::toJSON(self$metadata, pretty = TRUE, auto_unbox = TRUE)
      cat(json_meta, file = file.path(self$temp_directory,
                                      self$metadata_file))

      json_peaklist <- peak_list_2_json(peak_list_data$peak_list)
      cat(json_peaklist, file = file.path(self$temp_directory,
                                          "raw_peaklist.json"))

      self$peaks <- peak_list_data
    }
  ),
  private = list(
    generate_filename = function(out_file = NULL){

      is_zip_out <- regexpr("*.zip", self$zip_file)

      if (!is.null(out_file)) {

        out_file <- path.expand(out_file)
        #has_id <- regexpr(self$id, out_file)
        is_zip_out <- regexpr("*.zip", out_file)

        # if (has_id == -1) {
        #   out_file <- paste0(self$id, "_", out_file)
        # }

        if (is_zip_out == -1) {
          out_file <- paste0(tools::file_path_sans_ext(out_file), ".zip")
        }

      } else {
        out_file <- paste0(tools::file_path_sans_ext(self$zip_file), ".zip")
      }
      out_file
    },



    do_load_raw = NULL,
    do_load_peak_list = NULL,

    curr_md5 = list(metadata_file = numeric(0),
                           raw_metadata_file = numeric(0),
                           raw_data_file = numeric(0),
                           peaks_metadata_file = numeric(0),
                           peaks_data_file = numeric(0)),
    old_md5 = NULL,

    calc_md5_hashes = function(){

      if (!is.null(self$metadata_file)) {
        private$curr_md5$metadata_file <- tools::md5sum(file.path(self$temp_directory, self$metadata_file))
      }

      if (!is.null(self$raw_ms)) {
        private$curr_md5$raw_metadata_file <-
          tools::md5sum(file.path(self$temp_directory, self$metadata$raw$metadata))

        private$curr_md5$raw_data_file <-
          tools::md5sum(file.path(self$temp_directory, self$metadata$raw$raw_data))
      }

      private$old_md5 <- private$curr_md5
    }


  )
)

get_zip_raw_metdata <- function(zip_obj){
  zip_file_path <- dirname(zip_obj$zip_file)
  zip_file <- basename_no_file_ext(zip_obj$zip_file)
  json_file <- file.path(zip_file_path, paste0(zip_file, ".json"))

  # this first case should actually happen *almost* all the time, as the instantiation
  # of the zip container entails copying the meta-data (if present) into the raw
  # _metadata file and putting it into the temp directory that is the proxy of
  # our zip file
  if (file.exists(file.path(zip_obj$temp_directory, "raw_metadata.json"))) {
    file.path(zip_obj$temp_directory, "raw_metadata.json")
    raw_metadata <- import_json(file.path(zip_obj$temp_directory, "raw_metadata.json"))


    if (!is.null(raw_metadata$file)) {
      file_metadata <- raw_metadata$file
    } else {
      file_metadata <- list()
    }
  } else if (file.exists(json_file)) {
    json_metadata <- import_json(json_file)

    if (!is.null(json_metadata$file)) {
      file_metadata <- json_metadata$file
    } else {
      file_metadata <- list()
    }
  } else {
    file_metadata <- list()
  }

  zip_obj$zip_metadata <- file_metadata

  zip_obj

}

write_zip_file_metadata <- function(zip_obj){
  zip_metadata <- zip_obj$zip_metadata

  if (!is.null(zip_obj$raw_ms$ms_info)) {
    raw_ms_info <- zip_obj$raw_ms$ms_info
  } else {
    raw_ms_info <- NULL
  }

  if (!is.null(zip_obj$peak_finder$peak_meta)) {
    peak_info <- zip_obj$peak_finder$peak_meta()
  } else {
    peak_info <- NULL
  }

  if (file.exists(zip_obj$out_file)) {
    sha1 <- digest::digest(zip_obj$out_file, algo = "sha1", file = TRUE)

    zip_file_metadata <- list(file = basename(zip_obj$out_file),
                              saved_path = zip_obj$out_file,
                              sha1 = sha1)

    json_loc <- paste0(tools::file_path_sans_ext(zip_obj$out_file), ".json")

    zip_metadata$zip <- zip_file_metadata
    zip_metadata$raw <- raw_ms_info
    zip_metadata$peak <- peak_info
    cat(jsonlite::toJSON(zip_metadata, pretty = TRUE, auto_unbox = TRUE), file = json_loc)

  } else {
    stop("File path does not exist, cannot write JSON metadata!")
  }
}

#' plot raw peaks
#'
#' Given the raw_ms object, generate a plot of the peaks that are generated from
#' averaging the scans in the currently set scan-range.
#'
#' @param raw_ms either a RawMS object, or data.frame of mz / intensity
#' @param scan_range which scans to use
#' @param mz_range the range of M/z's to consider
#' @param transform apply a transform to the data
#'
#' @return ggplot2 object
#' @export
#'
plot_raw_peaks <- function(raw_ms, scan_range = NULL, mz_range = NULL, transform = NULL) {
  if (inherits(raw_ms, "RawMS")) {
    if (is.null(scan_range)) {
      scan_range <- raw_ms$scan_range
    }
    peaks <- raw_peak_intensity(raw_ms, scanrange = scan_range)

  } else if (inherits(raw_ms, "data.frame")) {
    peaks <- raw_ms

  }
  if (!is.null(mz_range)) {
    peaks <- peaks[(peaks$mz >= mz_range[1]) & (peaks$mz <= mz_range[2]), ]
  }

  if (!is.null(transform)) {
    peaks$intensity <- transform(peaks$intensity + 1)
  }

  ggplot(peaks, aes(x = mz, xend = mz, y = 0, yend = intensity)) + geom_segment() +
    labs(x = "M/Z", y = "Intensity")
}

#' raw peaks intensities
#'
#' generate raw peaks from just averaging scans via xcms
#'
#' @param raw_ms an RawMS object
#' @param scanrange the range of scans to use, default is derived from raw_ms
#'
#' @export
#' @importFrom xcms getSpec
#' @importFrom pracma findpeaks
#' @return data.frame of mz and intensity
#'
raw_peak_intensity <- function(raw_ms, scanrange = raw_ms$scan_range) {
  mean_scan <- xcms::getSpec(raw_ms$raw_data, scanrange = scanrange)
  mean_scan <- mean_scan[!is.na(mean_scan[, 2]), ]
  mean_peaks <- pracma::findpeaks(mean_scan[, 2], nups = 2)

  mean_peak_intensity <- data.frame(mz = mean_scan[mean_peaks[, 2], 1],
                             intensity = mean_scan[mean_peaks[, 2], 2])
  mean_peak_intensity
}


#' raw peaks
#'
#' generate raw peaks from just averaging scans via xcms
#'
#' @param raw_ms an RawMS object
#' @param scanrange the range of scans to use, default is derived from raw_ms
#'
#' @export
#' @importFrom xcms getSpec
#' @importFrom pracma findpeaks
#' @return numeric
#'
raw_peaks <- function(raw_ms, scanrange = raw_ms$scan_range) {
  mean_scan <- xcms::getSpec(raw_ms$raw_data, scanrange = scanrange)
  mean_scan <- mean_scan[!is.na(mean_scan[, 2]), ]
  mean_peaks <- pracma::findpeaks(mean_scan[, 2], nups = 2)

  mean_peak_mz <- mean_scan[mean_peaks[, 2], 1]
  mean_peak_mz
}

#' peak density
#'
#' calculates peak density in a sliding window of m/z
#'
#' @param peak_data numeric vector of m/z values representing peaks
#' @param use_range the range of m/z's to generate windows
#' @param window how big an m/z window to calculate density
#' @param delta how much to slide the window
#'
#' @return data.frame
#' @export
#'
calculate_density <- function(peak_data, use_range = NULL, window = 1, delta = 0.1){
  peak_data <- sort(peak_data, decreasing = FALSE)
  if (is.null(use_range)) {
    use_range <- range(peak_data)
  }

  window_locs <- data.frame(beg = round(seq(use_range[1], use_range[2] - window, delta), 2),
                            end = round(seq(use_range[1] + window, use_range[2], delta), 2))

  peak_density <- purrr::map_df(seq(1, nrow(window_locs)), function(in_window){
    out_val <- sum((peak_data >= window_locs[in_window, "beg"]) & (peak_data <= window_locs[in_window, "end"]))

    data.frame(window = (window_locs[in_window, "beg"] + window_locs[in_window, "end"]) / 2,
               density = out_val)
  })
  peak_density
}

#' determine sample run time
#'
#' @param zip the zip object you want to use
#' @param units what units should the run time be in? (s, m, h)
#'
#' @export
#' @return data.frame with sample, start and end time
sample_run_time = function(zip, units = "m"){
  if (inherits(zip, "character")) {
    zip = zip_ms(zip)
    cleanup = TRUE
  }
  if (is.null(zip$raw_ms)) {
    zip$load_raw()
    cleanup = FALSE
  } else {
    cleanup = FALSE
  }

  ms_data = get_ms_info(zip$raw_ms$raw_data, include_msn = TRUE, include_precursor = TRUE)
  ms_data = ms_data[order(ms_data$time), ]
  # assume that the last scan-scan time difference is how long the last scan should have taken as well
  last_diff = ms_data$time[nrow(ms_data)] - ms_data$time[nrow(ms_data) - 1]
  total_time = ms_data$time[nrow(ms_data)] + last_diff
  start_time = lubridate::as_datetime(zip$raw_ms$raw_metadata$run$startTimeStamp)
  end_time = start_time + total_time
  total_time_out = switch(units,
                          s = total_time,
                          m = total_time / 60,
                          h = total_time / 3600)
  if (cleanup) {
    zip$cleanup()
  }
  data.frame(sample = zip$id, start = start_time, run = total_time_out, end = end_time)
}

#' get frequencies for all data
#'
#' Given a zip object, gets frequency conversions for *all* of the scans, not just the MS1
#' or filtered scans.
#'
#' @param zip the zip object
#'
#' @export
#' @return data.frame
every_scan_frequency = function(zip){
  if (inherits(zip, "character")) {
    zip = zip_ms(zip)
    cleanup = TRUE
  }

  if (is.null(zip$raw_ms)) {
    zip$load_raw()
    cleanup = FALSE
  } else {
    cleanup = FALSE
  }

  scan_data = get_ms_info(zip$raw_ms$raw_data, include_msn = TRUE, include_precursor = TRUE)
  scan_data = scan_data[order(scan_data$time), ]

  mz_data = purrr::map(seq(1, nrow(scan_data)), function(scan_row){
    if (!is.na(scan_data[scan_row, "scan"])) {
      tmp = as.data.frame(xcms::getScan(zip$raw_ms$raw_data, scan_data[scan_row, "scan"]),
                          strinsAsFactors = FALSE)
    }  else {
      tmp = as.data.frame(xcms::getMsnScan(zip$raw_ms$raw_data, scan_data[scan_row, "scan_msn"]),
                          stringsAsFactors = FALSE)
    }
    tmp$scan = scan_data[scan_row, "acquisition"]
    tmp
  })

  if (cleanup){
    zip$cleanup()
  }

  ... = NULL
  frequency_fit_description = c(0, -1/2, -1/3)
  mz_fit_description = c(0, -1, -2, -3)
  mz_frequency = purrr::map(mz_data, function(in_scan){
    #message(scan_name)
    out_scan = FTMS.peakCharacterization::convert_mz_frequency(in_scan, ...)
    out_scan
  })
  frequency_fits = purrr::map(mz_frequency, function(in_freq){
    use_peaks = in_freq$convertable
    tmp_fit = FTMS.peakCharacterization:::fit_exponentials(in_freq$mean_mz[use_peaks], in_freq$mean_frequency[use_peaks], frequency_fit_description)
    tmp_fit$scan = in_freq[1, "scan"]
    tmp_fit
  })
  frequency_coefficients = purrr::map_df(frequency_fits, function(.x){
      tmp_df = as.data.frame(matrix(.x$coefficients, nrow = 1))
      names(tmp_df) = c("intercept", "sqrt", "cubert")
      tmp_df$scan = .x$scan
      tmp_df
  })

  coefficients = dplyr::left_join(scan_data, frequency_coefficients, by = c("acquisition" = "scan"))
  list(mz_data = mz_frequency, coefficients = coefficients)
}
