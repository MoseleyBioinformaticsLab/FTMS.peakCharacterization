# Functions related to converting to and from frequency space
#
# Note, frequency space is particularly useful b/c it has the property
# of having a constant difference between subsequent points. Of course
# we then have to worry about getting "to" frequency space, and possibly
# getting "back" to M/Z space.


#' Convert M/Z to Frequency
#'
#' Given a data.frame of m/z, generate frequency values for the data.
#'
#' The **M/Z** values from FTMS data do not have constant spacing between them.
#' This produces challenges in working with ranged intervals and windows. The
#' solution for FTMS data then is to convert them to **frequency** space. This
#' is done by:
#'   * taking subsequent M/Z points
#'   * averaging their M/Z
#'   * taking the difference to get an `offset` value
#'   * dividing averaged M/Z by offset to generate **frequency**
#'   * taking subsequent differences of frequency points
#'   * keep points with a difference in the supplied range as valid for modeling
#'
#' After deciding on the valid points for modeling, each point gets an interpolated
#' frequency value using the two averaged points to the left and right in M/Z.
#'
#' @param mz_data a data.frame with `mz`
#' @param keep_all keep **all** the variables generated, or just the original + **frequency**?
#'
#' @export
#'
#' @seealso mz_scans_to_frequency
#'
#' @return list
convert_mz_frequency = function(mz_data, keep_all = TRUE){
  original_vars = names(mz_data)
  stopifnot("mz" %in% names(mz_data))
  frequency_data = mz_data
  frequency_data = dplyr::mutate(frequency_data, mean_mz = NA,
                                 mean_offset = NA,
                                 mean_frequency = NA,
                                 mean_freq_diff = NA,
                                 convertable = FALSE)

  tmp_data = cbind(frequency_data$mz, dplyr::lag(frequency_data$mz))
  frequency_data$mean_mz = rowMeans(tmp_data)
  frequency_data$mean_offset = tmp_data[, 1] - tmp_data[, 2]
  frequency_data$mean_frequency = frequency_data$mean_mz / frequency_data$mean_offset
  valid_range = discover_frequency_offset(frequency_data$mean_frequency)
  frequency_data$mean_freq_diff = dplyr::lag(frequency_data$mean_frequency) - frequency_data$mean_frequency

  is_convertable = dplyr::between(frequency_data$mean_freq_diff, valid_range$range[1], valid_range$range[2])
  is_convertable[is.na(is_convertable)] = FALSE

  frequency_data[is_convertable, "convertable"] = TRUE

  rownames(frequency_data) = NULL
  frequency_data
}

fit_frequency_linear <- function(x, y, w = NULL){
  #center_x <- x - mean(x)

  X <- matrix(c(rep(1, length(x)), x), nrow = length(x), ncol = 2, byrow = FALSE)

  if (is.null(w)) {
    out_fit <- stats::lm.fit(X, y)
  } else {
    out_fit <- stats::lm.wfit(X, y, w)
  }
  names(out_fit$coefficients) <- NULL
  out_fit
}

predict_frequency_linear = function(data, fit_model){
  data2 = as.matrix(cbind(rep(1, length(data)), data))
  coef_matrix = matrix(fit_model$coefficients, nrow = 2, ncol = 1, byrow = TRUE)
  data2 %*% coef_matrix
}

fit_frequency_sr = function(x, y){
  sr_x = 1 / (x ^ 0.5)
  X = matrix(c(rep(1, length(x)), sr_x), nrow = length(x), ncol = 2, byrow = FALSE)

  sr_fit = stats::lm.fit(X, y)
  names(sr_fit$coefficients) = NULL

  sr_fit
}

predict_frequency_sr = function(x, coefficients){
  sr_x = 1 / (x ^ 0.5)
  X = matrix(c(rep(1, length(x)), sr_x), nrow = length(x), ncol = 2, byrow = FALSE)
  coef_matrix = matrix(coefficients, nrow = 2, ncol = 1, byrow = TRUE)
  out_frequency = X %*% coef_matrix
  out_frequency[, 1]
}

fit_mz_s2 = function(x, y){
  s_x = 1 / (x^2)
  X = matrix(c(rep(1, length(x)), s_x), nrow = length(x), ncol = 2, byrow = FALSE)

  s_fit = stats::lm.fit(X, y)
  names(s_fit$coefficients) = NULL
  s_fit$coefficients
}

predict_mz_s2 = function(x, coefficients){
  s_x = 1 / (x^2)
  X = matrix(c(rep(1, length(x)), s_x), nrow = length(x), ncol = 2, byrow = FALSE)
  coef_matrix = matrix(coefficients, nrow = 2, ncol = 1, byrow = TRUE)
  out_mz = X %*% coef_matrix
  out_mz[, 1]
}



#' convert mz to frequency across scans
#'
#' Given a multi-scan data.frame of m/z, generate frequency values for the data.
#'
#' @param mz_scan_df a data.frame with at least `mz` and `scan` columns
#' @param ... other parameters for `convert_mz_frequency`
#'
#' @seealso convert_mz_frequency
#'
#' @export
#' @return list
mz_scans_to_frequency = function(mz_scan_df, ...){
  split_mz = split(mz_scan_df, mz_scan_df$scan)

  mz_frequency = internal_map$map_function(split_mz, function(in_scan){
    #message(scan_name)
    out_scan = convert_mz_frequency(in_scan, ...)
    out_scan
  })

  frequency_coefficients = internal_map$map_function(mz_frequency, function(in_freq){
    use_peaks = in_freq$convertable
    tmp_fit = fit_frequency_sr(in_freq$mean_mz[use_peaks], in_freq$mean_frequency[use_peaks])$coefficients
    data.frame(intercept = tmp_fit[1], slope = tmp_fit[2], scan = in_freq[1, "scan"])
  })

  frequency_coefficients = do.call(rbind, frequency_coefficients)

  bad_coefficients = boxplot.stats(frequency_coefficients$intercept)$out

  frequency_coefficients = dplyr::filter(frequency_coefficients, !(intercept %in% bad_coefficients))

  model_coefficients = c(median(frequency_coefficients$intercept), median(frequency_coefficients$slope))

  # group them all together and remove any that were outliers based on statistics
  # on the coefficients
  mz_frequency_df = do.call(rbind, mz_frequency)
  mz_frequency_df = mz_frequency_df[mz_frequency_df$scan %in% frequency_coefficients$scan, ]

  convertable_df = mz_frequency_df[mz_frequency_df$convertable, ]
  all_model = fit_frequency_sr(convertable_df$mean_mz, convertable_df$mean_frequency)$coefficients


  mz_frequency_df$frequency = predict_frequency_sr(mz_frequency_df$mz, all_model)

  rownames(mz_frequency_df) = NULL
  valid_range = discover_frequency_offset(mz_frequency_df$frequency)
  mz_frequency_df$frequency_diff = dplyr::lag(mz_frequency_df$frequency) - mz_frequency_df$frequency
  mz_frequency_df$convertable = dplyr::between(mz_frequency_df$frequency_diff, valid_range$range[1], valid_range$range[2])
  mz_frequency_df[is.na(mz_frequency_df$convertable), "convertable"] = FALSE
  list(frequency = mz_frequency_df, all_coefficients = frequency_coefficients, coefficients = all_model, difference_range = valid_range)
}

#' convert mz to frequency using linear fit
#'
#' Given a query, and either two values of M/Z and two values of frequency or
#' a previously generated model, return a data.frame with the predicted value,
#' and the slope and the intercept so the model can be re-used later for
#' other points when needed.
#'
#' @param mz_query the M/Z value to fit
#' @param mz_values two M/Z values
#' @param frequency_values two frequency values
#' @param model a model to use instead of actual values
#'
#' @return data.frame with predicted_value, intercept, and slope
#'
#' @export
mz_frequency_interpolation = function(mz_query, mz_values = NULL, frequency_values = NULL, model = NULL){
  if (is.null(model) && is.null(mz_values) && is.null(frequency_values)) {
    stop("Please provide either a model or mz_values + frequency_values ...")
  }

  if ((!is.null(mz_values) && is.null(frequency_values)) || (is.null(mz_values) && !is.null(frequency_values))) {
    stop("Both mz_values and frequency_values must be provided ...")
  }

  if (!is.null(model) && is.null(mz_values) && is.null(frequency_values)) {
    slope = model$slope
    intercept = model$intercept
  } else {
    if ((length(mz_values) != 2) || (length(frequency_values) != 2)) {
      stop("mz_values and frequency_values MUST contain two entries ...")
    }
    diff_mz = mz_values[2] - mz_values[1]
    diff_frequency = frequency_values[2] - frequency_values[1]
    slope = diff_frequency / diff_mz
    intercept = frequency_values[1] - (slope * mz_values[1])
  }

  pred_frequency = intercept + (slope * mz_query)

  model_prediction = data.frame(predicted_frequency = -1 * pred_frequency,
                                intercept = intercept,
                                slope = slope)
  model_prediction
}



frequency_models_scans = function(scan_values){
  split_scans = split(scan_values, scan_values$scan)

  scan_models = internal_map$map_function(split_scans, function(x){
    frequency_model = create_frequency_model(x)
    frequency_model$scan = x$scan[1]
    frequency_model
  })

  all_scan_models = do.call(rbind, scan_models)
  all_scan_models$frequency = max(all_scan_models$frequency, na.rm = TRUE) - all_scan_models$frequency
  all_scan_models
}


convert_found_peaks = function(frequency_model, found_peaks){

  found_peaks$peak = seq_len(nrow(found_peaks))
  frequency_model = frequency_model[!is.na(frequency_model$frequency), ]
  split_model = split(frequency_model, frequency_model$scan)
  split_peaks = split(found_peaks, found_peaks$scan)

  split_model = split_model[names(split_peaks)]

  out_frequency = purrr::map2_df(split_model, split_peaks, function(in_model, in_peaks){
    tmp_res = purrr::map_df(seq_len(nrow(in_peaks)), function(in_speak){
      start_point = max(which(in_model$mean_mz <= in_peaks[in_speak, "ObservedCenter.mz"]))
      end_point = min(which(in_model$mean_mz >= in_peaks[in_speak, "ObservedCenter.mz"]))
      mz_frequency_interpolation(in_peaks$ObservedCenter.mz[in_speak], in_model$mean_mz[c(start_point, end_point)],
                                 in_model$mean_frequency[c(start_point, end_point)])

    })
    cbind(in_peaks, tmp_res)
  })

  out_frequency$frequency = out_frequency$predicted_frequency
  out_frequency$predicted_frequency = NULL
  cbind(found_peaks, out_frequency)
}

calculate_resolution_information = function(frequency_point_regions, use_mz = 400){
  frequency_points = S4Vectors::mcols(frequency_point_regions)
  valid_range = frequency_point_regions@metadata$difference_range
  med_frequency_difference = valid_range$most_common

  use_freq = predict_frequency_sr(use_mz, frequency_point_regions@metadata$mz_2_frequency)
  one_point_diff = use_freq + med_frequency_difference

  convertable = frequency_points$convertable
  convertable[is.na(convertable)] = FALSE
  mz_model = fit_mz_s2(frequency_points$frequency[convertable], frequency_points$mz[convertable])
  one_point_mz = predict_mz_s2(one_point_diff, mz_model)

  frequency_over_mz = med_frequency_difference / abs(one_point_mz - use_mz)

  list(frequency = list(point_point_differences = valid_range,
                        difference_mz = frequency_over_mz,
                        conversion = frequency_point_regions@metadata$mz_2_frequency),

       mz = list(value = use_mz,
                 ppm = abs(one_point_mz - use_mz) / use_mz * 1e6,
                 difference = abs(one_point_mz - use_mz))
       )
}

discover_frequency_offset = function(frequency_values, cutoff_range = 0.02){
  frequency_diffs = abs(frequency_values - dplyr::lead(frequency_values))

  round_diffs = round(frequency_diffs, digits = 4)
  round_diffs = round_diffs[!is.na(round_diffs)]
  rle_diffs = rle(sort(round_diffs))

  common_value = rle_diffs$values[which.max(rle_diffs$lengths)]

  cutoff_value = cutoff_range * common_value
  range_value = c(common_value - cutoff_value, common_value + cutoff_value)

  useful_points = frequency_diffs[dplyr::between(frequency_diffs, range_value[1], range_value[2])]
  list(most_common = median(useful_points, na.rm = TRUE),
       range = range_value)
}
