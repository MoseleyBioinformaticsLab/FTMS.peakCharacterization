#' log memory usage
#'
#' Logs the amount of memory being used to a log file if it is available, and
#' generating warnings if the amount of RAM hits zero.
#'
#' @export
#' @return NULL
log_memory = function(){
  if (get("memory", envir = has_logger)) {
    linux_memory = system("cat /proc/meminfo", intern = TRUE)
    linux_memory = grep("^MemTotal|^MemAvailable|^Active|^SwapTotal|^SwapFree", linux_memory, value = TRUE)
    linux_memory = grep("anon|active|file", linux_memory, value = TRUE, invert = TRUE)

    memory_values = stringr::str_extract(linux_memory, "[[:digit:]].*")
    memory_numbers = as.numeric(stringr::str_extract(memory_values, "[[:digit:]].* "))
    memory_ids = stringr::str_extract(linux_memory, "^[[:alpha:]]+")
    names(memory_numbers) = memory_ids
    memory_string = paste0("Memory: ", paste(paste(c("Total: ", "Available: ", "Active: ", "SwapTotal: ", "SwapFree: "), memory_values, sep = ""), collapse = ", "))

    active_to_total = memory_numbers["Active"] / memory_numbers["MemTotal"]
    swapfree_to_swap = memory_numbers["SwapFree"] / memory_numbers["SwapTotal"]
    if (is.nan(swapfree_to_swap)) {
      swapfree_to_swap = 1
    }

    swapfree_to_swap = memory_numbers["SwapFree"] / memory_numbers["SwapTotal"]
    if (is.nan(swapfree_to_swap)) {
      swapfree_to_swap = 1
    }


    if ((active_to_total >= 0.95) || (swapfree_to_swap <= 0.95)) {
      memory_string2 = paste0("HIGH MEMORY USAGE!!! ", memory_string)
      if (get("logger", envir = has_logger)) {
        logger::log_warn(memory_string2, namespace = "FTMS.peakCharacterization")
      } else {
        warning(memory_string2)
      }
    } else {
      if (get("logger", envir = has_logger)) {
        logger::log_info(memory_string, namespace = "FTMS.peakCharacterization")
      }
    }
  }

}

#' log messages
#'
#' If a log_appender is available, logs the given message at the `info` level.
#'
#' @param message_string the string to put in the message
#'
#' @export
#' @return NULL
log_message = function(message_string){
  if (get("logger", envir = has_logger)) {
    logger::log_info(message_string, namespace = "FTMS.peakCharacterization")
  }
  if (get("status", envir = pc_progress)) {
    message(message_string)
  }
}

#' turn logging off
#'
#' There may be good reasons to turn the logging off (running tests, etc). This
#' basically tells the package that the logger isn't available
#'
#' @export
#' @return NULL
disable_logging = function(){
  assign("logger", FALSE, envir = has_logger)
  message("Logging disabled.")
}
