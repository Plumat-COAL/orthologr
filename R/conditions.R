# Classed conditions. The errors and warnings of orthologr that a caller may want to handle have the class
# "orthologr_error" / "orthologr_warning" and one or two classes for their kind (see ?orthologr_options, section
# "Errors and warnings"), so that they can be caught without matching message text:
#   tryCatch(dNdS(...), orthologr_error_tool_missing = function(e) ...)
# The message and the absence of a call in it are the same as with stop(..., call. = FALSE).
#   ...   : the parts of the message, pasted together as by stop()
#   class : the kinds, without the prefix, most specific first, e.g. c("alignment_failed", "tool_failed")
orthologr_abort <- function(..., class = NULL, call. = FALSE) {
        cls <- c(if (length(class)) paste0("orthologr_error_", class), "orthologr_error")
        stop(errorCondition(paste(unlist(lapply(list(...), as.character)), collapse = ""), class = cls, call = NULL))
}

orthologr_warn <- function(..., class = NULL, call. = FALSE) {
        cls <- c(if (length(class)) paste0("orthologr_warning_", class), "orthologr_warning")
        warning(warningCondition(paste(unlist(lapply(list(...), as.character)), collapse = ""), class = cls, call = NULL))
}

# The kinds of condition orthologr signals (for the documentation and the tests).
orthologr_condition_classes <- function() {
        list(error   = c("tool_missing", "package_missing", "tool_failed", "alignment_failed", "no_hits",
                         "bad_input", "bad_argument"),
             warning = c("failed_pairs", "isoforms", "legacy_engine"))
}
