#' @title Filter dNdS values
#' @description This function takes the output data.table returned by \code{\link{dNdS}} and
#' filters the output by the following criteria:
#' 
#' 1) all dN values having an NA value are omitted
#' 
#' 2) all dS values having an NA value are omitted
#' 
#' 3) all dNdS values >= the specified \code{dnds.threshold} are omitted
#' 
#' @param dNdS_tbl a \code{data.table} returned by \code{\link{dNdS}}.
#' @param dnds.threshold a numeric value specifying the dnds threshold for genes that shall be retained.
#' Hence all genes having a dNdS value <= \code{dnds.threshold} are retained. Default is \code{dnds.threshold} = 2.
#' @param dS_zero what to do with gene pairs whose dS is 0, so that dN/dS is not defined (\code{dNdS = NA}): gene
#' pairs without any substitution (dN = dS = 0, e.g. identical coding sequences) and gene pairs with nonsynonymous
#' substitutions only (dN > 0, dS = 0).
#' \itemize{
#' \item \code{"drop"} (default): they are left out, as in every released version of orthologr (there because
#' \code{NA <= dnds.threshold} is \code{NA}); a message now says how many they are and how to keep them.
#' \item \code{"keep"}: they are kept, with \code{dNdS = NA}, and the logical column \code{dS_zero} marks
#' them (\code{FALSE} for all other rows): the result then has more rows and one more column.
#' }
#' The default can be set for the session with \code{options(orthologr.dS_zero = "keep")}. (During the development
#' of 0.5.0 the default was \code{"keep"} for a while; it was set back to \code{"drop"} before release so that
#' results and the shape of the table stay as they were for existing code.)
#' @return The rows of \code{dNdS_tbl} that pass the filter (with the column \code{dS_zero} when
#' \code{dS_zero = "keep"}). A message gives the number of rows left out for each reason.
#' @details
#' 
#' The dNdS ratio quantifies the selection pressure acting on a given protein sequence.
#' 
#' It is proposed that:
#' 
#'  1) dNdS values < 1 reflect stabilizing selection of a protein
#'  
#'  2) dNdS value = 1 reflect neutral selection of a protein
#'  
#'  3) dNdS values > 1 reflect variational selection of a protein
#' 
#' Now an assumption must be generated to allow for dNdS filtering.
#' Given a \code{dnds.threshold} all value above this threshold are omitted from the dataset.
#' 
#' @author Hajk-Georg Drost
#' @examples
#' # a dN/dS table with a pair without dS, a pair with dS = 0 and a pair above the threshold
#' dNdS_tbl <- data.frame(query_id   = paste0("q", 1:5), subject_id = paste0("s", 1:5),
#'                        dN   = c(0.010, 0.020, 0.000, 0.300, 0.050),
#'                        dS   = c(0.100, NA,    0.000, 0.100, 0.200),
#'                        dNdS = c(0.100, NA,    NA,    3.000, 0.250))
#' filter_dNdS(dNdS_tbl, dnds.threshold = 2)
#'
#' # keep the gene pairs with dS = 0 (dNdS = NA, column dS_zero = TRUE)
#' filter_dNdS(dNdS_tbl, dnds.threshold = 2, dS_zero = "keep")
#'
#' @examplesIf all(nzchar(Sys.which(c("diamond", "perl"))))
#' # on the result of dNdS()
#' res <- dNdS(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' filter_dNdS(res, dnds.threshold = 2)
#' @seealso \code{\link{dNdS}}, \code{\link{divergence_stratigraphy}}
#' @export
filter_dNdS <- function(dNdS_tbl, dnds.threshold = 2, dS_zero = getOption("orthologr.dS_zero", "drop")){
        dN <- dS <- NULL
        if (!is.character(dS_zero) || length(dS_zero) != 1 || !is.element(dS_zero, c("keep", "drop")))
                stop("dS_zero must be \"keep\" or \"drop\".", call. = FALSE)
        message("Filtering out NA values in dN or dS and all values with dNdS > ", dnds.threshold, " ...")
        message("Initial input contains ", nrow(dNdS_tbl), " rows.")
        
        # Gene pairs without synonymous divergence: dS = 0, so that dN/dS is not defined (NA).
        why <- dnds_filter_reason(dNdS_tbl, dnds.threshold)
        zero <- !is.na(why) & startsWith(why, "dS = 0")
        n_identical <- sum(zero & dNdS_tbl$dN == 0)
        n_dn_only <- sum(zero & dNdS_tbl$dN > 0)
        
        if (dS_zero == "drop") {
                # as up to 0.4.3: 'NA <= threshold' is NA, and such rows were left out
                res <- dplyr::filter(dNdS_tbl,!is.na(dN), !is.na(dS), dNdS <= dnds.threshold)
        } else {
                keep <- is.na(why) | zero
                res <- dNdS_tbl[keep, , drop = FALSE]
                res$dNdS[zero[keep]] <- NA_real_
                res$dS_zero <- zero[keep]
        }
        if (any(zero))
                message(sum(zero), " gene pairs have dS = 0, so that their dN/dS is not defined: ", n_identical,
                        " with dN = 0 as well (no substitutions between the aligned coding sequences) and ", n_dn_only,
                        " with dN > 0. ",
                        if (dS_zero == "keep")
                                "They are kept with dNdS = NA and dS_zero = TRUE (dS_zero = \"drop\", the default, leaves them out)."
                        else
                                paste0("They are left out, as in earlier versions of orthologr (dS_zero = \"drop\"); ",
                                       "dS_zero = \"keep\" or options(orthologr.dS_zero = \"keep\") keeps them with dNdS = NA."))
        other <- table(why[!is.na(why) & !zero])
        if (length(other) > 0)
                message("Left out: ", paste0(other, " gene pair", ifelse(other != 1, "s", ""), " with ", names(other), collapse = "; "), ".")
        message("Filtering done. New output table contains ", nrow(res), " rows.")
        return(res)
}

# Why filter_dNdS() leaves a row out (NA: it passes): one of
#   "dN or dS not available (NA)", "dS = 0 and dN = 0", "dS = 0 and dN > 0", "dNdS not available (NA)",
#   "dNdS above the threshold"
dnds_filter_reason <- function(dNdS_tbl, dnds.threshold) {
        dN <- dNdS_tbl$dN
        dS <- dNdS_tbl$dS
        r <- dNdS_tbl$dNdS
        why <- rep(NA_character_, nrow(dNdS_tbl))
        na_val <- is.na(dN) | is.na(dS)
        zero <- !na_val & dS == 0
        why[na_val] <- "dN or dS not available (NA)"
        why[zero & dN == 0] <- "dS = 0 and dN = 0"
        why[zero & dN > 0] <- "dS = 0 and dN > 0"
        rest <- !na_val & !zero
        why[rest & is.na(r)] <- "dNdS not available (NA)"
        why[rest & !is.na(r) & r > dnds.threshold] <- "dNdS above the threshold"
        why
}
