#' @title Helper function to select best BLAST hit based on minimum evalue
#' @details Extracts the best hit based on lowest e-value and longest splice variant 
#' when e-values are equal
#' @param x a tibble storing gene locus ids, splice variant ids, and blast output for filtering
#' @author Hajk-Georg Drost
filter_best_hits <- function(x) {
        
        if (nrow(x) == 0)
                stop("When trying to retrieve the best hit the column ",x," was empty and could not be processed. Please remove this line from the dataset, before extracting best hits.", call. = FALSE)
        min_val <- min(x$evalue)
        evalue <- alig_length <- NULL
        res <- dplyr::filter(x, evalue == min_val)
        if (nrow(res) > 1) {
                max_len <- max(res$alig_length)
                res <- dplyr::filter(res, alig_length == max_len)
        }
        
        if (nrow(res) > 1) 
                res <- dplyr::slice(res, 1)
        
        return(res)
}

# The rows filter_best_hits() selects for every group, for the whole table at once: the indices of the rows, one
# per group, in the order in which dplyr::do(dplyr::group_by(tbl, group), ...) returns the groups.
#   group : the grouping values (one per row)
#   keys  : list of vectors; within a group the row with the smallest values wins (first key first)
# Ties are broken as before: by the position of the row in the table (the first one wins).
first_row_per_group <- function(group, keys = list()) {
        idx <- seq_along(group)
        ord <- do.call(order, c(list(group), keys, list(idx), list(method = "radix", na.last = TRUE)))
        first <- ord[!duplicated(group[ord])]
        # the order of the groups is taken from dplyr itself
        g <- NULL
        group_order <- dplyr::group_keys(dplyr::group_by(tibble::tibble(g = group), g))$g
        first[match(group_order, group[first])]
}

# The best hit of every value of `group_col` (a grouped tibble, one row per group, groups in the order of
# dplyr::group_by()).
#
#   rule = "bitscore" (default since 0.5.0): the hit with the highest bit score; of several, the one with the
#          smallest e-value; then the longest alignment; then the highest percent identity; then the
#          alphabetically first subject_id (then query_id). The last keys make the choice deterministic: it does not
#          depend on the order of the rows of the hit table. Of several rows of the chosen subject that tie in all
#          of this, the one with the highest raw score, then the smallest values of the remaining columns.
#   rule = "evalue" (the rule of orthologr <= 0.4.3): the hit with the smallest e-value; of several, the one with
#          the longest alignment; of several, the first row of the hit table (which is sorted by query_id and
#          subject_id, i.e. the alphabetically first subject). The bit score is not used. DIAMOND and BLAST report
#          an e-value of 0 for every strong hit of a long protein, so that for such proteins the alignment length
#          decided - and a paralog with a longer but worse alignment could win. This is what
#          dplyr::do(dplyr::group_by(hit_tbl, <group_col>), filter_best_hits(.)) returned and is reproduced exactly.
best_hits <- function(hit_tbl, group_col = "query_id", rule = NULL) {
        rule <- orthologr_best_hit_rule(rule)
        # no hits: no best hits (a table without rows, grouped as usual)
        if (nrow(hit_tbl) == 0)
                return(dplyr::group_by(tibble::as_tibble(hit_tbl), .data[[group_col]]))
        if (rule == "bitscore") {
                need <- c("bit_score", "evalue", "alig_length", "perc_identity")
                if (!all(is.element(need, names(hit_tbl))))
                        stop("best_hit_rule = \"bitscore\" needs the columns '", paste(need, collapse = "', '"),
                             "' in the hit table (missing: '", paste(setdiff(need, names(hit_tbl)), collapse = "', '"),
                             "'). Use best_hit_rule = \"evalue\" for tables without them.", call. = FALSE)
                if (anyNA(hit_tbl$bit_score))
                        stop("best_hit_rule = \"bitscore\": the hit table has missing bit scores.", call. = FALSE)
                last <- lapply(intersect(c("subject_id", "query_id"), names(hit_tbl)), function(col) as.character(hit_tbl[[col]]))
                keys <- c(list(-as.numeric(hit_tbl$bit_score), as.numeric(hit_tbl$evalue), -as.numeric(hit_tbl$alig_length),
                               -as.numeric(hit_tbl$perc_identity)), last)
                # Several rows of one query and one subject (several HSPs; BLAST tables) can tie in all these keys:
                # which of them is kept (its q_start, s_start, ...) then depended on the order of the rows. The
                # raw score, which resolves ties hidden by the rounding of printed bit scores, and then every
                # other column of the table decide now. These keys come after subject_id, so they never change
                # which subject is the best hit - only which of its rows represents it.
                rest <- setdiff(names(hit_tbl), c("bit_score", "evalue", "alig_length", "perc_identity", "subject_id",
                                                  "query_id", "score_raw"))
                row_keys <- c(
                        if (is.element("score_raw", names(hit_tbl))) list(-as.numeric(hit_tbl$score_raw)),
                        lapply(rest, function(col) {
                                v <- hit_tbl[[col]]
                                if (is.numeric(v)) v else as.character(v)
                        })
                )
                rows <- first_row_per_group(hit_tbl[[group_col]], c(keys, row_keys))
                res <- tibble::as_tibble(hit_tbl)[rows, ]
                return(dplyr::group_by(res, .data[[group_col]]))
        }
        if (anyNA(hit_tbl$evalue) || anyNA(hit_tbl$alig_length) || !is.numeric(hit_tbl$alig_length))
                return(best_hits_do(hit_tbl, group_col))
        rows <- first_row_per_group(hit_tbl[[group_col]], list(hit_tbl$evalue, -hit_tbl$alig_length))
        res <- tibble::as_tibble(hit_tbl)[rows, ]
        dplyr::group_by(res, .data[[group_col]])
}

# The rule by which best hits are chosen: the argument when given, else options(orthologr.best_hit_rule = ),
# else "evalue" with options(orthologr.legacy = TRUE), else "bitscore".
orthologr_best_hit_rule <- function(rule = NULL) {
        if (is.null(rule))
                rule <- orthologr_option("orthologr.best_hit_rule", "bitscore")
        if (!is.character(rule) || length(rule) < 1 || is.na(rule[1]))
                orthologr_abort(class = "bad_argument", "best_hit_rule must be \"bitscore\" or \"evalue\".", call. = FALSE)
        # (the default of the exported functions is the vector of both choices)
        rule <- rule[1]
        if (!is.element(rule, c("bitscore", "evalue")))
                orthologr_abort(class = "bad_argument", "best_hit_rule = '", rule, "' is not supported. Please choose one of: 'bitscore', 'evalue'.", call. = FALSE)
        rule
}

# The implementation used until 0.4.3 (one call of filter_best_hits() per group); kept for tables best_hits()
# does not handle (missing e-values) and as the reference in the tests.
best_hits_do <- function(hit_tbl, group_col = "query_id") {
        '.' <- NULL
        dplyr::do(dplyr::group_by(hit_tbl, .data[[group_col]]), filter_best_hits(.))
}
