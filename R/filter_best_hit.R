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

# For every value of `group_col`: the hit with the smallest e-value; of several, the one with the longest
# alignment; of several, the first one in the table. Returns what
# dplyr::do(dplyr::group_by(hit_tbl, <group_col>), filter_best_hits(.)) returns (a grouped tibble, one row per
# group), without calling an R function for every group.
best_hits <- function(hit_tbl, group_col = "query_id") {
        if (nrow(hit_tbl) == 0 || anyNA(hit_tbl$evalue) || anyNA(hit_tbl$alig_length) || !is.numeric(hit_tbl$alig_length))
                return(best_hits_do(hit_tbl, group_col))
        rows <- first_row_per_group(hit_tbl[[group_col]], list(hit_tbl$evalue, -hit_tbl$alig_length))
        res <- tibble::as_tibble(hit_tbl)[rows, ]
        dplyr::group_by(res, .data[[group_col]])
}

# The implementation used until 0.4.3 (one call of filter_best_hits() per group); kept for tables best_hits()
# does not handle (missing e-values) and as the reference in the tests.
best_hits_do <- function(hit_tbl, group_col = "query_id") {
        '.' <- NULL
        dplyr::do(dplyr::group_by(hit_tbl, .data[[group_col]]), filter_best_hits(.))
}
