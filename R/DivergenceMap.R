#' @title Sort dNdS Values Into Divergence Strata
#' @description This function takes a data.table returned by dNdS
#' and sorts the corresponding dNdS value into divergence strata (10-quantile or deciles by default).
#' @param dNdS_tbl a data.table object returned by \code{\link{dNdS}}.
#' @param subject.id a logical value indicating whether \code{query_id} AND \code{subject_id} should be returned.
#' @param n_quantile a numeric value specifying the number of quantiles that should be returned (a whole number >= 1).
#' @details 
#' 
#' \strong{How values are assigned to strata.} Stratum \code{k} holds the dNdS values in \code{(q[k-1], q[k]]},
#' where \code{q[0], ..., q[n]} are the quantiles (\code{stats::quantile()}, type 7) for the probabilities
#' \code{0, 1/n, ..., 1}; the smallest value belongs to stratum 1. Rows with \code{dNdS = NA} get the stratum
#' \code{NA}.
#' 
#' \strong{Ties.} When many values are equal - typically \code{dNdS = 0} for more than \code{1/n} of the genes of
#' two closely related genomes - some of the thresholds coincide. Earlier versions stopped with an error in that
#' case ("at least one duplicate threshold"). The same rule is now applied through ranks: with \code{N} values and
#' \code{r} the smallest rank of a value (\code{rank(ties.method = "min")}), its stratum is
#' \code{max(1, ceiling((r - 1) * n / (N - 1)))}. Without ties this is exactly the quantile class above, so results
#' that could be computed before do not change. Equal values always share one stratum (the lowest one their ranks
#' reach); the strata such a block spans beyond that stay empty, so that fewer than \code{n_quantile} strata are
#' used and the strata differ in size. A warning names the empty strata. Example: 1,000 genes, 280 of them with
#' \code{dNdS = 0}, \code{n_quantile = 10}: the 280 genes are stratum 1, stratum 2 is empty, the next gene (rank
#' 281) is in stratum 3.
#' 
#' The input table is not modified.
#' 
#' 
#' Divergence Strata are typically decile (10-quantile) values of corresponding \code{\link{dNdS}} values.
#' The \code{\link{dNdS}} function returns dNdS values for orthologous genes
#' of a query species (versus subject species). These dNdS values are then
#' sorted into deciles and each orthologous protein coding gene of the
#' query species receives a corresponding decile value instead of the initial dNdS value.
#' 
#' This allows a better comparison between Phylostrata and Divergence Strata (for more details see package: \pkg{myTAI}).
#' 
#' It is also possible to specify other number of quantile (N-quantile) value such as quintile (5-quantile) 
#' rather than decile values.
#' 
#' @author Hajk-Georg Drost
#' @examples
#' # a dN/dS table as dNdS() returns it (first columns; 20 gene pairs)
#' set.seed(1)
#' dNdS_tbl <- data.frame(query_id   = sprintf("AT1G%05d.1", seq(1010, by = 10, length.out = 20)),
#'                        subject_id = sprintf("%d|PACid:%d", 470000 + 1:20, 16030000 + 1:20),
#'                        dN = round(stats::runif(20, 0.005, 0.1), 4),
#'                        dS = round(stats::runif(20, 0.1, 0.3), 4))
#' dNdS_tbl$dNdS <- round(dNdS_tbl$dN / dNdS_tbl$dS, 4)
#'
#' # the divergence map: the divergence stratum (decile of dN/dS) of every query gene
#' divMap <- divergence_map( dNdS_tbl )
#' head(divMap)
#'
#' # in case you want the subject_id as well, you can set
#' # the argument subject.id = TRUE
#' divMap <- divergence_map( dNdS_tbl = dNdS_tbl, subject.id = TRUE)
#'
#' # in case you want a divergence map with divergence stratum as quintile (5-quantile) values
#' # or any other N-quantile values rather than the default decile (10-quantile) values,
#' # you can specify this with n_quantile.
#' divMap <- divergence_map( dNdS_tbl = dNdS_tbl, subject.id = TRUE, n_quantile = 5)
#' table(divMap$divergence_strata)
#' @seealso \code{\link{divergence_stratigraphy}}
#' @return a data.table storing a standard divergence map.
#' @export

divergence_map <- function(dNdS_tbl, subject.id = FALSE, n_quantile = 10){
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- subject_id <- divergence_strata <- NULL
        
        if (!is.data.frame(dNdS_tbl) || !all(is.element(c("dNdS", "query_id", "subject_id"), names(dNdS_tbl))))
                stop("dNdS_tbl must be a table with the columns 'query_id', 'subject_id' and 'dNdS' (as returned by dNdS()).",
                     call. = FALSE)
        check_n_quantile(n_quantile)
        
        # a copy: the table of the caller is left as it is (setDT() / setkeyv() turned a tibble of the caller into
        # a re-sorted data.table)
        dNdS_tbl <- data.table::as.data.table(dNdS_tbl)
        data.table::setkeyv(dNdS_tbl, c("query_id","subject_id"))
        
        dNdS_tbl_divMap <- dNdS_tbl[, .(dNdS, query_id, subject_id)]
        
        strata <- divergence_strata_of(dNdS_tbl_divMap[ , dNdS], n_quantile)
        
        dNdS_tbl_divMap$dNdS <- strata
        data.table::setnames(dNdS_tbl_divMap, old = "dNdS", new = "divergence_strata")
        
        if (!subject.id)
                return(dNdS_tbl_divMap[ , list(divergence_strata, query_id)])
        if (subject.id)
                return(dNdS_tbl_divMap[ , list(divergence_strata, query_id, subject_id)])
}

# The divergence stratum (1 .. n) of every value of `x`; NA for NA.
#
# Strata are the n-quantile classes of the values, as before: stratum k holds the values in (q[k-1], q[k]], where
# q[0], ..., q[n] are the quantiles of stats::quantile() (type 7) for the probabilities 0, 1/n, ..., 1, and the
# smallest value belongs to stratum 1. When all n + 1 quantiles are different this is computed with cut(), as it
# always was, and gives the same strata as earlier versions.
#
# When quantiles coincide (many equal values: e.g. dN/dS = 0 for more than 1/n of the genes of two close genomes)
# cut() cannot be used; earlier versions stopped. The same rule is then applied through ranks: with N values and
# r the smallest rank of a value among them (ties.method = "min"), its stratum is max(1, ceiling((r - 1) n / (N - 1))).
# For values without ties this is exactly the quantile class. Equal values always get the same stratum, namely the
# lowest stratum their block of ranks starts in; the strata the block spans beyond that stay empty. So the strata
# are no longer of equal size and fewer than n strata are used - a warning says which are empty.
divergence_strata_of <- function(x, n) {
        n <- as.integer(n)
        res <- rep(NA_integer_, length(x))
        ok <- !is.na(x)
        N <- sum(ok)
        if (N == 0)
                return(res)
        if (n == 1L || N == 1L) {
                res[ok] <- 1L
                return(res)
        }
        # seq(0, 1, 1/n) does not always end at 1 (n = 49, 98: the largest value got NA)
        breaks <- stats::quantile(x[ok], probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE, names = FALSE)
        if (!any(duplicated(breaks))) {
                res[ok] <- base::cut(x[ok], breaks = breaks, include.lowest = TRUE, labels = FALSE)
                if (!anyNA(res[ok]))
                        return(res)
        }
        r <- rank(x[ok], ties.method = "min")
        # ceiling((r - 1) * n / (N - 1)) in whole numbers
        k <- as.integer(((r - 1) * n + (N - 2)) %/% (N - 1))
        k[k < 1L] <- 1L
        res[ok] <- k
        empty <- setdiff(seq_len(n), unique(k))
        if (length(empty) > 0) {
                top <- sort(table(x[ok]), decreasing = TRUE)[1]
                warning(
                        "divergence_map(): the ", n, " quantile thresholds are not all different because many values are equal (",
                        "the most frequent value, ", names(top), ", occurs ", top, " times among ", N, " values). ",
                        "Equal values are kept in one stratum (the lowest one their ranks reach), so ", length(empty),
                        " of the ", n, " strata ", if (length(empty) == 1) "is" else "are", " empty (",
                        paste(empty, collapse = ", "), ") and the strata differ in size. See ?divergence_map.",
                        call. = FALSE
                )
        }
        res
}

# n_quantile: one finite whole number >= 1 (Inf passed the earlier check and failed later with an error about NA).
check_n_quantile <- function(n_quantile) {
        if (!is.numeric(n_quantile) || length(n_quantile) != 1 || is.na(n_quantile) || !is.finite(n_quantile) ||
            n_quantile < 1 || n_quantile != round(n_quantile) || n_quantile > .Machine$integer.max)
                stop("n_quantile must be one whole number >= 1 (the number of divergence strata), not '",
                     paste(n_quantile, collapse = "', '"), "'.", call. = FALSE)
        invisible(n_quantile)
}
