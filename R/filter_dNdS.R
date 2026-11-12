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
#' \item \code{"keep"} (default): they are kept, with \code{dNdS = NA}, and the logical column \code{dS_zero} marks
#' them (\code{FALSE} for all other rows). \strong{This is a change of the default}: up to orthologr 0.4.3 these
#' rows were dropped without a message (because \code{NA <= dnds.threshold} is \code{NA}), so the result can now
#' hold more rows, and it has one more column.
#' \item \code{"drop"}: they are left out, as before; no column is added.
#' }
#' The default can be set for the session with \code{options(orthologr.dS_zero = "drop")}. A message gives the number
#' of such gene pairs in either case.
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
#' @examples \dontrun{
#' 
#' filter_dNdS( dNdS( query_file = 
#'                    system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                    subject_file = 
#'                    system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                    ortho_detection = "RBH", 
#'                    aa_aln_type = "pairwise",
#'                    aa_aln_tool = "NW",
#'                    codon_aln_tool = "pal2nal", 
#'                    dnds_est.method = "Comeron", 
#'                    comp_cores = 1), 
#'              dnds.threshold = 2)
#' 
#' # a small example using clustalw
#' filter_dNdS( dNdS( query_file      = 
#'                    system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                    subject_file    = 
#'                    system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                    ortho_detection = "RBH", 
#'                    aa_aln_type     = "pairwise",
#'                    aa_aln_tool     = "NW",
#'                    codon_aln_tool  = "pal2nal", 
#'                    dnds_est.method = "Comeron", 
#'                    comp_cores      = 1) , 
#'                dnds.threshold  = 2)
#'                   
#' }
#' @seealso \code{\link{dNdS}}, \code{\link{divergence_stratigraphy}}
#' @export
filter_dNdS <- function(dNdS_tbl, dnds.threshold = 2, dS_zero = getOption("orthologr.dS_zero", "keep")){
        dN <- dS <- NULL
        if (!is.character(dS_zero) || length(dS_zero) != 1 || !is.element(dS_zero, c("keep", "drop")))
                stop("dS_zero must be \"keep\" or \"drop\".", call. = FALSE)
        message("Filtering out NA values in dN or dS and all values with dNdS > ", dnds.threshold, " ...")
        message("Initial input contains ", nrow(dNdS_tbl), " rows.")
        
        # Gene pairs without synonymous divergence: dS = 0, so that dN/dS is not defined (NA).
        zero <- !is.na(dNdS_tbl$dN) & !is.na(dNdS_tbl$dS) & dNdS_tbl$dS == 0
        n_identical <- sum(zero & dNdS_tbl$dN == 0)
        n_dn_only <- sum(zero & dNdS_tbl$dN > 0)
        
        if (dS_zero == "drop") {
                # as up to 0.4.3: 'NA <= threshold' is NA, and such rows were left out
                res <- dplyr::filter(dNdS_tbl,!is.na(dN), !is.na(dS), dNdS <= dnds.threshold)
        } else {
                keep <- !is.na(dNdS_tbl$dN) & !is.na(dNdS_tbl$dS) &
                        (zero | (!is.na(dNdS_tbl$dNdS) & dNdS_tbl$dNdS <= dnds.threshold))
                res <- dNdS_tbl[keep, , drop = FALSE]
                res$dNdS[zero[keep]] <- NA_real_
                res$dS_zero <- zero[keep]
        }
        if (any(zero))
                message(sum(zero), " gene pairs have dS = 0, so that their dN/dS is not defined: ", n_identical,
                        " with dN = 0 as well (no substitutions between the aligned coding sequences) and ", n_dn_only,
                        " with dN > 0. ",
                        if (dS_zero == "keep")
                                "They are kept with dNdS = NA and dS_zero = TRUE (dS_zero = \"drop\" leaves them out, as orthologr <= 0.4.3 did)."
                        else
                                "They are left out (dS_zero = \"drop\").")
        message("Filtering done. New output table contains ", nrow(res), " rows.")
        return(res)
}
