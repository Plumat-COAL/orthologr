#' @title Realign sequences within \code{deepclust} clusters
#' @description Runs \code{diamond realign} on an existing set of
#' \code{\link{deepclust}} clusters to obtain per-pair alignment statistics
#' (approximate identity, e-value, bitscore, and query/subject coverages) for
#' every \code{(member, representative)} pair inside each cluster. The result
#' can be passed directly to \code{\link{deepclust_annotate}} via its
#' \code{realign} argument to enrich a cluster profile without re-running
#' \code{diamond deepclust}.
#' @param input_file a character string, or a character vector of file paths,
#' specifying the input sequence file(s) that were used for the original
#' \code{\link{deepclust}} run. When multiple paths are provided they are merged
#' into a single temporary FASTA before the DIAMOND2 database is built (same
#' behaviour as \code{\link{deepclust}}). Both plain and gzip-compressed
#' (\code{.gz}) files are supported.
#' @param clusters either
#' \itemize{
#'   \item a \code{character} string — path to a tab-separated
#'     \code{deepclust} output file whose first column is
#'     \code{representative_id} and second column is \code{member_id}, or
#'   \item a \code{data.frame} / \code{tibble} with columns
#'     \code{representative_id} and \code{member_id} (i.e. the direct return
#'     value of \code{\link{deepclust}} or \code{\link{deepclust_annotate}}).
#' }
#' @param seq_type a character string specifying the sequence type stored in the
#' input file(s). Options are: \code{"cds"}, \code{"protein"}, or \code{"dna"}.
#' In case of \code{"cds"}, sequences are translated to protein sequences before
#' building the database. Default is \code{seq_type = "protein"}.
#' @param format a character string specifying the file format of the sequence
#' file. Default is \code{format = "fasta"}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with
#' corrupt base triplets should be removed from the input file. Only relevant
#' when \code{seq_type = "cds"}. Default is \code{delete_corrupt_cds = TRUE}.
#' @param path a character string specifying the path to the DIAMOND2 executable
#' (in case it is not on the system \code{PATH}).
#' @param comp_cores a numeric value specifying the number of CPU threads to
#' use. Passed to \code{--threads}. Default is \code{comp_cores = 1}.
#' @param realign_params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string of additional \code{diamond realign}
#' parameters to append to the command. Default is \code{NULL}.
#' @param save.output a character string specifying the path where the output
#' TSV file should be saved. E.g. \code{save.output = getwd()} to save in the
#' current working directory. Default is \code{NULL} (file is only written to a
#' temporary directory).
#' @param quiet a logical value indicating whether DIAMOND2 should run in quiet
#' mode. Default is \code{quiet = TRUE}.
#' @details
#' This function wraps \code{diamond realign}, which re-aligns every cluster
#' member against its assigned representative using the cluster mapping produced
#' by \code{\link{deepclust}}. The output format is BLAST tabular
#' (\code{--outfmt 6}) with the fields:
#' \code{qseqid sseqid approx_pident evalue bitscore qcovhsp scovhsp}.
#'
#' The \code{--clusters} argument consumed by \code{diamond realign} expects a
#' two-column tab-separated file (\code{representative_id TAB member_id}).
#' When a tibble is supplied to \code{clusters}, it is written to a temporary
#' file automatically.
#'
#' Columns in the returned tibble are named to match those used by
#' \code{\link{deepclust}} (\code{representative_id}, \code{member_id}) so that
#' the result joins cleanly onto \code{\link{deepclust_annotate}} output.
#' Note that \code{diamond realign} uses the \emph{representative} as the
#' alignment query (\code{qseqid}) and the \emph{member} as the subject
#' (\code{sseqid}); the returned tibble is named accordingly.
#'
#' @author Jaruwatana Sodai Lotharukpong
#' @references
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments
#' at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki
#' @examples \dontrun{
#' # 1. Run deepclust first
#' clusters <- deepclust(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   )
#' )
#'
#' # 2. Realign each representative against its cluster members to get alignment stats
#' realign_result <- deepclust_realign(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   clusters = clusters
#' )
#'
#' # 2. Alternatively, supply the path to the deepclust TSV directly
#' realign_result <- deepclust_realign(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   clusters = "/path/to/deepclust_output.tsv"
#' )
#'
#' # 3. Incorporate into deepclust_annotate without re-running deepclust
#' profile <- deepclust_annotate(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   realign = realign_result
#' )
#' }
#' @return A \code{\link[tibble]{tibble}} with seven columns:
#' \itemize{
#'   \item \code{representative_id} — accession of the cluster representative
#'     (query in the realignment).
#'   \item \code{member_id} — accession of the cluster member (subject in the
#'     realignment).
#'   \item \code{approx_pident} — approximate percentage of identical matches.
#'   \item \code{evalue} — expect value.
#'   \item \code{bitscore} — bit score.
#'   \item \code{qcovhsp} — query coverage per HSP.
#'   \item \code{scovhsp} — subject (representative) coverage per HSP.
#' }
#' @seealso \code{\link{deepclust}}, \code{\link{deepclust_annotate}},
#'   \code{\link{diamond}}, \code{\link{set_diamond}}
#' @export
deepclust_realign <- function(
        input_file,
        clusters,
        seq_type           = "protein",
        format             = "fasta",
        delete_corrupt_cds = TRUE,
        path               = NULL,
        comp_cores         = 1,
        realign_params     = NULL,
        save.output        = NULL,
        quiet              = TRUE) {

        # validate that all supplied sequence files exist before doing any work
        missing_files <- input_file[!file.exists(input_file)]
        if (length(missing_files) > 0)
                stop(
                        "The following input file(s) were not found:\n",
                        paste(missing_files, collapse = "\n"),
                        call. = FALSE
                )

        path <- diamond_resolve_path(path)
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        
        message("Running ", diamond_version(path = path), " ...")
        
        # in case one tries to use more cores than are available
        comp_cores <- check_comp_cores(comp_cores)
        
        # every call works in a folder of its own, which is removed when the call ends
        work <- call_dir_start("deepclust_realign_")
        on.exit(call_dir_end(work), add = TRUE)
        
        # resolve clusters: accept either a file path or a data frame
        if (is.character(clusters) && length(clusters) == 1 && file.exists(clusters)) {
                clusters_file <- clusters
        } else if (is.data.frame(clusters)) {
                required_cols <- c("representative_id", "member_id")
                if (!all(required_cols %in% names(clusters)))
                        stop(
                                "When 'clusters' is a data frame it must have columns: ",
                                paste(required_cols, collapse = ", "),
                                call. = FALSE
                        )
                clusters_file <- file.path(work, "deepclust_realign_input.tsv")
                readr::write_tsv(
                        clusters[, required_cols],
                        file      = clusters_file,
                        col_names = FALSE
                )
        } else {
                stop(
                        "'clusters' must be either a path to a deepclust TSV file or a ",
                        "data frame with columns 'representative_id' and 'member_id'.",
                        call. = FALSE
                )
        }
        
        filename <- if (length(input_file) > 1)
                paste0("deepclust_merged_", length(input_file), "_files.fasta")
        else
                basename(input_file)
        output   <- file.path(work, paste0("deepclust_realign_", safe_file_name(filename), ".tsv"))
        
        message("Building DIAMOND2 database for realign ...")
        
        db_dt <- deepclust_protein_table(
                input_file         = input_file,
                seq_type           = seq_type,
                delete_corrupt_cds = delete_corrupt_cds,
                dir                = work
        )
        database <- make_diamond_db(db_dt, db_file = file.path(work, "deepclust_db.fasta"), path = path,
                                    comp_cores = comp_cores, quiet = quiet)
        
        message("Running diamond realign ...")
        t_realign_start   <- proc.time()["elapsed"]
        res <- run_tool(
                "diamond",
                args = list(
                        "realign", "--db", database, "--clusters", clusters_file, "--out", output,
                        "--threads", comp_cores,
                        "--outfmt", "6", "qseqid", "sseqid", "approx_pident", "evalue", "bitscore", "qcovhsp", "scovhsp",
                        if (!is.null(realign_params)) I(realign_params),
                        if (quiet) "--quiet"
                ),
                path = path, tool = "diamond realign", error = FALSE, echo = !quiet
        )
        t_realign_elapsed <- proc.time()["elapsed"] - t_realign_start
        message("diamond realign completed in ", format_elapsed(t_realign_elapsed))
        if (!identical(as.integer(res$status), 0L)) {
                said <- utils::tail(res$stderr[nzchar(trimws(res$stderr))], 5)
                stop(
                        "diamond realign exited with non-zero status: ", res$status,
                        "\nPlease check the path to the DIAMOND2 executable, the input_file, and clusters.",
                        if (length(said) > 0) paste0("\ndiamond said: ", paste(said, collapse = "\n")) else "",
                        call. = FALSE
                )
        }
        if (!file.exists(output))
                stop("diamond realign wrote no table for '", filename, "'.", call. = FALSE)
        
        col_Classes <- c(rep("character", 2),
                         rep("double", 5))
        
        tryCatch({
                recluster_table <- data.table::fread(
                        input      = output,
                        sep        = "\t",
                        header     = FALSE,
                        quote      = "",
                        na.strings = "",
                        colClasses = col_Classes
                )
                
                data.table::setnames(
                        recluster_table,
                        old = paste0("V", 1:7),
                        new = c("representative_id", "member_id",
                                "approx_pident", "evalue", "bitscore",
                                "qcovhsp", "scovhsp")
                )
                
                if (!is.null(save.output)) {
                        save_output_copy(output, save.output, paste0("deepclust_realign_", filename, ".tsv"))
                        message("Recluster TSV saved to: ", save.output)
                }
                
                return(tibble::as_tibble(recluster_table))
                
        }, error = function(e) {
                stop(
                        "The deepclust_realign output file '", output,
                        "' could not be read correctly.",
                        " Please verify that diamond realign produced output.",
                        "\n",
                        "Error: ", conditionMessage(e),
                        call. = FALSE
                )
        })
}
