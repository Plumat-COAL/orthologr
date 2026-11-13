#' @title Cluster sequences using \code{diamond deepclust}
#' @description This function clusters a set of protein (or CDS) sequences using the
#' DIAMOND2 \code{deepclust} algorithm. It prepares a DIAMOND2 database via
#' \code{\link{set_diamond}}, runs \code{diamond deepclust}, and returns a
#' two-column tibble mapping each representative sequence to its cluster members.
#' The result is also written to a TSV file.
#' Note, we recommend users to run \code{deepclust_annotate()} for more information  
#' about the cluster members and their representative sequences.
#' @param input_file a character string, or a character vector of file paths,
#' specifying the input sequence file(s). When multiple paths are provided the
#' sequences are merged into a single temporary FASTA before clustering, allowing
#' proteomes from different organisms to be clustered together. Both plain and
#' gzip-compressed (\code{.gz}) files are supported and may be mixed freely.
#' Supported sequence types are controlled by \code{seq_type}.
#' @param seq_type a character string specifying the sequence type stored in the
#' input file(s). Options are: \code{"cds"}, \code{"protein"}, or \code{"dna"}.
#' In case of \code{"cds"}, sequences are translated to protein sequences before
#' clustering. Default is \code{seq_type = "protein"}.
#' @param format a character string specifying the file format of the sequence file.
#' Default is \code{format = "fasta"}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with
#' corrupt base triplets should be removed from the input file. Only relevant when
#' \code{seq_type = "cds"}. Default is \code{delete_corrupt_cds = TRUE}.
#' @param path a character string specifying the path to the DIAMOND2 executable
#' (in case it is not on the system \code{PATH}).
#' @param comp_cores a numeric value specifying the number of CPU threads to use.
#' Passed to \code{--threads}. Default is \code{comp_cores = 1}.
#' @param mutual_cover a numeric value (percentage) specifying the minimum mutual
#' coverage of the cluster member and representative sequence. Passed to
#' \code{--mutual-cover}. Default is \code{NULL} (DIAMOND2 default applies).
#' @param approx_id a numeric value (percentage) specifying the minimum approximate
#' identity to report an alignment or to cluster sequences. Passed to
#' \code{--approx-id}. Default is \code{NULL} (DIAMOND2 default applies).
#' @param eval a numeric value specifying the maximum E-value to report alignments.
#' Passed to \code{--evalue}. Default is \code{NULL} (DIAMOND2 default of 0.001 applies).
#' @param deepclust_params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string of additional DIAMOND2 deepclust  
#' parameters to append to the \code{deepclust} call. Default is \code{NULL}.
#' @param save.output a character string specifying the path where the output TSV
#' file should be saved. E.g. \code{save.output = getwd()} to save in the current
#' working directory. Default is \code{NULL} (file is only written to a temporary
#' directory).
#' @param quiet a logical value indicating whether DIAMOND2 should run in quiet
#' mode. Default is \code{quiet = TRUE}.
#' @param direct_input a logical value. \code{FALSE} (Default): the sequences are read into R, sorted by their IDs
#' and written to one FASTA file from which the DIAMOND2 database is built (as in earlier versions). \code{TRUE}:
#' the input files are handed to \code{diamond makedb} as they are (several files are copied one after the other
#' into one file, without being parsed), so no sequence is held in R memory - meant for very large protein sets.
#' The sequence IDs are the same, but the sequences reach DIAMOND2 in file order instead of sorted by ID, and
#' \code{diamond deepclust} is sensitive to that order: \strong{cluster representatives and cluster membership
#' can differ from the default}. Only for \code{seq_type = "protein"}.
#' @details
#' This function wraps \code{diamond deepclust}, the graph-based sequence clustering
#' algorithm available in DIAMOND2 v2.1.0 and later. The output is a two-column
#' tab-separated file where the first column contains the representative (centroid)
#' sequence accession and the second column contains the cluster member accession.
#' A sequence that is its own representative will appear with itself in both columns.
#'
#' When multiple files are supplied via \code{input_file}, they are first merged
#' into a single temporary FASTA file (using \code{seqinr::read.fasta} and
#' \code{seqinr::write.fasta}) before the DIAMOND2 database is built. This allows
#' cross-species or multi-proteome clustering in a single run. The merged file is
#' written to the session temporary directory and is not retained after the R
#' session ends unless \code{save.output} is specified.
#'
#' The function uses \code{\link{set_diamond}} internally to handle CDS translation
#' and DIAMOND2 database creation.
#'
#' @author Jaruwatana Sodai Lotharukpong
#' @references
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments
#' at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki
#' @examplesIf nzchar(Sys.which("diamond"))
#' # Cluster a single proteome (protein sequences, the default)
#' deepclust(
#'   input_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr')
#' )
#'
#' # Cluster two proteomes together in a single run
#' deepclust(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   )
#' )
#'
#' \dontrun{
#' # Apply identity and coverage thresholds across multiple proteomes
#' deepclust(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   approx_id    = 50,
#'   mutual_cover = 80,
#'   eval         = 1e-5,
#'   comp_cores   = 4
#' )
#'
#' # Save the cluster TSV to the folder tempdir()
#' deepclust(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   save.output = tempdir()
#' )
#'
#' # Cluster CDS sequences (translated to protein internally)
#' deepclust(
#'   input_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   seq_type   = "cds"
#' )
#' }
#' @return A \code{\link[tibble]{tibble}} with two columns:
#' \itemize{
#'   \item \code{representative_id} — accession of the cluster representative sequence.
#'   \item \code{member_id} — accession of the cluster member sequence.
#' }
#' @seealso \code{\link{diamond}}, \code{\link{set_diamond}}, \code{\link{diamond_best}}, \code{\link{diamond_rec}}, \code{\link{deepclust_annotate}}, \code{\link{deepclust_realign}}
#' @export
deepclust <- function(
        input_file,
        seq_type           = "protein",
        format             = "fasta",
        delete_corrupt_cds = TRUE,
        path               = NULL,
        comp_cores         = 1,
        mutual_cover       = NULL,
        approx_id          = NULL,
        eval               = NULL,
        deepclust_params   = NULL,
        save.output        = NULL,
        quiet              = TRUE,
        direct_input       = FALSE) {

        # validate that all supplied files exist before doing any work
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
        work <- call_dir_start("deepclust_")
        on.exit(call_dir_end(work), add = TRUE)
        
        filename <- if (length(input_file) > 1)
                paste0("deepclust_merged_", length(input_file), "_files.fasta")
        else
                basename(input_file)
        output   <- file.path(work, paste0("deepclust_", safe_file_name(filename), ".tsv"))
        
        message("Building DIAMOND2 database for deepclust ...")
        
        if (isTRUE(direct_input)) {
                if (seq_type != "protein")
                        stop("direct_input = TRUE needs protein sequences (seq_type = 'protein').", call. = FALSE)
                # the files go to DIAMOND as they are
                db_in <- concatenate_files(input_file, file.path(work, "deepclust_input.fasta"))
                database <- diamond_makedb(in_file = db_in, db_file = file.path(work, "deepclust_db.fasta"),
                                           path = path, comp_cores = comp_cores, quiet = quiet)
        } else {
                db_dt <- deepclust_protein_table(
                        input_file         = input_file,
                        seq_type           = seq_type,
                        delete_corrupt_cds = delete_corrupt_cds,
                        dir                = work,
                        warn_duplicates    = TRUE
                )
                database <- make_diamond_db(db_dt, db_file = file.path(work, "deepclust_db.fasta"), path = path,
                                            comp_cores = comp_cores, quiet = quiet)
        }
        
        message("Running diamond deepclust ...")
        t_deepclust_start   <- proc.time()["elapsed"]
        res <- run_tool(
                "diamond",
                args = list(
                        "deepclust", "--db", database, "--out", output, "--threads", comp_cores,
                        if (!is.null(approx_id)) c("--approx-id", as.numeric(approx_id)),
                        if (!is.null(mutual_cover)) c("--mutual-cover", as.numeric(mutual_cover)),
                        if (!is.null(eval)) c("--evalue", as.numeric(eval)),
                        if (!is.null(deepclust_params)) I(deepclust_params),
                        if (quiet) "--quiet"
                ),
                path = path, tool = "diamond deepclust", error = FALSE, echo = !quiet
        )
        t_deepclust_elapsed <- proc.time()["elapsed"] - t_deepclust_start
        message("diamond deepclust completed in ", format_elapsed(t_deepclust_elapsed))
        if (!identical(as.integer(res$status), 0L)) {
                said <- utils::tail(res$stderr[nzchar(trimws(res$stderr))], 5)
                stop(
                        "diamond deepclust exited with non-zero status: ", res$status,
                        "\nPlease check the path to the DIAMOND2 executable and the input file.",
                        if (length(said) > 0) paste0("\ndiamond said: ", paste(said, collapse = "\n")) else "",
                        call. = FALSE
                )
        }
        if (!file.exists(output))
                stop("diamond deepclust wrote no cluster table for '", filename, "'.", call. = FALSE)
        
        tryCatch({
                cluster_table <- data.table::fread(
                        file    = output,
                        sep     = "\t",
                        header  = FALSE,
                        quote   = "",
                        na.strings = "",
                        col.names = c("representative_id", "member_id"),
                        colClasses = "character"
                )
                
                if (!is.null(save.output)) {
                        save_output_copy(output, save.output, paste0("deepclust_", filename, ".tsv"))
                        message("Cluster TSV saved to: ", save.output)
                }
                
                return(tibble::as_tibble(cluster_table))
                
        }, error = function(e) {
                stop(
                        "The deepclust output file '", output, "' could not be read correctly.",
                        " Please verify that diamond deepclust produced output.",
                        "\n",
                        "Error: ", conditionMessage(e),
                        call. = FALSE
                )
        })
}

# The proteins deepclust() and deepclust_realign() build their database from: a data.table (geneids, aa) sorted by
# ID - the table set_diamond() returned for the merged input file in earlier versions.
#   proteins : every file is read once with Biostrings and the tables are put together (several files used to be
#              read sequence by sequence with seqinr, written to one file and read again)
#   cds      : several files are merged into one FASTA file first, as before, and translated
deepclust_protein_table <- function(input_file, seq_type, delete_corrupt_cds = TRUE, dir, warn_duplicates = FALSE) {
        geneids <- NULL
        if (seq_type != "protein" || length(input_file) == 1) {
                merged <- merge_fasta_files(input_file, seq_type = seq_type, dir = dir, warn_duplicates = warn_duplicates)
                return(protein_table(file = merged, seq_type = seq_type, delete_corrupt_cds = delete_corrupt_cds))
        }
        
        message("Merging ", length(input_file), " input files ...")
        dt <- data.table::rbindlist(lapply(input_file, function(f) protein_table(file = f, seq_type = "protein")))
        message("Merged ", nrow(dt), " sequences into: ", paste0("deepclust_merged_", length(input_file), "_files.fasta"))
        
        if (warn_duplicates)
                warn_duplicate_ids(dt[ , geneids])
        
        data.table::setkey(dt, geneids)
        dt
}

warn_duplicate_ids <- function(all_names) {
        dup_ids <- all_names[duplicated(all_names)]
        if (length(dup_ids) > 0) {
                warning(
                        "Duplicate sequence IDs found across input files. ",
                        "DIAMOND2 may produce unexpected results.\n",
                        "Affected IDs include: ",
                        paste(utils::head(unique(dup_ids)), collapse = ", "),
                        call. = FALSE
                )
        }
        invisible(NULL)
}

# Copy files one after the other into `out` without parsing them (gzip-compressed files are decompressed).
concatenate_files <- function(files, out) {
        if (length(files) == 1 && !grepl("\\.gz$", files))
                return(files)
        con_out <- file(out, open = "w")
        on.exit(close(con_out), add = TRUE)
        for (f in files) {
                con_in <- gzfile(f, open = "r")
                repeat {
                        lines <- readLines(con_in, n = 100000L, warn = FALSE)
                        if (length(lines) == 0)
                                break
                        writeLines(lines, con_out)
                }
                close(con_in)
        }
        out
}

# Several FASTA files as one file in `dir` (one file is returned as it is). Used by deepclust() and
# deepclust_realign(), which hand one database to DIAMOND.
merge_fasta_files <- function(input_file, seq_type, dir, warn_duplicates = FALSE) {
        if (length(input_file) <= 1)
                return(input_file)
        
        n_files       <- length(input_file)
        seqtype_fasta <- if (seq_type == "protein") "AA" else "DNA"
        
        message("Merging ", n_files, " input files ...")
        
        all_seqs  <- list()
        all_names <- character(0)
        
        for (f in input_file) {
                seqs <- seqinr::read.fasta(
                        file            = f,
                        seqtype         = seqtype_fasta,
                        as.string       = TRUE,
                        forceDNAtolower = FALSE
                )
                all_seqs  <- c(all_seqs,  seqs)
                all_names <- c(all_names, names(seqs))
        }
        
        merged_filename <- paste0("deepclust_merged_", n_files, "_files.fasta")
        merged_file     <- file.path(dir, merged_filename)
        
        seqinr::write.fasta(
                sequences = all_seqs,
                names     = all_names,
                file.out  = merged_file,
                nbchar    = 80
        )
        
        message("Merged ", length(all_names), " sequences into: ", merged_filename)
        
        if (warn_duplicates)
                warn_duplicate_ids(all_names)
        
        merged_file
}
