#' @title Perform a DIAMOND2 search
#' @description This function performs a DIAMOND2 search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}.
#' Default is \code{format} = \code{"fasta"}.
#' @param diamond_algorithm a character string specifying the DIAMOND2 algorithm that shall be used, option is currently limited to: \code{diamond_algorithm} = \code{"blastp"}
#' @param sensitivity_mode specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' - sensitivity_mode = "faster" : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70
#' - sensitivity_mode = "default" : Default mode. Designed for finding hits of >70
#' - sensitivity_mode = "fast" : fast alignment mode, but least sensitive (default). Designed for finding hits of >70
#' - sensitivity_mode = "mid-sensitive" : fast alignments between the fast mode and the sensitive mode in sensitivity.
#' - sensitivity_mode = "sensitive" : fast alignments, but full sensitivity for hits >40
#' - sensitivity_mode = "more-sensitive" : more sensitive than the sensitive mode.
#' - sensitivity_mode = "very-sensitive" : sensitive alignment mode.
#' - sensitivity_mode = "ultra-sensitive" : most sensitive alignment mode (sensitivity as high as BLASTP).
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND2 hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param path a character string specifying the path to the DIAMOND2 program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run DIAMOND2 searches.
#' @param diamond_params a character string listing the input parameters that shall be passed to the executing DIAMOND2 program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running DIAMOND2.
#' @param clean_folders a boolean value specifying whether all internal folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the DIAMOND2 output shall be stored. E.g. \code{save.output} = \code{getwd()}
#' to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @param quiet a logical value indicating whether DIAMOND2 should be run with the quiet mode.
#' Default is \code{quiet} = \code{TRUE} (which adds \code{--quiet} to the diamond run).
#' @param database_maker a character string specifying whether the database should be made using diamond or blast.
#' Default is \code{database_maker} = \code{diamond}.
#' @details This function provides a fast communication between R and DIAMOND2. It is mainly used as internal functions
#' such as \code{\link{diamond_best}} and \code{\link{diamond_rec}} but can also be used to perform simple DIAMOND2 computations.
#' This function gives the same output as \code{\link{blast}} while being up to 10 000X faster in larger databases.
#'
#' @author Jaruwatana Sodai Lotharukpong
#' @references
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
#' @examples \dontrun{
#' # performing a DIAMOND2 search using diamond blastp (default)
#' diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'
#' # performing a DIAMOND2 search using diamond blastp (default) 
#' # using amino acid sequences as input file
#' diamond(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'         seq_type     = "protein")
#'
#'
#' # save the DIAMOND2 output table in your current working directory
#' diamond(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'         seq_type     = "protein",
#'         save.output  = getwd())
#'
#' # in case you are working with a multicore machine, you can also run parallel
#' # DIAMOND2 computations using the comp_cores parameter: here with 2 cores
#' diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'         comp_cores   = 2)
#'
#'
#'  # running diamond using additional parameters
#'  diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'          diamond_params = "--max-target-seqs 1")
#'
#'
#' # running diamond using additional parameters and an external diamond path
#'  diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'          diamond_params = "--max-target-seqs 1", path = "path/to/diamond/")
#' }
#'
#' @return A data.table storing the DIAMOND2 hit table returned by DIAMOND2. The format is the same as with BLAST.
#' @seealso \code{\link{diamond_best}}, \code{\link{diamond_rec}}, \code{\link{set_diamond}}, \code{\link{blast}}
#' @export
diamond <- function(
                query_file,
                subject_file,
                seq_type        = "cds",
                format          = "fasta",
                diamond_algorithm = "blastp",
                sensitivity_mode = "fast",
                eval            = "1E-5",
                max.target.seqs = 10000,
                delete_corrupt_cds = TRUE,
                path            = NULL,
                comp_cores      = 1,
                diamond_params    = NULL,
                clean_folders   = FALSE,
                save.output     = NULL,
                quiet           = TRUE,
                database_maker  = "diamond") {
        
        if (!is.element(diamond_algorithm, c("blastp")))
                stop(
                        "Please choose a valid DIAMOND mode. Only 'blastp' is available for this function.",
                        call. = FALSE
                )
        
        if (!is.element(database_maker, c("diamond", "blast")))
                stop("Please choose either: 'diamond' or 'blast' as database_maker.", call. = FALSE)
        
        path <- diamond_resolve_path(path)
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        
        message("Running ", diamond_version(path = path), " ...")
        
        message("sensitivity mode: ", sensitivity_mode)
        
        # in case one tries to use more cores than are available
        check_comp_cores(comp_cores)
        
        # every call works in a folder of its own, which is removed when the call ends: a table left by an
        # earlier or failed call cannot be returned as the result of this one
        work <- call_dir_start("diamond_")
        on.exit(call_dir_end(work), add = TRUE)
        
        if (database_maker == "diamond") {
                message("creating a diamond database")
                query.dt <- set_diamond(
                        file     = query_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds,
                        path     = path,
                        comp_cores = comp_cores
                )[[1]]
                
                subject.dt <- set_diamond(
                        file     = subject_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds,
                        path     = path,
                        comp_cores = comp_cores
                )[[1]]
                database <- make_diamond_db(subject.dt, db_file = file.path(work, "subject_db.fasta"),
                                            path = path, comp_cores = comp_cores)
        } else if (database_maker == "blast") {
                message("creating a blast database")
                query.dt <- set_blast(
                        file     = query_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds
                )[[1]]
                
                subject.dt <- set_blast(
                        file     = subject_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds
                )[[1]]
                database <- make_blast_db(subject.dt, db_file = file.path(work, "subject_db.fasta"))
        }
        
        filename <- basename(query_file)
        input    <- file.path(work, "query.fasta")
        output   <- file.path(work, "blastresult.tsv")
        
        # where the hit table is kept after the call (as before); a table of an earlier call is removed now
        db_dir <- file.path(orthologr_tmpdir(), "_blast_db")
        kept   <- file.path(db_dir, paste0("blastresult_", filename, ".csv"))
        if (file.exists(kept))
                unlink(kept)
        context  <- paste0("query '", filename, "'")
        
        tryCatch({
                write_protein_fasta(query.dt, input)
        }, error = function(e) {
                stop("The query sequences of '", query_file, "' could not be written to '", input, "': ",
                     conditionMessage(e), call. = FALSE)
        })
        
        message("Starting DIAMOND2 search ...")
        t_diamond_start   <- proc.time()["elapsed"]
        res <- run_tool(
                "diamond",
                args = c(
                        list("blastp",
                             "--db", database,
                             "--query", input,
                             "--evalue", as.numeric(eval),
                             paste0("--", sensitivity_mode),
                             "--max-target-seqs", format(max.target.seqs, scientific = FALSE),
                             "--out", output,
                             "--threads", comp_cores,
                             "--outfmt", "6"),
                        as.list(c("qseqid", "sseqid", "pident", "nident", "length", "mismatch", "gapopen", "gaps",
                                  "positive", "ppos", "qstart", "qend", "qlen", "qcovhsp", "qcovhsp", "sstart",
                                  "send", "slen", "evalue", "bitscore", "score")),
                        list(if (!is.null(diamond_params)) I(diamond_params),
                             if (quiet) "--quiet")
                ),
                path = path, tool = paste("diamond", diamond_algorithm), context = context, error = FALSE,
                echo = !quiet
        )
        t_diamond_elapsed <- proc.time()["elapsed"] - t_diamond_start
        message("DIAMOND2 search completed in ", format_elapsed(t_diamond_elapsed))
        
        if (!identical(as.integer(res$status), 0L)) {
                said <- utils::tail(res$stderr[nzchar(trimws(res$stderr))], 5)
                stop(
                        "diamond ", diamond_algorithm,
                        " exited with non-zero status: ", res$status,
                        "\nPlease check the path to the DIAMOND2 executable and the input file.",
                        if (length(said) > 0) paste0("\ndiamond said: ", paste(said, collapse = "\n")) else "",
                        call. = FALSE
                )
        }
        if (!file.exists(output))
                stop("diamond ", diamond_algorithm, " wrote no hit table for ", context, ".", call. = FALSE)
        
        diamond_table_names <-
                c(
                        "query_id",
                        "subject_id",
                        "perc_identity",
                        "num_ident_matches",
                        "alig_length",
                        "mismatches",
                        "gap_openings",
                        "n_gaps",
                        "pos_match",
                        "ppos",
                        "q_start",
                        "q_end",
                        "q_len",
                        "qcov",
                        "qcovhsp",
                        "s_start",
                        "s_end",
                        "s_len",
                        "evalue",
                        "bit_score",
                        "score_raw"
                )
        
        hit_table <- read_hit_table(output, diamond_table_names, tool = "DIAMOND")
        
        # save the DIAMOND output file to path save.output
        save_output_copy(output, save.output, paste0("blastresult_", filename, ".csv"))
        
        if (clean_folders) {
                clean_all_folders(db_dir)
        } else {
                if (!file.exists(db_dir))
                        dir.create(db_dir, recursive = TRUE)
                publish_file(output, kept)
        }
        
        return(hit_table)
}

# "diamond version 2.x.y" - the DIAMOND version line; stops when diamond cannot be started. Asked once per session.
diamond_version <- function(path = NULL) {
        tool_available(
                "diamond", version_args = "version", pattern = "diamond", path = path,
                not_found = paste0(
                        "It seems like you don't have DIAMOND2 installed locally on your machine or the PATH variable to the DIAMOND2 program is not set correctly.",
                        " Please use the 'diamond_exec_path' argument to specify a path to your diamond excecutable, e.g. diamond_exec_path = '/usr/local/bin' or diamond_exec_path = '/opt/miniconda3/bin'.")
        )
}
