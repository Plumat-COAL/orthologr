#' @title Compute dNdS Values For A Given Pairwise Alignment
#' @description This function takes a vector containing the
#' query amino acid sequence, subject amino acid sequence, query CDS sequence, and subject CDS sequence
#' and then runs the following pipieline:
#' 
#' \itemize{
#' 
#' \item 1) Multiple-Alignment of query amino acid sequence and subject amino acid sequence
#' 
#' \item 2) Codon-Alignment of the amino acid alignment returned by 1) and query CDS sequence + subject CDS sequence
#' 
#' \item 3) dNdS estimation of the codon alignment returned by 2)
#' 
#' }
#' @param complete_tbl a data.table object storing the query_id, subject_id, query_cds (sequence), 
#' subject_cds (sequence), query_aa (sequence), and subject_aa (sequence) of the organisms that shall be compared.
#' @param aa_aln_tool a character string specifying the multiple alignment tool that shall be used for pairwise protein alignments.
#' @param aa_aln_path a character string specifying the path to the corresponding multiple alignment tool.
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type} = "multiple" or \code{aa_aln_type} = \code{"pairwise"}.
#' Default is \code{aa_aln_type} = "multiple".
#' @param aa_aln_params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string specifying additional parameters that shall be passed to the multiple alignment system call.
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used for codon alignments. Default is \code{codon_aln_tool} = \code{"pal2nal"}.
#' @param dnds_est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' @param kaks_tool the program that computes the KaKs_Calculator methods: \code{"kaks_calculator"} (Default) or \code{"coal"}; see \code{\link{substitutionrate}}.
#' @param aa_sub_matrix the substitution matrix for pairwise amino acid alignments (\code{aa_aln_type = "pairwise"}); see \code{\link{pairwise_aln}}.
#' @param quiet \code{FALSE}: a progress message about every tenth of the gene pairs (the programs called for a
#' gene pair no longer report every file they write).
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform parallel computations on a multicore machine.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function takes the amino acid and CDS sequences two orthologous genes
#' and writes the corresponding amino acid and CDS sequences as fasta file into
#' the internal folder environment. The resulting fasta files (two files) store the 
#' amino acid sequence of the query_id and subject_id (file one) and the CDS sequence of
#' the query_id and subject_id (file two). These fasta files are then used to pass through the following pipeline:
#' 
#' 1) Multiple-Alignment or Pairwise-Alignment of query amino acid sequence and subject amino acid sequence
#' 
#' 2) Codon-Alignment of the amino acid alignment returned by 1) and query CDS sequence + subject CDS sequence
#' 
#' 3) dNdS estimation of the codon alignment returned by 2)
#' 
#' A gene pair for which one of the steps fails (for example because PAL2NAL or KaKs_Calculator exits with an error)
#' does not end the run: the pair is returned with \code{NA} values for dN, dS and dNdS, the reasons are collected in
#' the attribute \code{"failed_pairs"} of the result (a \code{data.frame} with \code{query_id}, \code{subject_id} and
#' \code{reason}), and a warning reports how many pairs failed.
#' 
#' @seealso \code{\link{multi_aln}}, \code{\link{substitutionrate}}, \code{\link{dNdS}}
compute_dnds <- function(complete_tbl,
                         aa_aln_type     = "multiple", 
                         aa_aln_tool     = "clustalw", 
                         aa_aln_path     = NULL,
                         aa_aln_params   = NULL, 
                         codon_aln_tool  = "pal2nal",
                         dnds_est.method = "YN",
                         store_locally   = FALSE,
                         kaks_calc_path  = NULL, 
                         quiet           = FALSE, 
                         comp_cores      = 1, 
                         clean_folders   = FALSE,
                         aa_sub_matrix   = "BLOSUM62",
                         kaks_tool       = "kaks_calculator"){
        
        # every argument is evaluated here: the gene pairs may be computed in other R processes (see force_args())
        force_args()
        comp_cores <- check_comp_cores(comp_cores)
        
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- subject_id <- query_cds <- subject_cds <- query_aa <- subject_aa <- i <- NULL
        
        multicore <- (comp_cores > 1)
        
        for (d in list(c("_alignment"), c("_alignment", "orthologs"), c("_alignment", "orthologs", "CDS"),
                       c("_alignment", "orthologs", "AA"))) {
                d <- do.call(file.path, as.list(c(orthologr_tmpdir(), d)))
                if (!file.exists(d))
                        dir.create(d, recursive = TRUE)
        }
        
        # In the files of a gene pair the two sequences are called 'query' and 'subject', whatever their IDs are:
        # PAL2NAL matches protein and coding sequence by name, so a pair whose two genes have the same ID (two
        # strains with the same locus tags, gene symbols or running numbers as IDs) could not be computed
        # ("inconsistency between pep and nuc"), and alignment programs change IDs they do not like (Clustal Omega:
        # '-1'). The IDs of the result come from the pair table. With store_locally = TRUE the table
        # 'pair_ids.tsv' next to the stored alignments says which pair a file 'q<i>_...' holds.
        pair_names <- dnds_pair_names()
        if (store_locally) {
                idx_dir <- file.path("orthologr_alignment_files", "_alignment")
                dir.create(idx_dir, recursive = TRUE, showWarnings = FALSE)
                utils::write.table(
                        data.frame(file_prefix = paste0("q", seq_len(nrow(complete_tbl))),
                                   query = as.character(complete_tbl[ , query_id]),
                                   subject = as.character(complete_tbl[ , subject_id]), stringsAsFactors = FALSE),
                        file.path(idx_dir, "pair_ids.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        }
        
        # The six columns of the pair table as plain vectors: one_pair() reads one element of each per gene pair,
        # and six look-ups in the data.table per pair took about 9 % of the run time.
        q_ids <- complete_tbl[["query_id"]]
        s_ids <- complete_tbl[["subject_id"]]
        q_cds <- complete_tbl[["query_cds"]]
        s_cds <- complete_tbl[["subject_cds"]]
        q_aa  <- complete_tbl[["query_aa"]]
        s_aa  <- complete_tbl[["subject_aa"]]
        
        # the complete pipeline for gene pair i: write the two FASTA files, align the proteins,
        # build the codon alignment, estimate dN/dS
        #   aln_rows : NULL, or the two rows of the pairwise protein alignment of the pair computed beforehand by
        #              pairwise_aln_rows() (they are written to the file pairwise_aln() would write)
        one_pair <- function(i, aln_rows = NULL) {
                # the IDs, coding sequences and amino acid sequences of the gene pair
                orthologs_names <- list(q_ids[i], s_ids[i])
                cds_seqs <- list(q_cds[i], s_cds[i])
                aa_seqs <- list(q_aa[i], s_aa[i])
                
                aa_aln_name <- paste0("q", i)
                
                # create cds fasta of orthologous gene pair having session name: 'aa_aln_name'
                cds_session_fasta <-
                        file.path(
                                orthologr_tmpdir(),
                                "_alignment",
                                "orthologs",
                                "CDS",
                                paste0(aa_aln_name, "_cds.fasta")
                        )
                
                write_pair_fasta(cds_seqs, pair_names, cds_session_fasta)
                
                # create aa fasta of orthologous gene pair having session name: 'aa_aln_name'
                aa_session_fasta <-
                        file.path(
                                orthologr_tmpdir(),
                                "_alignment",
                                "orthologs",
                                "AA",
                                paste0(aa_aln_name, "_aa.fasta")
                        )
                
                write_pair_fasta(aa_seqs, pair_names, aa_session_fasta)
                
                if (aa_aln_type == "multiple") {
                        # align aa -> <aa_aln_tool>.aln
                        multi_aln(
                                file           = aa_session_fasta,
                                tool           = aa_aln_tool,
                                get_aln        = FALSE,
                                multi_aln_name = aa_aln_name,
                                path           = aa_aln_path,
                                params         = aa_aln_params,
                                quiet          = TRUE
                        )
                        
                        aa_aln_file <-
                                file.path(
                                        orthologr_tmpdir(),
                                        "_alignment",
                                        "multi_aln",
                                        paste0(aa_aln_name, "_", aa_aln_tool, ".aln")
                                )
                }
                
                if (aa_aln_type == "pairwise" && !is.null(aln_rows)) {
                        # the alignment computed together with the other pairs of its chunk: the file pairwise_aln()
                        # writes (sequence names 'query' and 'subject', one line per sequence)
                        aa_aln_file <-
                                file.path(orthologr_tmpdir(), "_alignment", "pairwise_aln",
                                          paste0(aa_aln_name, "_", aa_aln_tool, "_AA.aln"))
                        write_pair_fasta(as.list(aln_rows), pair_names, aa_aln_file)
                } else if (aa_aln_type == "pairwise") {
                        # align aa -> <aa_aln_tool>.aln
                        pairwise_aln(
                                file              = aa_session_fasta,
                                tool              = aa_aln_tool,
                                seq_type          = "protein",
                                get_aln           = FALSE,
                                pairwise_aln_name = aa_aln_name,
                                path              = aa_aln_path,
                                quiet             = TRUE,
                                store_locally     = store_locally,
                                aa_sub_matrix     = aa_sub_matrix
                        )
                        
                        aa_aln_file <-
                                file.path(
                                        ifelse(store_locally, "orthologr_alignment_files", orthologr_tmpdir()),
                                        "_alignment",
                                        "pairwise_aln",
                                        paste0(aa_aln_name, "_", aa_aln_tool, "_AA.aln")
                                )
                }
                
                # align codon -> cds.aln
                codon_aln(
                        file_aln       = aa_aln_file,
                        file_nuc       = cds_session_fasta,
                        tool           = codon_aln_tool,
                        format         = "fasta",
                        codon_aln_name = aa_aln_name,
                        get_aln        = FALSE,
                        quiet          = TRUE
                )
                
                # compute kaks
                rates <- substitutionrate(
                        file           = file.path(
                                orthologr_tmpdir(),
                                "_alignment",
                                "codon_aln",
                                paste0(aa_aln_name, "_", codon_aln_tool, ".aln")
                        ),
                        subst_name     = aa_aln_name,
                        kaks_calc_path = kaks_calc_path,
                        est.method     = dnds_est.method,
                        quiet          = TRUE,
                        kaks_tool      = kaks_tool
                )
                
                # The IDs of the pair are those of the pair table - not what came back from the programs.
                # Alignment programs and file readers change IDs (MAFFT cuts them to 15 characters in CLUSTAL
                # output, fread() turned '007' into 7, KaKs_Calculator labels were split at '-'), and the result
                # was then joined with the hit table by these changed IDs: rows were lost or held IDs that are in
                # no input file.
                kaks_failed <- attr(rates, "kaks_failed")
                rates <- data.table::as.data.table(rates)
                data.table::set(rates, j = "query_id", value = as.character(orthologs_names[[1]]))
                data.table::set(rates, j = "subject_id", value = as.character(orthologs_names[[2]]))
                # an estimate of KaKs_Calculator that failed ('nan'): the row is a failed gene pair (with the method);
                # otherwise the column '.reason' is NA (set here with data.table::set(): ':=' in safe_pair() took
                # about 1 ms per pair)
                if (!is.null(kaks_failed) && nrow(kaks_failed) > 0 && is.element("method", names(rates))) {
                        data.table::set(rates, j = ".reason",
                                        value = dnds_pair_reason(kaks_failed$reason[match(rates$method, kaks_failed$method)],
                                                                 q_ids[i], s_ids[i]))
                } else if (!is.element(".reason", names(rates))) {
                        data.table::set(rates, j = ".reason", value = NA_character_)
                }
                rates
        }
        
        # one gene pair must not end the run: a pair whose pipeline fails is returned as a row of NA values
        # together with the reason (column '.reason', moved to the attribute "failed_pairs" below)
        safe_pair <- function(i, aln_rows = NULL) {
                res <- tryCatch(
                        one_pair(i, aln_rows),
                        error = function(e) {
                                dnds_na_row(
                                        query_id   = q_ids[i],
                                        subject_id = s_ids[i],
                                        est.method = dnds_est.method,
                                        # (the messages name the pair by the names in its files)
                                        reason     = dnds_pair_reason(conditionMessage(e), q_ids[i], s_ids[i])
                                )
                        }
                )
                # (both one_pair() and dnds_na_row() return a data.table with the column '.reason')
                return(res)
        }
        
        # The pairwise protein alignments (aa_aln_type = "pairwise", Needleman-Wunsch with a substitution matrix)
        # are computed for a chunk of gene pairs with one call of pwalign::pairwiseAlignment() (pairwise_aln_rows()):
        # the same alignments as one call per pair, without an S4 object and a FASTA file read per pair. Not with
        # store_locally = TRUE (pairwise_aln() also writes the alignment with its score there) or the pairwise
        # alignment of 0.4.3 (options(orthologr.legacy_pairwise_aln = TRUE), no substitution matrix).
        batch_aln <- aa_aln_type == "pairwise" && identical(aa_aln_tool, "NW") && !store_locally &&
                !is.null(aa_sub_matrix) && !isTRUE(orthologr_option("orthologr.legacy_pairwise_aln", FALSE))
        if (batch_aln) {
                d <- file.path(orthologr_tmpdir(), "_alignment", "pairwise_aln")
                if (!file.exists(d))
                        dir.create(d, recursive = TRUE)
        }
        # the gene pairs 'idx' (a chunk): a list with the result table of each pair
        run_pairs <- function(idx) {
                rows <- if (batch_aln) pairwise_aln_rows(q_aa[idx], s_aa[idx], aa_sub_matrix) else vector("list", length(idx))
                lapply(seq_along(idx), function(k) safe_pair(idx[k], rows[[k]]))
        }
        
        .reason <- NULL
        
        # Progress (only when not quiet, or when divergence_stratigraphy() asks for it in an interactive session):
        # a message about every tenth of the gene pairs. With several R processes the pairs are handed out in pieces
        # of that size, so that there is something to report in between; the result does not depend on it.
        show_progress <- !quiet || isTRUE(.orthologr_state$progress)
        n_pairs <- nrow(complete_tbl)
        step <- progress_step(n_pairs)
        t_start <- proc.time()[["elapsed"]]
        
        # the progress pieces (a message after each, every 'step' pairs), each split into chunks of at most
        # dnds_chunk_size() pairs, at least one chunk per R process
        pieces <- split(seq_len(n_pairs), ceiling(seq_len(n_pairs) / step))
        
        if (!multicore) {
                dNdS_values <- list()
                for (piece in pieces) {
                        for (chunk in dnds_chunks(piece, 1))
                                dNdS_values <- c(dNdS_values, run_pairs(chunk))
                        if (show_progress && length(dNdS_values) < n_pairs)
                                progress_message(length(dNdS_values), n_pairs, t_start)
                }
        }
        
        if (multicore) {
                # The gene pairs are spread over 'comp_cores' R processes of the 'parallel' package: forked
                # processes where the operating system can fork, otherwise (Windows) processes that are started
                # anew and get the package, the options and the working folder of this call.
                par_cores <- dnds_cluster(comp_cores)
                finished <- FALSE
                on.exit(dnds_cluster_stop(par_cores, finished = finished), add = TRUE)
                
                if (!show_progress || n_pairs <= step) {
                        dNdS_values <- do.call(c, parallel::parLapply(par_cores, dnds_chunks(seq_len(n_pairs), comp_cores), run_pairs))
                } else {
                        dNdS_values <- list()
                        for (piece in pieces) {
                                dNdS_values <- c(dNdS_values, do.call(c, parallel::parLapply(par_cores, dnds_chunks(piece, comp_cores), run_pairs)))
                                if (length(dNdS_values) < n_pairs)
                                        progress_message(length(dNdS_values), n_pairs, t_start)
                        }
                }
                finished <- TRUE
        }
        
        dNdS_tbl <- data.table::rbindlist(dNdS_values, use.names = TRUE, fill = TRUE)
        
        failed_pairs <- as.data.frame(dNdS_tbl[!is.na(.reason), list(query_id, subject_id, reason = .reason)])
        # several methods per gene pair ("ALL"): the method whose estimate failed (NA: the whole pair failed)
        if (identical(dnds_est.method, "ALL") && is.element("method", names(dNdS_tbl))) {
                m <- dNdS_tbl[!is.na(.reason)][["method"]]
                failed_pairs$method <- ifelse(m == "ALL", NA_character_, m)
        }
        dNdS_tbl[ , .reason := NULL]
        
        setkeyv(dNdS_tbl, c("query_id", "subject_id"))
        
        
        if (clean_folders)
                clean_all_folders(c(
                        file.path(orthologr_tmpdir(), "_alignment"),
                        file.path(orthologr_tmpdir(), "_blast_db"),
                        file.path(orthologr_tmpdir(), "_calculation")
                ))
        
        if (nrow(failed_pairs) > 0) {
                orthologr_warn(
                        class = "failed_pairs",
                        nrow(failed_pairs), " of ", nrow(complete_tbl),
                        " gene pairs could not be processed; their dN, dS and dNdS values are NA. ",
                        "The reasons are in attr(<result>, \"failed_pairs\"). First one: ",
                        "'", failed_pairs$query_id[1], "' / '", failed_pairs$subject_id[1], "': ",
                        strsplit(failed_pairs$reason[1], "\n", fixed = TRUE)[[1]][1],
                        call. = FALSE
                )
        }
        
        dNdS_tbl <- tibble::as_tibble(dNdS_tbl)
        attr(dNdS_tbl, "failed_pairs") <- failed_pairs
        # returning the dNdS table as data.table object
        return(dNdS_tbl)
        
}

# The gene pairs 'idx' in consecutive chunks for compute_dnds(): at most dnds_chunk_size() pairs each, and at least
# 'n_workers' chunks (when there are as many pairs), so that every R process gets work. The result does not depend
# on the chunks.
dnds_chunks <- function(idx, n_workers = 1) {
        if (length(idx) == 0)
                return(list())
        size <- max(1, min(dnds_chunk_size(), ceiling(length(idx) / n_workers)))
        unname(split(idx, ceiling(seq_along(idx) / size)))
}

# The number of gene pairs whose protein alignments are computed with one call of pwalign::pairwiseAlignment():
# 100 (from about 25 on the call costs less than 0.1 ms per pair).
dnds_chunk_size <- function() 100L

# What seqinr::write.fasta(sequences = seqs, names = as.list(names), file.out = file) writes for a list of strings -
# '>' and the name on one line, the whole sequence on the next - with one writeLines() call (seqinr writes line by
# line through a connection).
write_pair_fasta <- function(seqs, names, file) {
        seqs <- vapply(seqs, function(x) paste(as.character(x), collapse = ""), character(1), USE.NAMES = FALSE)
        writeLines(as.vector(rbind(paste0(">", unlist(names)), seqs)), file)
}

# The names of the two sequences in the files of a gene pair.
dnds_pair_names <- function() c("query", "subject")

# A message about "gene pair 'query' / 'subject'" with the IDs of the pair in place of the names in its files.
dnds_pair_reason <- function(reason, query_id, subject_id) {
        n <- dnds_pair_names()
        gsub(paste0("gene pair '", n[1], "' / '", n[2], "'"),
             paste0("gene pair '", query_id, "' / '", subject_id, "'"), reason, fixed = TRUE)
}

# A result row of NA values for a gene pair that could not be processed, with the columns substitutionrate()
# returns for 'est.method' plus the reason.
dnds_na_row <- function(query_id, subject_id, est.method, reason) {
        row <- data.table::data.table(
                query_id   = as.character(query_id),
                subject_id = as.character(subject_id),
                dN         = NA_real_,
                dS         = NA_real_,
                dNdS       = NA_real_
        )
        if (!is.element(est.method, c("Comeron", "Li")))
                row[["method"]] <- est.method
        row[[".reason"]] <- as.character(reason)
        return(row)
}

# A cluster of `n` R processes for compute_dnds().
#   "FORK"  : forked processes (not on Windows); they share the data and the state of this R session
#   "PSOCK" : new R processes; orthologr is loaded in each of them and the orthologr options of this session
#             (working folder of the call, legacy switches, ...) are set there
# options(orthologr.cluster_type = "PSOCK") chooses the second kind where both are possible.
dnds_cluster <- function(n) {
        type <- getOption("orthologr.cluster_type", if (.Platform$OS.type == "windows") "PSOCK" else "FORK")
        if (!is.element(type, c("FORK", "PSOCK")))
                orthologr_abort(class = "bad_argument", "options(orthologr.cluster_type = ) must be \"FORK\" or \"PSOCK\".", call. = FALSE)
        if (type == "FORK" && .Platform$OS.type == "windows")
                type <- "PSOCK"
        if (type == "FORK") {
                cl <- parallel::makeForkCluster(n)
                attr(cl, "orthologr_pids") <- tryCatch(unlist(parallel::clusterCall(cl, Sys.getpid)), error = function(e) integer(0))
                return(cl)
        }
        
        cl <- parallel::makePSOCKcluster(n)
        ok <- tryCatch({
                opts <- options()
                opts <- opts[startsWith(names(opts), "orthologr.")]
                # the new processes have a tempdir() of their own: they must work in the folder of this session
                # (it was not handed over when the caller had not set one, and every pair failed with
                # "cannot open the connection")
                opts[["orthologr.tmpdir"]] <- orthologr_tmpdir()
                parallel::clusterCall(cl, function(lib_paths, opts, path) {
                        .libPaths(lib_paths)
                        Sys.setenv(PATH = path)
                        loadNamespace("orthologr")
                        options(opts)
                        TRUE
                }, .libPaths(), opts, Sys.getenv("PATH"))
                TRUE
        }, error = function(e) {
                parallel::stopCluster(cl)
                stop("The R processes for comp_cores = ", n, " could not be started: ", conditionMessage(e), call. = FALSE)
        })
        attr(cl, "orthologr_pids") <- tryCatch(unlist(parallel::clusterCall(cl, Sys.getpid)), error = function(e) integer(0))
        cl
}

# Stop the processes of dnds_cluster(). When the work did not finish (error, interrupt) the workers may still be
# busy with a gene pair - a forked worker only reads the request to stop after its current task - and would go on
# writing into the working folder of the call after the caller has removed it (the folder 'orthologr_dNdS_*' was
# left in tempdir() after an interrupt). They are therefore ended, and given a moment to do so, before the
# caller's on.exit() removes the folder.
dnds_cluster_stop <- function(cl, finished = TRUE) {
        pids <- attr(cl, "orthologr_pids")
        try(parallel::stopCluster(cl), silent = TRUE)
        if (!isTRUE(finished) && length(pids) > 0) {
                try(tools::pskill(pids, tools::SIGTERM), silent = TRUE)
                Sys.sleep(0.3)
                try(tools::pskill(pids, tools::SIGKILL), silent = TRUE)
                Sys.sleep(0.1)
        }
        invisible(TRUE)
}
