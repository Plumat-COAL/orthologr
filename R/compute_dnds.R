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
        
        # The pipeline for gene pair i, in four steps: pair_files() writes the two FASTA files and aligns the
        # proteins, pair_codon() builds the codon alignment, pair_rates() estimates dN/dS, pair_finish() makes the
        # result rows. one_pair() runs them one after the other; run_pairs_staged() runs each step for a chunk of
        # pairs, with one PAL2NAL / KaKs_Calculator process for the chunk where it can.
        #   aln_rows : NULL, or the two rows of the pairwise protein alignment of the pair computed beforehand by
        #              pairwise_aln_rows() (they are written to the file pairwise_aln() would write)
        # (each step is evaluated before the next is called: a step handed over as an argument would be evaluated
        # inside the next one, and the error of PAL2NAL became one of gestimator)
        one_pair <- function(i, aln_rows = NULL) {
                f <- pair_files(i, aln_rows)
                codon_file <- pair_codon(i, f, aln_rows)
                rates <- pair_rates(i, codon_file)
                pair_finish(i, rates)
        }
        
        # step 1: the FASTA files of the pair and its protein alignment; list(name, cds = <CDS FASTA file>,
        # aa = <protein alignment file>)
        pair_files <- function(i, aln_rows = NULL) {
                # the coding sequences and amino acid sequences of the gene pair
                cds_seqs <- list(q_cds[i], s_cds[i])
                aa_seqs <- list(q_aa[i], s_aa[i])
                
                aa_aln_name <- paste0("q", i)
                
                # (P1 of REVIEW-OL-1) the codon alignment is built from the rows in memory (pair_codon(), same
                # condition): the two FASTA files and the protein alignment file would never be read, so they are
                # not written - unless the folder of the call is kept (options(orthologr.keep_tmp))
                if (skip_pair_files && !is.null(aln_rows) && plain_seqs(c(q_cds[i], s_cds[i])))
                        return(list(name = aa_aln_name, cds = NA_character_, aa = NA_character_))
                
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
                list(name = aa_aln_name, cds = cds_session_fasta, aa = aa_aln_file)
        }
        
        # step 2: the codon alignment of the pair (files 'f' of pair_files()); returns the file that holds it
        pair_codon <- function(i, f, aln_rows = NULL) {
                aa_aln_name <- f$name
                # align codon -> cds.aln
                if (codon_in_memory && !is.null(aln_rows) && plain_seqs(c(q_cds[i], s_cds[i]))) {
                        # COAL in process with the rows and coding sequences at hand: the file codon_aln() writes,
                        # without reading the two files of the pair back
                        codon <- coal_codon_rows(aln_rows, c(q_cds[i], s_cds[i]),
                                                 paste0("gene pair '", pair_names[1], "' / '", pair_names[2], "' (",
                                                        aa_aln_name, ")"))
                        names(codon) <- pair_names
                        codon_file <- file.path(orthologr_tmpdir(), "_alignment", "codon_aln",
                                                paste0(aa_aln_name, "_", codon_aln_tool, ".aln"))
                        if (file.exists(codon_file))
                                unlink(codon_file)
                        write_alignment(codon, codon_file, format = "fasta")
                } else {
                        codon_aln(
                                file_aln       = f$aa,
                                file_nuc       = f$cds,
                                tool           = codon_aln_tool,
                                format         = "fasta",
                                codon_aln_name = aa_aln_name,
                                get_aln        = FALSE,
                                quiet          = TRUE
                        )
                }
                file.path(orthologr_tmpdir(), "_alignment", "codon_aln", paste0(aa_aln_name, "_", codon_aln_tool, ".aln"))
        }
        
        # step 3: dN/dS of the codon alignment in 'codon_file'
        pair_rates <- function(i, codon_file) {
                # Comeron's method: what substitutionrate() computes, as a "dnds_values" record (no table per pair)
                if (comeron_values_ok)
                        return(comeron_rates(codon_file, subst_name = paste0("q", i), quiet = TRUE, values = TRUE))
                # compute kaks
                substitutionrate(
                        file           = codon_file,
                        subst_name     = paste0("q", i),
                        kaks_calc_path = kaks_calc_path,
                        est.method     = dnds_est.method,
                        quiet          = TRUE,
                        kaks_tool      = kaks_tool
                )
        }
        
        # step 4: the result rows of the pair from the table 'rates' of pair_rates()
        pair_finish <- function(i, rates) {
                # values of the pair (dnds_values()): they become rows of the table of the chunk (dnds_bind())
                if (inherits(rates, "dnds_values")) {
                        rates$query_id <- as.character(q_ids[i])
                        rates$subject_id <- as.character(s_ids[i])
                        rates$reason <- dnds_pair_reason(rates$reason, q_ids[i], s_ids[i])
                        return(rates)
                }
                # The IDs of the pair are those of the pair table - not what came back from the programs.
                # Alignment programs and file readers change IDs (MAFFT cuts them to 15 characters in CLUSTAL
                # output, fread() turned '007' into 7, KaKs_Calculator labels were split at '-'), and the result
                # was then joined with the hit table by these changed IDs: rows were lost or held IDs that are in
                # no input file.
                kaks_failed <- attr(rates, "kaks_failed")
                # (kaks_tool = "coal": COAL's status of the pair, attr "engine_status" of coal_kaks_file(), as the
                # column '.engine_status', which dNdS() takes out of the table for the column 'status')
                engine_status <- attr(rates, "engine_status")
                rates <- data.table::as.data.table(rates)
                if (!is.null(engine_status))
                        data.table::set(rates, j = ".engine_status", value = engine_status)
                data.table::set(rates, j = "query_id", value = as.character(q_ids[i]))
                data.table::set(rates, j = "subject_id", value = as.character(s_ids[i]))
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
                        error = function(e) pair_failed(i, e)
                )
                # (both one_pair() and dnds_na_row() return a data.table with the column '.reason')
                return(res)
        }
        # the result row of gene pair i whose pipeline stopped with the error 'e'
        pair_failed <- function(i, e) {
                dnds_na_row(
                        query_id   = q_ids[i],
                        subject_id = s_ids[i],
                        est.method = dnds_est.method,
                        # (the messages name the pair by the names in its files)
                        reason     = dnds_pair_reason(conditionMessage(e), q_ids[i], s_ids[i])
                )
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
        # With these alignments and codon_aln_tool = "coal" through coalr, the codon alignment is built from the
        # rows in memory (coal_codon_rows(), what codon_aln() does after reading the two files back).
        codon_in_memory <- batch_aln && identical(codon_aln_tool, "coal") &&
                isTRUE(tryCatch(identical(coal_backend()$type, "coalr"), error = function(e) FALSE))
        skip_pair_files <- codon_in_memory && !isTRUE(keep_tmp_setting()$keep)
        if (codon_in_memory) {
                d <- file.path(orthologr_tmpdir(), "_alignment", "codon_aln")
                if (!file.exists(d))
                        dir.create(d, recursive = TRUE)
        }
        # PAL2NAL and KaKs_Calculator for a chunk of gene pairs in one process each (pal2nal_batch(),
        # kaks_calc_batch()) instead of one process per pair: the same files and tables. A pair that the batch
        # does not compute, or that failed in it, goes through codon_aln() / substitutionrate() on its own, which
        # then also gives its error. Not with the settings of 0.4.3 (orthologr_legacy_options()), which handle one
        # pair at a time as 0.4.3 did, or with options(orthologr.batch_tools = FALSE).
        batch_tools <- isTRUE(getOption("orthologr.batch_tools", TRUE)) &&
                !isTRUE(orthologr_option("orthologr.legacy_pairwise_aln", FALSE))
        batch_p2n <- batch_tools && identical(codon_aln_tool, "pal2nal") && pal2nal_batch_possible()
        # (GMYN too: KaKs_Calculator 2.0 decides whether the header line is written for the gamma methods with a
        # member that is never initialised, 'tempt'; kaks_output_lines() supplies a missing header and drops repeated
        # ones, for one pair and for a chunk alike)
        batch_kc <- batch_tools && identical(kaks_tool, "kaks_calculator") &&
                is.element(dnds_est.method, setdiff(kaks_calculator_methods(), "kaks_calc"))
        
        # Comeron's method (substitutionrate()'s own checks of kaks_tool passed in dNdS()): pair_rates() gives
        # the values of a pair (comeron_rates(values = TRUE)); with the codon alignments in memory (the default
        # route: pairwise alignments, COAL through coalr, the files of the pair not kept) they are handed to
        # gestimator as strings (gestimator_pairs(), one call per chunk) instead of a file that is written and read
        # back - unless keep_tmp keeps the files, or a pair cannot go this way (then it takes the file route).
        comeron_values_ok <- identical(dnds_est.method, "Comeron") && is.element(kaks_tool, c("kaks_calculator", "coal"))
        comeron_in_memory <- comeron_values_ok && skip_pair_files
        
        # the gene pairs 'idx' (a chunk): a list with the result tables of the pairs (one table for the pairs
        # whose values were computed together, dnds_bind())
        run_pairs <- function(idx) {
                rows <- if (batch_aln) pairwise_aln_rows(q_aa[idx], s_aa[idx], aa_sub_matrix) else vector("list", length(idx))
                if (batch_p2n || batch_kc)
                        return(dnds_bind(run_pairs_staged(idx, rows)))
                if (comeron_in_memory)
                        return(dnds_bind(run_pairs_comeron(idx, rows)))
                dnds_bind(lapply(seq_along(idx), function(k) safe_pair(idx[k], rows[[k]])))
        }
        
        # run_pairs() for Comeron's method with the codon alignments in memory: the codon alignment of every pair of
        # the chunk (as pair_codon() builds it), then gestimator for all of them with one call. A pair whose
        # alignment is not a plain codon alignment, or whose numbers gestimator writes in another form than plain
        # decimals, goes through the files as before (safe_pair(): the same values, and its own error message).
        run_pairs_comeron <- function(idx, rows) {
                n <- length(idx)
                res <- vector("list", n)
                aln <- vector("list", n)
                for (k in seq_len(n)) {
                        i <- idx[k]
                        if (is.null(rows[[k]]) || !plain_seqs(c(q_cds[i], s_cds[i])))
                                next
                        aln[[k]] <- tryCatch(
                                coal_codon_rows(rows[[k]], c(q_cds[i], s_cds[i]),
                                                paste0("gene pair '", pair_names[1], "' / '", pair_names[2], "' (q", i, ")")),
                                error = function(e) e)
                        if (inherits(aln[[k]], "error")) {
                                res[[k]] <- pair_failed(i, aln[[k]])
                                aln[k] <- list(NULL)
                        } else if (!comeron_plain_alignment(aln[[k]])) {
                                aln[k] <- list(NULL)
                        }
                }
                todo <- which(!vapply(aln, is.null, TRUE))
                vals <- if (length(todo) > 0) tryCatch(
                        gestimator_pairs(vapply(aln[todo], `[`, "", 1), vapply(aln[todo], `[`, "", 2),
                                         legacy_saturation = comeron_legacy_saturation(),
                                         twosubs_libsequence = identical(comeron_twosubs(), "libsequence"),
                                         min_sites = comeron_min_sites()),
                        error = function(e) rep(NA_character_, length(todo))) else character(0)
                for (j in seq_along(todo)) {
                        k <- todo[j]
                        f <- if (is.na(vals[j])) NULL else strsplit(vals[j], "\t", fixed = TRUE)[[1]]
                        if (length(f) == 3 && gestimator_plain_numbers(f))
                                res[[k]] <- pair_finish(idx[k], comeron_values(list(dN = as.numeric(f[1]), dS = as.numeric(f[2]),
                                                                                    dNdS = as.numeric(f[3]))))
                }
                for (k in which(vapply(res, is.null, TRUE)))
                        res[[k]] <- safe_pair(idx[k], rows[[k]])
                res
        }
        
        # run_pairs() step by step: every step of one_pair() for all pairs of the chunk before the next step, the
        # codon alignments with PAL2NAL and the estimates of KaKs_Calculator with one process for the chunk. A pair
        # whose step stops with an error is the row of pair_failed() with that error, as in safe_pair().
        run_pairs_staged <- function(idx, rows) {
                n <- length(idx)
                res <- vector("list", n)
                files <- vector("list", n)
                for (k in seq_len(n)) {
                        files[[k]] <- tryCatch(pair_files(idx[k], rows[[k]]), error = function(e) e)
                        if (inherits(files[[k]], "error"))
                                res[[k]] <- pair_failed(idx[k], files[[k]])
                }
                todo <- which(vapply(res, is.null, TRUE))
                # the codon alignments
                codon <- rep(NA_character_, n)
                if (batch_p2n && length(todo) > 0) {
                        done <- pal2nal_batch(file_aln = vapply(files[todo], `[[`, "", "aa"),
                                              file_nuc = vapply(files[todo], `[[`, "", "cds"),
                                              codon_aln_name = vapply(files[todo], `[[`, "", "name"))
                        codon[todo[done]] <- file.path(orthologr_tmpdir(), "_alignment", "codon_aln",
                                                       paste0(paste0("q", idx[todo[done]]), "_", codon_aln_tool, ".aln"))
                }
                for (k in todo[is.na(codon[todo])]) {
                        x <- tryCatch(pair_codon(idx[k], files[[k]], rows[[k]]), error = function(e) e)
                        if (inherits(x, "error"))
                                res[[k]] <- pair_failed(idx[k], x)
                        else
                                codon[k] <- x
                }
                todo <- which(vapply(res, is.null, TRUE))
                # the estimates
                rates <- vector("list", n)
                if (batch_kc && length(todo) > 0)
                        rates[todo] <- kaks_calc_batch(files = codon[todo], subst_name = paste0("q", idx[todo]),
                                                       est.method = dnds_est.method, kaks_calc_path = kaks_calc_path)
                for (k in todo) {
                        res[[k]] <- tryCatch({
                                # (a table, a "dnds_values" record of kaks_calc_batch(), or NULL)
                                r <- rates[[k]]
                                if (is.null(r))
                                        r <- pair_rates(idx[k], codon[k])
                                pair_finish(idx[k], r)
                        }, error = function(e) pair_failed(idx[k], e))
                }
                res
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
        
        # (the tables of a chunk hold the rows of one or more pairs: the pairs done are counted with 'done')
        if (!multicore) {
                dNdS_values <- list()
                done <- 0L
                for (piece in pieces) {
                        for (chunk in dnds_chunks(piece, 1))
                                dNdS_values <- c(dNdS_values, run_pairs(chunk))
                        done <- done + length(piece)
                        if (show_progress && done < n_pairs)
                                progress_message(done, n_pairs, t_start)
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
                        done <- 0L
                        for (piece in pieces) {
                                dNdS_values <- c(dNdS_values, do.call(c, parallel::parLapply(par_cores, dnds_chunks(piece, comp_cores), run_pairs)))
                                done <- done + length(piece)
                                if (done < n_pairs)
                                        progress_message(done, n_pairs, t_start)
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
# 100. (A call costs about 3 - 5 ms besides the alignments: per pair 0.6 ms with 5 pairs per call, 0.13 ms with 25,
# 0.05 ms with 100.)
dnds_chunk_size <- function() 100L

# What seqinr::write.fasta(sequences = seqs, names = as.list(names), file.out = file) writes for a list of strings -
# '>' and the name on one line, the whole sequence on the next - with one writeLines() call (seqinr writes line by
# line through a connection).
write_pair_fasta <- function(seqs, names, file) {
        seqs <- vapply(seqs, function(x) paste(as.character(x), collapse = ""), character(1), USE.NAMES = FALSE)
        writeLines(as.vector(rbind(paste0(">", unlist(names)), seqs)), file)
}

# Are the sequences 'x' such that writing them with write_pair_fasta() and reading them back with seqinr::read.fasta()
# gives them unchanged (letters, '*' and '-' only; not NA, not empty)?
plain_seqs <- function(x) {
        !anyNA(x) && all(grepl("^[A-Za-z*-]+\\z", x, perl = TRUE, useBytes = TRUE))
}

# The names of the two sequences in the files of a gene pair.
dnds_pair_names <- function() c("query", "subject")

# A message about "gene pair 'query' / 'subject'" with the IDs of the pair in place of the names in its files.
dnds_pair_reason <- function(reason, query_id, subject_id) {
        n <- dnds_pair_names()
        gsub(paste0("gene pair '", n[1], "' / '", n[2], "'"),
             paste0("gene pair '", query_id, "' / '", subject_id, "'"), reason, fixed = TRUE)
}

# The result tables of a chunk of gene pairs from 'res' (per pair, in order: a table, or a "dnds_values" record of
# pair_finish()): consecutive records with the same columns become one table, with the rows, columns and types the
# tables of pair_finish() would have (query_id, subject_id, dN, dS, dNdS, [method,] .reason), so that rbindlist()
# of the list gives the same result. (A table per gene pair - data.table(), setkey(), as.data.table() and set() -
# took about a fifth of the run time of the YN and Comeron routes.)
dnds_bind <- function(res) {
        n <- length(res)
        if (n == 0)
                return(res)
        light <- vapply(res, inherits, TRUE, what = "dnds_values")
        if (!any(light))
                return(res)
        has_method <- vapply(res, function(r) inherits(r, "dnds_values") && !is.null(r$method), TRUE)
        # runs of consecutive records of the same kind
        grp <- cumsum(c(TRUE, light[-1] != light[-n] | has_method[-1] != has_method[-n] | !light[-1]))
        out <- vector("list", max(grp))
        for (g in unique(grp)) {
                at <- which(grp == g)
                if (!light[at[1]]) {
                        out[[g]] <- res[[at[1]]]
                        next
                }
                r <- res[at]
                len <- vapply(r, function(x) length(x$dN), 1L)
                col <- function(name) unlist(lapply(r, `[[`, name), use.names = FALSE)
                tbl <- list(query_id   = rep(vapply(r, `[[`, "", "query_id"), len),
                            subject_id = rep(vapply(r, `[[`, "", "subject_id"), len),
                            dN         = as.numeric(col("dN")),
                            dS         = as.numeric(col("dS")),
                            dNdS       = as.numeric(col("dNdS")))
                if (has_method[at[1]])
                        tbl$method <- col("method")
                tbl$.reason <- as.character(col("reason"))
                out[[g]] <- data.table::setDT(tbl)
        }
        out
}

# Is the codon alignment 'x' (two strings) one that gestimator_pairs() computes as gestimator() computes it after
# write_alignment() and check_codon_alignment_file() (letters, '*' and '-'; the same length, a multiple of 3 and at
# least 3)?
comeron_plain_alignment <- function(x) {
        length(x) == 2 && is.character(x) && plain_seqs(x) && nchar(x[1]) == nchar(x[2]) && nchar(x[1]) >= 3 &&
                nchar(x[1]) %% 3 == 0
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
