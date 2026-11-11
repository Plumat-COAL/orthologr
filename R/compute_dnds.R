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
#' @param aa_aln_params a character string specifying additional parameters that shall be passed to the multiple alignment system call.
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used for codon alignments. Default is \code{codon_aln_tool} = \code{"pal2nal"}.
#' @param dnds_est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' @param aa_sub_matrix the substitution matrix for pairwise amino acid alignments (\code{aa_aln_type = "pairwise"}); see \code{\link{pairwise_aln}}.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
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
#' @references \url{https://www.r-bloggers.com/2013/08/the-wonders-of-foreach/}
#' @seealso \code{\link{multi_aln}}, \code{\link{substitutionrate}}, \code{\link{dNdS}}
#' @import foreach 
#' @import data.table
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
                         aa_sub_matrix   = "BLOSUM62"){
        
        if (comp_cores > parallel::detectCores())
                stop("You assigned more cores to the comp_cores argument than are availible on your machine.")
        
        
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
        
        # the complete pipeline for gene pair i: write the two FASTA files, align the proteins,
        # build the codon alignment, estimate dN/dS
        one_pair <- function(i) {
                # storing the query gene id and subject gene id of the orthologous gene pair
                orthologs_names <-
                        list(complete_tbl[i , query_id], complete_tbl[i, subject_id])
                
                # storing the query CDS sequence and subject CDS sequence of the orthologous gene pair
                cds_seqs <-
                        list(complete_tbl[i, query_cds], complete_tbl[i, subject_cds])
                
                # storing the query amino acid sequence and subject amino acid sequence of the orthologous gene pair
                aa_seqs <-
                        list(complete_tbl[i, query_aa], complete_tbl[i, subject_aa])
                
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
                
                seqinr::write.fasta(
                        sequences = cds_seqs,
                        names     = orthologs_names,
                        file.out  = cds_session_fasta
                )
                
                # create aa fasta of orthologous gene pair having session name: 'aa_aln_name'
                aa_session_fasta <-
                        file.path(
                                orthologr_tmpdir(),
                                "_alignment",
                                "orthologs",
                                "AA",
                                paste0(aa_aln_name, "_aa.fasta")
                        )
                
                seqinr::write.fasta(
                        sequences = aa_seqs,
                        names     = orthologs_names,
                        file.out  = aa_session_fasta
                )
                
                if (aa_aln_type == "multiple") {
                        # align aa -> <aa_aln_tool>.aln
                        multi_aln(
                                file           = aa_session_fasta,
                                tool           = aa_aln_tool,
                                get_aln        = FALSE,
                                multi_aln_name = aa_aln_name,
                                path           = aa_aln_path,
                                quiet          = quiet
                        )
                        
                        aa_aln_file <-
                                file.path(
                                        orthologr_tmpdir(),
                                        "_alignment",
                                        "multi_aln",
                                        paste0(aa_aln_name, "_", aa_aln_tool, ".aln")
                                )
                }
                
                if (aa_aln_type == "pairwise") {
                        # align aa -> <aa_aln_tool>.aln
                        pairwise_aln(
                                file              = aa_session_fasta,
                                tool              = aa_aln_tool,
                                seq_type          = "protein",
                                get_aln           = FALSE,
                                pairwise_aln_name = aa_aln_name,
                                path              = aa_aln_path,
                                quiet             = quiet,
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
                        quiet          = quiet
                )
                
                # compute kaks
                substitutionrate(
                        file           = file.path(
                                orthologr_tmpdir(),
                                "_alignment",
                                "codon_aln",
                                paste0(aa_aln_name, "_", codon_aln_tool, ".aln")
                        ),
                        subst_name     = aa_aln_name,
                        kaks_calc_path = kaks_calc_path,
                        est.method     = dnds_est.method,
                        quiet          = quiet
                )
        }
        
        # one gene pair must not end the run: a pair whose pipeline fails is returned as a row of NA values
        # together with the reason (column '.reason', moved to the attribute "failed_pairs" below)
        safe_pair <- function(i) {
                res <- tryCatch(
                        one_pair(i),
                        error = function(e) {
                                dnds_na_row(
                                        query_id   = complete_tbl[i, query_id],
                                        subject_id = complete_tbl[i, subject_id],
                                        est.method = dnds_est.method,
                                        reason     = conditionMessage(e)
                                )
                        }
                )
                res <- data.table::as.data.table(res)
                if (!is.element(".reason", names(res)))
                        res[ , .reason := NA_character_]
                return(res)
        }
        
        .reason <- NULL
        
        if (!multicore) {
                dNdS_values <- vector(mode = "list", length = nrow(complete_tbl))
                for (i in seq_len(nrow(complete_tbl)))
                        dNdS_values[[i]] <- safe_pair(i)
        }
        
        if (multicore) {
                ### Parallellizing the sampling process using the 'doParallel' and 'parallel' package
                ### register all given cores for parallelization
                par_cores <- parallel::makeForkCluster(comp_cores)
                on.exit(parallel::stopCluster(par_cores), add = TRUE)
                doParallel::registerDoParallel(par_cores)
                
                dNdS_values <-
                        foreach::foreach(
                                i              = seq_len(nrow(complete_tbl)),
                                .packages      = c("seqinr", "data.table"),
                                .errorhandling = "stop"
                        ) %dopar% {
                                safe_pair(i)
                        }
        }
        
        dNdS_tbl <- data.table::rbindlist(dNdS_values, use.names = TRUE, fill = TRUE)
        
        failed_pairs <- as.data.frame(dNdS_tbl[!is.na(.reason), list(query_id, subject_id, reason = .reason)])
        dNdS_tbl[ , .reason := NULL]
        
        setkeyv(dNdS_tbl, c("query_id", "subject_id"))
        
        
        if (clean_folders)
                clean_all_folders(c(
                        file.path(orthologr_tmpdir(), "_alignment"),
                        file.path(orthologr_tmpdir(), "_blast_db"),
                        file.path(orthologr_tmpdir(), "_calculation")
                ))
        
        if (nrow(failed_pairs) > 0) {
                warning(
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
