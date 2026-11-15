#' @title Perform Divergence Stratigraphy
#' @description This function takes a query organism and performs
#' divergence stratigraphy (Quint et al.,2012 ; Drost et al. 2015) against a
#' closely related subject organism.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param aligner a character string specifying the sequence aligner. The options are \code{diamond} and \code{blast}.
#' The option \code{diamond} uses DIAMOND2 which is up to 10 000X folds faster than BLASTP while retaining the sensitivity of BLASTP.
#' Thus, the default is \code{aligner} = \code{diamond}.
#' @param sensitivity_mode specify the level of alignment sensitivity, when using DIAMOND2. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' \itemize{
#' \item \code{sensitivity_mode = "fast"} (Default): fast alignment mode, designed for finding hits of >90\% identity.
#' \item \code{sensitivity_mode = "faster"}: fastest alignment mode, least sensitive.
#' \item \code{sensitivity_mode = "default"}: the default mode of DIAMOND2 itself (no sensitivity option is passed), designed for finding hits of >70\% identity.
#' \item \code{sensitivity_mode = "mid-sensitive"}: between the default mode and the sensitive mode.
#' \item \code{sensitivity_mode = "sensitive"}: full sensitivity for hits of >40\% identity.
#' \item \code{sensitivity_mode = "more-sensitive"}: more sensitive than the sensitive mode.
#' \item \code{sensitivity_mode = "very-sensitive"}: sensitive alignment mode.
#' \item \code{sensitivity_mode = "ultra-sensitive"}: most sensitive alignment mode (sensitivity as high as BLASTP).
#' }
#' Any other value is an error.
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND or BLAST hit detection.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (DIAMOND or BLAST reciprocal best hit).
#' Available methods are: "BH" (DIAMOND or BLAST best hit), "RBH" (DIAMOND or BLAST reciprocal best hit).
#' @param dnds_est.method the dNdS estimation method that shall be used.
#' Options are:
#' \itemize{
#' \item \code{dnds_est.method = "Comeron"} (Default): Comeron's method (1995)
#' \item \code{dnds_est.method = "Li"}: Li's method (1993)
#' \item \code{dnds_est.method = "NG"}: Nei, M. and Gojobori, T. (1986)
#' \item \code{dnds_est.method = "LWL"}: Li, W.H., et al. (1985)
#' \item \code{dnds_est.method = "LPB"}: Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' \item \code{dnds_est.method = "MLWL"}: (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' \item \code{dnds_est.method = "YN"}: Yang, Z. and Nielsen, R. (2000)
#' \item \code{dnds_est.method ="MYN"} (Modified YN): Zhang, Z., et al. (2006)
#' }
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param aligner_path a character string specifying the path to the DIAMOND or BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform
#'  parallel computations on a multicore machine.
#' @param dnds.threshold a numeric value specifying the dnds threshold for genes that shall be retained.
#' Hence all genes having a dNdS value <= \code{dnds.threshold} are retained. Default is \code{dnds.threshold} = 2.
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out to the console.
#' @param clean_folders a logical value specifying whether the internal folder structure shall be deleted (cleaned) after
#'  processing this function. Default is \code{clean_folders} = \code{FALSE}.
#' @param ds.values a logical value specifying whether divergence stratum values (ds values) or dNdS values shall be returned
#' by \code{divergence_stratigraphy}. Default is \code{ds.values} = \code{TRUE}.
#' @param subject.id a logical value indicating whether \code{query_id} AND \code{subject_id} should be returned.
#' @param n_quantile a numeric value specifying the number of quantiles that should be returned. Many equal dN/dS
#' values (e.g. dN/dS = 0 for more than \code{1/n_quantile} of the genes of two closely related genomes) no longer
#' stop the function after the run: see the section on ties in \code{\link{divergence_map}}.
#' @param dS_zero what to do with gene pairs whose dS is 0 (dN/dS not defined): \code{"drop"} (default) leaves
#' them out of the divergence map, as every released version of orthologr did - now with a message, and listed in
#' \code{attr(result, "left_out")}; \code{"keep"} keeps them as rows with the stratum (or dNdS value) \code{NA}
#' (only choose this when the code that uses the map can deal with \code{NA} strata).
#' See \code{\link{filter_dNdS}}; \code{options(orthologr.dS_zero = "keep")} sets the default.
#' @return A \code{data.frame} with the divergence stratum (column \code{DS}; or \code{dNdS} with
#' \code{ds.values = FALSE}) and the \code{query_id} (and \code{subject_id}) of every gene pair that has a dN/dS
#' value of at most \code{dnds.threshold}. With the default \code{dS_zero = "drop"} no stratum is \code{NA}, as in
#' earlier versions; the function always returns this type of table.
#'
#' Two attributes come with it:
#' \itemize{
#' \item \code{"dNdS_tbl"}: the complete dN/dS table of the run as returned by \code{\link{dNdS}},
#' before filtering and binning, so that a finished run can be binned again (\code{\link{divergence_map}}) without
#' computing it again;
#' \item \code{"left_out"}: the gene pairs of that table that are not in the map (\code{query_id},
#' \code{subject_id}, \code{dN}, \code{dS}, \code{dNdS}, \code{reason}): pairs without dN or dS, pairs with
#' dS = 0 (identical coding sequences, or nonsynonymous changes only), pairs above \code{dnds.threshold}. A message
#' gives their numbers. \code{\link{divergence_left_out}} returns this table.
#' }
#' @details Introduced by Quint et al., 2012 and extended in Drost et al. 2015, divergence stratigraphy
#'  is the process of quantifying the selection pressure (in terms of amino acid sequence divergence) acting on
#'  orthologous genes between closely related species. The resulting sequence divergence map (short divergence map),
#'  stores the divergence stratum in the first column and the query_id of inferred orthologous genes in the second column.
#'  
#'  Following steps are performed to obtain a standard divergence map based on divergence_stratigraphy:
#'  
#'  Divergence Stratigraphy:
#'  
#'  \itemize{
#'  \item 1) Orthology Inference using DIAMOND or BLAST reciprocal best hit ("RBH") based on blastp
#'  
#'  \item 2) Pairwise global amino acid alignments of orthologous genes using the Needleman-Wunsch algorithm
#'  
#'  \item 3) Codon alignments of orthologous genes using PAL2NAL
#'  
#'  \item 4) dNdS estimation using Comeron's method (1995)
#'  
#'  \item 5) Assigning dNdS values to divergence strata (deciles)
#' }
#' @note Although this function has been heavily optimized and parallelized, performing
#' Divergence Stratigraphy using two genomes will take some computation time.
#' 
#' In our experience performing Divergence Stratigraphy using two genomes (one query and one subject genome)
#' on an 8 core machine can take up to 1,5 - 2 hours.
#'   
#' @author Hajk-Georg Drost and Jaruwatana Sodai Lotharukpong
#' @references
#'  
#'  Quint M et al. (2012). A transcriptomic hourglass in plant embryogenesis. Nature (490): 98-101.
#'  
#'  Drost HG et al. (2015). Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012
#'  
#' @examples \dontrun{
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE)
#'       
#'       
#'       
#'  # performing standard divergence stratigraphy using the aligner_path argument to specify
#'  # the path to diamond
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5", 
#'       ortho_detection = "RBH",
#'       aligner_path    = "path/to/diamond",
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE)
#'  
#'  
#'  
#'  # Divergence Stratigraphy can also be performed in parallel 
#'  # (on a multicore machine) using the 'comp_cores' argument
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5", 
#'       ortho_detection = "RBH", 
#'       comp_cores      = 2, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE)
#'  
#'  
#'  
#'  
#'  # in case you want a divergence map with DS values and query_id and subject_id
#'  # you can specify subject.id = TRUE
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE,
#'       subject.id      = TRUE)
#'       
#'       
#'       
#'       
#'  # in case you want a divergence map with KaKs values and query_id and subject_id
#'  # you can specify ds.values = FALSE and subject.id = TRUE
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE,
#'       subject.id      = TRUE)
#'  
#'       
#'       
#'  # in case you want a divergence map with divergence stratum as quintile (5-quantile) values 
#'  # or any other N-quantile values rather than the default decile (10-quantile) values,
#'  # you can specify this with n_quantile.
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE,
#'       n_quantile      = 5)
#'       
#'  }
#'  
#' @seealso 
#' \code{\link{dNdS}}, \code{\link{divergence_map}},\code{\link{substitutionrate}}, \code{\link{multi_aln}},
#'   \code{\link{codon_aln}}, \code{\link{diamond_best}}, \code{\link{diamond_rec}},
#'   \code{\link{codon_aln}}, \code{\link{blast_best}}, \code{\link{blast_rec}}, 
#'   \code{\link{map_generator_dnds}}
#' @param best_hit_rule how the best hit of a query is chosen among its hits:
#' \itemize{
#' \item \code{"bitscore"} (default since orthologr 0.5.0): the hit with the highest bit score; of several, the one
#' with the smallest e-value, then the longest alignment, then the highest percent identity, then the alphabetically
#' first \code{subject_id}. The choice does not depend on the order of the hit table.
#' \item \code{"evalue"}: the rule of orthologr <= 0.4.3 - the hit with the smallest e-value; of several, the one with
#' the longest alignment; of several, the alphabetically first \code{subject_id}. The bit score is not used.
#' }
#' \strong{This default changed in 0.5.0 and changes which orthologs are returned.} DIAMOND and BLAST report an
#' e-value of 0 for every strong hit of a long protein, so that under the old rule the alignment length decided
#' among them and a paralog could win over the true ortholog. Measured (DIAMOND, reciprocal best
#' hits): human against mouse, 19,746 x 21,584 genes, compared with the orthologs of Ensembl Compara - old rule
#' 16,458 pairs, 330 of them (2.0\%) pairing a human gene with another mouse gene than any Compara ortholog, 93.0\% of
#' the Compara one-to-one orthologs found; bit-score rule 16,827 pairs, 68 (0.4\%), 96.8\%; 3.2\% of the best hits
#' change. On the 1,000-gene example of the package 12 of 986 best hits change and there are 983 instead of 972
#' reciprocal best hits (970 in common).
#' \code{options(orthologr.best_hit_rule = "evalue")} restores the old rule for all functions of a session.
#' @export
divergence_stratigraphy <- function(query_file, 
                                    subject_file, 
                                    aligner         = "diamond",
                                    eval            = "1E-5",
                                    ortho_detection = "RBH",
                                    dnds_est.method = "Comeron",
                                    delete_corrupt_cds = FALSE,
                                    aligner_path      = NULL, 
                                    comp_cores      = 1,
                                    dnds.threshold  = 2,
                                    store_locally   = FALSE,
                                    quiet           = FALSE, 
                                    clean_folders   = FALSE, 
                                    ds.values       = TRUE,
                                    subject.id      = FALSE,
                                    n_quantile      = 10,
                                    sensitivity_mode = "fast",
                                    dS_zero         = getOption("orthologr.dS_zero", "drop"),
                                    best_hit_rule   = getOption("orthologr.best_hit_rule", "bitscore")){
        
        best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)
        
        if (!is.ortho_detection_method(ortho_detection))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        if (!is.element(aligner, c("blast", "diamond")))
                stop("Please choose an aligner that is supported by this function.")
        
        message("Proceeding with the aligner: ", aligner)
        
        # recommend users to choose DIAMOND2
        if (aligner == "blast")
                message("We recommend using DIAMOND2 for fast and sensitive alignment", "\n","P.S. the parameter 'sensitivity_mode' doesn't apply")
        
        message("Running Divergence Stratigraphy ...")
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- subject_id <- dNdS <- NULL
        
        # checked before the long run, not after it
        check_n_quantile(n_quantile)
        if (!is.character(dS_zero) || length(dS_zero) != 1 || !is.element(dS_zero, c("keep", "drop")))
                stop("dS_zero must be \"keep\" or \"drop\".", call. = FALSE)
        
        dNdS_full <- dNdS( query_file      = query_file,
                           subject_file    = subject_file,
                           aligner         = aligner,
                           sensitivity_mode = sensitivity_mode,
                           aligner_path    = aligner_path, 
                           ortho_detection = ortho_detection,
                           delete_corrupt_cds = delete_corrupt_cds,
                           aa_aln_type     = "pairwise", 
                           aa_aln_tool     = "NW",
                           codon_aln_tool  = "pal2nal", 
                           dnds_est.method = dnds_est.method,
                           eval            = eval,
                           comp_cores      = comp_cores,
                           store_locally   = store_locally,
                           quiet           = quiet,
                           best_hit_rule   = best_hit_rule )
        
        dNdS_tbl <- filter_dNdS(dNdS_full, dnds.threshold = dnds.threshold, dS_zero = dS_zero)
        
        # the gene pairs that are not in the map, and why
        left_out <- divergence_left_out_table(dNdS_full, dNdS_tbl, dnds.threshold)
        if (nrow(left_out) > 0) {
                tab <- table(left_out$reason)
                message(nrow(left_out), " of ", nrow(dNdS_full), " gene pairs are not in the divergence map: ",
                        paste0(tab, " with ", names(tab), collapse = "; "),
                        ". They are listed in attr(<result>, \"left_out\") (see divergence_left_out()).")
        }
        
        if (ds.values) {
                # divergence map: standard = col1: divergence stratum, col2: query_id
                # The dN/dS run is finished at this point and must not be lost: should the binning fail, a map
                # without rows is returned with a warning (always the documented type of table), and the dN/dS
                # table of the run is in its attribute "dNdS_tbl".
                dm_tbl <- tryCatch(
                        divergence_map( dNdS_tbl = dNdS_tbl , subject.id = subject.id , n_quantile = n_quantile),
                        error = function(e) e
                )
                if (inherits(dm_tbl, "error")) {
                        warning("The dN/dS values could not be sorted into divergence strata (", conditionMessage(dm_tbl),
                                "). A divergence map without rows is returned; the dN/dS table of the run is in ",
                                "attr(<result>, \"dNdS_tbl\"): use filter_dNdS() and divergence_map() on it.",
                                call. = FALSE)
                        dm_tbl <- data.frame(divergence_strata = integer(0), query_id = character(0), stringsAsFactors = FALSE)
                        if (subject.id)
                                dm_tbl$subject_id <- character(0)
                }
        }
        
        if (!ds.values) {
                # (with dS_zero = "keep" the rows with dNdS = NA are the gene pairs with dS = 0)
                if (!subject.id)
                        dm_tbl <- dplyr::select(dNdS_tbl, dNdS,query_id)
                if (subject.id)
                        dm_tbl <- dplyr::select(dNdS_tbl, dNdS,query_id,subject_id)
        }
        
        if (clean_folders)
                clean_all_folders(c(file.path(orthologr_tmpdir(),"_alignment"), file.path(orthologr_tmpdir(),"_blast_db"), file.path(orthologr_tmpdir(),"_calculation")))
        
        DivergenceMap.DF <- as.data.frame(dm_tbl)
        
        if (ds.values)
                colnames(DivergenceMap.DF)[1] <- "DS"
        if (!ds.values)
                colnames(DivergenceMap.DF)[1] <- "dNdS"
        
        # the complete dN/dS table of the run (before filtering and binning) comes with the result
        attr(DivergenceMap.DF, "dNdS_tbl") <- dNdS_full
        attr(DivergenceMap.DF, "left_out") <- left_out
        
        message("Divergence Stratigraphy completed successfully.")
        return( DivergenceMap.DF )
}

# The gene pairs of the dN/dS table `full` that are not in the filtered table `kept`, with the reason.
divergence_left_out_table <- function(full, kept, dnds.threshold) {
        full <- as.data.frame(full)
        why <- dnds_filter_reason(full, dnds.threshold)
        key <- function(x) paste(x$query_id, x$subject_id, sep = "\r")
        out <- !is.element(key(full), key(as.data.frame(kept)))
        why[out & is.na(why)] <- "another reason"
        res <- full[out, intersect(c("query_id", "subject_id", "dN", "dS", "dNdS"), names(full)), drop = FALSE]
        res$reason <- why[out]
        rownames(res) <- NULL
        res
}

#' @title Gene pairs that are not in a divergence map
#' @description \code{\link{divergence_stratigraphy}} leaves gene pairs without a usable dN/dS value out of the
#' divergence map: pairs for which dN or dS could not be estimated, pairs with dS = 0 (dN/dS is not defined: the
#' coding sequences are identical, or differ by nonsynonymous changes only) and pairs with dN/dS above
#' \code{dnds.threshold}. This function returns them.
#' @param x a divergence map returned by \code{\link{divergence_stratigraphy}}.
#' @return A \code{data.frame} with \code{query_id}, \code{subject_id}, \code{dN}, \code{dS}, \code{dNdS} and
#' \code{reason}, one row per gene pair that is in the dN/dS table of the run (\code{attr(x, "dNdS_tbl")}) but not in
#' the map. It has no rows when every pair is in the map.
#' @details Genes with identical coding sequences in both species are among the most conserved genes of a
#' comparison. They are not in the map because their dN/dS is 0/0. To give them a stratum anyway - for example the
#' lowest one - take them from this table and add them to the map; \code{divergence_stratigraphy(dS_zero = "keep")}
#' keeps them in the map with the stratum \code{NA}.
#' @examples \dontrun{
#' dm <- divergence_stratigraphy(
#'   query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' divergence_left_out(dm)
#' }
#' @seealso \code{\link{divergence_stratigraphy}}, \code{\link{filter_dNdS}}
#' @export
divergence_left_out <- function(x) {
        res <- attr(x, "left_out")
        if (is.null(res))
                stop("'x' has no attribute \"left_out\": it must be a divergence map returned by divergence_stratigraphy() ",
                     "of orthologr >= 0.5.0.", call. = FALSE)
        res
}
