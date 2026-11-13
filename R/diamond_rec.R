#' @title Perform a DIAMOND2 reciprocal best hit (RBH) search
#' @description This function performs a DIAMOND2 search (reciprocal best hit) of a given set of protein sequences against a second
#' set of protein sequences and vice versa.
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
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet. 
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND2 hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param path a character string specifying the path to the DIAMOND2 program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore DIAMOND2 computations.
#' @param diamond_params a character string listing the input paramters that shall be passed to the executing DIAMOND2 program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running DIAMOND2.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the DIAMOND2 output shall be stored. E.g. \code{save.output} = \code{getwd()}
#' to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @author Jaruwatana Sodai Lotharukpong
#' @details Given a set of protein sequences A and a different set of protein sequences B,
#' first a best hit diamond search is being performed from A to B: diamond(A,B) and afterwards
#' a best hit diamond search is being performed from B to A: diamond(B,A). Only protein sequences
#' that were found to be best hits in both directions are retained and returned.
#' 
#' This function gives the same output as \code{\link{blast_rec}} while being much much faster.
#'
#' This function can be used to perform orthology inference using DIAMOND2 best reciprocal hit methodology.
#' 
#' @references
#' 
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
#'
#' @examples \dontrun{
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' diamond_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'           
#'           
#'           
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' # starting with protein sequences
#' diamond_rec(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'             seq_type     = "protein")
#' 
#' 
#' 
#' # save the DIAMOND2 output file to the current working directory
#' diamond_rec(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'             seq_type     = "protein",
#'             save.output  = getwd())
#' 
#' 
#' 
#' # use multicore processing
#' diamond_rec(query_file    = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'             subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'             comp_cores   = 2)
#'            
#'            
#'            
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' # and external path to diamond
#' diamond_rec(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'             path         = "path/to/diamond")
#'           
#'           
#'           
#' }
#'
#' @return A data.table as returned by the \code{\link{diamond}} function, storing the geneids
#' of orthologous genes (reciprocal best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{diamond}}, \code{\link{diamond_best}}, \code{\link{set_diamond}}, \code{\link{blast_rec}}
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
#' among them and a paralog could win over the true ortholog. @@BHNUM@@
#' \code{options(orthologr.best_hit_rule = "evalue")} restores the old rule for all functions of a session.
#' @export
#' 
diamond_rec <- function(
                query_file, 
                subject_file, 
                seq_type        = "cds",
                format          = "fasta", 
                diamond_algorithm = "blastp",
                sensitivity_mode = "fast",
                delete_corrupt_cds = TRUE,
                eval            = "1E-5",
                max.target.seqs = 10000,
                path            = NULL, 
                comp_cores      = 1,
                diamond_params  = NULL, 
                clean_folders   = FALSE,
                save.output     = NULL,
                best_hit_rule   = getOption("orthologr.best_hit_rule", "bitscore")){
        
        best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)
        
        # both files are needed for both search directions: read and translate each of them once
        end_seq_cache <- seq_cache_start()
        on.exit(end_seq_cache(), add = TRUE)
        
        orthoA <- diamond_best(
                query_file      = query_file, 
                subject_file    = subject_file, 
                seq_type        = seq_type,
                format          = format, 
                diamond_algorithm = diamond_algorithm,
                sensitivity_mode = sensitivity_mode,
                delete_corrupt_cds = delete_corrupt_cds,
                eval            = eval,
                max.target.seqs = max.target.seqs,
                path            = path, 
                comp_cores      = comp_cores,
                diamond_params  = diamond_params, 
                clean_folders   = clean_folders,
                save.output     = save.output,
                best_hit_rule   = best_hit_rule)
        
        # now the reverse search
        orthoB <- diamond_best(
                query_file      = subject_file, 
                subject_file    = query_file, 
                seq_type        = seq_type,
                format          = format, 
                diamond_algorithm = diamond_algorithm,
                sensitivity_mode = sensitivity_mode,
                delete_corrupt_cds = delete_corrupt_cds,
                eval            = eval,
                max.target.seqs = max.target.seqs,
                path            = path, 
                comp_cores      = comp_cores,
                diamond_params  = diamond_params, 
                clean_folders   = clean_folders,
                save.output     = save.output,
                best_hit_rule   = best_hit_rule)
        
        # rename the colnames for the reverse search
        colnames(orthoB)[1:2] <- c("subject_id", "query_id")
        
        tryCatch({
                return(dplyr::semi_join(
                        orthoA,
                        orthoB,
                        by = c("query_id", "subject_id")
                ))
                
                
        }, error = function(e) {
                stop(
                        "The DIAMOND2 tables resulting from ",
                        query_file,
                        " and ",
                        subject_file,
                        " could not be joined properly to select only the reciprocal best hits."
                )
        })
}