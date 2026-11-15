#' @title Perform a DIAMOND2 best hit search
#' @description This function performs a DIAMOND2 search (best hit) of a given set of protein sequences against a given database.
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
#' @details Given a set of protein sequences (query sequences), a best hit diamond search (BH DIAMOND2) is being performed.
#' 
#' Internally to perform best hit searches, the DIAMOND2 parameter settings:
#' 
#' \code{"-best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1"}
#' 
#' are used to speed up best hit computations.
#' 
#' This function gives the same output as \code{\link{blast_best}} while being much much faster.
#' 
#' @references
#' 
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
#'
#' @examples \dontrun{
#' 
#' # performing gene orthology inference using the best hit (BH) method
#' diamond_best(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'            
#'            
#'            
#' # performing gene orthology inference using the best hit (BH) method starting with protein sequences
#' diamond_best(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'            seq_type     = "protein")
#' 
#' 
#' 
#' # save the DIAMOND2 output file to the current working directory
#' diamond_best(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'            seq_type     = "protein",
#'            save.output  = getwd())
#' 
#' 
#' 
#' 
#' 
#' # use multicore processing
#' diamond_best(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            comp_cores   = 2)
#'
#'
#'
#' # performing gene orthology inference using the best hit (BH) method and external
#' # diamond path
#' diamond_best(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            path         = "path/to/diamond/")
#'            
#'            
#' }
#'
#' @return A data.table as returned by the \code{diamond} function, storing the geneids
#' of orthologous genes (best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{diamond}}, \code{\link{diamond_rec}}, \code{\link{set_diamond}}, \code{\link{blast_best}}
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
diamond_best <- function(
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
        
        # default parameters for best hit filtering for blast
        # this is not needed for DIAMOND
        # default_pars <- "-best_hit_score_edge 0.05 -best_hit_overhang 0.25"
        
        
        # performing a DIAMOND search from query against subject: diamond(query,subject)
        
        hit_tbl <- diamond( 
                query_file      = query_file,
                subject_file    = subject_file,
                diamond_algorithm = diamond_algorithm,
                eval            = eval,
                max.target.seqs = max.target.seqs,
                delete_corrupt_cds = delete_corrupt_cds,
                seq_type        = seq_type,
                sensitivity_mode = sensitivity_mode,
                format          = format, 
                path            = path, 
                comp_cores      = comp_cores,
                diamond_params  = diamond_params,
                clean_folders   = clean_folders,
                save.output     = save.output)
        
        query_id <- '.' <- NULL
        hit_tbl <- best_hits(hit_tbl, "query_id", rule = best_hit_rule)
        return(hit_tbl)
}

