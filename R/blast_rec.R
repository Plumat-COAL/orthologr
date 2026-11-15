#' @title Perform a BLAST+ reciprocal best hit (RBH) search
#' @description This function performs a blast+ search (reciprocal best hit) of a given set of protein sequences against a second
#' set of protein sequences and vice versa.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, only "fasta" (Default) is supported.
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet. 
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the BLAST output shall be stored. E.g. \code{save.output} = \code{getwd()}
#'  to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @author Hajk-Georg Drost
#' @details Given a set of protein sequences A and a different set of protein sequences B,
#' first a best hit blast search is being performed from A to B: blast(A,B) and afterwards
#' a best hit blast search is being performed from B to A: blast(B,A). Only protein sequences
#' that were found to be best hits in both directions are retained and returned.
#'
#'
#' This function can be used to perform orthology inference using BLAST+ best reciprocal hit methodology.
#' 
#' @references
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schaeffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' @examples \dontrun{
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' blast_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'           
#'           
#'           
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' # starting with protein sequences
#' blast_rec(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type     = "protein")
#' 
#' 
#' 
#' # save the BLAST output file to the current working directory
#' blast_rec(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type     = "protein",
#'           save.output  = getwd())
#' 
#' 
#' 
#' # use multicore processing
#' blast_rec(query_file    = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            comp_cores   = 2)
#'            
#'            
#'            
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' # and external path to blastp
#' blast_rec(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           path         = "path/to/blastp")
#'           
#'           
#'           
#' }
#'
#' @return A data.table as returned by the \code{\link{blast}} function, storing the geneids
#' of orthologous genes (reciprocal best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{blast}}, \code{\link{blast_best}}, \code{\link{set_blast}}
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
blast_rec <- function(query_file, 
                      subject_file, 
                      seq_type        = "cds",
                      format          = "fasta", 
                      blast_algorithm = "blastp",
                      delete_corrupt_cds = TRUE,
                      eval            = "1E-5",
                      max.target.seqs = 10000,
                      path            = NULL, 
                      comp_cores      = 1, 
                      blast_params    = NULL, 
                      clean_folders   = FALSE,
                      save.output     = NULL,
                      best_hit_rule   = getOption("orthologr.best_hit_rule", "bitscore")){
        
        best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)
        
        # both files are needed for both search directions: read and translate each of them once
        end_seq_cache <- seq_cache_start()
        on.exit(end_seq_cache(), add = TRUE)
        
        orthoA <- blast_best( query_file      = query_file,
                              subject_file    = subject_file, 
                              eval            = eval,
                              max.target.seqs = max.target.seqs,
                              format          = format, 
                              seq_type        = seq_type,
                              blast_algorithm = blast_algorithm,
                              delete_corrupt_cds = delete_corrupt_cds,
                              path            = path, 
                              comp_cores      = comp_cores, 
                              blast_params    = blast_params,
                              save.output     = save.output,
                              best_hit_rule   = best_hit_rule )
        
        orthoB <- blast_best( query_file      = subject_file,
                              subject_file    = query_file,
                              seq_type        = seq_type,
                              eval            = eval,
                              max.target.seqs = max.target.seqs,
                              format          = format,
                              blast_algorithm = blast_algorithm,
                              delete_corrupt_cds = delete_corrupt_cds,
                              path            = path, 
                              comp_cores      = comp_cores, 
                              blast_params    = blast_params,
                              clean_folders   = clean_folders,
                              save.output     = save.output,
                              best_hit_rule   = best_hit_rule )
        
        colnames(orthoB)[1:2] <- c("subject_id", "query_id")
        
        tryCatch({
                return(dplyr::semi_join(
                        orthoA,
                        orthoB,
                        by = c("query_id", "subject_id")
                ))
                
                
        }, error = function(e) {
                stop(
                        "The BLAST tables resulting from ",
                        query_file,
                        " and ",
                        subject_file,
                        " could not be joined properly to select only the reciprocal best hits."
                )
        })
}

