#' @title Main Orthology Inference Function
#' @description This function takes nucleotide or protein sequences for a set of organisms
#' and performs orthology inference to detect orthologous genes within the given organisms
#' based on selected orthology inference programs.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_files a character string specifying the paths to the sequence files of interest (subject organisms).
#' @param outgroup_file deprecated and not used (a warning is given when it is set).
#' @param eval a numeric value specifying the E-Value cutoff for hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are: "cds", "protein", or "dna". In case of "cds", sequences are translated to protein sequences,
#' in case of "dna", CDS prediction is performed on the corresponding sequences which are subsequently
#' translated to protein sequences. Default is \code{seq_type} = "protein".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta". Default is "fasta".
#' @param ortho_detection a character string specifying the orthology inference method to use.
#' Default is \code{ortho_detection} = \code{"DIAMOND_RBH"} (DIAMOND2 reciprocal best hit).
#' Available methods:
#' \itemize{
#'   \item \code{"DIAMOND_RBH"}: DIAMOND2 reciprocal best hit (default; fast, recommended)
#'   \item \code{"DIAMOND_BH"}: DIAMOND2 best hit
#'   \item \code{"RBH"}: BLAST reciprocal best hit
#'   \item \code{"BH"}: BLAST best hit
#'   \item \code{"Orthofinder2"}: not implemented (an error). For orthogroups of two or more species use
#'   \code{\link{orthogroups}}, which needs DIAMOND only.
#' }
#' \code{orthologs()} compares one query with one subject; \code{\link{orthogroups}} infers groups of orthologs for
#' any number of species (reciprocal best hits of all pairs of species, or DIAMOND DeepClust clusters).
#' @param sensitivity_mode a character string specifying the DIAMOND2 sensitivity level.
#' Only used when \code{ortho_detection} is \code{"DIAMOND_RBH"} or \code{"DIAMOND_BH"}.
#' Options: \code{"faster"}, \code{"fast"} (default), \code{"mid-sensitive"},
#' \code{"sensitive"}, \code{"more-sensitive"}, \code{"very-sensitive"}, \code{"ultra-sensitive"}.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep (DIAMOND2 only).
#' Default is \code{10000}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should
#' be removed from the input \code{file}. Default is \code{FALSE}.
#' @param cdd.path deprecated and not used (a warning is given when it is set).
#' @param path a character string specifying the path to the orthology inference tool executable.
#' @param add_params a character string specifying additional parameters passed to the orthology inference tool.
#' For DIAMOND2 methods, use \code{diamond_params} instead. Default is \code{NULL}.
#' @param diamond_params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string specifying additional parameters passed to DIAMOND2.
#' Only used when \code{ortho_detection} is \code{"DIAMOND_RBH"} or \code{"DIAMOND_BH"}.
#' Default is \code{NULL}.
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore computations.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @param clean_folders a boolean value specifying whether all internal folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @details
#' This function takes sequence files of a query organism and a subject organism and performs orthology inference
#' using a defined orthology inference method to detect orthologous genes.
#'
#' The following interfaces are implemented in the \code{orthologs} function:
#'
#' DIAMOND2 based methods (recommended):
#' \itemize{
#'   \item DIAMOND2 reciprocal best hit (\code{"DIAMOND_RBH"}) — default
#'   \item DIAMOND2 best hit (\code{"DIAMOND_BH"})
#' }
#'
#' BLAST based methods:
#' \itemize{
#'   \item BLAST reciprocal best hit (\code{"RBH"})
#'   \item BLAST best hit (\code{"BH"})
#' }
#'
#' @author Hajk-Georg Drost
#' @return A data.table storing the query_ids of orthologous genes in the first column, the subject_ids of
#' orthologous genes in the second column and the amino acid sequences in the third column.
#' @examples \dontrun{
#'
#' ### DIAMOND2 Reciprocal Best Hit (default)
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein")
#'
#'
#' ### DIAMOND2 Best Hit
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein",
#'           ortho_detection = "DIAMOND_BH")
#'
#'
#' ### BLAST Reciprocal Best Hit
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein",
#'           ortho_detection = "RBH")
#'
#'
#' ### BLAST Best Hit
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein",
#'           ortho_detection = "BH")
#' }
#' @seealso \code{\link{orthogroups}} (two or more species), \code{\link{diamond_rec}}, \code{\link{diamond_best}}, \code{\link{blast_rec}}, \code{\link{dNdS}}
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

orthologs <- function(query_file,
                      subject_files,
                      seq_type           = "protein",
                      outgroup_file      = NULL,
                      eval               = "1E-5",
                      format             = "fasta",
                      ortho_detection    = "DIAMOND_RBH",
                      sensitivity_mode   = "fast",
                      max.target.seqs    = 10000,
                      delete_corrupt_cds = FALSE,
                      cdd.path           = NULL,
                      path               = NULL,
                      add_params         = NULL,
                      diamond_params     = NULL,
                      comp_cores         = 1,
                      quiet              = FALSE,
                      clean_folders      = FALSE,
                      best_hit_rule      = getOption("orthologr.best_hit_rule", "bitscore")) {
        
        best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)

        if (!is.element(ortho_detection,
                        c("DIAMOND_RBH", "DIAMOND_BH", "BH", "RBH", "Orthofinder2")))
                stop(
                        paste0(
                                "Unsupported orthology detection method: '", ortho_detection, "'. ",
                                "Choose one of: 'DIAMOND_RBH', 'DIAMOND_BH', 'RBH', 'BH', 'Orthofinder2'."
                        ),
                        call. = FALSE
                )

        for (unused in c("outgroup_file", "cdd.path")) {
                if (!is.null(get(unused)))
                        warning("The argument '", unused, "' is deprecated and not used.", call. = FALSE)
        }

        # 'Orthofinder2' was announced but never implemented (the call returned NULL with a message)
        if (ortho_detection == "Orthofinder2")
                stop("ortho_detection = 'Orthofinder2' is not implemented in orthologs(). ",
                     "Please choose one of: 'DIAMOND_RBH', 'DIAMOND_BH', 'RBH', 'BH'. For orthogroups of two or more ",
                     "species use orthogroups() (DIAMOND only), or run OrthoFinder yourself ",
                     "(see ?orthofinder2_retrieve_core_orthologs to import its single-copy orthogroups).", call. = FALSE)

        i <- NULL

        # (For seq_type = "cds" the coding sequences used to be written to '<file>_translated.fasta' first - as
        # DNA, not translated - and these copies were searched. The searches read and translate the input files
        # themselves; the copies are gone.)
        
        # several transcripts per gene in an input file? (a warning; see ?dNdS, section 'Splice variants')
        if (is.element(seq_type, c("cds", "protein"))) {
                isoform_warning(query_file, "query_file")
                for (f in subject_files)
                        isoform_warning(f, "subject_files")
        }

        # ---- DIAMOND2 reciprocal best hit (default) ----
        if (ortho_detection == "DIAMOND_RBH") {

                if (length(subject_files) > 1)
                        stop("The DIAMOND2 reciprocal best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        diamond_rec(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                sensitivity_mode   = sensitivity_mode,
                                max.target.seqs    = max.target.seqs,
                                diamond_params     = diamond_params,
                                seq_type           = seq_type,
                                format             = format,
                                clean_folders      = clean_folders,
                                best_hit_rule      = best_hit_rule
                        )
                )
        }

        # ---- DIAMOND2 best hit ----
        if (ortho_detection == "DIAMOND_BH") {

                if (length(subject_files) > 1)
                        stop("The DIAMOND2 best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        diamond_best(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                sensitivity_mode   = sensitivity_mode,
                                max.target.seqs    = max.target.seqs,
                                diamond_params     = diamond_params,
                                seq_type           = seq_type,
                                format             = format,
                                clean_folders      = clean_folders,
                                best_hit_rule      = best_hit_rule
                        )
                )
        }

        # ---- BLAST best hit ----
        if (ortho_detection == "BH") {

                if (length(subject_files) > 1)
                        stop("The BLAST best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        blast_best(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                blast_params       = add_params,
                                seq_type           = seq_type,
                                format             = format,
                                best_hit_rule      = best_hit_rule
                        )
                )

                if (clean_folders)
                        clean_all_folders(file.path(orthologr_tmpdir(), "_blast_db"))
        }

        # ---- BLAST reciprocal best hit ----
        if (ortho_detection == "RBH") {

                if (length(subject_files) > 1)
                        stop("The BLAST reciprocal best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        blast_rec(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                blast_params       = add_params,
                                seq_type           = seq_type,
                                format             = format,
                                best_hit_rule      = best_hit_rule
                        )
                )

                if (clean_folders)
                        clean_all_folders(file.path(orthologr_tmpdir(), "_blast_db"))
        }

        # ---- OrthoFinder2 (under development) ----
        return(ortho_tbl)
}
