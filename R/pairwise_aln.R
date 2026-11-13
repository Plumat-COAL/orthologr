#' @title Compute Pairwise Alignments
#' @description This function takes a FASTA file containing two DNA or amino acid sequences
#' that shall be aligned and computes a paiwise alignment using a defined alignment method.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program/algorithm that should be used: \code{tool} = \code{"NW"}.
#' @param seq_type a character string specifying the sequence type stored within the given FASTA file. Options are
#' \code{seq_type} = \code{"protein"}, \code{seq_type} = \code{"cds"}, \code{seq_type} = \code{"dna"}.
#'  Default is \code{seq_type} = \code{"protein"}. 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param pairwise_aln_name a character string specifying the name of the stored alignment file. 
#' Default is \code{pairwise_aln_name} = \code{NULL} denoting a default name: 'toolname_seq_type.aln' .
#' @param path a character string specifying the path to the pairwise alignment program (in case you don't use the default path).
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' If \code{store_locally = TRUE}, the alignment is also written with its score to
#' \code{orthologr_alignment_files/_pairwise_alignment_with_score} (\code{\link[pwalign]{writePairwiseAlignments}}).
#' This score file is no longer written when \code{store_locally = FALSE}.
#' @param aa_sub_matrix the substitution matrix used to score amino acid alignments (\code{seq_type = "protein"}):
#' the name of a matrix shipped with \pkg{pwalign} (\code{"BLOSUM62"} (Default), \code{"BLOSUM45"}, \code{"BLOSUM50"},
#' \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"}, \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, \code{"PAM250"}),
#' or a numeric substitution matrix. \code{aa_sub_matrix = NULL} gives the scoring of \pkg{orthologr} <= 0.4.3
#' (\code{\link[pwalign]{pairwiseAlignment}} without a substitution matrix, i.e. a match/mismatch score that ignores which
#' amino acids are exchanged). The argument is ignored for DNA sequences.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details This function provides an interface between R and common pairwise alignment computation methods.
#' 
#' The current version of this function computes pairwise alignments based on the \code{\link[pwalign]{pairwiseAlignment}}
#' function implemented in the \pkg{pwalign} package.
#' 
#' The default pairwise alignment method is based on the \emph{Needleman-Wunsch Algorithm}.
#' 
#' Since \pkg{orthologr} 0.4.4 (see \code{NEWS.md}):
#' \itemize{
#' \item amino acid sequences are scored with \code{BLOSUM62} by default (argument \code{aa_sub_matrix});
#' \item the alignment that is written is full length: every residue of both sequences is in it, and the overhanging
#' ends of the longer sequence are aligned against gaps. Earlier versions wrote only the aligned range and thereby cut
#' the overhangs off;
#' \item a stop codon symbol (\code{*}) at the end of an amino acid sequence is removed before aligning, so that the
#' stop codon is not part of the alignment (and of the codon alignment and dN/dS estimate derived from it).
#' }
#' All three change the alignments and therefore, slightly, the dN/dS values computed from them.
#' The behaviour of \pkg{orthologr} <= 0.4.3 (old scoring, truncated overhangs, stop codon aligned, score file always
#' written) can be restored for reproducing earlier results with \code{options(orthologr.legacy_pairwise_aln = TRUE)}.
#' 
#' @examples \dontrun{        
#'                                             
#' # Needleman-Wunsch Example:  
#' 
#' # in case Biostrings works properly
#' pairwise_aln( file     = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'               tool     = "NW", 
#'               get_aln  = TRUE, 
#'               seq_type = "protein")
#' 
#'                                                    
#'                                                                                                                                                          
#' }
#' @return In case the argument \code{get_aln} is set \code{TRUE}, an object of class alignment of the seqinr package is returned.
#' @export
pairwise_aln <- function(file, 
                         tool              = "NW", 
                         seq_type          = "protein", 
                         get_aln           = FALSE, 
                         pairwise_aln_name = NULL,
                         path              = NULL, 
                         quiet             = FALSE,
                         store_locally     = FALSE,
                         aa_sub_matrix     = "BLOSUM62"){
        
        
        if (!is.pairwise_aln_tool(tool))
                stop("Please choose a tool that is supported by this function.")
        
        if (!is.element(seq_type, c("dna", "cds", "protein")))
                stop("Please choose a seq_type that is supported by this function.")
        
        if (is.element(seq_type, c("dna", "cds")))
                seqtype <- "DNA"
        
        if (seq_type == "protein")
                seqtype <- "AA"
        
        # options(orthologr.legacy_pairwise_aln = TRUE): the behaviour of orthologr <= 0.4.3
        legacy <- isTRUE(getOption("orthologr.legacy_pairwise_aln", FALSE))
        
        if (legacy)
                aa_sub_matrix <- NULL
        
        check_aa_sub_matrix(aa_sub_matrix)
        
        if (!file.exists(file.path(ifelse(store_locally, "orthologr_alignment_files", orthologr_tmpdir()), "_alignment", "pairwise_aln"))) {
                dir.create(file.path(ifelse(store_locally, "orthologr_alignment_files", orthologr_tmpdir()), "_alignment", "pairwise_aln"),
                           recursive = TRUE)
        }
        
        if (is.null(pairwise_aln_name)) {
                file.out <-
                        file.path(ifelse(store_locally, "orthologr_alignment_files", orthologr_tmpdir()),
                                  "_alignment",
                                  "pairwise_aln",
                                  paste0(tool, "_", seqtype, ".aln"))
        }
        
        if (!is.null(pairwise_aln_name)) {
                file.out <-
                        file.path(
                                ifelse(store_locally, "orthologr_alignment_files", orthologr_tmpdir()),
                                "_alignment",
                                "pairwise_aln",
                                paste0(pairwise_aln_name, "_", tool, "_", seqtype, ".aln")
                        )
        }
        
        # Needleman-Wunsch using pwalign::pairwiseAlignment( type=global)
        if (tool == "NW") {
                #read file
                if (is.element(seq_type, c("dna", "cds"))) {
                        #                         not possible as is sorts the input by name
                        #                         seqs <- read.cds(file=file, format=format)
                        #                         names <- dt[,geneids]
                        #                         seqs <- sapply(dt[,seqs] , Biostrings::DNAString)
                        
                        input <- seqinr::read.fasta(file    = file,
                                                    seqtype = seqtype)
                        
                        names <- names(input)
                        seqs <-
                                lapply(input, function(x) {
                                        return(Biostrings::DNAString(seqinr::c2s(x)))
                                })
                }
                
                if (seq_type == "protein") {
                        #                         dt <- read.proteome(file=file, format=format)
                        #                         names <- dt[,geneids]
                        #                         seqs <- sapply(dt[,seqs] , Biostrings::AAString)
                        
                        input <-
                                seqinr::read.fasta(file    = file,
                                                   seqtype = seqtype)
                        
                        names <- names(input)
                        seqs <-
                                lapply(input, function(x) {
                                        return(Biostrings::AAString(seqinr::c2s(x)))
                                })
                }
                
                if (length(seqs) < 2)
                        stop("The file '", file, "' must store two sequences to compute a pairwise alignment.", call. = FALSE)
                
                # a stop codon symbol at the end of a protein is not aligned (it would otherwise become a column
                # of the codon alignment and be scored in the dN/dS estimate)
                if (seq_type == "protein" && !legacy) {
                        seqs <- lapply(seqs, function(x) {
                                n <- length(x)
                                if (n > 1 && as.character(Biostrings::subseq(x, n, n)) == "*")
                                        x <- Biostrings::subseq(x, 1, n - 1)
                                return(x)
                        })
                }
                
                # align
                if (seq_type == "protein" && !is.null(aa_sub_matrix)) {
                        aln  <-
                                pwalign::pairwiseAlignment(pattern = seqs[[1]],
                                                           subject = seqs[[2]],
                                                           type    = "global",
                                                           substitutionMatrix = aa_sub_matrix)
                } else {
                        aln  <-
                                pwalign::pairwiseAlignment(pattern = seqs[[1]],
                                                           subject = seqs[[2]],
                                                           type    = "global")
                }
                
                # write file
                # pwalign::pattern(aln) / pwalign::subject(aln) only hold the aligned range of each sequence:
                # for a global alignment the overhanging ends of the longer sequence are cut off.
                # pwalign::alignedPattern(aln) / pwalign::alignedSubject(aln) are the full-length rows
                # (overhangs aligned against gaps).
                if (legacy) {
                        aln_seqs <- list(pwalign::pattern(aln), pwalign::subject(aln))
                } else {
                        aln_seqs <- list(as.character(pwalign::alignedPattern(aln)),
                                         as.character(pwalign::alignedSubject(aln)))
                }
                
                seqinr::write.fasta(
                        sequences = aln_seqs,
                        names     = names,
                        file.out  = file.out
                )
                
                # the alignment with its score: only when alignment files are stored locally
                if (store_locally) {
                        local_store_path <- file.path("orthologr_alignment_files", "_pairwise_alignment_with_score")
                        if (!file.exists(local_store_path))
                                dir.create(local_store_path, recursive = TRUE)
                        pwalign::writePairwiseAlignments(aln, file.path(local_store_path, paste0(pairwise_aln_name, "_", tool, "_", seqtype, ".aln")))
                } else if (legacy) {
                        local_store_path <- file.path(orthologr_tmpdir(), "_pairwise_alignment_with_score")
                        if (!file.exists(local_store_path))
                                dir.create(local_store_path, recursive = TRUE)
                        pwalign::writePairwiseAlignments(aln, file.path(local_store_path, paste0(pairwise_aln_name, "_", tool, "_", seqtype, ".aln")))
                }
                
                if (!quiet) {
                        print(paste0("File successfully written to ", file.out, "."))
                }
                
                if (get_aln)
                        return(aln)
        }
}

# 'aa_sub_matrix' must be the name of a substitution matrix, a numeric matrix, or NULL.
check_aa_sub_matrix <- function(aa_sub_matrix) {
        if (!is.null(aa_sub_matrix) &&
            !((is.character(aa_sub_matrix) && length(aa_sub_matrix) == 1 && !is.na(aa_sub_matrix)) ||
              (is.matrix(aa_sub_matrix) && is.numeric(aa_sub_matrix))))
                stop("'aa_sub_matrix' must be the name of a substitution matrix (e.g. \"BLOSUM62\"), a numeric substitution matrix, or NULL.", call. = FALSE)
        invisible(TRUE)
}
