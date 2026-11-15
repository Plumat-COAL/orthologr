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
#' amino acids are exchanged). The argument is ignored for DNA sequences. \strong{The scoring of 0.4.3 is slow:}
#' \pkg{pwalign} aligns without a matrix about 15 times slower than with one, about 0.6 - 0.8 s per gene pair of
#' typical length (1,000-gene example, 1 core: 632 s against 45 s for the whole \code{dNdS()} run); a genome-wide
#' run with \code{aa_sub_matrix = NULL} or \code{options(orthologr.legacy_pairwise_aln = TRUE)} (also set by
#' \code{\link{orthologr_legacy_options}}) takes hours per core.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details This function provides an interface between R and common pairwise alignment computation methods.
#' 
#' The current version of this function computes pairwise alignments based on the \code{\link[pwalign]{pairwiseAlignment}}
#' function implemented in the \pkg{pwalign} package.
#' 
#' The default pairwise alignment method is based on the \emph{Needleman-Wunsch Algorithm}.
#' 
#' Since \pkg{orthologr} 0.5.0 (see \code{NEWS.md}):
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
#' \strong{Memory.} \code{\link[pwalign]{pairwiseAlignment}} keeps the whole dynamic-programming matrix: about 3 bytes
#' per cell, i.e. per pair of residues (two proteins of 10,000 residues: 0.3 GB; a titin pair of 36,000 x 35,400
#' residues: 3.9 GB). An alignment of more than \code{options(orthologr.max_aln_cells)} cells (default 5e8, about
#' 1.5 GB) is not computed by \pkg{pwalign}. With the R package \pkg{coalr} installed (and the default
#' \code{aa_sub_matrix = "BLOSUM62"}) such a pair is aligned by COAL's aligner with its \code{"pwalign"} preset instead
#' (a message names the pair), which returns the alignment \pkg{pwalign} would return - the same scoring and
#' tie-break; titin: identical rows and score - in a fraction of the memory (titin, 35,991 x 35,390 residues: 424 MB
#' for the R process and 49 s, against 3.9 GB and 18 s with \pkg{pwalign}). Without
#' \pkg{coalr} (or with \code{options(orthologr.large_pair_aligner = "none")}, or with \code{get_aln = TRUE})
#' \code{pairwise_aln()} stops with an error of class \code{orthologr_error_alignment_too_large}, and \code{\link{dNdS}}
#' returns the pair with status \code{"failed"}. \code{options(orthologr.max_aln_cells = Inf)} removes the limit:
#' \pkg{pwalign} aligns every pair (see \code{\link{orthologr_options}}).
#' 
#' @examples
#' # Needleman-Wunsch alignment of the two proteins of a file
#' # (the alignment file is written to tempdir())
#' pairwise_aln( file     = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'               tool     = "NW",
#'               get_aln  = TRUE,
#'               seq_type = "protein",
#'               quiet    = TRUE)
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
                orthologr_abort(class = "bad_argument", "Please choose a tool that is supported by this function.")
        
        if (!is.element(seq_type, c("dna", "cds", "protein")))
                orthologr_abort(class = "bad_argument", "Please choose a seq_type that is supported by this function.")
        
        if (is.element(seq_type, c("dna", "cds")))
                seqtype <- "DNA"
        
        if (seq_type == "protein")
                seqtype <- "AA"
        
        # options(orthologr.legacy_pairwise_aln = TRUE): the behaviour of orthologr <= 0.4.3
        legacy <- isTRUE(orthologr_option("orthologr.legacy_pairwise_aln", FALSE))
        
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
                
                # The dynamic programming of pwalign::pairwiseAlignment() needs about 3 bytes per cell of the
                # n x m matrix (measured: 1e8 cells 0.31 GB; a titin pair of 36,000 x 35,400 residues 3.9 GB), with
                # comp_cores > 1 once per running pair. Above options(orthologr.max_aln_cells) the pair is not
                # aligned (an error, which dNdS() reports as a failed pair) - unless the R package coalr is installed:
                # then COAL's aligner with its pwalign preset aligns it, which returns pwalign's alignment (same
                # scoring and tie-break; titin: identical rows and score) in a ninth of the memory.
                if (seq_type == "protein" && !legacy && !get_aln && aln_too_large(length(seqs[[1]]), length(seqs[[2]]))) {
                        rows <- large_pair_coal_rows(seqs, names, aa_sub_matrix)
                        if (!is.null(rows)) {
                                seqinr::write.fasta(sequences = as.list(rows), names = names, file.out = file.out)
                                if (!quiet)
                                        message("File successfully written to ", file.out, ".")
                                return(invisible(NULL))
                        }
                }
                check_aln_cells(length(seqs[[1]]), length(seqs[[2]]))
                
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
                # The full-length rows (overhangs aligned against gaps) are what
                # pwalign::alignedPattern(aln) / pwalign::alignedSubject(aln) return; aligned_rows() builds the
                # same two strings from the aligned ranges, which takes about a hundredth of the time.
                if (legacy) {
                        aln_seqs <- list(pwalign::pattern(aln), pwalign::subject(aln))
                } else {
                        aln_seqs <- as.list(aligned_rows(aln, as.character(seqs[[1]]), as.character(seqs[[2]])))
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
                        message("File successfully written to ", file.out, ".")
                }
                
                if (get_aln)
                        return(aln)
        }
}






# The two full-length rows of a global pairwise alignment as character strings: what
# as.character(pwalign::alignedPattern(aln)) and as.character(pwalign::alignedSubject(aln)) return, built from the
# aligned ranges (pwalign::pattern(aln), pwalign::subject(aln)) and the ends of the two sequences that lie outside
# them. An unaligned end of the pattern stands against gaps in the subject row and the other way round. (Only
# for type = "global": for local and overlap alignments pwalign's functions leave the unaligned ends out.)
#   aln         : a PairwiseAlignments object holding one alignment
#   pattern_seq : the complete (unaligned) pattern sequence, one string
#   subject_seq : the complete (unaligned) subject sequence, one string
aligned_rows <- function(aln, pattern_seq, subject_seq) {
        rows <- aligned_rows_list(aln, pattern_seq, subject_seq)
        c(rows[[1]], rows[[2]])
}

# aligned_rows() for a PairwiseAlignments object of any length (alignment i of pattern_seq[i] and subject_seq[i]):
# list(<pattern rows>, <subject rows>), two character vectors of length(aln).
aligned_rows_list <- function(aln, pattern_seq, subject_seq) {
        p_range <- pwalign::pattern(aln)
        s_range <- pwalign::subject(aln)
        p_start <- pwalign::start(p_range)
        p_end   <- pwalign::end(p_range)
        s_start <- pwalign::start(s_range)
        s_end   <- pwalign::end(s_range)
        n_p <- nchar(pattern_seq)
        n_s <- nchar(subject_seq)
        gaps <- function(n) strrep("-", pmax(n, 0))
        list(paste0(gaps(s_start - 1), substr(pattern_seq, 1, p_start - 1), as.character(p_range),
                    substr(pattern_seq, p_end + 1, n_p), gaps(n_s - s_end)),
             paste0(substr(subject_seq, 1, s_start - 1), gaps(p_start - 1), as.character(s_range),
                    gaps(n_p - p_end), substr(subject_seq, s_end + 1, n_s)))
}

# The full-length alignment rows of many protein pairs (proteins aa1[i] and aa2[i]), with one call of
# pwalign::pairwiseAlignment() for all of them: the rows pairwise_aln() writes for each pair (BLOSUM62 or the
# matrix given, stop codon symbol at the end removed, overhangs against gaps), without building an AAString object
# and reading and writing a file per pair (about a third of the run time of dNdS() went into that). pwalign aligns
# the i-th pattern with the i-th subject, each pair on its own (the same alignment and tie-break as a call for the
# pair alone; checked on every pair of the 20- and 1,000-gene examples and of S. cerevisiae / C. glabrata).
# Returns a list with one element per pair: c(<row 1>, <row 2>), or NULL where the pair is left to pairwise_aln()
# itself: no substitution matrix (the scoring of 0.4.3), a sequence that is empty, NA or holds a character other
# than A - Z and '*', a pair above options(orthologr.max_aln_cells) (aligned by COAL or refused there), or when
# anything in the joint call fails (the pairs are then aligned one by one, and a pair that fails does so with its
# own error message).
pairwise_aln_rows <- function(aa1, aa2, aa_sub_matrix = "BLOSUM62") {
        n <- length(aa1)
        out <- vector("list", n)
        if (n == 0 || is.null(aa_sub_matrix))
                return(out)
        limit <- tryCatch({
                check_aa_sub_matrix(aa_sub_matrix)
                max_aln_cells()
        }, error = function(e) NULL)
        if (is.null(limit))
                return(out)
        strip_stop <- function(x) {
                k <- nchar(x)
                cut <- !is.na(x) & k > 1 & substr(x, k, k) == "*"
                x[cut] <- substr(x[cut], 1, k[cut] - 1)
                x
        }
        a <- strip_stop(as.character(aa1))
        b <- strip_stop(as.character(aa2))
        use <- !is.na(a) & !is.na(b) & grepl("^[A-Z*]+\\z", a, perl = TRUE, useBytes = TRUE) &
                grepl("^[A-Z*]+\\z", b, perl = TRUE, useBytes = TRUE)
        use[use] <- !(is.finite(limit) & as.numeric(nchar(a[use])) * as.numeric(nchar(b[use])) > limit)
        idx <- which(use)
        if (length(idx) == 0)
                return(out)
        rows <- tryCatch({
                aln <- pwalign::pairwiseAlignment(pattern = Biostrings::AAStringSet(a[idx]),
                                                  subject = Biostrings::AAStringSet(b[idx]),
                                                  type    = "global",
                                                  substitutionMatrix = aa_sub_matrix)
                aligned_rows_list(aln, a[idx], b[idx])
        }, error = function(e) NULL)
        if (is.null(rows) || length(rows[[1]]) != length(idx))
                return(out)
        out[idx] <- mapply(c, rows[[1]], rows[[2]], SIMPLIFY = FALSE, USE.NAMES = FALSE)
        out
}

# The largest pairwise alignment (n x m cells) that pairwise_aln() computes: options(orthologr.max_aln_cells), by
# default 5e8 cells (about 1.5 GB of memory in pwalign::pairwiseAlignment(), e.g. two proteins of 22,000 residues).
# Inf (the value of 0.4.3, also under the legacy switch) has no limit.
max_aln_cells <- function() {
        v <- orthologr_option("orthologr.max_aln_cells", 5e8)
        if (!is.numeric(v) || length(v) != 1 || is.na(v) || v <= 0)
                orthologr_abort(class = "bad_argument", "options(orthologr.max_aln_cells) must be one positive number (Inf: no limit), not ",
                                if (is.character(v)) paste0("'", paste(v, collapse = "', '"), "' (text)") else paste(format(v), collapse = ", "),
                                ".", call. = FALSE)
        v
}

# Has the alignment of two sequences of n and m residues more cells than max_aln_cells() allows?
aln_too_large <- function(n, m) {
        limit <- max_aln_cells()
        is.finite(limit) && as.numeric(n) * as.numeric(m) > limit
}

# A pair above options(orthologr.max_aln_cells) aligned with COAL (decision D3 of 2026-11-14): the two full-length
# rows of the alignment coalr::coal_align_pairs(preset = "pwalign") returns - pwalign's own scoring (BLOSUM62 as
# pwalign has it, gap opening 10 + extension 4) and tie-break, so the rows are those pwalign::pairwiseAlignment()
# would return (titin, 35,991 x 35,390 residues: identical rows and score 166,805; through coalr 49 s and 424 MB for
# the whole R process, pwalign 18 s and 3.9 GB: slower, in a ninth of the memory). NULL when that is not possible - coalr is not installed (or options(orthologr.coal_backend) is
# "none" / "binary"), the matrix is not BLOSUM62 (the only one checked to be pwalign's), COAL fails, or
# options(orthologr.large_pair_aligner = "none") - and then the pair is not aligned (the error of check_aln_cells()).
large_pair_coal_rows <- function(seqs, names, aa_sub_matrix) {
        aligner <- orthologr_option("orthologr.large_pair_aligner", "coal")
        if (!is.character(aligner) || length(aligner) != 1 || !is.element(aligner, c("coal", "none")))
                orthologr_abort(class = "bad_argument", "options(orthologr.large_pair_aligner = ) must be \"coal\" or \"none\".",
                                call. = FALSE)
        backend <- getOption("orthologr.coal_backend", NULL)
        if (aligner == "none" || isTRUE(backend %in% c("none", "binary")) || !identical(aa_sub_matrix, "BLOSUM62") ||
            !requireNamespace("coalr", quietly = TRUE))
                return(NULL)
        matrix <- coal_matrix_arg(aa_sub_matrix, "pwalign")
        if (!identical(matrix, "BLOSUM62-PWALIGN"))
                return(NULL)
        a <- as.character(seqs[[1]])
        b <- as.character(seqs[[2]])
        cells <- as.numeric(nchar(a)) * as.numeric(nchar(b))
        message("The pair '", names[1], "' / '", names[2], "' (", nchar(a), " x ", nchar(b), " residues, ",
                format(cells, big.mark = ",", scientific = FALSE), " cells) is above options(orthologr.max_aln_cells) = ",
                format(max_aln_cells(), big.mark = ",", scientific = FALSE), ": aligned with COAL (coalr, preset ",
                "\"pwalign\"), which returns pwalign's alignment in a fraction of the memory.")
        r <- tryCatch(coalr::coal_align_pairs(stats::setNames(a, names[1]), stats::setNames(b, names[2]),
                                              preset = "pwalign", matrix = matrix, threads = 1L),
                      error = function(e) NULL)
        if (is.null(r) || nrow(r) != 1 || !is.na(r$error[1]) || is.na(r$aligned_a[1]))
                return(NULL)
        c(r$aligned_a[1], r$aligned_b[1])
}

# Stop when the alignment of two sequences of n and m residues has more cells than max_aln_cells() allows.
check_aln_cells <- function(n, m) {
        limit <- max_aln_cells()
        cells <- as.numeric(n) * as.numeric(m)
        if (is.finite(limit) && cells > limit)
                orthologr_abort(class = "alignment_too_large",
                     "the pairwise alignment of two sequences of ", n, " and ", m, " residues has ",
                     format(cells, big.mark = ",", scientific = FALSE), " cells and would need about ",
                     signif(3 * cells / 1e9, 2), " GB of memory; options(orthologr.max_aln_cells) is ",
                     format(limit, big.mark = ",", scientific = FALSE), ". The pair is not aligned; ",
                     "options(orthologr.max_aln_cells = Inf) aligns it anyway (about 3 bytes per cell, with comp_cores > 1 for every pair aligned at the same time); ",
                     "with the R package coalr installed such pairs are aligned by COAL with pwalign's scoring (the same alignment ",
                     "in a fraction of the memory: titin, 36,000 x 35,000 residues, 0.4 GB), and dNdS(engine = \"coal\") ",
                     "aligns every pair with COAL.",
                     call. = FALSE)
        invisible(cells)
}
