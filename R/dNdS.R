#' @title Compute dNdS values for two organisms
#' @description This function takes the CDS files of two organisms of interest (query_file and subject_file)
#' and computes the dNdS estimation values for orthologous gene pairs between these organisms. 
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param aligner a character string specifying the sequence aligner. Options are: 
#' \itemize{
#' \item \code{aligner = "diamond"} (Default): sequence alignment is performed using DIAMOND2
#' \item \code{aligner = "blast"}: sequence alignment is performed using BLASTP.
#' }
#' The option diamond uses DIAMOND2 which is up to 10 000X folds faster than BLASTP while retaining the sensitivity of BLASTP.
#' @param seq_type the sequence type stored in the input files. Only \code{seq_type = "cds"} (Default) is possible:
#' dN/dS is computed from codon alignments, so both files must hold coding sequences. Any other value is an error.
#' @param format the file format of the sequence files. Only \code{format = "fasta"} (Default) is supported.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed to detect orthologous genes. Options are:
#' \itemize{
#' \item \code{ortho_detection = "RBH"}: DIAMOND or BLAST reciprocal/bidirectional best hit (Default).
#' \item \code{ortho_detection = "BH"}: DIAMOND or BLAST unidirectional best hit.
#' }
#' Any other value is an error (\code{"Orthofinder2"}, \code{"PO"} and \code{"DELTA"} are not available).
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet. 
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' @param cdd.path deprecated and not used (a warning is given when it is set).
#' @param aligner_params a character string specifying additional parameters that shall be passed to DIAMOND or BLAST. Default is \code{aligner_params} = \code{NULL}. 
#' @param aligner_path a character string specifying the path to the DIAMOND or BLAST program (in case you don't use the default path).
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
#' \strong{Choosing the mode.} \code{"fast"} is tuned for closely related species. For distant species it misses
#' orthologs at almost no gain in time: \emph{S. cerevisiae} against \emph{C. glabrata} (reciprocal best hits, 2
#' cores) gives 4,164 pairs with \code{"fast"} (7 s), 4,384 with \code{"default"} (+5.3 \%, 8 s) and 4,502 with
#' \code{"sensitive"} (+8.1 \%, 24 s); the pairs found only by \code{"default"} have a median identity of 37 \%
#' and a median length of 298 amino acids. On the 1,000-gene \emph{A. thaliana} / \emph{A. lyrata} example the
#' three modes give 983, 985 and 985 pairs. Use \code{"default"} or \code{"sensitive"} when the proteins of the
#' two species are less than about 70 \% identical.
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND or BLAST hit detection.
#' @param ortho_path deprecated and not used (a warning is given when it is set).
#' @param aa_aln_type a character string specifying the amino acid alignment type:
#' \itemize{ 
#' \item \code{aa_aln_type = "pairwise"} (Default)
#' \item \code{aa_aln_type = "multiple"}
#' }.
#' @param aa_aln_tool a character string specifying the program that aligns the two proteins of a gene pair:
#' \code{"NW"} (Default; Needleman-Wunsch in R, with \code{aa_aln_type = "pairwise"}), or with
#' \code{aa_aln_type = "multiple"} one of \code{"clustalw"}, \code{"clustalo"}, \code{"mafft"}, \code{"muscle"},
#' \code{"t_coffee"} (see \code{\link{multi_aln}}).
#' @param aa_aln_path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @param aa_aln_params a character string with command-line parameters for the multiple alignment program selected via
#' \code{aa_aln_tool}. It is handed to \code{\link{multi_aln}} as \code{params} and \emph{replaces} the default parameters
#' of that program (see there; the alignment must still be written in CLUSTAL or FASTA format, e.g.
#' \code{aa_aln_tool = "mafft", aa_aln_params = "--maxiterate 1 --clustalout"}). It is not used (with a warning)
#' with \code{aa_aln_type = "pairwise"}. Until version 0.4.3 the argument was accepted but never passed on.
#' Note: this controls the amino acid alignment step only, not the initial DIAMOND or BLAST homology search
#' (use \code{aligner_params} for that). Default is \code{aa_aln_params} = \code{NULL} (the default parameters).
#' @param aa_sub_matrix the substitution matrix used to score pairwise amino acid alignments (\code{aa_aln_type = "pairwise"}):
#' the name of a matrix shipped with \pkg{pwalign}, e.g. \code{"BLOSUM62"} (Default), a numeric substitution matrix, or \code{NULL}
#' for the scoring of \pkg{orthologr} <= 0.4.3 (no substitution matrix). See \code{\link{pairwise_aln}} for details and for how to
#' reproduce dN/dS values computed with earlier versions.
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used. Default is \code{codon_aln_tool} = \code{"pal2nal"}.
#' Right now only "pal2nal" can be selected as codon alignment tool.
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
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
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform parallel computations on a multicore machine.
#' When it is larger than the number of available cores, the number of available cores is used and a warning is given.
#' @param quiet \code{TRUE} (default): one line at the start and one summary at the end of the run (as messages,
#' i.e. on \code{stderr}; nothing is written to \code{stdout}). \code{FALSE}: also the parameters, the messages of
#' the search program and a progress message about every tenth of the gene pairs.
#' Default is \code{quiet} = \code{FALSE}.
#' @param clean_folders kept for compatibility. All intermediate files of a call are written to a folder of its own inside
#' \code{tempdir()}, which is always removed when the call ends (see Details), so there is nothing left to clean.
#' @param print_citation a logical value indicating whether the citation note is shown (as a message, once per R session).
#' @author Hajk-Georg Drost and Jaruwatana Sodai Lotharukpong
#' @details 
#' 
#' The dN/dS ratio quantifies the mode and strength of selection acting on a pair of orthologous genes.
#' This selection pressure can be quantified by comparing synonymous substitution rates (dS) that are assumed to be neutral
#' with nonsynonymous substitution rates (dN), which are exposed to selection as they 
#' change the amino acid composition of a protein (Mugal et al., 2013 http://mbe.oxfordjournals.org/content/31/1/212).
#' 
#' The \pkg{orthologr} package provides the \code{\link{dNdS}} function to perform dNdS estimation on pairs of orthologous genes.
#' This function takes the CDS files of two organisms of interest (\code{query_file} and \code{subject_file}) 
#' and computes the dNdS estimation values for orthologous gene pairs between these organisms.
#' 
#' The following pipieline resembles the dNdS estimation process:
#' 
#' \itemize{
#'                 
#' \item 1) Orthology Inference: e.g. DIAMOND or BLAST reciprocal best hit (RBH)
#' 
#' \item 2) Pairwise sequence alignment: e.g. clustalw for pairwise amino acid sequence alignments
#' 
#' \item 3) Codon Alignment: e.g. pal2nal program
#' 
#' \item 4) dNdS estimation: e.g. Yang, Z. and Nielsen, R. (2000) \url{http://mbe.oxfordjournals.org/content/17/1/32.short}
#' 
#' }
#' Note: it is assumed that when using \code{dNdS()} all corresponding multiple sequence alignment programs you
#' want to use are already installed on your machine and are executable via either
#' the default execution \code{PATH} or you specifically define the location of the executable file
#' via the \code{aa_aln_path} or \code{aligner_path} argument that can be passed to \code{dNdS()}.
#' 
#' The \code{dNdS()} function can be used choosing the folllowing options:
#' 
#' \itemize{
#' \item \code{ortho_detection} : 
#'    \itemize{  
#'    \item \code{"RBH"} (DIAMOND or BLAST best reciprocal hit)
#'    \item \code{"BH"} (DIAMOND or BLAST best hit)
#'    \item \code{"Orthofinder2"}
#'      }
#'      
#' \item \code{aa_aln_type} : 
#'  \itemize{
#'   \item \code{"multiple"}
#'   \item \code{"pairwise"}
#'  }
#' 
#' \item \code{aa_aln_tool} :
#' \itemize{
#'  \item \code{"clustalw"}
#'  \item \code{"t_coffee"}
#'  \item \code{"muscle"}
#'  \item \code{"clustalo"}
#'  \item \code{"mafft"}
#'  \item \code{"NW"} (in case \code{aa_aln_type = "pairwise"})
#' }
#' 
#' \item \code{codon_aln_tool} :
#' \itemize{
#'  \item \code{"pal2nal"}
#'  }
#' \item \code{dnds_est.method} : 
#' \itemize{
#' \item "Li" : Li's method (1993)
#' \item "Comeron" : Comeron's method (1995)
#' \item "NG": Nei, M. and Gojobori, T. (1986)
#' \item "LWL": Li, W.H., et al. (1985)
#' \item "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' \item "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' \item "YN": Yang, Z. and Nielsen, R. (2000)
#' \item "MYN" (Modified YN): Zhang, Z., et al. (2006)
#' 
#' }
#' 
#' }
#' 
#' @references
#' 
#' seqinr: \url{https://cran.r-project.org/package=seqinr}
#' 
#' Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator: Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics 2006 , 4:259-263. 
#' 
#' KaKs_Calculator: \url{https://code.google.com/archive/p/kaks-calculator} [GNU GPL-3 license]
#' 
#' Paradis, E. (2012) Analysis of Phylogenetics and Evolution with R (Second Edition). New York: Springer.
#'
#' Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289-290.
#' 
#' More information on \pkg{ape} can be found at \url{https://emmanuelparadis.github.io/} or \url{https://github.com/emmanuelparadis/ape}.
#' 
#' Pages H, Aboyoun P, Gentleman R and DebRoy S. Biostrings: String objects representing biological sequences, and
#' matching algorithms. R package version 2.32.1.
#' 
#' @return A tibble with one row per gene pair: \code{query_id}, \code{subject_id}, \code{dN}, \code{dS},
#' \code{dNdS}, the columns of the hit table of the search, and, as the last column, \code{status} (see the
#' section "Status of a gene pair").
#' 
#' A gene pair for which a step of the pipeline fails (e.g. PAL2NAL or KaKs_Calculator exit with an error) does not end
#' the run: it is returned with \code{NA} for dN, dS and dNdS, a warning reports how many pairs failed, and
#' \code{attr(result, "failed_pairs")} is a \code{data.frame} (\code{query_id}, \code{subject_id}, \code{reason})
#' naming the pair, the program and its message.
#' 
#' Each call keeps its intermediate files (search database, alignments, KaKs_Calculator files) in a folder of its own
#' inside \code{tempdir()} and removes it when the call ends; \code{options(orthologr.keep_tmp = TRUE)} keeps it.
#' 
#' \code{attr(result, "dropped_sequences")} is a \code{data.frame} (\code{file}, \code{id}, \code{reason}) with the
#' input sequences that were left out before the search: coding sequences with characters other than A, C, G, T (they
#' cannot be translated; e.g. 5'-incomplete Ensembl CDS that are padded with N), empty sequences and, with
#' \code{delete_corrupt_cds = TRUE}, sequences whose length is not a multiple of 3. The summary message gives their
#' number. A sequence with a single \code{N} is left out as a whole. How much that costs depends on the source of
#' the files: 0 to 1 sequence in each of eight RefSeq sets of yeasts and bacteria and 2 of 1,000 in the
#' \emph{A. lyrata} example, but 51 of 2,048 (2.5 \%) of the Ensembl transcripts of human chromosome 21 and 365 of
#' 7,728 (4.7 \%) of those of mouse chromosomes 10, 16 and 17 - there typically one leading \code{N} of a transcript
#' whose start is not annotated (\code{cds_start_NF}). Sequences whose length is not a multiple of 3 (97 of 2,048
#' and 934 of 7,728 in the same files) are trimmed by default (\code{delete_corrupt_cds = FALSE}).
#' 
#' @section Status of a gene pair:
#' The last column of the result, \code{status}, says how each gene pair came out. The same words are used by
#' \code{dNdS()} and in the \code{pairs} tables of \code{\link{dNdS_orthogroups}}:
#' \describe{
#' \item{\code{"ok"}}{dN, dS and dN/dS were estimated;}
#' \item{\code{"dS_zero"}}{dN and dS were estimated and dS is 0 (no synonymous change), so that the ratio is not
#' defined: \code{dNdS} is \code{NA}. \code{\link{filter_dNdS}} leaves such pairs out or keeps them (\code{dS_zero});}
#' \item{\code{"saturated"}}{dN or dS is \code{NA} because the correction for multiple substitutions at a site is
#' not defined: the two sequences differ at so many sites of a class that the number of substitutions cannot be
#' estimated. A pair that is saturated at synonymous sites usually still has its dN;}
#' \item{\code{"no_data"}}{the pair has no codons that can be compared (the programs of \code{dNdS()} do not tell
#' this case from \code{"saturated"});}
#' \item{\code{"failed"}}{a step of the pipeline failed for the pair; the reason is in
#' \code{attr(result, "failed_pairs")}.}
#' }
#' At the end of a run one message gives the numbers by status and the median dS, and says so when a large part of
#' the pairs is saturated or the median dS is above 1.5: the two species are then too distant for pairwise dN/dS
#' (S. cerevisiae against C. glabrata: 72 \% of 4,164 pairs saturated, median dS of the others 1.9).
#' The counts are also in \code{attr(result, "orthologr_run")$status_counts}.
#' \code{options(orthologr.status_column = FALSE)} returns the table without the column (as 0.4.3 did; also set by
#' \code{\link{orthologr_legacy_options}}).
#'
#' @section Record of the run:
#' \code{attr(result, "orthologr_run")} is a list with the version of the package, the arguments of the call, every
#' \code{orthologr.*} option in force, the engine, the versions of the programs, the date, size and MD5 checksum of
#' the input files and the numbers of gene pairs by status. \code{\link{write_dNdS}} writes the table together with
#' this record and \code{\link{read_dNdS}} reads both back.
#'
#' @section Sequence IDs:
#' The ID of a sequence is its FASTA header up to the first white space. IDs must be unique within a file (a
#' duplicated ID is an error that names it) and are returned exactly as they are in the input files, whatever they
#' look like (numbers, \code{|}, \code{-}, more than 15 characters): the IDs of a gene pair are carried through the
#' alignment and estimation steps and are not read back from the output of the external programs.
#' 
#' @section Splice variants:
#' Every sequence of an input file is treated as a gene. A file that holds several transcripts per gene (for example
#' the \code{cds/all} files of Ensembl or files retrieved with \code{biomartr::getCDS()}) gives several rows per gene
#' with different partners and different dN/dS values, and the splice variants of a gene compete as best hits.
#' Please keep one transcript per gene first, e.g. with \code{\link{retrieve_longest_isoforms}}.
#' \code{dNdS()} warns when a file seems to hold several transcripts per gene. The check is a conservative heuristic
#' that never changes a result: (1) when at least 90\% of the FASTA headers name their gene (Ensembl
#' \code{gene:<ID>}, NCBI \code{[locus_tag=<ID>]} or \code{[gene=<symbol>]}) the genes are counted from these fields;
#' (2) otherwise a transcript suffix is removed from the IDs (\code{.<number>}, \code{.t<number>}, \code{-T<number>},
#' \code{-R<letter>}, \code{-P<letter>}, \code{-mRNA-<number>}, \code{_i<number>}, \code{.p<number>}) and a warning
#' is given only when at least 5 of the remaining stems, and at least 1\% of them, occur more than once.
#' \code{options(orthologr.isoform_warning = FALSE)} switches the check off.
#' @examplesIf all(nzchar(Sys.which(c("diamond", "perl"))))
#' # get a dNdS table using:
#' # 1) reciprocal best hit for orthology inference (RBH)
#' # 2) Needleman-Wunsch for pairwise amino acid alignments
#' # 3) pal2nal for codon alignments
#' # 4) Comeron for dNdS estimation
#' # 5) single core processing 'comp_cores = 1'
#' dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      ortho_detection = "RBH",
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW",
#'      codon_aln_tool  = "pal2nal",
#'      dnds_est.method = "Comeron",
#'      comp_cores      = 1 )
#'
#' \dontrun{
#' # running dNdS using the 'aa_aln_path' argument to specify the path to
#' # the corresponding alignment tool
#' dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      ortho_detection = "RBH",
#'      aa_aln_type     = "multiple",
#'      aa_aln_tool     = "clustalw",
#'      aa_aln_path     = "/path/to/clustalw/",
#'      codon_aln_tool  = "pal2nal",
#'      dnds_est.method = "Comeron",
#'      comp_cores      = 1 )
#'
#' # The same result can be obtained using multicore processing using: comp_cores = 2 or 3 or more ...
#' dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      ortho_detection = "RBH",
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW",
#'      codon_aln_tool  = "pal2nal",
#'      dnds_est.method = "Comeron",
#'      comp_cores      = 2 )
#' }
#' @seealso \code{\link{divergence_stratigraphy}}, \code{\link{orthologs}}, 
#' \code{\link{substitutionrate}}, \code{\link{multi_aln}}, \code{\link{codon_aln}}, 
#' \code{\link{diamond_best}}, \code{\link{diamond_rec}},
#' \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{read.cds}} 
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

dNdS <- function(query_file, 
                 subject_file, 
                 aligner         = "diamond",
                 sensitivity_mode = "fast",
                 aligner_path    = NULL, 
                 seq_type        = "cds",
                 format          = "fasta", 
                 ortho_detection = "RBH",
                 delete_corrupt_cds = FALSE,
                 store_locally   = FALSE,
                 cdd.path        = NULL,
                 aligner_params  = NULL, 
                 eval            = "1E-5", 
                 ortho_path      = NULL, 
                 aa_aln_type     = "pairwise", 
                 aa_aln_tool     = "NW", 
                 aa_aln_path     = NULL, 
                 aa_aln_params   = NULL, 
                 codon_aln_tool  = "pal2nal", 
                 kaks_calc_path  = NULL, 
                 dnds_est.method = "Comeron", 
                 comp_cores      = 1, 
                 quiet           = TRUE, 
                 clean_folders   = FALSE,
                 print_citation = TRUE,
                 aa_sub_matrix   = "BLOSUM62",
                 best_hit_rule   = getOption("orthologr.best_hit_rule", "bitscore")
                 ){
        
        best_hit_rule <- orthologr_best_hit_rule(if (missing(best_hit_rule)) NULL else best_hit_rule)
        
        # what is not supported is rejected now, before any program is started
        dnds_check_args(
                aligner          = aligner,
                sensitivity_mode = sensitivity_mode,
                seq_type         = seq_type,
                format           = format,
                ortho_detection  = ortho_detection,
                aa_aln_type      = aa_aln_type,
                aa_aln_tool      = aa_aln_tool,
                aa_aln_params    = aa_aln_params
        )
        
        if (!is.codon_aln_tool(codon_aln_tool))
                stop("Please choose a codon alignment tool that is supported by this function.", call. = FALSE)
        
        if (!is.dnds_est_method(dnds_est.method))
                stop("Please choose a dNdS estimation method that is supported by this function.", call. = FALSE)
        
        for (unused in c("cdd.path", "ortho_path")) {
                if (!is.null(get(unused)))
                        warning("The argument '", unused, "' is deprecated and not used.", call. = FALSE)
        }
        
        comp_cores <- check_comp_cores(comp_cores)
        
        check_aa_sub_matrix(aa_sub_matrix)
        
        # the programs of the chosen pipeline must be there before the search starts
        dnds_require_programs(aa_aln_tool = aa_aln_tool, aa_aln_path = aa_aln_path, codon_aln_tool = codon_aln_tool,
                              dnds_est.method = dnds_est.method, kaks_calc_path = kaks_calc_path)
        
        # the arguments of this call, for the record of the run (attr(<result>, "orthologr_run"))
        run_args <- run_arguments()
        
        # one line that says what is computed; the parameters in detail only when not quiet
        message("dNdS(): ", if (ortho_detection == "RBH") "reciprocal best hits" else "best hits", " of '",
                basename(query_file), "' against '", basename(subject_file), "' (", aligner,
                if (aligner == "diamond") paste0(", sensitivity_mode '", sensitivity_mode, "'"),
                ", e-value ", eval, ", best_hit_rule '", best_hit_rule, "'), then dN/dS by method '", dnds_est.method,
                "' (", aa_aln_type, " alignments with '", aa_aln_tool, "', comp_cores = ", comp_cores, ").")
        detail_message(quiet, "seq_type = '", seq_type, "'; codon_aln_tool = '", codon_aln_tool, "'; aa_sub_matrix = ",
                       if (is.character(aa_sub_matrix)) paste0("'", aa_sub_matrix, "'") else if (is.null(aa_sub_matrix)) "NULL" else "<matrix>",
                       "; delete_corrupt_cds = ", delete_corrupt_cds)
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # every call works in a folder of its own, which is removed when the call ends (also after an error or an
        # interrupt): files of an earlier or failed run can never be picked up as results of this one
        # (keep the folder with options(orthologr.keep_tmp = TRUE))
        end_tmpdir <- orthologr_tmpdir_start("orthologr_dNdS_")
        on.exit(end_tmpdir(), add = TRUE)
        
        # the two input files are read and translated once for all steps of this call
        end_seq_cache <- seq_cache_start()
        on.exit(end_seq_cache(), add = TRUE)
        
        # sequences that are left out are collected for attr(<result>, "dropped_sequences")
        end_dropped <- dropped_start()
        on.exit(end_dropped(), add = TRUE)
        
        # several transcripts per gene in an input file? (a warning; see ?dNdS, section 'Splice variants')
        isoforms_query <- isoform_warning(query_file, "query_file")
        isoform_warning(subject_file, "subject_file")
        
        if (store_locally) {
                if (file.exists("orthologr_alignment_files")) {
                        message("The folder 'orthologr_alignment_files' seems to exist already.",
                                "Please make sure to delete this folder or store the previous results in a different location if you don't want files to be overwritten.",
                                "The existing folder 'orthologr_alignment_files' is used to store alignment files ...")
                } else {
                        message("Creating folder 'orthologr_alignment_files' to store alignment files ...")
                        dir.create("orthologr_alignment_files")
                }
                if (!file.exists(file.path("orthologr_alignment_files", "_pairwise_alignment_with_score")))
                        dir.create(file.path("orthologr_alignment_files", "_pairwise_alignment_with_score"))
        }
        
                
        aa <- geneids <- NULL
        
        # blast/diamond each translated aminoacid sequence against the related database to get a 
        # hit table with pairs of geneids  
        
        detail_message(quiet, "Starting Orthology Inference ...")
        
        # choosing the aligner between diamond and blast
        # (aligner and ortho_detection were checked above)
        {
                hit.table <- dnds_search(query_file, subject_file, aligner = aligner, ortho_detection = ortho_detection,
                                         aligner_params = aligner_params, sensitivity_mode = sensitivity_mode,
                                         aligner_path = aligner_path, delete_corrupt_cds = delete_corrupt_cds,
                                         comp_cores = comp_cores, seq_type = seq_type, eval = eval, format = format,
                                         best_hit_rule = best_hit_rule, quiet = quiet)
                
                # (both files were read and translated for the search: they are not read again)
                q_cds <- read_cds_once(query_file, delete_corrupt_cds = delete_corrupt_cds)
                
                s_cds <- read_cds_once(subject_file, delete_corrupt_cds = delete_corrupt_cds)
                
                # The proteins of both organisms, translated from the coding sequences read above (the same
                # translation set_diamond() / set_blast() perform for the search).
                # They used to be read back from the files 'query_<basename(file)>.fasta' that the search step leaves
                # in '_blast_db': when query_file and subject_file have the same file name in different folders
                # (e.g. 'a/cds.fa' and 'b/cds.fa') both names are the same file, and the query proteins were
                # silently replaced by the subject's.
                q_aa <- cds_protein_table_file(query_file, delete_corrupt_cds = delete_corrupt_cds)
                s_aa <- cds_protein_table_file(subject_file, delete_corrupt_cds = delete_corrupt_cds)
        }

        data.table::setnames(q_cds, old = c("geneids", "seqs"), new = c("query_id","query_cds"))
        data.table::setnames(s_cds, old = c("geneids", "seqs"), new = c("subject_id","subject_cds"))
        data.table::setnames(q_aa, old = c("geneids", "seqs"), new = c("query_id","query_aa"))
        data.table::setnames(s_aa, old = c("geneids", "seqs"), new = c("subject_id","subject_aa"))
        
        # joining all tables to a final table containing: query_id, subject_id, query_aa_seq, subject_aa_seq, query_cds_seq, and subject_cds_seq
        query_tbl <- dplyr::inner_join(tibble::as_tibble(q_cds), tibble::as_tibble(q_aa), by = "query_id")
        subject_tbl <- dplyr::inner_join(tibble::as_tibble(s_cds), tibble::as_tibble(s_aa), by = "subject_id")
        joint_query_tbl <- dplyr::inner_join(tibble::as_tibble(hit.table), tibble::as_tibble(query_tbl), by = "query_id")
        joint_subject_tbl <- dplyr::inner_join(tibble::as_tibble(hit.table), tibble::as_tibble(subject_tbl), by = "subject_id")
        final_tbl <- dplyr::inner_join(tibble::as_tibble(joint_query_tbl), tibble::as_tibble(joint_subject_tbl), by = c("query_id","subject_id"))
        
        sensitivity_hint(hit.table, aligner = aligner, sensitivity_mode = sensitivity_mode)
        
        if (nrow(hit.table) == 0)
                orthologr_abort(class = "no_hits", "No orthologs could be found: the ", if (aligner == "diamond") "DIAMOND" else "BLAST",
                     " search of '", basename(query_file), "' against '", basename(subject_file),
                     "' returned no ", if (ortho_detection == "RBH") "reciprocal best hits" else "hits",
                     " at e-value ", eval, ". Please check your input files (are both sets of coding sequences related?) ",
                     "or relax 'eval'.", call. = FALSE)
        
        if (nrow(final_tbl) == 0)
                orthologr_abort("No orthologs could be found! Please check your input files!", class = "no_hits")
        
        
        detail_message(quiet, "Orthology Inference Completed: ", nrow(final_tbl), " gene pairs.")
        detail_message(quiet, "Starting dN/dS Estimation ...")
        dNdS_tbl <- compute_dnds( complete_tbl    = data.table::as.data.table(final_tbl),
                                  aa_aln_type     = aa_aln_type,
                                  aa_aln_tool     = aa_aln_tool,
                                  aa_aln_path     = aa_aln_path,
                                  aa_aln_params   = aa_aln_params,
                                  codon_aln_tool  = codon_aln_tool, 
                                  kaks_calc_path  = kaks_calc_path, 
                                  store_locally   = store_locally,
                                  dnds_est.method = dnds_est.method, 
                                  quiet           = quiet,
                                  comp_cores      = comp_cores, 
                                  clean_folders   = clean_folders,
                                  aa_sub_matrix   = aa_sub_matrix )
        
        # the dN/dS rows are joined with the hit table by the pair (query_id and subject_id; it was query_id only)
        res <- dplyr::inner_join(dNdS_tbl, tibble::as_tibble(hit.table), by = c("query_id", "subject_id"))
        
        if (nrow(res) != nrow(dNdS_tbl))
                warning(nrow(dNdS_tbl) - nrow(res), " of ", nrow(dNdS_tbl), " gene pairs with dN/dS values have no row ",
                        "in the hit table and are missing in the result. Please report this as a bug of orthologr.",
                        call. = FALSE)
        
        # The end of the run: the status of every gene pair (last column), the sequences that were left out before
        # the search (attr "dropped_sequences": file, id, reason), the gene pairs that could not be processed
        # (attr "failed_pairs"), the record of the run (attr "orthologr_run") and one summary message.
        kaks_used <- is.element(dnds_est.method, kaks_calculator_methods())
        res <- dnds_finish(
                res, failed_pairs = attr(dNdS_tbl, "failed_pairs"), dropped = end_dropped(),
                record = list(fun = "dNdS", args = run_args, engine = "orthologr",
                              tools = c(run_tool_versions(aligner = aligner, aligner_path = aligner_path,
                                                          aa_aln_tool = aa_aln_tool, aa_aln_path = aa_aln_path,
                                                          codon_aln_tool = codon_aln_tool, kaks = kaks_used,
                                                          kaks_calc_path = kaks_calc_path)),
                              files = c(query_file, subject_file)),
                print_citation = print_citation
        )
        isoform_result_note(res, isoforms_query)
        
       # return the dNdS table for all query_ids and subject_ids
           return(res)   
}

# The search step shared by dNdS() and dNdS_coal(): best hits (ortho_detection = "BH") or reciprocal best hits
# ("RBH") with DIAMOND or BLAST+, as a keyed data.table. (When quiet, the banners and timings of the search program
# are not shown.)
dnds_search <- function(query_file, subject_file, aligner, ortho_detection, aligner_params, sensitivity_mode,
                        aligner_path, delete_corrupt_cds, comp_cores, seq_type, eval, format, best_hit_rule, quiet) {
        align <- switch(paste(aligner, ortho_detection),
                "diamond BH"  = function(...) diamond_best(..., diamond_params = aligner_params, sensitivity_mode = sensitivity_mode),
                "diamond RBH" = function(...) diamond_rec(..., diamond_params = aligner_params, sensitivity_mode = sensitivity_mode),
                "blast BH"    = function(...) blast_best(..., blast_params = aligner_params),
                "blast RBH"   = function(...) blast_rec(..., blast_params = aligner_params),
                orthologr_abort("aligner = '", aligner, "' with ortho_detection = '", ortho_detection, "' is not supported.",
                                class = "bad_argument"))
        hit.table <- quiet_search(
                align(
                        query_file         = query_file,
                        subject_file       = subject_file,
                        delete_corrupt_cds = delete_corrupt_cds,
                        path               = aligner_path,
                        comp_cores         = comp_cores,
                        seq_type           = seq_type,
                        eval               = eval,
                        format             = format,
                        best_hit_rule      = best_hit_rule
                ), quiet = quiet)
        data.table::setDT(hit.table)
        data.table::setkeyv(hit.table, c("query_id", "subject_id"))
        hit.table
}
