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
#' @param aligner_params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string specifying additional parameters that shall be passed to DIAMOND or BLAST. Default is \code{aligner_params} = \code{NULL}. 
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
#' @param aa_aln_params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string with command-line parameters for the multiple alignment program selected via
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
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used: \code{"pal2nal"} (Default;
#' the PAL2NAL Perl script shipped with this package) or \code{"coal"} (COAL: the same codon alignments without Perl, computed
#' inside R when the package \pkg{coalr} is installed; see \code{\link{codon_aln}}).
#' @param kaks_tool the program that computes the KaKs_Calculator methods (\code{dnds_est.method} = \code{"NG"}, \code{"LWL"},
#' \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"}, \code{"YN"}, \code{"MYN"}): \code{"kaks_calculator"} (Default; the external
#' program) or \code{"coal"} (COAL's implementation of these methods, no external program; see \code{\link{substitutionrate}}).
#' With a COAL that computes them, \code{kaks_tool = "coal"} also computes \code{"Comeron"} and \code{"Li"} (same
#' values as the built-in code for the same codon alignment).
#' @param engine \code{"orthologr"} (Default): each gene pair is processed on its own (protein alignment, codon alignment and
#' dN/dS estimation with the tools chosen above, with intermediate files per pair). \code{"coal"}: the batched fast path
#' \code{\link{dNdS_coal}}: orthologs are detected as usual, then all gene pairs are aligned and estimated by COAL in one
#' call, in memory and on \code{comp_cores} threads. \code{engine = "coal"} needs the package \pkg{coalr}, a
#' \code{dnds_est.method} that COAL computes (\code{"NG"}, \code{"LWL"}, \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"},
#' \code{"YN"}, \code{"MYN"}; and \code{"Comeron"} - the default method - and \code{"Li"} with a \pkg{coalr} that has
#' them, see \code{\link{dNdS_coal}}), \code{ortho_detection} \code{"RBH"} or \code{"BH"}, and ignores
#' \code{aa_aln_type}, \code{aa_aln_tool}, \code{codon_aln_tool}, \code{kaks_tool} and \code{store_locally}. It
#' returns the same table, with the same attributes.
#' @param aa_aln_scoring with \code{engine = "coal"}: how COAL scores the pairwise protein alignments -
#' \code{"pwalign"} (Default; the scoring of \code{\link{pairwise_aln}}, so that the alignments are those of
#' \code{engine = "orthologr"} up to co-optimal alignments) or \code{"needle"} (COAL's own defaults). The same
#' argument as in \code{\link{dNdS_coal}}, \code{\link{dNdS_clusters}} and \code{\link{dNdS_orthogroups}}.
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
#' @param quiet a logical value specifying whether the output of the corresponding alignment tool shall be printed out to the console.
#' Default is \code{quiet} = \code{FALSE}.
#' @param clean_folders kept for compatibility. All intermediate files of a call are written to a folder of its own inside
#' \code{tempdir()}, which is always removed when the call ends (see Details), so there is nothing left to clean.
#' @param print_citation a logical value indicating whether or not the citation message shall be printed.
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
#'  \item \code{"coal"}
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
#' @return A data.table storing the dNdS values of the correspnding genes.
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
#' \code{delete_corrupt_cds = TRUE}, sequences whose length is not a multiple of 3. A message gives their number.
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
#' @examples \dontrun{
#' 
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
#' 
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
#' 
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
                 kaks_tool       = "kaks_calculator",
                 engine          = "orthologr",
                 aa_aln_scoring  = "pwalign",
                 best_hit_rule   = getOption("orthologr.best_hit_rule", "bitscore")
                 ){
        
        best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)
        
        if (!is.element(engine, c("orthologr", "coal")))
                stop("Please choose engine = \"orthologr\" or engine = \"coal\".", call. = FALSE)
        
        # what is not supported is rejected now, before any program is started
        dnds_check_args(
                aligner          = aligner,
                sensitivity_mode = sensitivity_mode,
                seq_type         = seq_type,
                format           = format,
                ortho_detection  = ortho_detection,
                # engine = "coal" does not use the alignment arguments
                aa_aln_type      = if (engine == "coal") "pairwise" else aa_aln_type,
                aa_aln_tool      = if (engine == "coal") "NW" else aa_aln_tool,
                aa_aln_params    = if (engine == "coal") NULL else aa_aln_params
        )
        
        # (engine = "coal" checks its methods itself - it takes several at once - and has no codon alignment tool)
        if (engine != "coal") {
                if (!is.codon_aln_tool(codon_aln_tool))
                        stop("Please choose a codon alignment tool that is supported by this function.", call. = FALSE)
                
                if (length(dnds_est.method) != 1 || !is.dnds_est_method(dnds_est.method))
                        stop("Please choose a dNdS estimation method that is supported by this function.", call. = FALSE)
        }
        
        for (unused in c("cdd.path", "ortho_path")) {
                if (!is.null(get(unused)))
                        warning("The argument '", unused, "' is deprecated and not used.", call. = FALSE)
        }
        
        comp_cores <- check_comp_cores(comp_cores)
        
        check_aa_sub_matrix(aa_sub_matrix)
        
        if (!is.element(kaks_tool, c("kaks_calculator", "coal")))
                stop("Please choose kaks_tool = \"kaks_calculator\" or kaks_tool = \"coal\".", call. = FALSE)
        
        if (kaks_tool == "coal" && !is.element(dnds_est.method, c("Comeron", "Li", coal_methods())))
                stop("kaks_tool = \"coal\" computes the methods ", paste(coal_methods(), collapse = ", "),
                     "; dnds_est.method = \"", dnds_est.method, "\" needs kaks_tool = \"kaks_calculator\".", call. = FALSE)
        
        # COAL must be there before the search starts, not when the first gene pair needs it
        if (engine == "orthologr" && (identical(codon_aln_tool, "coal") ||
                                      (kaks_tool == "coal" && is.element(dnds_est.method, coal_methods()))))
                coal_backend()
        
        # the programs of the chosen pipeline must be there before the search starts (engine = "coal" needs none)
        if (engine == "orthologr")
                dnds_require_programs(aa_aln_tool = aa_aln_tool, aa_aln_path = aa_aln_path, codon_aln_tool = codon_aln_tool,
                                      dnds_est.method = dnds_est.method, kaks_calc_path = kaks_calc_path,
                                      need_kaks = !(kaks_tool == "coal" && is.element(dnds_est.method, coal_methods())))
        
        check_aa_aln_scoring(aa_aln_scoring)
        if (engine == "orthologr" && aa_aln_scoring != "pwalign")
                stop("aa_aln_scoring = \"", aa_aln_scoring, "\" needs engine = \"coal\" (engine = \"orthologr\" aligns with ",
                     "pairwise_aln(), i.e. \"pwalign\").", call. = FALSE)
        
        # the batched fast path: all gene pairs in one call to COAL
        if (engine == "coal") {
                if (is.null(aa_sub_matrix) || !is.character(aa_sub_matrix))
                        stop("engine = \"coal\" needs the name of a substitution matrix in aa_sub_matrix (e.g. \"BLOSUM62\").", call. = FALSE)
                return(dNdS_coal(
                        query_file         = query_file,
                        subject_file       = subject_file,
                        aligner            = aligner,
                        sensitivity_mode   = sensitivity_mode,
                        aligner_path       = aligner_path,
                        seq_type           = seq_type,
                        format             = format,
                        ortho_detection    = ortho_detection,
                        delete_corrupt_cds = delete_corrupt_cds,
                        aligner_params     = aligner_params,
                        eval               = eval,
                        dnds_est.method    = dnds_est.method,
                        aa_aln_scoring     = aa_aln_scoring,
                        aa_sub_matrix      = aa_sub_matrix,
                        comp_cores         = comp_cores,
                        quiet              = quiet,
                        print_citation     = print_citation,
                        best_hit_rule      = best_hit_rule
                ))
        }
        
        message("\n")
        message("Starting orthology inference (",ortho_detection, ") and dNdS estimation (",dnds_est.method,") using the follwing parameters:")
        message("query = '", basename(query_file), "'")
        message("subject = '", basename(subject_file), "'")
        message("aligner = ", "'", basename(aligner), "'")
        if (aligner == "diamond")
                message("sensitivity_mode = ", "'", basename(sensitivity_mode), "'")
        message("seq_type = '", seq_type, "'")
        message("e-value: ", eval)
        message("best_hit_rule = '", best_hit_rule, "'")
        message("aa_aln_type = '", aa_aln_type, "'")
        message("aa_aln_tool = '", aa_aln_tool, "'")
        message("comp_cores = '", comp_cores, "'")
        message("\n")
        
        
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
        
        message("Starting Orthology Inference ...")
        
        # choosing the aligner between diamond and blast
        # (aligner and ortho_detection were checked above)
        {
                if (aligner == "diamond" & ortho_detection == "BH")
                        align <- function(...) 
                                diamond_best(...,
                                             diamond_params = aligner_params,
                                             sensitivity_mode = sensitivity_mode)
                if (aligner == "diamond" & ortho_detection == "RBH")
                        align <- function(...) 
                                diamond_rec(...,
                                            diamond_params = aligner_params,
                                            sensitivity_mode = sensitivity_mode)
                if (aligner == "blast" & ortho_detection == "BH")
                        align <- function(...) 
                                blast_best(...,
                                           blast_params = aligner_params)
                if (aligner == "blast" & ortho_detection == "RBH")
                        align <- function(...) 
                                blast_rec(...,
                                          blast_params = aligner_params)
                
                # run orthology inference.
                hit.table <- 
                        align(
                                query_file   = query_file,
                                subject_file = subject_file,
                                delete_corrupt_cds = delete_corrupt_cds,
                                path         = aligner_path,
                                comp_cores   = comp_cores,
                                seq_type     = seq_type,
                                eval         = eval,
                                format       = format,
                                best_hit_rule = best_hit_rule
                        )
                
                data.table::setDT(hit.table)
                data.table::setkeyv(hit.table, c("query_id","subject_id"))
                
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
        
        if (nrow(hit.table) == 0)
                stop("No orthologs could be found: the ", if (aligner == "diamond") "DIAMOND" else "BLAST",
                     " search of '", basename(query_file), "' against '", basename(subject_file),
                     "' returned no ", if (ortho_detection == "RBH") "reciprocal best hits" else "hits",
                     " at e-value ", eval, ". Please check your input files (are both sets of coding sequences related?) ",
                     "or relax 'eval'.", call. = FALSE)
        
        if (nrow(final_tbl) == 0)
                stop("No orthologs could be found! Please check your input files!", call. = FALSE)
        
        
        message("Orthology Inference Completed.")
        message("Starting dN/dS Estimation ...")
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
                                  aa_sub_matrix   = aa_sub_matrix,
                                  kaks_tool       = kaks_tool )
        
        # the dN/dS rows are joined with the hit table by the pair (query_id and subject_id; it was query_id only)
        res <- dplyr::inner_join(dNdS_tbl, tibble::as_tibble(hit.table), by = c("query_id", "subject_id"))
        
        if (nrow(res) != nrow(dNdS_tbl))
                warning(nrow(dNdS_tbl) - nrow(res), " of ", nrow(dNdS_tbl), " gene pairs with dN/dS values have no row ",
                        "in the hit table and are missing in the result. Please report this as a bug of orthologr.",
                        call. = FALSE)
        
        # sequences that were left out before the search (file, id, reason), and a note on splice variants
        dropped <- end_dropped()
        if (!is.null(dropped)) {
                attr(res, "dropped_sequences") <- dropped
                if (nrow(dropped) > 0)
                        message(nrow(dropped), " input sequence", if (nrow(dropped) > 1) "s were" else " was",
                                " left out before the search (", paste0(table(dropped$reason), " x ", names(table(dropped$reason)), collapse = "; "),
                                "). They are listed in attr(<result>, \"dropped_sequences\").")
        }
        isoform_result_note(res, isoforms_query)
        
        # gene pairs that could not be processed (their dN, dS and dNdS are NA) and why
        attr(res, "failed_pairs") <- attr(dNdS_tbl, "failed_pairs")
        
        message("dN/dS Estimation Completed.")
        
        if (print_citation) {
                message("\n")
                message("Please cite the following paper when using orthologr for your own research:")
                cat("Drost et al. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. 2015. Mol. Biol. Evol. 32 (5): 1221-1231.")
                message("\n")
                message("\n")
        }
        
       # return the dNdS table for all query_ids and subject_ids
           return(res)   
}






