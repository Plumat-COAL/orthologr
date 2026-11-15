#' @title Internal function for dNdS computations
#' @description This function takes a pairwise alignment as input file and estimates the
#' dNdS ratio of the corresponding alignment. Nevertheless, this function is a helper function for
#' \code{\link{dNdS}}. For dNdS computations you should use the function: \code{\link{dNdS}}.
#' @param file a character string specifying the path to a codon alignment file
#' @param est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li" .
#' Note, that when using "Comeron" as dNdS estimation method, the program 'gestimator' is used to compute the
#' corresponding dNdS values from a given alignment. The program 'gestimator' can only read "fasta" files,
#' hence it is important to use format = "fasta" when choosing est.method = "Comeron".
#' @param format a character string specifying the file format in which the alignment is stored:  
#' "mase", "clustal", "phylip", "fasta" , "msf"
#' @param quiet a logical value specifying whether the output of the coresponding interface shall be printed out.
#' @param kaks_calc.params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string storing additional parameters for KaKs_Calculator. Default is \code{NULL}.
#' The parameters are added to the method chosen with \code{est.method} (e.g. \code{est.method = "YN"},
#' \code{kaks_calc.params = "-c 11"} runs \code{-m YN -c 11}; up to orthologr 0.4.3 they replaced the method).
#' With \code{est.method = "kaks_calc"} the methods themselves are given here, e.g.
#' \code{kaks_calc.params} = "-m NG -m YN". 
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
#' @param subst_name a character string specifying the substitution name that shall be added to the internal folder path naming. 
#' Default is \code{subst_name} = \code{NULL}.
#' @param kaks_tool the program that computes the KaKs_Calculator methods (\code{est.method} = \code{"NG"}, \code{"LWL"},
#' \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"}, \code{"YN"}, \code{"MYN"}): \code{"kaks_calculator"} (Default; the
#' external program KaKs_Calculator) or \code{"coal"} (COAL's implementation of the same methods, which gives
#' KaKs_Calculator 2.0's values for the same alignment and runs inside R when the package \pkg{coalr} is installed;
#' see \code{\link{codon_aln}} for how to install COAL). With \code{"coal"} a dN or dS of exactly 0 is returned as
#' 0 (KaKs_Calculator prints \code{NA} for values below 1e-6) and a saturated estimate is \code{NA}.
#' The methods \code{"MA"}, \code{"MS"}, \code{"GY"}, \code{"GMYN"} and \code{"ALL"} need KaKs_Calculator.
#' \code{est.method = "Comeron"} and \code{"Li"} are always computed by \pkg{orthologr}'s own code (gestimator,
#' \code{seqinr::kaks()}), whatever \code{kaks_tool} is: \code{kaks_tool} only chooses the program for the
#' KaKs_Calculator methods. (COAL's ports of the two methods, which \code{dNdS(engine = "coal")} uses, give the same
#' table for the same alignment.)
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function takes a pairwise alignments file as input and estimates dNdS ratios
#' of the corresponding alignments using predefined dNdS estimation methods.
#' 
#' The dNdS estimation methods available in this function are:
#' 
#' \itemize{
#' \item "Li" : Li's method (1993) -> provided by the seqinr package (\code{seqinr::kaks()})
#' 
#' \item "Comeron" : Comeron's method (1995)
#' }
#' dNdS estimation methods provided by KaKs_Calculator 1.2 :
#' 
#' Approximate Methods:
#' 
#' \itemize{
#' \item "NG": Nei, M. and Gojobori, T. (1986)
#' 
#' \item "LWL": Li, W.H., et al. (1985)
#' 
#' \item "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' 
#' \item "MLWL": (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' 
#' \item "YN": Yang, Z. and Nielsen, R. (2000)
#' 
#' \item "MYN": (Modified YN): Zhang, Z., et al. (2006)
#' 
#' \item "GMYN": Wang, D.P., et al. Biology Direct. (2009)
#' \item "GY": Goldman, N. and Yang, Z. (1994)
#' \item "MS": (Model Selection)
#' \item "MA": (Model Averaging): based on a set of candidate models, Posada, D. (2003)
#' }
#' 
#' Maximum-Likelihood Methods:
#' \itemize{
#' \item GY: Goldman, N. and Yang, Z. (1994)
#' }
#' 
#' @examples
#' # estimate the dNdS rate using Li's method
#' substitutionrate(
#'    file       = system.file("seqs/pal2nal.aln", package = "orthologr"),
#'    est.method = "Li",
#'    format     = "fasta")
#'
#' # estimate the dNdS rate using Comeron's method
#' substitutionrate(
#'    file       = system.file("seqs/pal2nal.aln", package = "orthologr"),
#'    est.method = "Comeron",
#'    format     = "fasta")
#'
#' @examplesIf nzchar(Sys.which("KaKs_Calculator"))
#' # estimate the dNdS rate using model averaging provided by
#' # the KaKs_Calculator 1.2 program
#'  substitutionrate(
#'     file       = system.file("seqs/pal2nal.aln", package = "orthologr"),
#'     est.method = "MA",
#'     format     = "fasta")
#'
#'  # estimate the dNdS rate using Nei and Gojobori's method provided by the
#'  # KaKs_Calculator 1.2 program
#'  substitutionrate(
#'     file       = system.file("seqs/pal2nal.aln", package = "orthologr"),
#'     est.method = "NG",
#'     format     = "fasta")
#'
#'   # estimate the dNdS rate using Nei and Gojobori's method AND Yang and Nielsen's
#'   # method provided by the KaKs_Calculator 1.2 program
#'   # for this purpose we choose:
#'   # est.method = "kaks_calc" and kaks_calc.params = "-m NG -m YN"
#'  substitutionrate(
#'    file             = system.file("seqs/pal2nal.aln", package = "orthologr"),
#'    est.method       = "kaks_calc",
#'    format           = "fasta",
#'    kaks_calc.params = "-m NG -m YN")
#' @references 
#' Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol., 36:96-99.
#' 
#' Charif, D. and Lobry, J.R. (2007) SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.
#' 
#' Thornton, K. (2003) libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics 19(17): 2325-2327
#' 
#' Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator:
#' Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics 2006 , 4:259-263.
#' 
#' https://code.google.com/archive/p/kaks-calculator
#' 
#' @section Gene pairs without nonsynonymous change (KaKs_Calculator methods):
#' KaKs_Calculator prints \code{NA} for every Ka or Ks below 1e-6, also when the estimate is exactly 0. Up to
#' orthologr 0.4.3 these \code{NA} were returned as they are: a pair without nonsynonymous change had
#' \code{dN = NA} and \code{dNdS = 0} (and \code{\link{filter_dNdS}} / \code{\link{divergence_stratigraphy}}
#' removed it - the most conserved genes of a comparison), and a pair of identical sequences had \code{NA} for dN,
#' dS and dNdS (status \code{"saturated"}). Since 0.5.0:
#' \itemize{
#' \item \code{dN = 0} when Ka is \code{NA}, Ks is a positive number and Ka/Ks is 0;
#' \item \code{dN = 0} and \code{dS = 0} (\code{dNdS = NA}, status \code{"dS_zero"}) when no nucleotide differs
#' (column "Substitutions" of KaKs_Calculator 0);
#' \item \code{dS = 0} (\code{dNdS = NA}, status \code{"dS_zero"}) when Ks is \code{NA}, Ka is a number and no
#' difference is synonymous (methods NG, YN, MYN, GMYN: "S-Substitutions"; LWL, MLWL, LPB, MLPB: no transition at
#' 2-fold and no difference at 4-fold sites);
#' \item \code{dN = 0} (and \code{dNdS = 0} when dS is a positive number) when Ka is \code{NA} and no difference is
#' nonsynonymous (NG, YN, MYN, GMYN: "N-Substitutions"; LWL family: no difference at 0-fold and no transversion at
#' 2-fold sites). (GMYN: since 0.6.0 exactly 0; it was Ka/Ks * Ks, a rounding residue of about 1e-15.)
#' }
#' These are the values and the status the COAL engine (\code{kaks_tool = "coal"}, \code{engine = "coal"}) gives
#' for the same pairs. When Ka/Ks is \code{NA} for another reason (a model that failed) the values stay \code{NA}.
#' \code{options(orthologr.legacy_kaks_zero_dn = TRUE)} restores the old table.
#'
#' @section Saturation (KaKs_Calculator methods):
#' Where the correction for multiple substitutions is not defined, KaKs_Calculator prints a capped number instead
#' of \code{NA}: Ka or Ks = 99 (YN family, when the proportion of differences P + Q exceeds 1); for YN and MYN
#' JC69's distance at P = 0.75 (n - 1) / n, i.e. 0.75 log(n) for n synonymous or nonsynonymous sites (Ks of 3 - 5
#' for typical genes); for LWL, MLWL, LPB and MLPB 0 substitutions for a class of sites whose correction is not
#' defined (an underestimate). Since 0.5.0 these estimates are \code{NA} (and \code{dNdS} with them; status
#' \code{"saturated"}), as with the COAL engine; they are recognised in the columns KaKs_Calculator prints (sites
#' and differences by class), so the other values are unchanged. The methods GY, MA, MS and GMYN are not checked.
#' (GMYN: KaKs_Calculator 2.0 decides whether it writes the header line of its output with memory it never
#' initialises, so the line can be missing; orthologr then supplies it, so the values read do not depend on it.
#' Since 0.6.0 GMYN is computed in chunks as the other methods; a run that came back without GMYN values for many
#' pairs with an earlier version had lost these header lines.)
#' Up to 0.4.3 such a value was returned as an estimate: e.g. dS = 99, which makes dN/dS look like strong
#' purifying selection, put the pair into the most conserved stratum of \code{\link{divergence_map}} and passed
#' \code{\link{filter_dNdS}}. \code{options(orthologr.legacy_kaks_saturation = TRUE)} restores the old values.
#'
#' For the LWL family the printed counts (6 significant digits) cannot tell the sign of \code{1 - 2P - Q} when it
#' is within about 1e-5 of 0 - typically exactly 0, where KaKs_Calculator's own rounding decides. For such pairs
#' (an argument within 1e-3 of 0 by the printed counts) the counts are computed again from the codon alignment, with
#' the operations of KaKs_Calculator 2.0 in the same order, and its own test is applied to them; so the status
#' follows what KaKs_Calculator did (checked on all LWL, MLWL, LPB and MLPB pairs of S. cerevisiae / C. glabrata,
#' S. cerevisiae / K. lactis and E. coli K-12 / O157 against a KaKs_Calculator that reports its decision: no pair
#' differs; with the printed counts alone 21 pairs differed, among them 2 MLPB pairs whose undefined 4-fold class
#' gave underestimated Ks of 0.97 and 2.26 that passed as estimates). Where KaKs_Calculator computed such a class
#' its Ks is large and unreliable (9 - 19 on these pairs). The recount needs the genetic code: it is read from
#' \code{-c} in \code{kaks_calc.params} (default 1); with a code KaKs_Calculator 2.0 does not have, or a
#' KaKs_Calculator whose counts differ from the recount, the printed counts decide. A build of KaKs_Calculator
#' that rounds differently (e.g. with fused multiply-add on arm64) may decide an exact 0 the other way.
#'
#' \strong{Iterations that do not converge (YN, MYN).} KaKs_Calculator iterates the estimates of YN and MYN for at
#' most 10 rounds and prints the values of the last round without saying whether they converged. Such a pair
#' cannot be told from an estimate and has status \code{"ok"}: S. cerevisiae / C. glabrata 18 of 605 YN estimates
#' (3\%) and 7 of 352 MYN estimates, S. cerevisiae / K. lactis 17 of 597 YN estimates, E. coli K-12 / O157 3 of
#' 3,412 YN estimates; all of them with dS above 1.4. \code{kaks_tool = "coal"} (and \code{dNdS(engine = "coal")})
#' computes the same methods and reports these pairs (\code{NA}; counted as "not converged" in the summary message
#' of \code{dNdS(engine = "coal")}); for the other pairs its values agree with KaKs_Calculator's within 5e-6. Treat
#' dS above about 1.5 from KaKs_Calculator with caution, or use COAL where the status matters.
#'
#' \strong{Estimates that failed (\code{nan}).} When KaKs_Calculator prints \code{nan} for Ka or Ks (MYN:
#' S. cerevisiae / C. glabrata 50 of 4,164 pairs, all of them iterations that did not converge), the row is
#' \code{NA} and is a failed estimate: \code{substitutionrate()} names it in \code{attr(result, "kaks_failed")}
#' (method, reason), and \code{\link{dNdS}} gives the pair the status \code{"failed"} with the reason in
#' \code{attr(result, "failed_pairs")}. It is not saturated (COAL estimates 43 of these 50 pairs). Up to this
#' version the row was \code{NaN} (in one of the 50 pairs dN a number) with status \code{"saturated"};
#' \code{options(orthologr.legacy_kaks_saturation = TRUE)} returns it as printed.
#'
#' @section Saturation (Comeron's method):
#' The corrections of Comeron (1995) take logarithms of \code{1 - 2P - Q} and \code{1 - 2Q} (\code{P}, \code{Q}:
#' proportions of transitional and transversional differences in a class of sites). When such an argument is not
#' positive, the observed differences exceed what the model can correct for and the distance is not defined. In
#' that case \code{dS} (synonymous classes) or \code{dN} (nonsynonymous classes) and \code{dNdS} are returned as
#' \code{NA}. Up to orthologr 0.4.3 the logarithm was replaced by 0 and a finite, too small value was returned
#' (e.g. \code{dS = 0} for \code{ATGGCTAAA} / \code{ATGGCCAAA}, which differ by one synonymous substitution; now
#' \code{NA}); \code{options(orthologr.legacy_comeron_saturation = TRUE)} reproduces those values. Values of gene pairs
#' without such a term are the same in both modes.
#'
#' A class of sites that has no sites (weight 0 in the estimate) never makes an estimate \code{NA}; this is also
#' the rule of COAL's Comeron method, so \code{engine = "coal"} gives the same \code{NA} pattern.
#' \strong{Short alignments:} with few codons a class with very few sites (typically the twofold sites whose
#' synonymous change is a transversion: a median of 1.75 sites in 30-codon alignments) is easily "saturated" by one
#' or two differences. Since 0.5.0 a class with fewer than \code{options(orthologr.comeron_min_sites)} sites
#' (default 3) cannot make \code{dN} or \code{dS} \code{NA}: its undefined logarithm is 0, as in 0.4.3 - but only
#' where the value that results is greater than 0 (since 0.6.0): a 0 from such a class would read as "no change"
#' (e.g. \code{dS = 0} for \code{ATGGCTAAA} / \code{ATGGCCAAA}, whose fourfold class has a single site) and stays
#' \code{NA}. A value that is not \code{NA} otherwise never changes. In a simulation without selection this keeps \code{dN} for the
#' 3.7 \% of 30-codon pairs that lost it (their mean absolute error 8.8 \%, against 11.2 \% for all 30-codon pairs);
#' \code{options(orthologr.comeron_min_sites = 0)} gives \code{NA} for them again. Comeron's method (like the
#' others) is unreliable below about 50 - 100 codons.
#'
#' @section Codons that differ at two positions (Comeron's method):
#' For such codons Comeron (1995) counts the changes along the two pathways codon 1 - intermediate 1 - codon 2 and
#' codon 1 - intermediate 2 - codon 2. libsequence's \code{TwoSubs}, from which this code comes, passes its
#' arguments in another order and counts codon 1 - intermediate 1, intermediate 1 - intermediate 2, codon 1 -
#' codon 2 and codon 2 - intermediate 2 instead, with which the estimate depends on which sequence is first
#' (swapping the two sequences of 300 random pairs changed Ka in 279 of them). Since 0.5.0 the two pathways are
#' counted (\code{options(orthologr.comeron_twosubs = "comeron")}, the default); \code{"libsequence"} gives the
#' numbers of 0.4.3 (set by \code{\link{orthologr_legacy_options}}). The COAL engine
#' (\code{engine = "coal"}, \code{kaks_tool = "coal"}) follows both options.
#' Effect on dS (same alignments): 1,000-gene example median 0 \%, 95th percentile 2.9 \% (\code{dNdS()}); human / mouse
#' median 1.0 \%, 95th percentile 5.0 \% (measured with the saturation rule of 0.4.3).
#' 
#' The alignment must be a codon alignment: at least two sequences of the same length, which is a multiple of 3.
#' Anything else is an error.
#' @return A data.table storing the query_id, subject_id, dN, dS, and dNdS values or 
#' a data.table storing the query_id, method, dN, dS, and dNdS values when using KaKs_Calculator.
#' If the dNdS value cannot be calculated NA is returned. This can happen because of constraints
#' of the used model. As each program throws different exception values we set all of them to NA instead. 
#' @seealso \code{\link{dNdS}}, \code{\link{multi_aln}}, \code{\link{codon_aln}}, \code{\link{blast_best}},
#' \code{\link{blast_rec}}, \code{\link{read.cds}}
#' @export
substitutionrate <- function(file, 
                             est.method, 
                             format           = "fasta",
                             quiet            = FALSE, 
                             kaks_calc.params = NULL,
                             kaks_calc_path   = NULL, 
                             subst_name       = NULL,
                             kaks_tool        = "kaks_calculator"){
        
        # dNdS estimation methods provided by the KaKs_Calculator program
        if (!is.null(kaks_calc_path))
                kaks_calc_path <- normalizePath(kaks_calc_path, mustWork = FALSE)
        
        kaks_calc_methods <- kaks_calculator_methods()
        
        if (length(est.method) != 1 || !is.dnds_est_method(est.method))
                dnds_method_abort(est.method, "substitutionrate()", kaks_calc = TRUE)
        
        if (!is.element(kaks_tool, c("kaks_calculator", "coal")))
                orthologr_abort(class = "bad_argument", "Please choose kaks_tool = \"kaks_calculator\" or kaks_tool = \"coal\".", call. = FALSE)
        
        if (kaks_tool == "coal" && !is.element(est.method, c(coal_native_methods(), coal_kaks_calculator_methods())))
                kaks_tool_coal_abort("est.method", est.method)
        
        if (!is.element(format, c("mase", "clustal", "phylip", "fasta" , "msf")))
                orthologr_abort(class = "bad_argument", "Please choose a format that is supported by seqinr::read.alignment.", call. = FALSE)
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- dN <- dS <- NULL
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if (!file.exists(file.path(orthologr_tmpdir(), "_calculation"))) {
                dir.create(file.path(orthologr_tmpdir(), "_calculation"), recursive = TRUE)
        }
        
        
        if (est.method == "Comeron") {
                # file in fasta required
                if (format != "fasta")
                        stop("To use gestimator an alignment file in fasta format is required.")
                return(comeron_rates(file, subst_name = subst_name, quiet = quiet))
        }
        
        if (est.method == "Li") {
                aln <- seqinr::read.alignment(file   = file,
                                              format = format)
                
                res_list <- seqinr::kaks(aln)
                
                # seqinr::kaks documentation:
                # When the alignment does not contain enough information (i.e we approach saturation),
                # the Ka and Ks values take the value 10. Negative values indicate that Ka and Ks can
                # not be computed.
                # we set both of them to NA
                #
                # dN, dS and dNdS are numeric columns. (They used to be built with t(c(names, values)), which made
                # every column a character column: the comparisons with 0 and 9.9 were then comparisons of strings
                # ("9.95e-05" > "9.9" is TRUE, so a tiny value was taken for saturation and set to NA), and the
                # result held text where every other method returns numbers.)
                li <- li_values(res_list$ka, res_list$ks)
                
                res <-
                        data.table::data.table(
                                query_id   = as.character(unlist(aln$nam))[1],
                                subject_id = as.character(unlist(aln$nam))[2],
                                dN         = li$dN,
                                dS         = li$dS,
                                dNdS       = li$dNdS
                        )
                
                data.table::setkey(res, query_id)
                
                if (!quiet) {
                        message("Substitutionrate successfully calculated using Li's method.")
                }
                
                return(res)
        }
        
        # the KaKs_Calculator methods computed by COAL (in process, or by the coal-cpp program): no KaKs_Calculator
        # (Comeron and Li were computed above, by orthologr's own code, whatever kaks_tool is)
        if (kaks_tool == "coal" && is.element(est.method, coal_kaks_calculator_methods())) {
                if (format != "fasta")
                        orthologr_abort(class = "bad_argument", "kaks_tool = \"coal\" needs the codon alignment in fasta format.", call. = FALSE)
                res <- coal_kaks_file(file = file, est.method = est.method, subst_name = subst_name)
                if (!quiet)
                        message("Substitutionrate successfully calculated by COAL.")
                return(res)
        }
        
        if (is.element(est.method, kaks_calc_methods)) {
                
                # (with kaks_calc_path, a program that is not there is reported below, together with the gene pair)
                if (is.null(kaks_calc_path))
                        is_installed_kaks_calculator()
                
                # the gene pair, for messages
                context <- pair_context(file, subst_name)
                
                # The AXT input of KaKs_Calculator: a name line, the two aligned sequences, an empty line. It is
                # written here (not by the Perl script parseFastaIntoAXT.pl, which labelled the pair with
                # '<query>-<subject>' and so made IDs with a hyphen impossible to recover).
                aln_seqs <- seqinr::read.fasta(file, as.string = TRUE, set.attributes = FALSE, forceDNAtolower = FALSE)
                if (length(aln_seqs) != 2)
                        stop("The codon alignment of ", context, " must hold two sequences, but '", file, "' holds ",
                             length(aln_seqs), ".", call. = FALSE)
                pair_ids <- fasta_ids(file)[1:2]
                
                # input, output and log of this call have names of their own and are removed afterwards: a result
                # of an earlier (or failed) call can never be read in place of the result of this one
                calc_dir <- file.path(orthologr_tmpdir(), "_calculation")
                axt_file <- tempfile(paste0(safe_file_name(basename(file)), "_"), tmpdir = calc_dir, fileext = ".axt")
                kaks_file <- paste0(axt_file, ".kaks")
                on.exit(unlink(c(axt_file, kaks_file)), add = TRUE)
                
                writeLines(c(if (is.null(subst_name)) "pair" else subst_name,
                             as.character(aln_seqs[[1]]), as.character(aln_seqs[[2]]), ""), axt_file)
                
                run_tool(
                        "KaKs_Calculator",
                        args = c(list("-i", axt_file, "-o", kaks_file),
                                 kaks_calc_args(est.method, kaks_calc.params)),
                        path    = kaks_calc_path,
                        tool    = "KaKs_Calculator",
                        context = context,
                        output  = kaks_file
                )
                
                kaks_tbl <- tryCatch({
                        utils::read.csv(text = kaks_output_lines(readLines(kaks_file, warn = FALSE)), sep = "\t", header = TRUE)
                }, error = function(e) {
                        stop("The output of KaKs_Calculator for ", context, " ('", kaks_file, "') could not be read: ",
                             conditionMessage(e), call. = FALSE)
                })
                
                return(kaks_result_table(kaks_tbl, pair_ids, paste0(context, " (method ", est.method, ")"),
                                         aln = c(as.character(aln_seqs[[1]]), as.character(aln_seqs[[2]])),
                                         genetic_code = kaks_genetic_code(kaks_calc.params)))
        }
        
        # every method accepted above is handled by one of the branches; a method that is not must never give NULL
        orthologr_abort(class = "bad_argument", "est.method = \"", est.method, "\" is not handled by substitutionrate().", call. = FALSE)
}

# Comeron's method (gestimator) for the codon alignment in 'file' (FASTA): the table of substitutionrate()
# (query_id, subject_id, dN, dS, dNdS; key query_id), or with values = TRUE its values as a "dnds_values" record
# (dnds_values(); compute_dnds() makes one table of the records of a chunk of gene pairs instead of a table per pair).
comeron_rates <- function(file, subst_name = NULL, quiet = TRUE, values = FALSE) {
        query_id <- NULL
        if (!file.exists(file.path(orthologr_tmpdir(), "_calculation")))
                dir.create(file.path(orthologr_tmpdir(), "_calculation"), recursive = TRUE)
        # gestimator writes its table to a file that only this call uses and that is removed afterwards
        file.out <- tempfile(
                paste0(if (is.null(subst_name)) "" else paste0(subst_name, "_"), "gestimout_"),
                tmpdir = file.path(orthologr_tmpdir(), "_calculation")
        )
        on.exit(unlink(file.out), add = TRUE)
        
        # the options are checked here, so that a wrong value is an error of its own (not a failed pair)
        twosubs_libsequence <- identical(comeron_twosubs(), "libsequence")
        min_sites <- comeron_min_sites()
        
        tryCatch({
                # the alignment is checked here as well as in the C++ code: at least two sequences of the
                # same length, which is a multiple of 3 and not 0
                check_codon_alignment_file(file)
                gestimator(file     = file,
                           file_out = file.out,
                           legacy_saturation = comeron_legacy_saturation(),
                           twosubs_libsequence = twosubs_libsequence,
                           min_sites = min_sites)
                
                # The columns are read as what they are: without colClasses fread() guesses the types,
                # and IDs that look like numbers ("1001", "007", "1e5", "1.10") or like logical values
                # ("TRUE", "NA") came back as another ID or as NA.
                if (values) {
                        # (one line: the table of a pair, which setkey() below leaves in its order)
                        v <- read_gestimator_values(readLines(file.out, warn = FALSE))
                        if (!is.null(v) && length(v$dN) == 1)
                                return(comeron_values(v))
                }
                hit.table <- read_gestimator_table(file.out)
                
                data.table::setkey(hit.table, query_id)
                
                # for consistent output format we set the default value of gestimator 999
                # to NA for output
                # Therefor the user cannot make a mistake if using the dNdS results without
                # filtering for some value.
                # (data.table::set() with the rows: three ':=' sub-assignments took about 3 ms per gene pair)
                for (col in c("dS", "dN", "dNdS"))
                        data.table::set(hit.table, i = which(hit.table[[col]] == 999), j = col, value = NA_real_)
                
                if (!quiet) {
                        message("Substitutionrate successfully calculated by gestimator.")
                }
                
                if (values)
                        return(comeron_values(hit.table))
                return(hit.table)
        }, error = function(e) {
                stop(
                        "dN/dS estimation with Comeron's method (gestimator) failed for ",
                        pair_context(file, subst_name), ": ", conditionMessage(e),
                        call. = FALSE
                )
        })
}

# The "dnds_values" record of the gestimator table 'v' (list or data.table with dN, dS, dNdS): gestimator's 999
# is NA, as in comeron_rates().
comeron_values <- function(v) {
        na999 <- function(x) { x[which(x == 999)] <- NA_real_; x }
        dnds_values(dN = na999(v$dN), dS = na999(v$dS), dNdS = na999(v$dNdS))
}

# options(orthologr.legacy_comeron_saturation = TRUE): Comeron's method as up to 0.4.3, where an undefined
# logarithm of the correction (saturation) was replaced by 0 and a finite dS or dN was returned.
comeron_legacy_saturation <- function() {
        isTRUE(orthologr_option("orthologr.legacy_comeron_saturation", FALSE))
}

# options(orthologr.comeron_twosubs): how Comeron's method counts codons that differ at two positions. "comeron"
# (default): the two pathways of Comeron (1995). "libsequence": the branches libsequence's TwoSubs counts
# (codon1 - int1, int1 - int2, codon1 - codon2, codon2 - int2), with which the estimate depends on which sequence
# is first; the numbers of orthologr 0.4.3 (set by orthologr_legacy_options()).
comeron_twosubs <- function() {
        v <- orthologr_option("orthologr.comeron_twosubs", "comeron")
        if (!is.character(v) || length(v) != 1 || !is.element(v, c("comeron", "libsequence")))
                orthologr_abort(class = "bad_argument",
                                "options(orthologr.comeron_twosubs = ) must be \"comeron\" or \"libsequence\".",
                                call. = FALSE)
        v
}

# options(orthologr.comeron_min_sites): a site class of Comeron's method with fewer sites than this cannot make dN or
# dS NA (its undefined logarithm is 0, as in 0.4.3); default 3; 0: every class with sites can.
comeron_min_sites <- function() {
        v <- orthologr_option("orthologr.comeron_min_sites", 3)
        if (!is.numeric(v) || length(v) != 1 || is.na(v) || !is.finite(v) || v < 0)
                orthologr_abort(class = "bad_argument",
                                "options(orthologr.comeron_min_sites = ) must be one number, 0 or more.", call. = FALSE)
        as.numeric(v)
}

# Stop unless `file` (FASTA) is a codon alignment: at least two sequences, all of the same length, the length a
# multiple of 3 and at least 3. Line ends may be LF or CRLF.
check_codon_alignment_file <- function(file) {
        if (!file.exists(file))
                orthologr_abort(class = "bad_input", "The alignment file '", file, "' does not exist.", call. = FALSE)
        x <- sub("\r$", "", readLines(file, warn = FALSE))
        x <- x[nzchar(x)]
        head_at <- startsWith(x, ">")
        if (sum(head_at) < 2)
                stop("the alignment file holds ", sum(head_at), " sequence", if (sum(head_at) != 1) "s",
                     "; a codon alignment of at least two sequences is needed.", call. = FALSE)
        len <- as.vector(tapply(nchar(x) * !head_at, cumsum(head_at), sum))
        if (length(unique(len)) != 1)
                stop("the aligned sequences have unequal lengths (", paste(utils::head(len, 5), collapse = ", "),
                     "): this is not a codon alignment.", call. = FALSE)
        if (len[1] < 3)
                stop("the aligned sequences are shorter than one codon (length ", len[1], ").", call. = FALSE)
        if (len[1] %% 3 != 0)
                stop("the length of the aligned sequences (", len[1], ") is not a multiple of 3: this is not a codon alignment.",
                     call. = FALSE)
        invisible(TRUE)
}

# The table gestimator writes (tab-separated, no header: query id, subject id, dN, dS, dN/dS), with the IDs read
# as character whatever they look like.
# read_gestimator_values() reads the lines where that is possible (see below; NULL otherwise).
# (P2 of REVIEW-OL-1) The table of one gene pair is a single line: it is split with strsplit() and the numbers are
# converted with as.numeric(), which is much faster than starting fread() for it. Only where fread() reads the
# fields the same way (5 fields per line, no empty line; plain decimal numbers; IDs of printable ASCII without white
# space - fread() strips white space at their ends - and not empty); anything else is read by fread() as before.
read_gestimator_table <- function(file) {
        lines <- tryCatch(readLines(file, warn = FALSE), error = function(e) NULL)
        v <- read_gestimator_values(lines)
        if (!is.null(v))
                return(data.table::data.table(query_id = v$query_id, subject_id = v$subject_id,
                                              dN = v$dN, dS = v$dS, dNdS = v$dNdS))
        res <- data.table::fread(
                file       = file,
                sep        = "\t",
                header     = FALSE,
                quote      = "",
                colClasses = c("character", "character", "numeric", "numeric", "numeric"),
                col.names  = c("query_id", "subject_id", "dN", "dS", "dNdS"),
                na.strings = ""
        )
        res
}

# The columns of the lines of a gestimator table (query_id, subject_id, dN, dS, dNdS) where they can be read without
# fread() and give what fread() reads (see read_gestimator_table()); NULL otherwise.
read_gestimator_values <- function(lines) {
        if (length(lines) == 0 || !all(nzchar(lines)))
                return(NULL)
        parts <- strsplit(lines, "\t", fixed = TRUE)
        if (!all(lengths(parts) == 5L))
                return(NULL)
        m <- matrix(unlist(parts, use.names = FALSE), nrow = 5)
        ids <- m[1:2, , drop = FALSE]
        if (!(all(nzchar(ids)) && !any(grepl("[^!-~]", ids, perl = TRUE)) && gestimator_plain_numbers(m[3:5, ])))
                return(NULL)
        list(query_id = m[1, ], subject_id = m[2, ], dN = as.numeric(m[3, ]), dS = as.numeric(m[4, ]),
             dNdS = as.numeric(m[5, ]))
}

# Are the numbers 'x' (text) plain decimal numbers, which as.numeric() and fread() read the same way?
gestimator_plain_numbers <- function(x) {
        all(grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", x, perl = TRUE))
}

# The method arguments of a KaKs_Calculator call. The method asked for with 'est.method' is always handed over
# ('-m <est.method>'); 'kaks_calc.params' is added to it. (It used to replace it: est.method = "YN" with
# kaks_calc.params = "-c 1" ran KaKs_Calculator's default model MA.) est.method = "kaks_calc" stands for "the
# methods are given in kaks_calc.params" and adds nothing.
kaks_calc_args <- function(est.method, kaks_calc.params = NULL) {
        has_params <- !is.null(kaks_calc.params) && any(nzchar(trimws(kaks_calc.params)))
        if (identical(est.method, "kaks_calc")) {
                if (!has_params)
                        stop("est.method = \"kaks_calc\" needs the methods in 'kaks_calc.params', e.g. kaks_calc.params = \"-m NG -m YN\".",
                             call. = FALSE)
                return(list(I(kaks_calc.params)))
        }
        c(list("-m", est.method), if (has_params) list(I(kaks_calc.params)))
}

# The result table of substitutionrate() from KaKs_Calculator's output table 'kaks_tbl' (columns Sequence, Method,
# Ka, Ks, Ka/Ks, ...) for the gene pair 'pair_ids' (query id, subject id).
#   aln          : NULL, or the two aligned coding sequences handed to KaKs_Calculator (for the saturation test of
#                  the LWL family, see kaks_saturated())
#   genetic_code : KaKs_Calculator's genetic code ('-c', 1 by default)
kaks_result_table <- function(kaks_tbl, pair_ids, context = "the gene pair", aln = NULL, genetic_code = 1L) {
        query_id <- NULL
        v <- kaks_result_values(kaks_tbl, context = context, aln = aln, genetic_code = genetic_code)
        kaks_tbl_res <- data.frame(
                query_id   = rep(as.character(pair_ids[1]), nrow(kaks_tbl)),
                subject_id = rep(as.character(pair_ids[2]), nrow(kaks_tbl)),
                dN         = v$dN,
                dS         = v$dS,
                dNdS       = v$dNdS,
                method     = v$method,
                stringsAsFactors = FALSE
        )
        kaks_tbl_res <- data.table::as.data.table(kaks_tbl_res)
        data.table::setkey(kaks_tbl_res, query_id)
        if (any(v$failed))
                data.table::setattr(kaks_tbl_res, "kaks_failed",
                                    data.frame(method = v$method[v$failed], reason = v$reason[v$failed],
                                               stringsAsFactors = FALSE))
        return(kaks_tbl_res)
}

# The values of kaks_result_table() for the rows of 'kaks_tbl', which may belong to several gene pairs (the rows of a
# chunk computed by one KaKs_Calculator process, kaks_calc_batch()): every step works row by row, so the values of
# a row are those of a table of its gene pair alone. Returns list(dN, dS, dNdS, method, failed, reason): 'failed'
# marks the rows KaKs_Calculator printed 'nan' for (attr(, "kaks_failed") of kaks_result_table(): their method and
# 'reason'; NA for the other rows).
#   context : the gene pair, for messages - one string, one per row, or a function of row numbers that returns
#             them (called only for the rows that need one)
#   aln     : NULL, the two aligned sequences of the pair (one pair), or a list with them for every row
kaks_result_values <- function(kaks_tbl, context = "the gene pair", aln = NULL, genetic_code = 1L) {
        ctx_of <- function(rows) if (is.function(context)) context(rows) else if (length(context) == 1) context else context[rows]
        if (nrow(kaks_tbl) == 0 || ncol(kaks_tbl) < 5)
                stop("KaKs_Calculator returned no result for ", ctx_of(1L)[1], ".", call. = FALSE)
        num <- function(x) suppressWarnings(as.numeric(as.character(x)))
        dN   <- num(kaks_tbl[ , 3])
        dS   <- num(kaks_tbl[ , 4])
        dNdS <- num(kaks_tbl[ , 5])
        method <- kaks_tbl[ , 2]
        if (!kaks_legacy_zero_dn()) {
                fixed <- kaks_zero_dn(dN, dS, dNdS)
                dN <- fixed$dN
                dNdS <- fixed$dNdS
                zero <- kaks_zero_ds(kaks_tbl)
                if (!is.null(zero)) {
                        dN[zero$identical | zero$no_nonsynonymous] <- 0
                        dS[zero$identical | zero$no_synonymous] <- 0
                        dNdS[zero$identical | zero$no_synonymous] <- NA_real_
                        # (Ka/Ks * Ks of kaks_zero_dn() can be a rounding residue such as 2.5e-12 there)
                        no_non <- zero$no_nonsynonymous & is.finite(dS) & dS > 0
                        dNdS[no_non] <- 0
                }
        }
        failed <- rep(FALSE, length(dN))
        reason <- rep(NA_character_, length(dN))
        if (!kaks_legacy_saturation()) {
                sat <- kaks_saturated(kaks_tbl, aln = aln, genetic_code = genetic_code)
                if (!is.null(sat)) {
                        dS[sat$Ks] <- NA_real_
                        dN[sat$Ka] <- NA_real_
                        dNdS[sat$Ks | sat$Ka] <- NA_real_
                }
                # Ka or Ks printed as 'nan' (or '-nan'): the estimate failed (MYN of S. cerevisiae / C. glabrata: 50 of
                # 4,164 pairs, all of them iterations that did not converge). This is not saturation: COAL estimates
                # most of these pairs (43 of the 50). The row is NA and named in attr(, "kaks_failed"); dNdS() reports
                # it as a failed pair (status "failed", the reason in attr(<result>, "failed_pairs")). It was NaN
                # (the dN of one pair a number) and the status "saturated".
                nan <- is.nan(num(kaks_tbl[ , 3])) | is.nan(num(kaks_tbl[ , 4]))
                if (any(nan)) {
                        dN[nan] <- NA_real_
                        dS[nan] <- NA_real_
                        dNdS[nan] <- NA_real_
                        ctx <- ctx_of(which(nan))
                        failed <- nan
                        reason[nan] <- paste0("KaKs_Calculator printed 'nan' for ",
                                              ifelse(is.nan(num(kaks_tbl[nan, 3])) & is.nan(num(kaks_tbl[nan, 4])), "Ka and Ks",
                                                     ifelse(is.nan(num(kaks_tbl[nan, 4])), "Ks", "Ka")),
                                              " (method ", method[nan], ") of ",
                                              sub(" \\(method [^)]*\\)$", "", ctx),
                                              ": its estimate failed (for YN / MYN typically an iteration that did ",
                                              "not converge); the pair is not saturated")
                }
        }
        list(dN = dN, dS = dS, dNdS = dNdS, method = method, failed = failed, reason = reason)
}

# The columns of a KaKs_Calculator table by position (1.2 and 2.0 write the same 21 columns; the names are read by
# read.csv() with 'check.names', so they are not used): NULL when the table does not have them.
kaks_columns <- function(kaks_tbl) {
        if (ncol(kaks_tbl) < 15)
                return(NULL)
        num <- function(x) suppressWarnings(as.numeric(as.character(x)))
        # "a:b:c" (the three fold classes 0, 2, 4) -> a matrix of three columns; "NA" means all three are 0
        fold <- function(x, na_is_zero) {
                x <- as.character(x)
                m <- matrix(if (na_is_zero) 0 else NA_real_, length(x), 3)
                has <- !is.na(x) & grepl(":", x, fixed = TRUE)
                if (any(has))
                        m[has, ] <- matrix(num(unlist(strsplit(x[has], ":", fixed = TRUE))), ncol = 3, byrow = TRUE)
                m
        }
        list(
                method = as.character(kaks_tbl[, 2]),
                Ka     = num(kaks_tbl[, 3]),
                Ks     = num(kaks_tbl[, 4]),
                S      = num(kaks_tbl[, 8]),
                N      = num(kaks_tbl[, 9]),
                L      = fold(kaks_tbl[, 10], na_is_zero = FALSE),
                snp    = num(kaks_tbl[, 11]),
                Sd     = num(kaks_tbl[, 12]),
                Nd     = num(kaks_tbl[, 13]),
                Si     = fold(kaks_tbl[, 14], na_is_zero = TRUE),
                Vi     = fold(kaks_tbl[, 15], na_is_zero = TRUE)
        )
}

# The methods of the LWL family (Li, Wu and Luo 1985; Li 1993; Tzeng et al. 2004) and of Yang and Nielsen (2000).
kaks_lwl_methods <- function() c("LWL", "MLWL", "LPB", "MLPB")
kaks_yn_methods <- function() c("YN", "MYN")
# The methods whose columns "S-Substitutions" and "N-Substitutions" count the synonymous and nonsynonymous
# differences (KaKs_Calculator prints NA for a count of 0). GMYN counts as MYN does (its Ka/Ks for a pair without
# nonsynonymous difference is a residue such as 1.3e-15, its dN Ka/Ks * Ks; GY, MA, MS print expected numbers).
kaks_counted_methods <- function() c("NG", kaks_yn_methods(), "GMYN")

# Gene pairs for which KaKs_Calculator printed 'NA' for a dS that is 0 (it prints NA for every value below 1e-6):
#   identical     : no nucleotide differs ("Substitutions" 0) - dN = dS = 0 (KaKs_Calculator: NA NA NA)
#   no_synonymous : Ks NA, Ka a number, and no synonymous difference - dS = 0 (all differences are
#                   nonsynonymous). A synonymous difference is counted in "S-Substitutions" (NG, YN, MYN, GMYN) or, for
#                   the LWL family, as a transition at 2-fold sites or any difference at 4-fold sites (columns
#                   "Fold-S-Substitutions" / "Fold-N-Substitutions": transitions / transversions per class). The
#                   other methods (GY, MA, MS, ...) are left as they are.
#   no_nonsynonymous : Ka NA and no nonsynonymous difference ("N-Substitutions"; LWL family: differences at 0-fold
#                   sites and transversions at 2-fold sites) - dN = 0 (also when Ks is saturated, where
#                   kaks_zero_dn() cannot use Ka/Ks).
# The COAL engine reports the same pairs with dS = 0 (status "dS_zero").
kaks_zero_ds <- function(kaks_tbl) {
        k <- kaks_columns(kaks_tbl)
        if (is.null(k))
                return(NULL)
        identical <- is.na(k$Ka) & is.na(k$Ks) & !is.na(k$snp) & k$snp == 0
        no_syn <- rep(FALSE, length(k$Ka))
        counted <- is.element(k$method, kaks_counted_methods())
        no_syn[counted] <- is.na(k$Sd[counted]) | k$Sd[counted] <= 1e-6
        lwl <- is.element(k$method, kaks_lwl_methods())
        no_syn[lwl] <- (k$Si[lwl, 2, drop = FALSE] + k$Si[lwl, 3, drop = FALSE] + k$Vi[lwl, 3, drop = FALSE]) == 0
        no_syn <- !identical & no_syn & is.na(k$Ks) & is.finite(k$Ka) & k$Ka > 0 & !is.na(k$snp) & k$snp > 0
        # the same for dN when Ks is not a number either (kaks_zero_dn() needs Ka/Ks): no nonsynonymous difference
        no_non <- rep(FALSE, length(k$Ka))
        no_non[counted] <- is.na(k$Nd[counted]) | k$Nd[counted] <= 1e-6
        no_non[lwl] <- (k$Si[lwl, 1, drop = FALSE] + k$Vi[lwl, 1, drop = FALSE] + k$Vi[lwl, 2, drop = FALSE]) == 0
        no_non <- !identical & no_non & is.na(k$Ka) & !is.na(k$snp) & k$snp > 0
        list(identical = identical, no_synonymous = no_syn %in% TRUE, no_nonsynonymous = no_non %in% TRUE)
}

# options(orthologr.legacy_kaks_saturation = TRUE): KaKs_Calculator's values for saturated estimates are returned as
# printed, as up to 0.4.3.
kaks_legacy_saturation <- function() {
        isTRUE(orthologr_option("orthologr.legacy_kaks_saturation", FALSE))
}

# Estimates that KaKs_Calculator replaced with a capped value because the correction for multiple substitutions is
# not defined (saturation); they are NA, as with the COAL engine (status "saturated"). Returns list(Ks, Ka) of
# logicals, or NULL for a table without the columns needed.
#   - every method: Ka or Ks of exactly 99, the largest distance of the YN family ('maxt = 99' when P + Q > 1);
#   - YN, MYN: JC69's distance at its cap. When F84 and K80 fail and the proportion of differences P is at least
#     0.75, KaKs_Calculator (as PAML yn00) sets P to 0.75 (n - 1) / n, which gives t = 0.75 log(n) for n sites:
#     Ks = 0.75 log(S-Sites) or Ka = 0.75 log(N-Sites) to the printed precision;
#   - LWL, MLWL, LPB, MLPB: a class of sites (0-, 2-, 4-fold) with 1 - 2P - Q <= 0 or 1 - 2Q <= 0 (P, Q:
#     transitions and transversions per site of the class); KaKs_Calculator then counts 0 substitutions for that
#     class. 2- or 4-fold classes saturate Ks, 0- or 2-fold classes saturate Ka (the rule of COAL's LWL family).
#     P and Q are taken from the counts KaKs_Calculator prints (6 significant digits). Where that gives an
#     argument within 1e-3 of 0, the printed digits cannot tell the sign (S. cerevisiae / K. lactis: a 4-fold class
#     with 1 - 2P - Q = 0 exactly passed as 'ok', with Ks 0.97 instead of saturated; others were taken for saturated
#     although KaKs_Calculator gave a value): the counts are then computed again from the codon alignment 'aln', with
#     the same operations in the same order as KaKs_Calculator 2.0 (kaks_lwl_counts()), and its own test is applied
#     to them. When 'aln' is not given, or the counts computed do not agree with the printed ones (another version
#     of KaKs_Calculator), the printed counts decide.
#   aln          : NULL, the two aligned sequences of the pair (one pair for all rows of 'kaks_tbl'), or a list with
#                  the two aligned sequences (or NULL) for every row
#   genetic_code : the genetic code of KaKs_Calculator ('-c')
kaks_saturated <- function(kaks_tbl, aln = NULL, genetic_code = 1L) {
        k <- kaks_columns(kaks_tbl)
        if (is.null(k))
                return(NULL)
        n <- length(k$Ka)
        ks <- k$Ks %in% 99
        ka <- k$Ka %in% 99
        yn <- is.element(k$method, kaks_yn_methods())
        at_cap <- function(d, sites) yn & is.finite(d) & is.finite(sites) & sites > 1 & abs(d - 0.75 * log(sites)) <= 1e-5 * d
        ks <- ks | at_cap(k$Ks, k$S)
        ka <- ka | at_cap(k$Ka, k$N)
        lwl <- which(is.element(k$method, kaks_lwl_methods()))
        if (length(lwl) > 0) {
                L <- k$L[lwl, , drop = FALSE]
                P <- k$Si[lwl, , drop = FALSE] / L
                Q <- k$Vi[lwl, , drop = FALSE] / L
                a1 <- 1 - 2 * P - Q
                a2 <- 1 - 2 * Q
                sat <- !is.na(L) & L > 0 & !(a1 > 0 & a2 > 0)
                sat[is.na(sat)] <- FALSE
                # near 0 the printed digits do not tell the sign: KaKs_Calculator's own numbers decide
                near <- which(rowSums(!is.na(L) & L > 0 & ((!is.na(a1) & abs(a1) < 1e-3) | (!is.na(a2) & abs(a2) < 1e-3))) > 0)
                if (length(near) > 0 && (is.list(aln) || length(aln) == 2)) {
                        for (r in near) {
                                # (a list: the alignment of the gene pair of each row)
                                a <- if (is.list(aln)) aln[[lwl[r]]] else aln
                                if (length(a) != 2)
                                        next
                                exact <- kaks_lwl_counts(a[1], a[2], method = k$method[lwl[r]], genetic_code = genetic_code)
                                if (is.null(exact) || !kaks_counts_agree(exact, L[r, ], k$Si[lwl[r], ], k$Vi[lwl[r], ]))
                                        next
                                sat[r, ] <- exact$L > 0 & !((1 - 2 * exact$P - exact$Q) > 0 & (1 - 2 * exact$Q) > 0)
                        }
                }
                ks[lwl] <- ks[lwl] | sat[, 2] | sat[, 3]
                ka[lwl] <- ka[lwl] | sat[, 1] | sat[, 2]
        }
        list(Ks = ks, Ka = ka)
}

# options(orthologr.legacy_kaks_zero_dn = TRUE): the table of KaKs_Calculator as up to 0.4.3, where a gene pair
# without nonsynonymous change has dN = NA and dNdS = 0.
kaks_legacy_zero_dn <- function() {
        isTRUE(orthologr_option("orthologr.legacy_kaks_zero_dn", FALSE))
}

# dN of gene pairs for which KaKs_Calculator printed 'NA' although the value is known.
#
# KaKs_Calculator (1.2 and 2.0) prints "NA" in place of every Ka or Ks below 1e-6, and prints Ka/Ks whenever Ks is
# at least 1e-6 and Ka was computed. A pair without nonsynonymous change therefore comes back as
# Ka = NA, Ks = <value>, Ka/Ks = 0 (or -0): the estimate is dN = 0, and only its printed form is 'NA'. (Found for
# NG, LWL, LPB, MLWL, MLPB, YN; the maximum-likelihood methods print small positive numbers instead.)
#   Ka NA, Ks a number > 0, Ka/Ks a number with |Ka/Ks * Ks| < 1e-6   ->  dN = Ka/Ks * Ks (0 when Ka/Ks is 0)
# Everything else is left as it is: when Ka itself could not be computed KaKs_Calculator prints NA (or nan) for
# Ka/Ks as well - identical sequences (NA NA NA), a failed model (nan nan nan), Ks below 1e-6 (Ka/Ks NA).
kaks_zero_dn <- function(dN, dS, dNdS) {
        known <- is.na(dN) & is.finite(dS) & dS > 0 & is.finite(dNdS)
        known[known] <- abs(dNdS[known] * dS[known]) < 1e-6
        if (any(known)) {
                dN[known] <- pmax(0, dNdS[known] * dS[known])
                dNdS[known] <- pmax(0, dNdS[known])          # '-0' is 0
        }
        list(dN = dN, dNdS = dNdS, changed = known)
}

# KaKs_Calculator's genetic code: the value of '-c' in kaks_calc.params (1, the standard code, when not given; NA
# when it cannot be read).
kaks_genetic_code <- function(kaks_calc.params = NULL) {
        p <- paste(unlist(kaks_calc.params), collapse = " ")
        m <- regmatches(p, regexec("(^|\\s)-c\\s*([0-9]+)(\\s|$)", p))[[1]]
        if (length(m) == 0)
                return(if (grepl("(^|\\s)-c", p)) NA_integer_ else 1L)
        as.integer(m[3])
}

# The genetic codes of KaKs_Calculator 2.0 (base.cpp, 'transl_table'; codons in the order TTT, TTC, TTA, TTG, TCT,
# ... with T, C, A, G = 0, 1, 2, 3; '!' a stop codon). NULL for a code it does not have.
kaks_transl_table <- function(genetic_code) {
        tables <- c(
                "1"  = "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "2"  = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG",
                "3"  = "FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "4"  = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "5"  = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                "6"  = "FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "9"  = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "10" = "FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "11" = "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "12" = "FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "13" = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                "14" = "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "15" = "FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "16" = "FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "21" = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "22" = "FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "23" = "FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "24" = "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                "25" = "FFLLSSSSYY!!CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "26" = "FFLLSSSSYY!!CC!WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "27" = "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "28" = "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "29" = "FFLLSSSSYYYYCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "30" = "FFLLSSSSYYEECC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "31" = "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "33" = "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
        )
        if (length(genetic_code) != 1 || is.na(genetic_code) || !is.element(as.character(genetic_code), names(tables)))
                return(NULL)
        strsplit(tables[[as.character(genetic_code)]], "")[[1]]
}

# The per-codon tables of the LWL family for one genetic code and one rule of transitions (LWL85 for LWL and LPB;
# Tzeng et al. 2004 for MLWL and MLPB), as KaKs_Calculator 2.0 computes them (LWL85.cpp, LPB93.cpp):
#   L  : 64 x 3, the number of positions of a codon that are 0-, 2- and 4-fold degenerate
#   S, V : 4096 x 3, what a pair of codons (row 64 * i + j + 1 for codons i, j) adds to the transitions (S) and
#        transversions (V) of the three classes: the differences of all evolutionary pathways that do not pass a
#        stop codon, 0.5 to the class of either codon at the position, divided by the number of these pathways
# Computed once per code and rule.
kaks_lwl_tables <- function(genetic_code, modified) {
        key <- paste0("lwl_", genetic_code, "_", modified)
        if (!is.null(kaks_lwl_cache[[key]]))
                return(kaks_lwl_cache[[key]])
        aa <- kaks_transl_table(genetic_code)
        if (is.null(aa))
                return(NULL)
        # codon i (1..64) as bases 0..3 (T, C, A, G) and back
        base <- cbind((0:63) %/% 16L, ((0:63) %/% 4L) %% 4L, (0:63) %% 4L)
        id <- function(b) b[1] * 16L + b[2] * 4L + b[3] + 1L
        named <- function(x) id(match(strsplit(x, "")[[1]], c("T", "C", "A", "G")) - 1L)
        # the class (0, 2, 4) of every position of every codon: how many of the other three bases give the same
        # amino acid (a stop codon is never the same)
        cls <- matrix(0L, 64, 3)
        for (i in 1:64)
                for (pos in 1:3) {
                        n <- 0L
                        for (b in setdiff(0:3, base[i, pos])) {
                                t <- base[i, ]
                                t[pos] <- b
                                a <- aa[id(t)]
                                if (a != "!" && a == aa[i])
                                        n <- n + 1L
                        }
                        cls[i, pos] <- if (n == 3L) 4L else if (n > 0L) 2L else 0L
                }
        L <- cbind(rowSums(cls == 0L), rowSums(cls == 2L), rowSums(cls == 4L))
        slot <- cls %/% 2L + 1L
        ATA <- named("ATA"); ATG <- named("ATG"); ATC <- named("ATC"); ATT <- named("ATT")
        arg_pairs <- rbind(c(named("CGA"), named("AGA")), c(named("CGG"), named("AGG")))
        is_arg <- function(c1, c2) any((arg_pairs[, 1] == c1 & arg_pairs[, 2] == c2) | (arg_pairs[, 1] == c2 & arg_pairs[, 2] == c1))
        # one difference at position 'pos' between codons c1 and c2: what it adds to s / v (class slots 1..3)
        tt <- function(c1, c2, pos, sv) {
                k1 <- slot[c1, pos]
                k2 <- slot[c2, pos]
                syn <- -1L
                if (modified && pos == 3L) {
                        if ((c1 == ATA && c2 == ATG) || (c1 == ATG && c2 == ATA))
                                syn <- 0L
                        if ((c1 == ATA && (c2 == ATC || c2 == ATT)) || ((c1 == ATC || c1 == ATT) && c2 == ATA))
                                syn <- 1L
                }
                if (pos == 1L && is_arg(c1, c2)) {
                        if (!modified) {
                                # LWL85: half a transition to the class of the first codon, half a transversion to the second
                                sv[k1] <- sv[k1] + 0.5
                                sv[3L + k2] <- sv[3L + k2] + 0.5
                                return(sv)
                        }
                        syn <- 1L
                }
                if (syn == -1L) {
                        sm <- base[c1, pos] + base[c2, pos]
                        syn <- if (sm == 5L || sm == 1L) 1L else 0L
                }
                off <- if (syn == 1L) 0L else 3L
                sv[off + k1] <- sv[off + k1] + 0.5
                sv[off + k2] <- sv[off + k2] + 0.5
                sv
        }
        perms <- list(list(1L), list(1:2, 2:1), list(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L), c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)))
        SV <- matrix(0, 4096, 6)
        for (i in which(aa != "!")) {
                for (j in which(aa != "!")) {
                        if (i == j)
                                next
                        diff <- which(base[i, ] != base[j, ])
                        num <- length(diff)
                        sv <- numeric(6)
                        stop <- 0L
                        for (pp in perms[[num]]) {
                                path <- diff[pp]
                                cur <- base[i, ]
                                steps <- i
                                for (p in path[-num]) {
                                        cur[p] <- base[j, p]
                                        steps <- c(steps, id(cur))
                                }
                                if (any(aa[steps[-1]] == "!")) {
                                        stop <- stop + 1L
                                        next
                                }
                                steps <- c(steps, j)
                                for (h in seq_len(num))
                                        sv <- tt(steps[h], steps[h + 1L], path[h], sv)
                        }
                        n_path <- length(perms[[num]]) - stop
                        if (n_path > 0)
                                SV[64L * (i - 1L) + j, ] <- sv / n_path
                }
        }
        out <- list(aa = aa, L = L, S = SV[, 1:3], V = SV[, 4:6])
        assign(key, out, envir = kaks_lwl_cache)
        out
}

kaks_lwl_cache <- new.env(parent = emptyenv())

# The sites (L), transitions (Si) and transversions (Vi) of the 0-, 2- and 4-fold classes of a codon alignment as
# KaKs_Calculator 2.0 computes them for 'method' (LWL, MLWL, LPB, MLPB), and P = Si / L, Q = Vi / L. Codons with a
# gap, a character other than A, C, G, T or a stop codon in either sequence are left out, as by KaKs_Calculator;
# the sums are built codon by codon in the order of the alignment, so that the numbers are KaKs_Calculator's to
# the last bit. NULL when the alignment or the genetic code cannot be handled here (e.g. 'U' in the sequences).
kaks_lwl_counts <- function(seq1, seq2, method, genetic_code = 1L) {
        if (!is.element(method, kaks_lwl_methods()) || length(seq1) != 1 || length(seq2) != 1 || is.na(seq1) || is.na(seq2))
                return(NULL)
        tab <- kaks_lwl_tables(genetic_code, modified = is.element(method, c("MLWL", "MLPB")))
        if (is.null(tab))
                return(NULL)
        s1 <- toupper(seq1)
        s2 <- toupper(seq2)
        if (nchar(s1) != nchar(s2) || nchar(s1) %% 3 != 0 || grepl("U", s1, fixed = TRUE) || grepl("U", s2, fixed = TRUE))
                return(NULL)
        n <- nchar(s1) %/% 3
        if (n == 0)
                return(NULL)
        code <- function(s) {
                b <- match(strsplit(s, "")[[1]], c("T", "C", "A", "G")) - 1L
                m <- matrix(b, nrow = 3)
                m[1, ] * 16L + m[2, ] * 4L + m[3, ] + 1L
        }
        i1 <- code(s1)
        i2 <- code(s2)
        keep <- !is.na(i1) & !is.na(i2)
        keep[keep] <- tab$aa[i1[keep]] != "!" & tab$aa[i2[keep]] != "!"
        i1 <- i1[keep]
        i2 <- i2[keep]
        L <- (colSums(tab$L[i1, , drop = FALSE]) + colSums(tab$L[i2, , drop = FALSE])) / 2
        row <- 64L * (i1 - 1L) + i2
        # sums in the order of the alignment, in double precision (sum() and cumsum() add in extended precision)
        add <- function(x) Reduce(`+`, x, 0)
        Si <- c(add(tab$S[row, 1]), add(tab$S[row, 2]), add(tab$S[row, 3]))
        Vi <- c(add(tab$V[row, 1]), add(tab$V[row, 2]), add(tab$V[row, 3]))
        list(L = L, Si = Si, Vi = Vi, P = Si / L, Q = Vi / L)
}

# Do the counts of kaks_lwl_counts() agree with those KaKs_Calculator printed (6 significant digits)?
kaks_counts_agree <- function(exact, L, Si, Vi) {
        close <- function(a, b) all(is.finite(a) & is.finite(b) & abs(a - b) <= 1e-5 * pmax(1, abs(a)))
        close(exact$L, L) && close(exact$Si, Si) && close(exact$Vi, Vi)
}

# The header line of KaKs_Calculator 2.0's output (all methods).
kaks_calc_header <- function() {
        paste(c("Sequence", "Method", "Ka", "Ks", "Ka/Ks", "P-Value(Fisher)", "Length", "S-Sites", "N-Sites",
                "Fold-Sites(0:2:4)", "Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)",
                "Fold-N-Substitutions(0:2:4)", "Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)",
                "GC(1:2:3)", "ML-Score", "AICc", "Akaike-Weight", "Model"), collapse = "\t")
}

# The lines of an output file of KaKs_Calculator with exactly one header line, the first. KaKs_Calculator 2.0
# (bioconda 2.0.1) decides whether it writes the header for the gamma methods (GMYN; also in "ALL") with a class
# member it never initialises ('tempt'), so the line can be missing - read.csv() then took the first row of values
# for the header and the pair came back without a result - or be written again further down. When the first line
# is not a header, the known header is put first; later lines that are the header are dropped (AUDIT-OLP2 M1).
# Output that has the header once, first (every other method), is returned as it is.
# A line is the header when its second field is "Method" (a row has the name of a method there): a row of a pair
# whose name is "Sequence" - the first field of the header - is a row (it was taken for the header, and the pair
# came back without a result).
kaks_output_lines <- function(lines) {
        if (length(lines) == 0)
                return(lines)
        header <- kaks_calc_header()
        if (!startsWith(lines[1], "Sequence\tMethod\t"))
                lines <- c(header, lines)
        rest <- lines[-1]
        drop <- rest == header | rest == lines[1]
        if (any(drop))
                lines <- c(lines[1], rest[!drop])
        lines
}

# KaKs_Calculator for the codon alignments of several gene pairs with one process (compute_dnds()): one AXT file
# with a record per pair, named 'subst_name', and one output table. Returns a list with, for every pair, the table
# substitutionrate(file = files[k], est.method = est.method, quiet = TRUE, kaks_calc_path = kaks_calc_path,
# subst_name = subst_name[k]) returns - or NULL where this function cannot be sure of it; the caller then calls
# substitutionrate() for that pair, which also gives its error.
# KaKs_Calculator (2.0) computes every record of an AXT file on its own: the objects of the methods are made anew for
# each record, and its output is the header line and then the rows of each record, which start with the name of
# the record. The rows of a pair are read as the table of a call for that pair alone (the header line and these
# rows, read by read.csv() as substitutionrate() reads its file) and go through kaks_result_table() as there.
# Every pair is NULL when the program fails or its rows do not belong to the records one to one and in their order
# (a record that KaKs_Calculator rejects is left out of its output), when a record has another number of rows than
# the method gives (1; 17 for "ALL"), or when a line has another number of fields than the header (AUDIT-OLP2 L1);
# a pair whose alignment does not have two
# sequences of letters, '*' and '-' (an empty sequence or a line break would shift the records) is not handed over.
kaks_calc_batch <- function(files, subst_name, est.method, kaks_calc_path = NULL) {
        n <- length(files)
        out <- vector("list", n)
        if (n == 0)
                return(out)
        tryCatch({
                if (!is.null(kaks_calc_path))
                        kaks_calc_path <- normalizePath(kaks_calc_path, mustWork = FALSE)
                if (is.null(kaks_calc_path))
                        is_installed_kaks_calculator()
                calc_dir <- file.path(orthologr_tmpdir(), "_calculation")
                if (!file.exists(calc_dir))
                        dir.create(calc_dir, recursive = TRUE)
                seqs <- vector("list", n)
                for (k in seq_len(n)) {
                        s <- tryCatch(seqinr::read.fasta(files[k], as.string = TRUE, set.attributes = FALSE,
                                                         forceDNAtolower = FALSE), error = function(e) NULL)
                        if (length(s) != 2)
                                next
                        s <- c(as.character(s[[1]]), as.character(s[[2]]))
                        if (!all(grepl("^[A-Za-z*-]+\\z", s, perl = TRUE)) || !grepl("^[A-Za-z0-9_.-]+\\z", subst_name[k], perl = TRUE))
                                next
                        seqs[[k]] <- s
                }
                use <- which(!vapply(seqs, is.null, TRUE))
                if (length(use) == 0 || anyDuplicated(subst_name[use]))
                        return(out)
                axt_file <- tempfile("batch_", tmpdir = calc_dir, fileext = ".axt")
                kaks_file <- paste0(axt_file, ".kaks")
                on.exit(unlink(c(axt_file, kaks_file)), add = TRUE)
                writeLines(as.vector(rbind(subst_name[use], vapply(seqs[use], `[`, "", 1), vapply(seqs[use], `[`, "", 2), "")),
                           axt_file)
                run_tool(
                        "KaKs_Calculator",
                        args    = c(list("-i", axt_file, "-o", kaks_file), kaks_calc_args(est.method, NULL)),
                        path    = kaks_calc_path,
                        tool    = "KaKs_Calculator",
                        context = paste0(length(use), " gene pairs"),
                        output  = kaks_file
                )
                lines <- kaks_output_lines(readLines(kaks_file, warn = FALSE))
                if (length(lines) < 2)
                        return(out)
                key <- sub("\t.*$", "", lines[-1])
                blocks <- rle(key)
                # A record KaKs_Calculator writes no row for (GMYN of a single codon: only the header line, which
                # kaks_output_lines() dropped) is left out of its output: the other records are still those of the
                # chunk, in their order and with names that occur once, and are read; the pair without rows gets
                # NULL (the caller computes it on its own, which gives its error). Before, one such pair sent the
                # whole chunk to one KaKs_Calculator process per pair.
                if (!identical(blocks$values, subst_name[use])) {
                        at <- match(blocks$values, subst_name[use])
                        if (anyNA(at) || anyDuplicated(blocks$values) || is.unsorted(at, strictly = TRUE))
                                return(out)
                        use <- use[at]
                }
                rows_per_pair <- if (identical(est.method, "ALL")) 17L else 1L
                if (any(blocks$lengths != rows_per_pair))
                        return(out)
                n_tabs <- lengths(regmatches(lines, gregexpr("\t", lines, fixed = TRUE)))
                if (any(n_tabs != n_tabs[1]))
                        return(out)
                ends <- cumsum(blocks$lengths)
                starts <- ends - blocks$lengths + 1L
                genetic_code <- kaks_genetic_code(NULL)
                # The table is read once when there is no quote character (every line has the fields of the
                # header, checked above); the rows of a pair are then a table with the values of read.csv() for its
                # rows alone:
                # kaks_result_table() takes every column as character or numbers (as.numeric(as.character())),
                # so a column that read.csv() types differently for all rows than for the rows of one pair (e.g.
                # numbers instead of logical NA, double instead of integer) gives the same values. Otherwise the rows
                # of each pair are read on their own.
                whole <- if (!any(grepl("\"", lines, fixed = TRUE)))
                        utils::read.csv(text = lines, sep = "\t", header = TRUE) else NULL
                if (!is.null(whole) && nrow(whole) != length(lines) - 1L)
                        whole <- NULL
                # The values of all rows at once (kaks_result_values() works row by row: the values of a pair are
                # those of its rows alone), handed back as a "dnds_values" record per pair, which compute_dnds()
                # turns into rows of the result together with the other pairs of the chunk (one table per chunk
                # instead of one per pair). Should that fail, every pair is read on its own as below.
                if (!is.null(whole)) {
                        vals <- tryCatch({
                                pair_of_row <- rep(seq_along(use), blocks$lengths)
                                # (the gene pair named in messages, read from its file only for a row that needs it)
                                ctx <- function(rows) vapply(use[pair_of_row[rows]], function(k)
                                        paste0(pair_context(files[k], subst_name[k]), " (method ", est.method, ")"), "")
                                v <- kaks_result_values(whole, context = ctx, aln = seqs[use][pair_of_row],
                                                        genetic_code = genetic_code)
                                lapply(seq_along(use), function(j) {
                                        r <- starts[j]:ends[j]
                                        dnds_values(dN = v$dN[r], dS = v$dS[r], dNdS = v$dNdS[r], method = v$method[r],
                                                    reason = kaks_failed_reason(v$method[r], v$failed[r], v$reason[r]))
                                })
                        }, error = function(e) NULL)
                        if (!is.null(vals)) {
                                out[use] <- vals
                                return(out)
                        }
                }
                for (j in seq_along(use)) {
                        k <- use[j]
                        out[k] <- list(tryCatch({
                                kaks_tbl <- if (!is.null(whole)) whole[starts[j]:ends[j], , drop = FALSE] else
                                        utils::read.csv(text = c(lines[1], lines[-1][starts[j]:ends[j]]), sep = "\t", header = TRUE)
                                context <- pair_context(files[k], subst_name[k])
                                kaks_result_table(kaks_tbl, fasta_ids(files[k])[1:2], paste0(context, " (method ", est.method, ")"),
                                                  aln = seqs[[k]], genetic_code = genetic_code)
                        }, error = function(e) NULL))
                }
                out
        }, error = function(e) vector("list", n))
}

# The reason of each row of one gene pair as compute_dnds() gives it from attr(, "kaks_failed") of
# kaks_result_table(): the reason of the first failed row with the row's method, NA where there is none.
#   method, failed, reason : per row of the pair (kaks_result_values())
kaks_failed_reason <- function(method, failed, reason) {
        if (!any(failed))
                return(rep(NA_character_, length(method)))
        reason[failed][match(method, method[failed])]
}

# The values of one gene pair that compute_dnds() puts into its result table together with those of other pairs
# (rows of the same columns as the table of substitutionrate() plus '.reason'): numeric dN, dS, dNdS (one element
# per row), 'method' (NULL for Comeron and Li, which have no such column) and 'reason' (NA: no failure).
dnds_values <- function(dN, dS, dNdS, method = NULL, reason = rep(NA_character_, length(dN))) {
        structure(list(dN = dN, dS = dS, dNdS = dNdS, method = method, reason = reason), class = "dnds_values")
}
