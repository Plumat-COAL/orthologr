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
#' With a COAL that computes them (\pkg{coalr} >= 0.5.0 listing them in \code{coalr::coal_info()$methods}),
#' \code{kaks_tool = "coal"} also computes \code{est.method = "Comeron"} and \code{"Li"}: COAL's ports of gestimator
#' and of \code{seqinr::kaks()}, with the saturation rule of this page, rounded to the 4 significant digits
#' gestimator writes (Comeron) - the same table as the built-in code gives for the alignment. With an older COAL
#' these two methods are computed by the built-in code.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function takes a pairwise alignments file as input and estimates dNdS ratios
#' of the corresponding alignments using predefined dNdS estimation methods.
#' 
#' The dNdS estimation methods available in this function are:
#' 
#' \itemize{
#' \item "Li" : Li's method (1993) -> provided by the ape package
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
#' @examples \dontrun{
#' 
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
#'                   
#'                                           
#' }
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
#' KaKs_Calculator prints \code{NA} for every Ka below 1e-6, also when the estimate is exactly 0, and still prints
#' \code{Ka/Ks = 0}. Up to orthologr 0.4.3 such a pair was returned with \code{dN = NA} and \code{dNdS = 0}, and
#' \code{\link{filter_dNdS}} / \code{\link{divergence_stratigraphy}} then removed it - the most conserved genes of
#' a comparison. Since 0.5.0 \code{dN = 0} is returned when Ka is \code{NA}, Ks is a positive number and Ka/Ks is
#' 0. When Ka/Ks is \code{NA} as well (identical sequences, a model that failed) \code{dN} stays \code{NA}.
#' \code{options(orthologr.legacy_kaks_zero_dn = TRUE)} restores the old table.
#'
#' @section Saturation (Comeron's method):
#' The corrections of Comeron (1995) take logarithms of \code{1 - 2P - Q} and \code{1 - 2Q} (\code{P}, \code{Q}:
#' proportions of transitional and transversional differences in a class of sites). When such an argument is not
#' positive, the observed differences exceed what the model can correct for and the distance is not defined. In
#' that case \code{dS} (synonymous classes) or \code{dN} (nonsynonymous classes) and \code{dNdS} are returned as
#' \code{NA}. Up to orthologr 0.4.3 the logarithm was replaced by 0 and a finite, too small value was returned
#' (e.g. \code{dS = 0} for \code{ATGGCTAAA} / \code{ATGGCCAAA}, which differ by one synonymous substitution);
#' \code{options(orthologr.legacy_comeron_saturation = TRUE)} reproduces those values. Values of gene pairs
#' without such a term are the same in both modes.
#'
#' A class of sites that has no sites (weight 0 in the estimate) never makes an estimate \code{NA}; this is also
#' the rule of COAL's Comeron method, so \code{engine = "coal"} gives the same \code{NA} pattern.
#' \strong{Short alignments:} with few codons a class with very few sites (typically the twofold sites whose
#' synonymous change is a transversion) is easily "saturated" by one or two differences, and \code{dN} is then
#' \code{NA} although the legacy value was a usable estimate: in a simulation of random sequences without
#' selection, 4 to 13\% of the pairs of 30 codons lost \code{dN} at 10 to 35\% differing sites, 3 of 600 pairs of
#' 100 codons, none of 300 codons or more. For short sequences consider the legacy option, or a method with
#' pooled site classes (\code{"NG"}, \code{"YN"}).
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
        
        kaks_calc_methods <-
                c("MA",
                  "MS",
                  "NG",
                  "LWL",
                  "LPB",
                  "MLWL",
                  "YN",
                  "MYN",
                  "GY",
                  "GMYN",
                  "kaks_calc",
                  "ALL")
        
        if (!is.dnds_est_method(est.method))
                stop("Please choose a dNdS estimation method that is supported by this function.", call. = FALSE)
        
        if (!is.element(kaks_tool, c("kaks_calculator", "coal")))
                stop("Please choose kaks_tool = \"kaks_calculator\" or kaks_tool = \"coal\".", call. = FALSE)
        
        if (kaks_tool == "coal" && !is.element(est.method, c("Comeron", "Li", coal_methods())))
                stop("kaks_tool = \"coal\" computes the methods ", paste(coal_methods(), collapse = ", "),
                     "; est.method = \"", est.method, "\" needs kaks_tool = \"kaks_calculator\".", call. = FALSE)
        
        if (!is.element(format, c("mase", "clustal", "phylip", "fasta" , "msf")))
                stop("Please choose a format that is supported by seqinr::read.alignment.", call. = FALSE)
        
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
                # gestimator writes its table to a file that only this call uses and that is removed afterwards
                file.out <- tempfile(
                        paste0(if (is.null(subst_name)) "" else paste0(subst_name, "_"), "gestimout_"),
                        tmpdir = file.path(orthologr_tmpdir(), "_calculation")
                )
                on.exit(unlink(file.out), add = TRUE)
                
                # file in fasta required
                if (format != "fasta")
                        stop("To use gestimator an alignment file in fasta format is required.")
                
                tryCatch({
                        # To use gestimator a file in fasta format is required
                        
                        # include this line instead of the following, to use internal gestimator
                        # the alignment is checked here as well as in the C++ code: at least two sequences of the
                        # same length, which is a multiple of 3 and not 0
                        check_codon_alignment_file(file)
                        gestimator(file     = file,
                                   file_out = file.out,
                                   legacy_saturation = comeron_legacy_saturation())
                        
                        #system(paste0("gestimator -i ",file," -o ",file.out))
                        
                        # The columns are read as what they are: without colClasses fread() guesses the types,
                        # and IDs that look like numbers ("1001", "007", "1e5", "1.10") or like logical values
                        # ("TRUE", "NA") came back as another ID or as NA.
                        hit.table <- read_gestimator_table(file.out)
                        
                        data.table::setkey(hit.table, query_id)
                        
                        # for consistent output format we set the default value of gestimator 999
                        # to NA for output
                        # Therefor the user cannot make a mistake if using the dNdS results without
                        # filtering for some value.
                        hit.table[which(hit.table[, dS == 999]), dS := NA]
                        hit.table[which(hit.table[, dN == 999]), dN := NA]
                        hit.table[which(hit.table[, dNdS == 999]), dNdS := NA]
                        
                        if (!quiet) {
                                message("Substitutionrate successfully calculated by gestimator.")
                        }
                        
                        return(hit.table)
                }, error = function(e) {
                        stop(
                                "dN/dS estimation with Comeron's method (gestimator) failed for ",
                                pair_context(file, subst_name), ": ", conditionMessage(e),
                                call. = FALSE
                        )
                })
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
        if (kaks_tool == "coal" && is.element(est.method, coal_methods())) {
                if (format != "fasta")
                        stop("kaks_tool = \"coal\" needs the codon alignment in fasta format.", call. = FALSE)
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
                        utils::read.csv(kaks_file, sep = "\t", header = TRUE)
                }, error = function(e) {
                        stop("The output of KaKs_Calculator for ", context, " ('", kaks_file, "') could not be read: ",
                             conditionMessage(e), call. = FALSE)
                })
                
                return(kaks_result_table(kaks_tbl, pair_ids, context))
        }
}

# options(orthologr.legacy_comeron_saturation = TRUE): Comeron's method as up to 0.4.3, where an undefined
# logarithm of the correction (saturation) was replaced by 0 and a finite dS or dN was returned.
comeron_legacy_saturation <- function() {
        isTRUE(orthologr_option("orthologr.legacy_comeron_saturation", FALSE))
}

# Stop unless `file` (FASTA) is a codon alignment: at least two sequences, all of the same length, the length a
# multiple of 3 and at least 3. Line ends may be LF or CRLF.
check_codon_alignment_file <- function(file) {
        if (!file.exists(file))
                stop("The alignment file '", file, "' does not exist.", call. = FALSE)
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
read_gestimator_table <- function(file) {
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
kaks_result_table <- function(kaks_tbl, pair_ids, context = "the gene pair") {
        query_id <- NULL
        if (nrow(kaks_tbl) == 0 || ncol(kaks_tbl) < 5)
                stop("KaKs_Calculator returned no result for ", context, ".", call. = FALSE)
        num <- function(x) suppressWarnings(as.numeric(as.character(x)))
        kaks_tbl_res <- data.frame(
                query_id   = rep(as.character(pair_ids[1]), nrow(kaks_tbl)),
                subject_id = rep(as.character(pair_ids[2]), nrow(kaks_tbl)),
                dN         = num(kaks_tbl[ , 3]),
                dS         = num(kaks_tbl[ , 4]),
                dNdS       = num(kaks_tbl[ , 5]),
                method     = kaks_tbl[ , 2],
                stringsAsFactors = FALSE
        )
        if (!kaks_legacy_zero_dn()) {
                fixed <- kaks_zero_dn(kaks_tbl_res$dN, kaks_tbl_res$dS, kaks_tbl_res$dNdS)
                kaks_tbl_res$dN <- fixed$dN
                kaks_tbl_res$dNdS <- fixed$dNdS
        }
        kaks_tbl_res <- data.table::as.data.table(kaks_tbl_res)
        data.table::setkey(kaks_tbl_res, query_id)
        return(kaks_tbl_res)
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
