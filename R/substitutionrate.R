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
#' @param kaks_calc.params a character string storing additional parameters for KaKs_Calculator. Default is \code{NULL}.
#' The parameters are added to the method chosen with \code{est.method} (e.g. \code{est.method = "YN"},
#' \code{kaks_calc.params = "-c 11"} runs \code{-m YN -c 11}; up to orthologr 0.4.3 they replaced the method).
#' With \code{est.method = "kaks_calc"} the methods themselves are given here, e.g.
#' \code{kaks_calc.params} = "-m NG -m YN". 
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
#' @param subst_name a character string specifying the substitution name that shall be added to the internal folder path naming. 
#' Default is \code{subst_name} = \code{NULL}.
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
                             subst_name       = NULL){
        
        # dNdS estimation methods provided by the KaKs_Calculator program
        if (!is.null(kaks_calc_path))
                kaks_calc_path <- normalizePath(kaks_calc_path, mustWork = FALSE)
        
        kaks_calc_methods <- kaks_calculator_methods()
        
        if (!is.dnds_est_method(est.method))
                stop("Please choose a dNdS estimation method that is supported by this function.", call. = FALSE)
        
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
                        gestimator(file     = file,
                                   file_out = file.out)
                        
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
                        print("Substitutionrate successfully calculated using Li's method.")
                }
                
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
                
                return(kaks_result_table(kaks_tbl, pair_ids, paste0(context, " (method ", est.method, ")")))
        }
        
        # every method accepted above is handled by one of the branches; a method that is not must never give NULL
        stop("est.method = \"", est.method, "\" is not handled by substitutionrate().", call. = FALSE)
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
        kaks_tbl_res <- data.table::as.data.table(kaks_tbl_res)
        data.table::setkey(kaks_tbl_res, query_id)
        return(kaks_tbl_res)
}
