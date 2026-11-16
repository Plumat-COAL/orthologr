# Internal checks of the arguments of dNdS(): what is not supported is rejected before any work starts, with a
# message that names the argument, the value and the valid values.

# The DIAMOND sensitivity modes. "default" is DIAMOND's own default mode (no sensitivity option on the command line).
diamond_sensitivity_modes <- function() {
        c("fast", "faster", "default", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive")
}

check_sensitivity_mode <- function(sensitivity_mode) {
        if (!is.character(sensitivity_mode) || length(sensitivity_mode) != 1 ||
            !is.element(sensitivity_mode, diamond_sensitivity_modes()))
                orthologr_abort(class = "bad_argument", "sensitivity_mode = '", paste(sensitivity_mode, collapse = "', '"), "' is not a DIAMOND2 sensitivity mode. ",
                     "Please choose one of: '", paste(diamond_sensitivity_modes(), collapse = "', '"), "'.", call. = FALSE)
        invisible(sensitivity_mode)
}

# The DIAMOND sensitivity mode: the argument when given, else options(orthologr.sensitivity_mode = ), else "fast"
# with options(orthologr.legacy = TRUE) (the default of 0.4.3), else "default" (DIAMOND's own default mode).
orthologr_sensitivity_mode <- function(mode = NULL) {
        if (is.null(mode))
                mode <- orthologr_option("orthologr.sensitivity_mode", "default")
        check_sensitivity_mode(mode)
        mode
}

# One value out of a fixed set, with an error in the style of match.arg() that also names the argument.
check_choice <- function(value, choices, arg, hint = NULL) {
        if (!is.character(value) || length(value) != 1 || is.na(value) || !is.element(value, choices))
                orthologr_abort(class = "bad_argument", arg, " = ", choice_text(value), " is not supported. Please choose ",
                     if (length(choices) > 1) "one of: " else "",
                     paste0("'", choices, "'", collapse = ", "), ".",
                     if (!is.null(hint)) paste0(" ", hint) else "", call. = FALSE)
        invisible(value)
}

# A value as the messages of check_choice() show it: one text in single quotes ('coal'), anything else as R code
# (NULL, NA, c("coal", "orthologr"), 1).
choice_text <- function(value) {
        if (is.character(value) && length(value) == 1 && !is.na(value))
                return(paste0("'", value, "'"))
        paste(deparse(value, width.cutoff = 500L, nlines = 1L), collapse = "")
}

dnds_check_args <- function(aligner, sensitivity_mode, seq_type, format, ortho_detection, aa_aln_type, aa_aln_tool,
                            aa_aln_params = NULL) {
        
        check_choice(aligner, c("diamond", "blast"), "aligner")
        
        if (identical(ortho_detection, "Orthofinder2"))
                orthologr_abort(class = "bad_argument", "ortho_detection = 'Orthofinder2' is not available in dNdS(): OrthoFinder is not run by this function. ",
                     "Please choose 'RBH' (reciprocal best hit, default) or 'BH' (best hit).", call. = FALSE)
        if (is.element(ortho_detection, c("PO", "DELTA")))
                orthologr_abort(class = "bad_argument", "ortho_detection = '", ortho_detection, "' (", if (ortho_detection == "PO") "ProteinOrtho" else "DELTA-BLAST",
                     ") is no longer available in dNdS(). ",
                     "Please choose 'RBH' (reciprocal best hit, default) or 'BH' (best hit).", call. = FALSE)
        check_choice(ortho_detection, c("RBH", "BH"), "ortho_detection")
        
        if (identical(seq_type, "protein"))
                orthologr_abort(class = "bad_argument", "seq_type = 'protein' is not possible in dNdS(): dN/dS is computed from codon alignments, ",
                     "so both files must hold coding sequences (seq_type = 'cds').", call. = FALSE)
        check_choice(seq_type, "cds", "seq_type")
        
        check_choice(format, "fasta", "format", hint = "dNdS() reads coding sequences from FASTA files.")
        
        if (aligner == "diamond")
                check_sensitivity_mode(sensitivity_mode)
        
        check_choice(aa_aln_type, c("pairwise", "multiple"), "aa_aln_type")
        
        # options(orthologr.max_aln_cells) is checked now: an invalid value (-1, "5e8") ran the search and then
        # failed every gene pair with the message about the option; the same for the options of Comeron's method and
        # options(orthologr.ambiguous_cds)
        max_aln_cells()
        comeron_twosubs()
        comeron_min_sites()
        ambiguous_cds_mode()
        
        if (aa_aln_type == "pairwise") {
                check_choice(aa_aln_tool, "NW", "aa_aln_tool",
                             hint = "(with aa_aln_type = 'pairwise'; the tools 'clustalw', 'clustalo', 'mafft', 'muscle' and 't_coffee' need aa_aln_type = 'multiple').")
                if (!is.null(aa_aln_params))
                        warning("aa_aln_params is not used with aa_aln_type = 'pairwise' (aa_aln_tool = 'NW' has no command line).",
                                call. = FALSE)
        } else {
                check_choice(aa_aln_tool, c("clustalw", "clustalo", "mafft", "muscle", "t_coffee"), "aa_aln_tool",
                             hint = "(with aa_aln_type = 'multiple'; 'NW' needs aa_aln_type = 'pairwise').")
        }
        
        invisible(TRUE)
}

# The e-value threshold of a search: one number >= 0, given as a number or as text ("1E-5"). Anything else (e.g.
# "1e-5; echo hi") is rejected; it used to be handed to the program (DIAMOND got 'NA'). The value is returned as
# it was given, so that messages and file names do not change.
check_evalue <- function(eval) {
        value <- suppressWarnings(as.numeric(eval))
        if (length(eval) != 1 || is.na(value) || value < 0)
                orthologr_abort(class = "bad_argument", "eval = '", paste(eval, collapse = "', '"), "' is not an e-value threshold. Please give one number >= 0, e.g. eval = \"1E-5\".",
                     call. = FALSE)
        eval
}

# ---- the dN/dS methods: one table ------------------------------------------------------------------------------

# Every dN/dS method orthologr knows, one entry each; every list of methods, every check of a method name and every
# message that names methods is derived from this table (they were seven separate lists).
#   kind    : "own"    - computed by orthologr's own code whatever kaks_tool is (Comeron: gestimator; Li:
#                        seqinr::kaks())
#             "kaks"   - a method of KaKs_Calculator (engine = "orthologr", kaks_tool = "kaks_calculator")
#             "pseudo" - not a method but a choice of KaKs_Calculator: "ALL" (every method) and "kaks_calc" (the
#                        methods are in kaks_calc.params; substitutionrate() only)
#   coal    : whether COAL computes it: "always" (every COAL), "native" (a COAL whose coalr::coal_info()$methods
#             lists "Comeron" and "Li", coalr >= 0.5.0), "gy" (a COAL with Goldman & Yang's maximum likelihood
#             estimate, coalr::coal_info()$features "gy94_ml", coalr >= 0.6.0), NA (never)
#   family  : the LWL family (Li, Wu and Luo 1985; Li 1993; Tzeng et al. 2004) or Yang and Nielsen (2000), for
#             reading KaKs_Calculator's output (kaks_saturated(), kaks_zero_ds())
#   counted : KaKs_Calculator's columns "S-Substitutions" and "N-Substitutions" count the synonymous and
#             nonsynonymous differences (NA for a count of 0). GMYN counts as MYN does; GY, MA, MS print
#             expected numbers.
dnds_method_registry <- function() {
        # (a list of columns, not a data.frame: the views below are called for every gene pair)
        list(
                method  = c("Comeron", "Li", "NG",     "LWL",    "MLWL",   "LPB",    "MLPB",   "YN",     "MYN",
                            "GY",   "MA",   "MS",   "GMYN", "kaks_calc", "ALL"),
                kind    = c("own",     "own", "kaks",  "kaks",   "kaks",   "kaks",   "kaks",   "kaks",   "kaks",
                            "kaks", "kaks", "kaks", "kaks", "pseudo",    "pseudo"),
                coal    = c("native",  "native", "always", "always", "always", "always", "always", "always", "always",
                            "gy",   NA,     NA,     NA,     NA,          NA),
                family  = c(NA,        NA,    NA,      "LWL",    "LWL",    "LWL",    "LWL",    "YN",     "YN",
                            NA,     NA,     NA,     NA,     NA,          NA),
                counted = c(FALSE,     FALSE, TRUE,    FALSE,    FALSE,    FALSE,    FALSE,    TRUE,     TRUE,
                            FALSE,  FALSE,  FALSE,  TRUE,   FALSE,       FALSE)
        )
}

# The methods of the table with the given properties (NULL: any), in the order of the table.
dnds_methods <- function(kind = NULL, coal = NULL, family = NULL, counted = NULL) {
        r <- dnds_method_registry()
        keep <- rep(TRUE, length(r$method))
        if (!is.null(kind))
                keep <- keep & is.element(r$kind, kind)
        if (!is.null(coal))
                keep <- keep & is.element(r$coal, coal)
        if (!is.null(family))
                keep <- keep & is.element(r$family, family)
        if (!is.null(counted))
                keep <- keep & r$counted == counted
        r$method[keep]
}

# Views of the table that are used in many places.
#   kaks_calculator_methods() : what KaKs_Calculator is asked for, the pseudo methods included
#   coal_native_methods()     : orthologr's own methods that newer versions of COAL also compute
#   coal_kaks_methods()       : the KaKs_Calculator methods COAL can compute (kaks_tool = "coal"); "GY" only with a
#                               COAL that has it (coal_method_missing())
#   is.dnds_est_method()      : is it a name of the table?
kaks_calculator_methods <- function() dnds_methods(kind = c("kaks", "pseudo"))
coal_native_methods <- function() dnds_methods(kind = "own", coal = "native")
coal_kaks_methods <- function() dnds_methods(kind = "kaks", coal = c("always", "gy"))
is.dnds_est_method <- function(method = NULL) {
        # (without an argument it lists the methods, as the other is.*() helpers do)
        if (is.null(method)) {
                message("The following methods are provided: ", paste(dnds_methods(), collapse = ", "))
                return(FALSE)
        }
        is.element(method, dnds_methods())
}

# Messages for a dnds_est.method that is not a method (U09): the methods of each engine, and the method a mistyped
# name probably meant - the one that matches it ignoring case ("comeron" -> "Comeron"), else the only one at an edit
# distance of 1 ("YNN" -> "YN"). Method names stay case-sensitive: the hint only names the method.
dnds_method_hint <- function(method, candidates = dnds_methods()) {
        if (!is.character(method) || length(method) != 1 || is.na(method) || !nzchar(method))
                return("")
        same <- candidates[tolower(candidates) == tolower(method)]
        if (length(same) == 0) {
                d <- utils::adist(tolower(method), tolower(candidates))[1, ]
                same <- if (sum(d == 1) == 1) candidates[d == 1] else character(0)
        }
        if (length(same) == 0) "" else paste0(" Did you mean \"", same[1], "\"?")
}

# "Methods: ..." of engine = "orthologr" (both engines, and those of KaKs_Calculator only), for the messages of
# dNdS(), substitutionrate() and dNdS_orthogroups(); `kaks_calc`: whether "kaks_calc" is a choice (substitutionrate())
dnds_method_list <- function(kaks_calc = FALSE) {
        q <- function(x) paste0("\"", x, "\"", collapse = ", ")
        paste0("Methods: ", q(coal_native_methods()), " (both engines), ", q(dnds_methods(kind = "kaks", coal = "always")),
               " (both engines; KaKs_Calculator or COAL), ", q(dnds_methods(kind = "kaks", coal = "gy")),
               " (KaKs_Calculator, or a COAL that has it - other numbers, see ?dNdS), ",
               q(dnds_methods(kind = "kaks", coal = NA)), ", \"ALL\"",
               if (kaks_calc) ", \"kaks_calc\"", " (KaKs_Calculator only, engine = \"orthologr\").")
}

# The error for a dnds_est.method that engine = "orthologr" does not compute: `what` names the function.
dnds_method_abort <- function(method, what, kaks_calc = FALSE, arg = "dnds_est.method") {
        shown <- if (is.character(method) && length(method) == 1) paste0("\"", method, "\"") else
                paste(utils::capture.output(dput(method)), collapse = " ")
        orthologr_abort(class = "bad_argument", arg, " = ", shown, " is not a method of ", what, ". ",
                        dnds_method_list(kaks_calc = kaks_calc), dnds_method_hint(method), call. = FALSE)
}

# Stop now when a program the chosen pipeline needs is not there: perl (PAL2NAL), the alignment program, and
# KaKs_Calculator. Called before the search starts. (A missing program used to show only after the complete
# search, as a table in which every gene pair was NA and one warning.)
#   aa_aln_tool : "NW" / "pairwise" (no program), or the multiple alignment program
dnds_require_programs <- function(aa_aln_tool = "NW", aa_aln_path = NULL, codon_aln_tool = "pal2nal",
                                  dnds_est.method = "Comeron", kaks_calc_path = NULL, need_kaks = TRUE) {
        if (identical(codon_aln_tool, "pal2nal") && !nzchar(Sys.which("perl")))
                orthologr_abort(class = "tool_missing", "The program 'perl' was not found on the PATH. It is needed for the codon alignments (PAL2NAL). ",
                     "Please install Perl (Windows: e.g. Strawberry Perl) and make sure 'perl' can be started from a terminal.",
                     call. = FALSE)
        switch(as.character(aa_aln_tool)[1],
               clustalw = is_installed_clustalw(path = aa_aln_path),
               clustalo = is_installed_clustalo(path = aa_aln_path),
               mafft    = is_installed_mafft(path = aa_aln_path),
               muscle   = is_installed_muscle(path = aa_aln_path),
               t_coffee = is_installed_tcoffee(path = aa_aln_path),
               TRUE)
        if (need_kaks && length(dnds_est.method) == 1 && is.element(dnds_est.method, kaks_calculator_methods())) {
                kp <- if (is.null(kaks_calc_path)) NULL else normalizePath(kaks_calc_path, mustWork = FALSE)
                tryCatch(is_installed_kaks_calculator(path = kp), error = function(e) {
                        orthologr_abort(class = "tool_missing", conditionMessage(e), if (!is.null(kp)) paste0(" (kaks_calc_path = '", kaks_calc_path, "')"),
                             " It is needed for dnds_est.method = \"", dnds_est.method,
                             "\"; the methods \"Comeron\" and \"Li\" need no external program.", call. = FALSE)
                })
        }
        invisible(TRUE)
}

# A hint when the best hits are far below the identity the DIAMOND mode is designed for (DIAMOND's documentation:
# "fast" > 90 %, "default" > 60 %, "sensitive" > 40 %): "sensitive" finds more reciprocal best hits there.
# Measured (reciprocal best hits as dNdS() makes them, e-value 1e-5; DECISIONS of 2026-11-14, item 1): S. cerevisiae
# against C. glabrata / K. lactis / A. gossypii 4,165 / 3,846 / 3,695 pairs with "fast", 4,385 / 4,206 / 4,044 with
# "default", 4,503 / 4,366 / 4,253 with "sensitive" (+2.7 / +3.8 / +5.2 % over "default"; recall of OrthoFinder's
# 1:1 orthologs 93.8 / 94.3 / 93.2 % -> 96.7 / 98.1 / 98.1 %), at about 3 times the search time of "default".
sensitivity_hint <- function(hits, aligner, sensitivity_mode) {
        threshold <- c(faster = 70, fast = 70, default = 60, `mid-sensitive` = 60)[sensitivity_mode]
        if (!identical(aligner, "diamond") || length(threshold) != 1 || is.na(threshold) || is.null(hits) ||
            nrow(hits) < 20 || !is.element("perc_identity", names(hits)))
                return(invisible(FALSE))
        med <- stats::median(as.numeric(hits$perc_identity), na.rm = TRUE)
        if (!is.finite(med) || med >= threshold)
                return(invisible(FALSE))
        message("The best hits are ", round(med), " % identical (median). sensitivity_mode = \"", sensitivity_mode,
                "\" is designed for hits of more than ", if (threshold == 70) 90 else 60, " % identity; ",
                "sensitivity_mode = \"sensitive\" finds about 3 - 5 % more reciprocal best hits than \"default\" ",
                "between species of this distance (e.g. two yeasts), at about 3 times the search time.")
        invisible(TRUE)
}
