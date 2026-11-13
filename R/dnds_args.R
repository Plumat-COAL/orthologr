# Internal checks of the arguments of dNdS(): what is not supported is rejected before any work starts, with a
# message that names the argument, the value and the valid values.

# The DIAMOND sensitivity modes. "default" is DIAMOND's own default mode (no sensitivity option on the command line).
diamond_sensitivity_modes <- function() {
        c("fast", "faster", "default", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive")
}

check_sensitivity_mode <- function(sensitivity_mode) {
        if (!is.character(sensitivity_mode) || length(sensitivity_mode) != 1 ||
            !is.element(sensitivity_mode, diamond_sensitivity_modes()))
                stop("sensitivity_mode = '", paste(sensitivity_mode, collapse = "', '"), "' is not a DIAMOND2 sensitivity mode. ",
                     "Please choose one of: '", paste(diamond_sensitivity_modes(), collapse = "', '"), "'.", call. = FALSE)
        invisible(sensitivity_mode)
}

# One value out of a fixed set, with an error in the style of match.arg() that also names the argument.
check_choice <- function(value, choices, arg, hint = NULL) {
        if (!is.character(value) || length(value) != 1 || is.na(value) || !is.element(value, choices))
                stop(arg, " = '", paste(value, collapse = "', '"), "' is not supported. Please choose ",
                     if (length(choices) > 1) "one of: " else "",
                     paste0("'", choices, "'", collapse = ", "), ".",
                     if (!is.null(hint)) paste0(" ", hint) else "", call. = FALSE)
        invisible(value)
}

dnds_check_args <- function(aligner, sensitivity_mode, seq_type, format, ortho_detection, aa_aln_type, aa_aln_tool,
                            aa_aln_params = NULL) {
        
        check_choice(aligner, c("diamond", "blast"), "aligner")
        
        if (identical(ortho_detection, "Orthofinder2"))
                stop("ortho_detection = 'Orthofinder2' is not available in dNdS(): OrthoFinder is not run by this function. ",
                     "Please choose 'RBH' (reciprocal best hit, default) or 'BH' (best hit).", call. = FALSE)
        if (is.element(ortho_detection, c("PO", "DELTA")))
                stop("ortho_detection = '", ortho_detection, "' (", if (ortho_detection == "PO") "ProteinOrtho" else "DELTA-BLAST",
                     ") is no longer available in dNdS(). ",
                     "Please choose 'RBH' (reciprocal best hit, default) or 'BH' (best hit).", call. = FALSE)
        check_choice(ortho_detection, c("RBH", "BH"), "ortho_detection")
        
        if (identical(seq_type, "protein"))
                stop("seq_type = 'protein' is not possible in dNdS(): dN/dS is computed from codon alignments, ",
                     "so both files must hold coding sequences (seq_type = 'cds').", call. = FALSE)
        check_choice(seq_type, "cds", "seq_type")
        
        check_choice(format, "fasta", "format", hint = "dNdS() reads coding sequences from FASTA files.")
        
        if (aligner == "diamond")
                check_sensitivity_mode(sensitivity_mode)
        
        check_choice(aa_aln_type, c("pairwise", "multiple"), "aa_aln_type")
        
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

# The dN/dS methods that are computed by KaKs_Calculator.
kaks_calculator_methods <- function() {
        c("MA", "MS", "NG", "LWL", "LPB", "MLWL", "YN", "MYN", "GY", "GMYN", "kaks_calc", "ALL")
}

# Stop now when a program the chosen pipeline needs is not there: perl (PAL2NAL), the alignment program, and
# KaKs_Calculator. Called before the search starts. (A missing program used to show only after the complete
# search, as a table in which every gene pair was NA and one warning.)
#   aa_aln_tool : "NW" / "pairwise" (no program), or the multiple alignment program
dnds_require_programs <- function(aa_aln_tool = "NW", aa_aln_path = NULL, codon_aln_tool = "pal2nal",
                                  dnds_est.method = "Comeron", kaks_calc_path = NULL, need_kaks = TRUE) {
        if (identical(codon_aln_tool, "pal2nal") && !nzchar(Sys.which("perl")))
                stop("The program 'perl' was not found on the PATH. It is needed for the codon alignments (PAL2NAL). ",
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
                        stop(conditionMessage(e), if (!is.null(kp)) paste0(" (kaks_calc_path = '", kaks_calc_path, "')"),
                             " It is needed for dnds_est.method = \"", dnds_est.method,
                             "\"; the methods \"Comeron\" and \"Li\" need no external program.", call. = FALSE)
                })
        }
        invisible(TRUE)
}
