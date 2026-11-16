# Integration of COAL (https://github.com/Plumat-COAL/COAL): codon alignments and Ka/Ks estimation without Perl
# or KaKs_Calculator, and a batched dN/dS path that hands every gene pair to COAL in one call.
#
# COAL is optional. It is used in process through the R package 'coalr' when that is installed, otherwise through
# the 'coal-cpp' program (part of the Python package 'coal-align').

# Which methods COAL computes is in dnds_method_registry() (column 'coal'): seven of KaKs_Calculator's methods
# always; "Comeron" and "Li" (orthologr's own methods) and "GY" only with a COAL that has them - asked of the COAL at
# hand, never by its version number (a pinned development coalr reports 0.6.0.9000 without GY).

# The coalr that has everything orthologr uses (a release tag of COAL), for the install hints.
coalr_install_hint <- function() {
        "remotes::install_github(\"Plumat-COAL/COAL/packaging/r/coalr@v0.6.0\")"
}

# Does the COAL at hand have a feature? `what`: "native" (Comeron's and Li's method: coalr::coal_info()$methods
# lists both; the coal-cpp program mentions --comeron-flavour in `coal dnds --help`) or "gy" (Goldman & Yang's
# maximum likelihood estimate: coalr::coal_info()$features lists "gy94_ml"; the coal-cpp program lists GY among the
# methods of `coal dnds --help`). FALSE when there is no COAL at all. The answer is remembered for the session, per
# backend and program.
coal_has <- function(what = c("native", "gy")) {
        what <- match.arg(what)
        backend <- tryCatch(coal_backend(), error = function(e) NULL)
        if (is.null(backend))
                return(FALSE)
        key <- paste0("coal_", what, "\r", backend$type, "\r", if (is.null(backend$path)) "" else backend$path)
        hit <- .orthologr_state[[key]]
        if (!is.null(hit))
                return(hit)
        res <- if (backend$type == "coalr") {
                if (what == "native")
                        isTRUE(tryCatch(all(is.element(coal_native_methods(), coalr::coal_info()$methods)), error = function(e) FALSE))
                else
                        is.element("gy94_ml", coal_features())
        } else {
                help <- tryCatch(run_tool(basename(backend$path), args = list("dnds", "--help"), path = dirname(backend$path),
                                          error = FALSE), error = function(e) NULL)
                txt <- c(help$stdout, help$stderr)
                if (what == "native")
                        any(grepl("comeron-flavour", txt, fixed = TRUE))
                else
                        # (the line of -m: "-m METHOD  YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, GY, Comeron ...")
                        any(grepl("-m METHOD", txt, fixed = TRUE) & grepl("(^|[^A-Za-z])GY([^A-Za-z]|$)", txt))
        }
        .orthologr_state[[key]] <- res
        res
}
coal_has_native_methods <- function() coal_has("native")
coal_has_gy <- function() coal_has("gy")

# The methods the COAL at hand computes: the seven KaKs_Calculator methods every COAL has and, with a COAL that has
# them, "GY", "Comeron" and "Li" (in the order of dnds_method_registry()).
coal_methods <- function() {
        dnds_methods(coal = c("always", if (coal_has_native_methods()) "native", if (coal_has_gy()) "gy"))
}

# The methods among `methods` that COAL could compute but the COAL at hand does not, with what is missing:
# a named character vector, method -> requirement (empty when nothing is missing).
coal_method_missing <- function(methods) {
        need <- c(native = "Comeron's and Li's method (coalr >= 0.5.0 whose coalr::coal_info()$methods lists \"Comeron\" and \"Li\")",
                  gy = "Goldman & Yang's maximum likelihood estimate (coalr >= 0.6.0 whose coalr::coal_info()$features lists \"gy94_ml\", or a coal-cpp whose 'coal dnds --help' lists GY)")
        r <- dnds_method_registry()
        kind <- r$coal[match(methods, r$method)]
        miss <- (kind %in% "native" & !coal_has_native_methods()) | (kind %in% "gy" & !coal_has_gy())
        stats::setNames(unname(need[kind[miss]]), methods[miss])
}

# Stop unless COAL computes all `methods`; `what` names the caller for the message ("engine = \"coal\"").
#   other   : what computes the methods COAL never computes
#   methods_of_what : the methods `what` takes at all (default: every method COAL can compute)
coal_check_methods <- function(methods, what, other = "engine = \"orthologr\"",
                               methods_of_what = dnds_methods(coal = c("always", "native", "gy"))) {
        if (!is.character(methods) || length(methods) < 1 || anyNA(methods))
                orthologr_abort(class = "bad_argument", "dnds_est.method must be the name of a method.", call. = FALSE)
        known <- methods_of_what
        unknown <- setdiff(methods, known)
        q <- function(x) paste0("\"", x, "\"", collapse = ", ")
        if (length(unknown) > 0) {
                # (U09) "only available with <other>" only for the methods the other engine / tool has; a name that is
                # no method at all is said to be one, with the method it probably meant
                elsewhere <- unknown[is.element(unknown, dnds_methods())]
                none <- setdiff(unknown, elsewhere)
                orthologr_abort(class = "bad_argument", what, " computes dnds_est.method = ", q(known), ".",
                     if (length(elsewhere) > 0) paste0(" ", q(elsewhere), if (length(elsewhere) > 1) " are" else " is",
                                                       " only available with ", other, "."),
                     if (length(none) > 0) paste0(" ", q(none), if (length(none) > 1) " are not methods" else " is not a method",
                                                  " of either.", if (length(none) == 1) dnds_method_hint(none) else ""),
                     call. = FALSE)
        }
        miss <- coal_method_missing(methods)
        if (length(miss) > 0) {
                have <- tryCatch(as.character(utils::packageVersion("coalr")), error = function(e) NULL)
                backend <- tryCatch(coal_backend()$type, error = function(e) NULL)
                orthologr_abort(class = "package_missing", what, ": dnds_est.method = ", q(names(miss)),
                     " needs a COAL that computes ", paste(unique(miss), collapse = " and "), "; ",
                     if (identical(backend, "binary")) "the coal-cpp program in use does not" else
                             if (is.null(have)) "coalr is not installed" else paste0("the installed coalr ", have, " does not"),
                     ". Update it with ",
                     if (identical(backend, "binary")) paste0("'pip install -U coal-align' (the coal-cpp program), or install coalr with ",
                                                              coalr_install_hint()) else coalr_install_hint(),
                     ", choose one of ",
                     q(intersect(known, coal_methods())), ", or use ", other, ".", call. = FALSE)
        }
        invisible(methods)
}

# The versions of COAL in use, for the record of a run.
# The backend that is used (coal_backend(): "coalr" or "binary", also when options(orthologr.coal_backend) forces
# the program although coalr is installed) is named in `coal_backend`.
coal_versions <- function() {
        v <- list()
        backend <- tryCatch(coal_backend(), error = function(e) NULL)
        if (is.null(backend))
                return(v)
        v$coal_backend <- backend$type
        if (backend$type == "coalr") {
                v$coalr <- package_version_text("coalr")
                # (coalr names the version of its engine coal_version; there is no field 'version', so the engine's
                # version was never recorded - M3 of REVIEW-OL-1)
                engine <- tryCatch(coalr::coal_info()$coal_version, error = function(e) NULL)
                if (!is.null(engine) && length(engine) == 1)
                        v$coal <- as.character(engine)
        } else {
                v$`coal-cpp` <- tool_version_line(backend$path, "--version")
        }
        v
}

# The arguments that make COAL's Comeron and Li give orthologr's numbers, for a function `fun` of coalr (only the
# arguments that function has - older versions of coalr have none of them):
#   comeron_flavour = "gestimator", li_flavour = "seqinr" : the pair is read as gestimator / seqinr::kaks() read it
#   terminal_stops = "drop"    : orthologr 0.5.0 translates and aligns without the stop codon;
#                    "compare" : with options(orthologr.legacy_pairwise_aln = TRUE), i.e. the alignments of 0.4.3,
#                                in which the two stop codons are an aligned column that gestimator compares
#   comeron_twosubs, comeron_min_sites : Comeron's pathways for codons that differ at two positions and the
#                    smallest site class that can make dN or dS NA, as the native engine computes them
#                    (options(orthologr.comeron_twosubs), options(orthologr.comeron_min_sites))
coal_flavour_args <- function(fun, methods = NULL) {
        args <- list(
                comeron_flavour   = "gestimator",
                comeron_twosubs   = comeron_twosubs(),
                comeron_min_sites = comeron_min_sites(),
                li_flavour        = "seqinr",
                terminal_stops    = if (isTRUE(orthologr_option("orthologr.legacy_pairwise_aln", FALSE))) "compare" else "drop"
        )
        have <- names(formals(fun))
        # A coalr whose coal_info()$features lacks "comeron_min_sites_positive" (COAL before #76) computes Comeron
        # otherwise than engine = "orthologr": without comeron_twosubs / comeron_min_sites (older still) along
        # libsequence's branches and with NA for small site classes, and with them a 0 where the rule of the
        # smallest site class gives 0 (e.g. dS = 0 for ATGGCTAAA / ATGGCCAAA, NA with engine = "orthologr").
        # The numbers are still computed (as by every coalr since 0.5.0), so this is a warning, of the class that
        # documents an old coalr, not an error.
        if (is.element("Comeron", methods) && !is.element("comeron_min_sites_positive", coal_features()))
                orthologr_warn(class = "coal_version",
                               "This version of coalr predates COAL's Comeron refinements (coalr::coal_info()$features ",
                               "does not list \"comeron_min_sites_positive\"): its Comeron can differ from engine = ",
                               "\"orthologr\" (codons that differ at two positions, small site classes; e.g. dS = 0 ",
                               "instead of NA). Install a current coalr: ",
                               coalr_install_hint(), ".", call. = FALSE)
        args[intersect(names(args), have)]
}

# The features the installed coalr reports (coalr::coal_info()$features; character(0) for a coalr without that
# field or without coalr).
coal_features <- function() {
        f <- tryCatch(coalr::coal_info()$features, error = function(e) NULL)
        if (is.character(f)) f else character(0)
}

# The matrix to give COAL: with aa_aln_scoring = "pwalign" and BLOSUM62, pwalign's copy of BLOSUM62 (it scores B, Z
# and X otherwise than NCBI's current file, which COAL's "BLOSUM62" is) when the installed coalr has it, so that the
# protein alignments are those of pwalign::pairwiseAlignment() also for proteins with X.
coal_matrix_arg <- function(aa_sub_matrix, aa_aln_scoring) {
        if (identical(aa_aln_scoring, "pwalign") && identical(aa_sub_matrix, "BLOSUM62") &&
            is.element("BLOSUM62-PWALIGN", tryCatch(coalr::coal_info()$matrices, error = function(e) character(0))))
                return("BLOSUM62-PWALIGN")
        aa_sub_matrix
}

# What orthologr gives COAL, for the record of the run ("engine_settings"), and the `coal pipeline` command line
# that computes the same pairs with the same settings (C05 of the integration pass; ?dNdS_coal, "The same run with
# the COAL command line").
#   fun : the coalr function;  args : the arguments orthologr passes (without the sequences)
coal_engine_settings <- function(fun_name, fun, args, mode = c("pairs", "clusters")) {
        mode <- match.arg(mode)
        args <- c(args, coal_flavour_args(fun))
        has <- names(formals(fun))
        tie <- if (is.element("tie_break", has)) {
                if (identical(args$preset, "pwalign")) "pwalign" else "leftmost"
        } else "leftmost (coalr without tie_break)"
        matrix <- args$matrix
        if (is.null(matrix))
                matrix <- if (identical(args$preset, "pwalign") && is.element("tie_break", has)) "BLOSUM62-PWALIGN" else "BLOSUM62"
        flags <- c(paste("-m", args$methods), "--preset", args$preset, "--matrix", matrix,
                   if (is.element("tie_break", has)) c("--tie-break", tie),
                   if (!is.null(args$comeron_flavour)) c("--comeron-flavour", args$comeron_flavour),
                   # (M1 of REVIEW-OL-1) with --comeron-flavour gestimator COAL's command line counts the codons
                   # that differ at two positions along libsequence's branches unless it is told otherwise, and
                   # orthologr's default is Comeron's two pathways: without these two flags the command gave other
                   # numbers than the run it records
                   if (!is.null(args$comeron_twosubs)) c("--comeron-twosubs", args$comeron_twosubs),
                   if (!is.null(args$comeron_min_sites)) c("--comeron-min-sites", format(args$comeron_min_sites)),
                   if (!is.null(args$li_flavour)) c("--li-flavour", args$li_flavour),
                   if (!is.null(args$terminal_stops)) c("--terminal-stops", args$terminal_stops))
        cmd <- if (mode == "pairs") {
                c("coal pipeline pairs --query QUERY_CDS --subject SUBJECT_CDS --pairs PAIRS.tsv", flags,
                  "-o OUT_DIR")
        } else {
                c("coal pipeline clusters --cds CDS --clusters CLUSTERS.tsv", flags, "-o OUT_DIR")
        }
        list(coalr_function = fun_name,
             coalr = package_version_text("coalr"),
             preset = args$preset, matrix = matrix, tie_break = tie,
             comeron_flavour = args$comeron_flavour, comeron_twosubs = args$comeron_twosubs,
             comeron_min_sites = args$comeron_min_sites, li_flavour = args$li_flavour,
             terminal_stops = args$terminal_stops,
             coal_command = paste(cmd, collapse = " "))
}

# options(orthologr.legacy_pairwise_aln = TRUE) - set on its own, by orthologr_legacy_options() or by
# options(orthologr.legacy = TRUE) - asks for the pairwise protein alignments of 0.4.3 (no substitution matrix, cut
# to the aligned range). COAL aligns itself and cannot make them, so the numbers are not those of 0.4.3: say so
# instead of returning a table that looks like a legacy run.
coal_legacy_warning <- function(what) {
        if (isTRUE(orthologr_option("orthologr.legacy_pairwise_aln", FALSE)))
                orthologr_warn(class = "legacy_engine", what, ": the COAL engine cannot reproduce the pairwise ",
                               "protein alignments of orthologr 0.4.3 (options(orthologr.legacy_pairwise_aln = TRUE), ",
                               "set by orthologr_legacy_options()), so its dN and dS differ from those of 0.4.3. ",
                               "Use engine = \"orthologr\" for the numbers of 0.4.3.")
        invisible(NULL)
}

# The substitution matrices COAL has (aa_sub_matrix of dNdS(engine = "coal") / dNdS_coal()).
coal_sub_matrix_names <- function() c("BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250")

# The alignment scoring of the COAL entry points: one argument `aa_aln_scoring` everywhere.
check_aa_aln_scoring <- function(aa_aln_scoring) {
        check_choice(aa_aln_scoring, c("pwalign", "needle"), "aa_aln_scoring")
}

# dN, dS and dN/dS of COAL as orthologr's own methods return them, so that engine = "coal" / kaks_tool = "coal"
# and the native path give the same table for the same codon alignment. `tbl`: a table of coalr (columns method,
# dN, dS, dNdS, note). The KaKs_Calculator methods are left as they are.
#   Comeron : COAL's gestimator flavour returns gestimator's value also where a correction is undefined (the
#             finite value of orthologr 0.4.3) and says so in the note ("Ka: correction undefined ...",
#             "Ks: correction undefined ..."). orthologr returns NA for such a dN / dS and for dN/dS, unless
#             options(orthologr.legacy_comeron_saturation = TRUE). COAL applies the same rule as
#             src/Comeron95.h for classes without sites, so the NA pattern is that of the native path.
#             gestimator writes 4 significant digits, COAL returns doubles: the values are rounded to 4
#             significant digits, which makes them identical to the native table.
#   Li      : li_values() as for seqinr::kaks(): negative values and the saturation marker are NA.
coal_orthologr_values <- function(tbl) {
        tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
        if (nrow(tbl) == 0 || !is.element("method", names(tbl)))
                return(tbl)
        note <- if (is.element("note", names(tbl))) as.character(tbl$note) else rep("", nrow(tbl))
        note[is.na(note)] <- ""
        com <- !is.na(tbl$method) & tbl$method == "Comeron"
        if (any(com)) {
                if (!comeron_legacy_saturation()) {
                        ka <- com & grepl("Ka: correction undefined", note, fixed = TRUE)
                        ks <- com & grepl("Ks: correction undefined", note, fixed = TRUE)
                        tbl$dN[ka] <- NA_real_
                        tbl$dS[ks] <- NA_real_
                        tbl$dNdS[ka | ks] <- NA_real_
                }
                for (col in c("dN", "dS", "dNdS"))
                        tbl[[col]][com] <- signif(as.numeric(tbl[[col]][com]), 4)
        }
        li <- !is.na(tbl$method) & tbl$method == "Li"
        if (any(li)) {
                dN <- as.numeric(tbl$dN[li])
                dS <- as.numeric(tbl$dS[li])
                bad_n <- is.na(dN) | dN < 0 | dN > 9.9
                bad_s <- is.na(dS) | dS < 0 | dS > 9.9
                ratio <- dN / dS
                dN[bad_n] <- NA_real_
                dS[bad_s] <- NA_real_
                ratio[bad_n | bad_s | !is.finite(ratio)] <- NA_real_
                tbl$dN[li] <- dN
                tbl$dS[li] <- dS
                tbl$dNdS[li] <- ratio
        }
        tbl
}

coal_install_message <- function() {
        paste0(
                "COAL is not available. Install its R package (recommended; runs inside R, no other program needed):\n",
                "    ", coalr_install_hint(), "\n",
                "or put COAL's 'coal-cpp' program on the PATH (pip install coal-align; ",
                "or set options(orthologr.coal_path = \"/path/to/coal-cpp\"))."
        )
}

# Which COAL to use: "coalr" (the R package, in process) or "binary" (the coal-cpp program).
# options(orthologr.coal_backend = "coalr" | "binary") forces one; 'need' = "coalr" asks for the package.
# Returns list(type = , path = ) or stops with the installation message.
coal_backend <- function(need = NULL) {
        forced <- getOption("orthologr.coal_backend", NULL)
        if (!is.null(forced) && !is.element(forced, c("coalr", "binary", "none")))
                orthologr_abort(class = "bad_argument", "options(orthologr.coal_backend) must be \"coalr\", \"binary\" or \"none\" (as if COAL were not installed; for tests).", call. = FALSE)
        if (identical(forced, "none"))
                orthologr_abort(class = c("package_missing", "tool_missing"), coal_install_message())

        if (!identical(forced, "binary") && requireNamespace("coalr", quietly = TRUE))
                return(list(type = "coalr", path = NULL))

        if (identical(forced, "coalr") || identical(need, "coalr"))
                orthologr_abort(class = "package_missing",
                        "This needs the R package 'coalr' (COAL's C++ engine inside R). Install it with:\n",
                        "    ", coalr_install_hint(),
                        call. = FALSE
                )

        exe <- coal_binary()
        if (is.null(exe))
                orthologr_abort(class = c("package_missing", "tool_missing"), coal_install_message())
        return(list(type = "binary", path = exe))
}

# The coal-cpp program: options(orthologr.coal_path), the environment variable COAL_CPP, then the PATH.
# NULL when there is none.
coal_binary <- function() {
        candidates <- c(getOption("orthologr.coal_path", ""), Sys.getenv("COAL_CPP", ""),
                        unname(Sys.which("coal-cpp")))
        for (exe in candidates) {
                if (is.character(exe) && length(exe) == 1 && nzchar(exe) && file.exists(exe) && !dir.exists(exe))
                        return(normalizePath(exe))
        }
        return(NULL)
}

# Is COAL available (as 'coalr' or as the coal-cpp program)?
is_installed_coal <- function() {
        isTRUE(tryCatch({
                coal_backend()
                TRUE
        }, error = function(e) FALSE))
}

# Run coal-cpp with 'args' (a list; quoted by run_tool()), stdout to 'stdout' (a file) when given. Stops when it
# fails or when 'output' was not written. Started through run_tool() like every other external program.
coal_run_binary <- function(exe, args, context, output, stdout = NULL) {
        run_tool(
                basename(exe),
                args    = args,
                path    = dirname(exe),
                stdout  = stdout,
                tool    = "COAL (coal-cpp)",
                context = context,
                output  = output
        )
        invisible(TRUE)
}

# Read a FASTA file as a named character vector (ids up to the first white space; sequences as written).
read_fasta_strings <- function(file) {
        x <- seqinr::read.fasta(file, as.string = TRUE, set.attributes = FALSE, forceDNAtolower = FALSE,
                                seqtype = "AA", strip.desc = TRUE)
        res <- vapply(x, function(s) as.character(s)[1], character(1))
        names(res) <- sub("\\s.*$", "", names(x))
        return(res)
}

# Write a named character vector as FASTA (one line per sequence) or in CLUSTAL format.
write_alignment <- function(seqs, file, format = "fasta") {
        if (format == "fasta") {
                writeLines(as.vector(rbind(paste0(">", names(seqs)), unname(seqs))), file)
                return(invisible(file))
        }
        width <- 60
        n <- max(nchar(seqs))
        pad <- max(nchar(names(seqs))) + 8   # as PAL2NAL writes it
        lines <- c("CLUSTAL W multiple sequence alignment", "")
        for (start in seq(1, max(n, 1), by = width)) {
                lines <- c(lines,
                           paste0(formatC(names(seqs), width = -pad), substr(seqs, start, start + width - 1)),
                           "")
        }
        writeLines(lines, file)
        invisible(file)
}

# The codon alignment of a protein alignment file and a CDS file with COAL, written to 'file.out'.
coal_codon_aln_file <- function(file_aln, file_nuc, file.out, format = "fasta", aln_format = NULL, context = NULL) {
        backend <- coal_backend()
        if (is.null(context))
                context <- pair_context(file_nuc)
        # as for PAL2NAL: an alignment left by an earlier call is removed, the new one is written under a name of
        # its own and becomes 'file.out' only when it is complete
        if (file.exists(file.out))
                unlink(file.out)
        tmp.out <- tempfile("aln_", tmpdir = dirname(file.out), fileext = ".tmp")
        on.exit(unlink(tmp.out), add = TRUE)

        if (backend$type == "binary") {
                # coal-cpp takes PAL2NAL's command line
                coal_run_binary(backend$path,
                                list(file_aln, file_nuc, "-output", format, "-o", tmp.out),
                                context = context, output = tmp.out)
                publish_file(tmp.out, file.out)
                return(invisible(file.out))
        }

        first <- readLines(file_aln, n = 1, warn = FALSE)
        if (length(first) > 0 && startsWith(first, ">")) {
                rows <- read_fasta_strings(file_aln)
        } else {
                aln <- tryCatch(suppressWarnings(seqinr::read.alignment(file_aln, format = "clustal", forceToLower = FALSE)),
                                error = function(e) NULL)
                if (!is.null(aln) && isTRUE(aln$nb > 0)) {
                        rows <- structure(toupper(unlist(aln$seq)), names = sub("\\s.*$", "", aln$nam))
                } else {
                        # seqinr reads no sequence from a CLUSTAL file without a conservation line (the format
                        # orthologr's write_alignment(), PAL2NAL, T-Coffee and MUSCLE write), as in codon_aln()
                        rows <- toupper(read_clustal(file_aln))
                }
        }
        cds <- read_fasta_strings(file_nuc)
        # PAL2NAL pairs the two files by position; so does orthologr (both files hold the pair in the same order)
        # (a failed codon alignment has the class of PAL2NAL's, orthologr_error_alignment_failed, so that callers
        # can handle both tools alike; it was a plain error)
        if (length(rows) != length(cds))
                orthologr_abort(class = c("alignment_failed", "tool_failed"),
                                "COAL: the protein alignment holds ", length(rows), " sequences and the CDS file ",
                                length(cds), " for ", context, ".", call. = FALSE)
        codon <- coal_codon_rows(rows, cds, context)
        names(codon) <- names(cds)
        write_alignment(codon, tmp.out, format = format)
        publish_file(tmp.out, file.out)
        invisible(file.out)
}

# The codon alignment rows of the protein alignment rows 'rows' and the coding sequences 'cds' (same order) with
# coalr (in process); a failure is an error of class orthologr_error_alignment_failed naming 'context'.
coal_codon_rows <- function(rows, cds, context) {
        tryCatch(
                coalr::coal_codon_align(unname(rows), unname(cds)),
                error = function(e) orthologr_abort(class = c("alignment_failed", "tool_failed"),
                                                    "COAL could not build the codon alignment for ", context, ": ",
                                                    conditionMessage(e), call. = FALSE)
        )
}

# kaks_tool = "coal" with a method that is neither one of the KaKs_Calculator methods COAL computes nor one of
# orthologr's own (Comeron, Li, which are computed by orthologr whatever kaks_tool is).
kaks_tool_coal_abort <- function(arg, method) {
        orthologr_abort(class = "bad_argument", "kaks_tool = \"coal\" takes ", arg, " = ",
                        paste0("\"", kaks_tool_coal_methods(), "\"", collapse = ", "),
                        " (GY with a COAL that has it; Comeron and Li are computed by orthologr itself, with either ",
                        "kaks_tool); ", arg, " = \"", paste(method, collapse = "\", \""),
                        "\" needs kaks_tool = \"kaks_calculator\".", call. = FALSE)
}

# The methods dNdS() and substitutionrate() take with kaks_tool = "coal": orthologr's own (Comeron, Li) and those
# COAL computes in place of KaKs_Calculator.
kaks_tool_coal_methods <- function() c(dnds_methods(kind = "own"), coal_kaks_methods())

# Ka/Ks of a two-sequence codon alignment file with COAL: the table substitutionrate() returns for the
# KaKs_Calculator methods (query_id, subject_id, dN, dS, dNdS, method). Only for the KaKs_Calculator methods COAL
# computes (coal_kaks_methods()): substitutionrate() computes Comeron and Li with orthologr's own code whatever
# kaks_tool is.
# COAL's status of the pair ("ok", "warning", "not_converged", "saturated", "no_data") is the attribute
# "engine_status" of the table: compute_dnds() carries it to dNdS(), which reports an estimate that did not
# converge as status "not_converged" (up to 0.6.0 such a pair was reported as "saturated").
coal_kaks_file <- function(file, est.method, subst_name = NULL) {
        query_id <- NULL
        coal_check_methods(est.method, "kaks_tool = \"coal\"", other = "kaks_tool = \"kaks_calculator\"",
                           methods_of_what = c(coal_native_methods(), coal_kaks_methods()))
        if (!is.element(est.method, coal_kaks_methods()))
                orthologr_abort(class = "bad_argument", "coal_kaks_file() computes the KaKs_Calculator methods ",
                                paste(coal_kaks_methods(), collapse = ", "), "; \"", est.method,
                                "\" is computed by substitutionrate() itself.", call. = FALSE)
        backend <- coal_backend()
        context <- pair_context(file, subst_name)
        seqs <- read_fasta_strings(file)
        if (length(seqs) != 2)
                orthologr_abort(class = "bad_input", "The codon alignment of ", context, " must hold two sequences, but '",
                                file, "' holds ", length(seqs), ".", call. = FALSE)

        if (backend$type == "coalr") {
                k <- do.call(coalr::coal_kaks, c(list(seqs[[1]], seqs[[2]], method = est.method),
                                                 coal_flavour_args(coalr::coal_kaks, est.method)))
        } else {
                calc_dir <- file.path(orthologr_tmpdir(), "_calculation")
                if (!file.exists(calc_dir))
                        dir.create(calc_dir, recursive = TRUE)
                # input and output of this call have names of their own and are removed afterwards
                axt_file <- tempfile(paste0(safe_file_name(basename(file)), "_"), tmpdir = calc_dir, fileext = ".axt")
                out_file <- paste0(axt_file, ".coal.tsv")
                on.exit(unlink(c(axt_file, out_file)), add = TRUE)
                writeLines(c(if (is.null(subst_name)) "pair" else subst_name, seqs[[1]], seqs[[2]], ""), axt_file)
                coal_run_binary(backend$path,
                                list("dnds", axt_file, "-o", out_file, "-m", est.method, "-j", "1", "-q",
                                     "--digits", "17"),
                                context = context, output = out_file)
                k <- utils::read.delim(out_file, sep = "\t", quote = "", na.strings = "NA", stringsAsFactors = FALSE)
                if (nrow(k) < 1)
                        orthologr_abort(class = "tool_failed", "COAL returned no result for ", context, ".", call. = FALSE)
        }
        if (identical(k$status[1], "error"))
                orthologr_abort(class = "tool_failed", "COAL could not estimate dN/dS for ", context, ": ", k$note[1], call. = FALSE)
        k <- coal_orthologr_values(k)
        res <- data.table::data.table(query_id = names(seqs)[1], subject_id = names(seqs)[2],
                                      dN = as.numeric(k$dN[1]), dS = as.numeric(k$dS[1]),
                                      dNdS = as.numeric(k$dNdS[1]), method = as.character(k$method[1]))
        data.table::setkey(res, query_id)
        data.table::setattr(res, "engine_status", as.character(k$status[1]))
        return(res)
}

# A CDS table of read.cds() as the named character vector COAL takes.
cds_vector <- function(cds_dt) {
        geneids <- seqs <- NULL
        structure(cds_dt[ , seqs], names = cds_dt[ , geneids])
}

#' @title Compute dNdS values for two organisms with the COAL engine (batched)
#' @description The fast path of \code{\link{dNdS}}: orthologs are detected with \pkg{orthologr}'s own DIAMOND or
#' BLAST best (reciprocal) hit search, and all orthologous gene pairs are then handed to COAL
#' (\url{https://github.com/Plumat-COAL/COAL}) in one call. COAL translates the coding sequences, aligns each pair
#' of proteins globally (Needleman-Wunsch, BLOSUM62), threads the codons onto the protein alignment and estimates
#' dN, dS and dN/dS, inside the R process and on \code{comp_cores} threads. No file is written and no program is
#' started per gene pair (no Perl / PAL2NAL, no KaKs_Calculator).
#'
#' The result is the table \code{\link{dNdS}} returns: the same columns and types, joined with the hit table.
#' The default method differs: \code{dNdS_coal()} computes \code{"YN"}, \code{dNdS()} (also with
#' \code{engine = "coal"}) \code{"Comeron"}.
#' \code{dNdS(..., engine = "coal")} calls this function.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param aligner a character string specifying the sequence aligner used to detect orthologs:
#' \code{"diamond"} (Default) or \code{"blast"}.
#' @param sensitivity_mode the DIAMOND sensitivity mode, see \code{\link{dNdS}}.
#' @param aligner_path a character string specifying the path to the DIAMOND or BLAST program (in case you don't use the default path).
#' @param seq_type a character string specifying the sequence type stored in the input files. Only \code{"cds"} is supported.
#' @param format a character string specifying the file format of the sequence files. Default is \code{format = "fasta"}.
#' @param ortho_detection the orthology inference method: \code{"RBH"} (Default; reciprocal best hit) or \code{"BH"} (best hit).
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input files.
#' @param aligner_params a character string specifying additional parameters that shall be passed to DIAMOND or BLAST.
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND or BLAST hit detection.
#' @param dnds_est.method the dN/dS estimation method(s) computed by COAL: \code{"YN"} (Default; Yang and Nielsen 2000),
#' \code{"NG"}, \code{"LWL"}, \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"} or \code{"MYN"}. The values agree with
#' KaKs_Calculator 2.0's for the same codon alignment. Several methods can be given; the result then has one row
#' per gene pair and method.
#' \code{"Comeron"} and \code{"Li"} - the two methods \pkg{orthologr} computes itself - are available when the
#' installed \pkg{coalr} lists them in \code{coalr::coal_info()$methods} (coalr >= 0.5.0), \code{"GY"} (Goldman and
#' Yang 1994, maximum likelihood) when \code{coalr::coal_info()$features} lists \code{"gy94_ml"} (coalr >= 0.6.0);
#' otherwise an error says so. COAL's GY uses F3x4 codon frequencies and is fitted as PAML codeml's pairwise
#' estimate; it is not KaKs_Calculator's GY (GY-HKY), so its values differ from those of
#' \code{dNdS(dnds_est.method = "GY")} with KaKs_Calculator. See Details.
#' @param aa_aln_scoring how the pairwise protein alignment is scored:
#' \itemize{
#' \item \code{"pwalign"} (Default): the scoring of \code{\link{pairwise_aln}} (BLOSUM62, gap opening 10, gap
#' extension 4, end gaps penalised), so that the alignments are those \code{dNdS()} computes itself;
#' \item \code{"needle"}: the defaults of EMBOSS needle and of COAL (BLOSUM62, gap opening 10, gap extension 0.5,
#' end gaps free).
#' }
#' \code{dNdS(engine = "coal")}, \code{\link{dNdS_clusters}} and \code{dNdS_orthogroups(engine = "coal")} have the
#' same argument with the same default, so that the same gene pair gets the same alignment and values from each
#' of them.
#' @param aa_sub_matrix the substitution matrix: \code{"BLOSUM62"} (Default), \code{"BLOSUM45"}, \code{"BLOSUM80"} or \code{"PAM250"}.
#' @param comp_cores the number of threads COAL (and the DIAMOND or BLAST search) uses.
#' @param quiet \code{TRUE} (default): one line at the start and one summary at the end of the run (messages).
#' \code{FALSE}: also the messages of the search program.
#' @param print_citation a logical value indicating whether the citation note is shown (as a message, once per R session).
#' @param best_hit_rule how best hits are chosen: \code{"bitscore"} (default since orthologr 0.5.0) or \code{"evalue"}
#' (the rule of orthologr <= 0.4.3); see \code{\link{dNdS}} and \code{\link{diamond_best}}.
#' \code{options(orthologr.best_hit_rule = "evalue")} sets the old rule for a session.
#' @details
#' Differences from \code{dNdS()} with \code{engine = "orthologr"}:
#' \itemize{
#' \item \code{"Comeron"} and \code{"Li"} are COAL's ports of the programs \pkg{orthologr} uses (libsequence's
#' gestimator, \code{seqinr::kaks()}), asked to read a pair as those programs do (\code{comeron_flavour =
#' "gestimator"}, \code{li_flavour = "seqinr"}; the stop codon at the end is not compared, as in the alignments of
#' \code{dNdS()} - it is with \code{options(orthologr.legacy_pairwise_aln = TRUE)}). \pkg{orthologr} then applies
#' its own rules to the values, so that both engines return the same numbers for the same codon alignment:
#' Comeron's saturation rule (\code{NA} where a correction is undefined in a class that has sites; finite values
#' with \code{options(orthologr.legacy_comeron_saturation = TRUE)}; see \code{\link{substitutionrate}}), rounding
#' to the 4 significant digits gestimator writes (Comeron only), and \code{NA} for negative or saturated values of
#' Li's method. \code{attr(result, "coal")} holds COAL's unrounded values and notes. With
#' \code{aa_aln_scoring = "pwalign"} and a \pkg{coalr} that has COAL's pwalign tie-break (\code{tie_break}),
#' COAL's protein alignments are those of \code{pwalign::pairwiseAlignment()}, row for row, so both engines
#' compute the same codon alignments (see NEWS for the measured agreement); with an older \pkg{coalr} they may
#' place a gap differently where several alignments score the same.
#' \item A dN or dS of exactly 0 is returned as 0 (KaKs_Calculator prints values below 1e-6 as \code{NA}); a
#' saturated estimate is \code{NA} (KaKs_Calculator prints a capped value).
#' \item Codon columns with a gap, an ambiguous base or an internal stop codon are skipped; a terminal stop codon is
#' not scored.
#' \item A gene pair that cannot be computed is returned with \code{NA} values and is listed in
#' \code{attr(result, "failed_pairs")}, as in \code{dNdS()}.
#' }
#' COAL's complete table for the gene pairs (standard errors, numbers of sites and substitutions, codons compared,
#' protein identity, status and note) is attached as \code{attr(result, "coal")}.
#'
#' This function needs the R package \pkg{coalr}:
#' \code{remotes::install_github("Plumat-COAL/COAL/packaging/r/coalr@v0.6.0")} (COAL's release 0.6.0; a newer
#' release works as well).
#'
#' \strong{The same run with the COAL command line, and the other way.} The record of the run,
#' \code{attr(result, "orthologr_run")}, has in \code{engine_settings} what orthologr gave COAL (\code{preset},
#' \code{matrix}, \code{tie_break}, \code{comeron_flavour}, \code{comeron_twosubs}, \code{comeron_min_sites},
#' \code{li_flavour}, \code{terminal_stops}) and in
#' \code{engine_settings$coal_command} the \code{coal pipeline pairs} command line with these settings; its
#' \code{commands} are the DIAMOND (or BLAST+) searches of the run. With the gene pairs of the result as
#' \code{--pairs}, the command line computes the same pairs with the same alignments:
#' \preformatted{
#' res <- dNdS_coal("query.fa", "subject.fa", dnds_est.method = "YN")
#' pairs <- unique(as.data.frame(res)[, c("query_id", "subject_id")])
#' write.table(pairs, "pairs.tsv", sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#' attr(res, "orthologr_run")$engine_settings$coal_command
#' # coal pipeline pairs --query QUERY_CDS --subject SUBJECT_CDS --pairs PAIRS.tsv -m YN --preset pwalign
#' #   --matrix BLOSUM62-PWALIGN --tie-break pwalign --comeron-flavour gestimator --comeron-twosubs comeron
#' #   --comeron-min-sites 3 --li-flavour seqinr --terminal-stops drop -o OUT_DIR
#' }
#' Its \code{dnds.tsv} holds COAL's own values, which \code{attr(result, "coal")} keeps too; orthologr's
#' result applies its rules to Comeron and Li (above). The pairs must come from the result: COAL's own search
#' takes DIAMOND's first target as the best hit, orthologr the hit with the highest bit score, and orthologr
#' leaves out coding sequences with ambiguous bases (\code{options(orthologr.ambiguous_cds)}), which COAL keeps.
#' The other way, a \code{coal pipeline pairs} run in R: \code{coalr::coal_dnds_pairs()} with the run's
#' \code{orthologs.tsv} as \code{pairs} and the run's options (\code{report.json}, \code{parameters}) as arguments of
#' the same names; a \code{coal pipeline clusters} run: \code{\link{dNdS_clusters}} with its \code{clusters.tsv} as
#' \code{clusters}.
#' @author Chloe Newton and Keanu Murati
#' @return A tibble with the columns of \code{\link{dNdS}}: \code{query_id}, \code{subject_id}, \code{dN}, \code{dS},
#' \code{dNdS}, \code{method} and the columns of the hit table.
#' @examplesIf nzchar(system.file(package = "coalr"))
#' if (nzchar(Sys.which("diamond"))) {
#'   dNdS_coal(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'             dnds_est.method = "YN", comp_cores = 2)
#'
#'   # the same through dNdS()
#'   dNdS(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'        dnds_est.method = "YN", engine = "coal")
#' }
#' @seealso \code{\link{dNdS}}, \code{\link{dNdS_clusters}}, \code{\link{diamond_rec}}
#' @export
dNdS_coal <- function(query_file,
                      subject_file,
                      aligner            = "diamond",
                      sensitivity_mode   = getOption("orthologr.sensitivity_mode", "default"),
                      aligner_path       = NULL,
                      seq_type           = "cds",
                      format             = "fasta",
                      ortho_detection    = "RBH",
                      delete_corrupt_cds = FALSE,
                      aligner_params     = NULL,
                      eval               = "1E-5",
                      dnds_est.method    = "YN",
                      aa_aln_scoring     = "pwalign",
                      aa_sub_matrix      = "BLOSUM62",
                      comp_cores         = 1,
                      quiet              = TRUE,
                      print_citation     = TRUE,
                      best_hit_rule      = getOption("orthologr.best_hit_rule", "bitscore")) {
        sensitivity_mode <- orthologr_sensitivity_mode(if (missing(sensitivity_mode)) NULL else sensitivity_mode)

        best_hit_rule <- orthologr_best_hit_rule(if (missing(best_hit_rule)) NULL else best_hit_rule)

        if (!is.element(aligner, c("diamond", "blast")))
                orthologr_abort(class = "bad_argument", "engine = \"coal\" detects orthologs with aligner = \"diamond\" or \"blast\".", call. = FALSE)
        if (!is.element(ortho_detection, c("BH", "RBH")))
                orthologr_abort(class = "bad_argument", "engine = \"coal\" supports ortho_detection = \"RBH\" or \"BH\".", call. = FALSE)
        if (!identical(seq_type, "cds"))
                orthologr_abort(class = "bad_argument", "dN/dS estimation needs coding sequences: seq_type = \"cds\".", call. = FALSE)
        check_choice(format, "fasta", "format", hint = "dNdS_coal() reads coding sequences from FASTA files.")
        check_input_files(query_file, subject_file)
        if (aligner == "diamond")
                check_sensitivity_mode(sensitivity_mode)
        comp_cores <- check_comp_cores(comp_cores)
        check_aa_aln_scoring(aa_aln_scoring)
        check_choice(aa_sub_matrix, coal_sub_matrix_names(), "aa_sub_matrix", hint = "(the matrices of the COAL engine)")
        coal_backend(need = "coalr")
        coal_check_methods(dnds_est.method, "engine = \"coal\"")
        coal_legacy_warning("dNdS(engine = \"coal\") / dNdS_coal()")

        # the arguments of this call, for the record of the run (attr(<result>, "orthologr_run"))
        run_args <- run_arguments()

        message("dNdS_coal(): ", if (ortho_detection == "RBH") "reciprocal best hits" else "best hits", " of '",
                basename(query_file), "' against '", basename(subject_file), "' (", aligner,
                if (aligner == "diamond") paste0(", sensitivity_mode '", sensitivity_mode, "'"),
                ", e-value ", eval, ", best_hit_rule '", best_hit_rule, "'), then dN/dS by method '",
                paste(dnds_est.method, collapse = "', '"), "' with the COAL engine (aa_aln_scoring '", aa_aln_scoring,
                "', ", aa_sub_matrix, ", comp_cores = ", comp_cores, ").")

        # the search step works in a folder of its own, removed when the call ends
        end_tmpdir <- orthologr_tmpdir_start("orthologr_dNdS_")
        on.exit(end_tmpdir(), add = TRUE)

        # the two input files are read and translated once for all steps of this call
        end_seq_cache <- seq_cache_start()
        on.exit(end_seq_cache(), add = TRUE)

        # as dNdS(): sequences that are left out are collected for attr(<result>, "dropped_sequences"), and
        # several transcripts per gene in an input file give a warning
        end_dropped <- dropped_start(quiet = quiet)
        on.exit(end_dropped(), add = TRUE)
        isoforms_query <- isoform_warning(query_file, "query_file")
        isoform_warning(subject_file, "subject_file")

        detail_message(quiet, "Starting Orthology Inference ...")

        hit.table <- dnds_search(query_file, subject_file, aligner = aligner, ortho_detection = ortho_detection,
                                 aligner_params = aligner_params, sensitivity_mode = sensitivity_mode,
                                 aligner_path = aligner_path, delete_corrupt_cds = delete_corrupt_cds,
                                 comp_cores = comp_cores, seq_type = seq_type, eval = eval, format = format,
                                 best_hit_rule = best_hit_rule, quiet = quiet)

        # (both files were read and translated for the search: they are not read again)
        q_cds <- read_cds_once(query_file, delete_corrupt_cds = delete_corrupt_cds)
        s_cds <- read_cds_once(subject_file, delete_corrupt_cds = delete_corrupt_cds)

        # the gene pairs dNdS() would process: hits whose two genes have a coding sequence that passes the
        # filters of the search step (cds2aa())
        q_ok <- cds_protein_table_file(query_file, delete_corrupt_cds = delete_corrupt_cds)$geneids
        s_ok <- cds_protein_table_file(subject_file, delete_corrupt_cds = delete_corrupt_cds)$geneids
        pairs <- unique(as.data.frame(hit.table)[ , c("query_id", "subject_id")])
        pairs <- pairs[is.element(pairs$query_id, q_ok) & is.element(pairs$subject_id, s_ok), , drop = FALSE]

        if (nrow(pairs) == 0)
                orthologr_abort("No orthologs could be found! Please check your input files!", class = "no_hits")

        detail_message(quiet, "Orthology Inference Completed: ", nrow(pairs), " gene pairs.")
        detail_message(quiet, "Starting dN/dS Estimation (COAL, all gene pairs in one call) ...")

        coal_tbl <- do.call(coalr::coal_dnds_pairs, c(list(
                query_cds   = cds_vector(q_cds),
                subject_cds = cds_vector(s_cds),
                pairs       = pairs,
                methods     = dnds_est.method,
                threads     = comp_cores,
                matrix      = coal_matrix_arg(aa_sub_matrix, aa_aln_scoring),
                preset      = aa_aln_scoring
        ), coal_flavour_args(coalr::coal_dnds_pairs, dnds_est.method)))
        # (the attribute "coal" keeps COAL's own values; the result holds them as orthologr's methods define them)
        coal_raw <- coal_tbl
        coal_tbl <- coal_orthologr_values(coal_tbl)

        # (.coal_status: what COAL said about the row; used for the column 'status' and removed)
        dNdS_tbl <- data.table::data.table(
                query_id   = coal_tbl$seq1,
                subject_id = coal_tbl$seq2,
                dN         = coal_tbl$dN,
                dS         = coal_tbl$dS,
                dNdS       = coal_tbl$dNdS,
                method     = coal_tbl$method,
                .coal_status = as.character(coal_tbl$status)
        )
        data.table::setkeyv(dNdS_tbl, c("query_id", "subject_id"))
        dNdS_tbl <- tibble::as_tibble(dNdS_tbl)

        failed <- coal_tbl[coal_tbl$status == "error", c("seq1", "seq2", "note"), drop = FALSE]
        failed <- unique(data.frame(query_id = failed$seq1, subject_id = failed$seq2, reason = failed$note,
                                    stringsAsFactors = FALSE))
        if (nrow(failed) > 0) {
                orthologr_warn(
                        class = "failed_pairs",
                        nrow(failed), " of ", nrow(pairs),
                        " gene pairs could not be processed; their dN, dS and dNdS values are NA. ",
                        "The reasons are in attr(<result>, \"failed_pairs\"). First one: ",
                        "'", failed$query_id[1], "' / '", failed$subject_id[1], "': ", failed$reason[1],
                        call. = FALSE
                )
        }

        # as dNdS(): Comeron's and Li's method have no 'method' column (one method only: several methods keep it)
        if (length(dnds_est.method) == 1 && is.element(dnds_est.method, dnds_methods(kind = "own")))
                dNdS_tbl$method <- NULL

        # joined with the hit table by the pair (query_id and subject_id), as dNdS() does
        res <- dplyr::inner_join(dNdS_tbl, tibble::as_tibble(hit.table), by = c("query_id", "subject_id"))

        # The end of the run, as in dNdS(): the status of every row (last column), the attributes
        # "dropped_sequences", "failed_pairs", "coal" (COAL's own table) and "orthologr_run", one summary message.
        engine_status <- res[[".coal_status"]]
        res[[".coal_status"]] <- NULL
        trimmed <- trimmed_table()
        ambiguous <- ambiguous_table()
        stops <- stop_codon_table()
        res <- dnds_finish(
                res, failed_pairs = failed, dropped = end_dropped(), ambiguous = ambiguous, engine_status = engine_status, trimmed = trimmed,
                stops = stops,
                record = list(fun = "dNdS_coal", args = run_args, engine = "coal",
                              tools = c(run_tool_versions(aligner = aligner, aligner_path = aligner_path), coal_versions()),
                              files = c(query_file, subject_file),
                              settings = coal_engine_settings("coal_dnds_pairs", coalr::coal_dnds_pairs, list(
                                      methods = dnds_est.method, preset = aa_aln_scoring,
                                      matrix = coal_matrix_arg(aa_sub_matrix, aa_aln_scoring)))),
                print_citation = print_citation, extra_attr = list(coal = coal_raw)
        )
        isoform_result_note(res, isoforms_query)

        return(res)
}

#' @title Compute dNdS values for all gene pairs within sequence clusters
#' @description Estimates dN, dS and dN/dS for the pairs of genes in each cluster of homologous sequences, e.g. the
#' clusters of \code{\link{deepclust}} or \code{\link{deepclust_realign}}, and summarises each cluster. For every
#' cluster the proteins are aligned (one multiple alignment per cluster with Clustal Omega or MAFFT, or pairwise),
#' the codons are threaded onto the alignment and Ka/Ks is estimated for the pairs of the cluster by COAL
#' (\url{https://github.com/Plumat-COAL/COAL}), in process and on \code{comp_cores} threads.
#' @param cds_file a character vector with the path(s) to the CDS file(s) (fasta) holding the coding sequences of all cluster members.
#' @param clusters the clusters: the tibble returned by \code{\link{deepclust}}, \code{\link{deepclust_realign}} or
#' \code{\link{deepclust_annotate}} (first two columns: \code{representative_id}, \code{member_id}), any data.frame
#' with these two columns first, or the path to such a tab-separated file. If \code{clusters = NULL} (Default),
#' \code{\link{deepclust}} is run on \code{cds_file}.
#' @param dnds_est.method the dN/dS estimation method(s): \code{"YN"} (Default), \code{"NG"}, \code{"LWL"},
#' \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"} or \code{"MYN"}. \code{"Comeron"}, \code{"Li"} and \code{"GY"} with a
#' \pkg{coalr} that has them (see \code{\link{dNdS_coal}}; Comeron and Li as \code{\link{substitutionrate}} returns
#' them).
#' @param aa_aln_tool how the proteins of clusters with three or more members are aligned:
#' \itemize{
#' \item \code{"clustalo"} (Default): one Clustal Omega multiple alignment per cluster;
#' \item \code{"mafft"}: one MAFFT multiple alignment per cluster;
#' \item \code{"pairwise"}: every pair is aligned on its own (global, BLOSUM62) by COAL; no external program is needed.
#' }
#' Clusters of two members are always aligned pairwise by COAL.
#' @param aa_aln_scoring how the pairwise protein alignments COAL makes itself are scored (\code{aa_aln_tool =
#' "pairwise"}, clusters of two members, large clusters): \code{"pwalign"} (Default; the scoring of
#' \code{\link{pairwise_aln}} and of \code{dNdS(engine = "coal")}: gap opening 10, extension 4, end gaps
#' penalised) or \code{"needle"} (COAL's own defaults: extension 0.5, end gaps free; the scoring this function
#' used before the argument existed). The same argument as in \code{\link{dNdS_coal}}.
#' @param aa_aln_path a character string specifying the path to the Clustal Omega or MAFFT executable (in case it is not on the PATH).
#' @param max_cluster_size clusters with more usable members than this are not computed for all pairs (the number of
#' pairs grows with the square of the cluster size), but as \code{large_clusters} says. Default is \code{200}.
#' @param large_clusters what to do with clusters above \code{max_cluster_size}: \code{"representative"} (Default; the
#' representative against every other member), \code{"sample"} (\code{pair_sample} pairs drawn reproducibly with
#' \code{seed}) or \code{"skip"}. Such clusters are aligned pairwise.
#' @param pair_sample the number of pairs drawn per large cluster for \code{large_clusters = "sample"}.
#' @param seed the seed of the draw for \code{large_clusters = "sample"}.
#' @param delete_corrupt_cds a logical value indicating whether coding sequences whose length is not a multiple of 3 shall be left out.
#' @param comp_cores the number of threads (COAL, DIAMOND, and the multiple alignment program).
#' @param path a character string specifying the path to the DIAMOND executable (only used when \code{clusters = NULL}).
#' @param ... additional arguments passed to \code{\link{deepclust}} when \code{clusters = NULL}
#' (e.g. \code{approx_id}, \code{mutual_cover}).
#' @details This function needs the R package \pkg{coalr}:
#' \code{remotes::install_github("Plumat-COAL/COAL/packaging/r/coalr@v0.6.0")} (COAL's release 0.6.0; a newer
#' release works as well).
#' A gene, pair or cluster that cannot be used never ends the run: see the \code{status} and \code{note} columns and
#' the \code{failed} table.
#' @author Chloe Newton and Keanu Murati
#' @return A list of three tibbles:
#' \itemize{
#' \item \code{pairs}: one row per gene pair and method: \code{cluster_id}, \code{query_id}, \code{subject_id},
#' \code{method}, \code{dN}, \code{dS}, \code{dNdS}, standard errors, numbers of sites and differences,
#' \code{codons_compared}, \code{prot_aln_length}, \code{prot_identity}, \code{status} and \code{note};
#' \item \code{summary}: one row per cluster and method: \code{cluster_id}, \code{n_members}, \code{n_used},
#' \code{method}, \code{pair_selection}, \code{aligner}, \code{n_pairs}, \code{n_ok}, \code{n_not_converged},
#' \code{n_saturated}, \code{n_no_data}, \code{n_failed} (the pairs by their \code{status} in the table \code{pairs},
#' see \code{\link{dNdS}}, section "Status of a gene pair"; up to 0.6.0 COAL's own counts), \code{n_dN}, \code{n_dS}, \code{n_dNdS} (the numbers of pairs with a dN, a dS
#' and a dN/dS value), and median, mean, quartiles (\code{_q25}, \code{_q75}), minimum and maximum of dN, dS and
#' dN/dS. \strong{Which pairs enter the statistics:} every statistic of dN (dS, dN/dS) is taken over all pairs of the
#' cluster that have a dN (dS, dN/dS) value, whatever their status - a pair with status \code{"saturated"} that
#' has a dN counts with it. \code{n_dN}, \code{n_dS} and \code{n_dNdS} are the numbers behind the statistics
#' (\code{n_ok} counts the pairs with status \code{"ok"}). \code{\link{dNdS_orthogroups}} uses the same rules and the
#' same column names (the names of \pkg{coalr}'s summary).
#' The column \code{status} of the summary: \code{"ok"} (at least one pair has an estimate), \code{"failed"} (no
#' estimate, and a pair failed or the cluster could not be aligned), \code{"saturated"} / \code{"no_data"} /
#' \code{"not_converged"} (no estimate because the pairs are saturated, have no codons to compare or did not
#' converge; \code{note} gives the numbers), or \code{"skipped"} (a cluster above \code{max_cluster_size}).
#' (\pkg{coalr}'s own summary takes the statistics over the pairs with status \code{"ok"} only.)
#' \item \code{failed}: genes and clusters that could not be used, and why.
#' }
#' @examplesIf nzchar(system.file(package = "coalr"))
#' cds <- system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr')
#' if (nzchar(Sys.which("diamond"))) {
#'   # cluster with DIAMOND DeepClust, then dN/dS within each cluster
#'   res <- dNdS_clusters(cds, aa_aln_tool = "pairwise")
#'   res$pairs
#'   res$summary
#' }
#'
#' \dontrun{
#' # with clusters computed before, one Clustal Omega alignment per cluster (needs Clustal Omega)
#' cl <- deepclust(cds, seq_type = "cds")
#' res <- dNdS_clusters(cds, clusters = cl, aa_aln_tool = "clustalo", comp_cores = 2)
#' }
#' @seealso \code{\link{deepclust}}, \code{\link{deepclust_realign}}, \code{\link{dNdS_coal}}, \code{\link{dNdS}}
#' @export
dNdS_clusters <- function(cds_file,
                          clusters           = NULL,
                          dnds_est.method    = "YN",
                          aa_aln_tool        = "clustalo",
                          aa_aln_path        = NULL,
                          max_cluster_size   = 200,
                          large_clusters     = "representative",
                          pair_sample        = 1000,
                          seed               = 1,
                          delete_corrupt_cds = FALSE,
                          comp_cores         = 1,
                          path               = NULL,
                          aa_aln_scoring     = "pwalign",
                          ...) {

        missing_files <- cds_file[!file.exists(cds_file)]
        if (length(missing_files) > 0)
                orthologr_abort(class = "bad_input", "The following input file(s) were not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
        if (!is.element(aa_aln_tool, c("clustalo", "mafft", "pairwise")))
                orthologr_abort(class = "bad_argument", "aa_aln_tool must be \"clustalo\", \"mafft\" or \"pairwise\".", call. = FALSE)
        if (!is.element(large_clusters, c("representative", "sample", "skip")))
                orthologr_abort(class = "bad_argument", "large_clusters must be \"representative\", \"sample\" or \"skip\".", call. = FALSE)
        check_aa_aln_scoring(aa_aln_scoring)
        coal_backend(need = "coalr")
        coal_check_methods(dnds_est.method, "dNdS_clusters()", other = "dNdS() or dNdS_orthogroups(engine = \"orthologr\")")
        coal_legacy_warning("dNdS_clusters()")

        # the arguments of this call, for the record of the run (attr(<result>, "orthologr_run"))
        run_args <- run_arguments()

        if (is.null(clusters)) {
                message("Clustering the sequences with DIAMOND DeepClust ...")
                clusters <- deepclust(input_file = cds_file, seq_type = "cds", delete_corrupt_cds = delete_corrupt_cds,
                                      path = path, comp_cores = comp_cores, ...)
        } else if (is.character(clusters) && length(clusters) == 1) {
                if (!file.exists(clusters))
                        orthologr_abort(class = "bad_input", "The clusters file '", clusters, "' was not found.", call. = FALSE)
        } else if (!is.data.frame(clusters) || ncol(clusters) < 2) {
                orthologr_abort(class = "bad_argument", "'clusters' must be the table returned by deepclust() / deepclust_realign() (representative_id, member_id), ",
                     "a file of that table, or NULL.", call. = FALSE)
        }
        if (is.data.frame(clusters))
                clusters <- as.data.frame(clusters)[ , 1:2]

        geneids <- seqs <- NULL
        # (unname: c() on a named list would prefix every gene id with the name of its file, e.g.
        # cds_file = c(A = "a.fasta", B = "b.fasta"), and no member would be found)
        cds <- do.call(c, lapply(unname(cds_file), function(f) {
                cds_vector(read.cds(file = f, format = "fasta", delete_corrupt_cds = delete_corrupt_cds))
        }))

        msa_aligner <- if (aa_aln_tool != "pairwise") coal_msa_aligner_for(aa_aln_tool, aa_aln_path) else NULL

        message("Estimating dN/dS within ", length(unique(clusters[[1]])), " clusters with COAL ...")
        res <- do.call(coalr::coal_dnds_clusters, c(list(
                cds              = cds,
                clusters         = clusters,
                aligner          = msa_aligner,
                methods          = dnds_est.method,
                threads          = comp_cores,
                max_cluster_size = max_cluster_size,
                large_clusters   = large_clusters,
                pair_sample      = pair_sample,
                seed             = seed,
                # pairs that COAL aligns itself (aa_aln_tool = "pairwise", clusters of two, large clusters)
                preset           = aa_aln_scoring
        ), coal_flavour_args(coalr::coal_dnds_clusters, dnds_est.method)))

        pairs <- coal_orthologr_values(res$pairs)
        names(pairs)[names(pairs) == "seq1"] <- "query_id"
        names(pairs)[names(pairs) == "seq2"] <- "subject_id"
        # the column 'status' in the words used by all dN/dS functions (see ?dNdS, "Status of a gene pair"; COAL's
        # own words are "ok", "warning", "not_converged", "saturated", "no_data", "error")
        has_status <- nrow(pairs) > 0 && all(is.element(c("dN", "dS", "dNdS", "status"), names(pairs)))
        status <- if (has_status) dnds_status(pairs$dN, pairs$dS, pairs$dNdS, engine_status = pairs$status) else NULL
        # the summary by the rules of dNdS_orthogroups(): the statistics over all pairs that have the value, the
        # counts by the status of the pairs as the table 'pairs' gives it (coalr counts COAL's own statuses)
        res$summary <- dnds_cluster_summary_stats(res$summary, pairs)
        if (has_status)
                res$summary <- dnds_cluster_status_counts(res$summary, pairs, status)
        res$summary <- dnds_cluster_status(res$summary)

        if (nrow(res$failed) > 0)
                message(nrow(res$failed), " gene(s) or cluster(s) could not be used; see the 'failed' table of the result.")

        # one summary message, the record of the run
        if (has_status) {
                if (dnds_status_column())
                        pairs$status <- status
                message(dnds_summary_text(status, dS = pairs$dS, fun = "dNdS_clusters()",
                                          what = paste0("gene pairs in ", length(unique(pairs$cluster_id)), " clusters"),
                                          failed_where = "the column 'note' of the table 'pairs'"))
        }

        out <- list(
                pairs   = tibble::as_tibble(pairs),
                summary = tibble::as_tibble(res$summary),
                failed  = tibble::as_tibble(res$failed)
        )
        attr(out, "orthologr_run") <- run_record(
                fun = "dNdS_clusters", args = run_args, engine = "coal",
                tools = c(if (is.null(run_args$clusters)) run_tool_versions(aligner = "diamond", aligner_path = path),
                          if (aa_aln_tool != "pairwise") run_tool_versions(aa_aln_tool = aa_aln_tool, aa_aln_path = aa_aln_path),
                          coal_versions()),
                files = cds_file, status = status,
                settings = coal_engine_settings("coal_dnds_clusters", coalr::coal_dnds_clusters, list(
                        methods = dnds_est.method, preset = aa_aln_scoring), mode = "clusters"))
        out
}

# The multiple alignment program of coalr (coalr::coal_msa_aligner()) for dNdS_clusters() and
# dNdS_orthogroups(engine = "coal"). `aa_aln_path`: NULL (the program on the PATH), the program, or the folder that
# holds it - there it is found as every other program of orthologr is (tool_command(): on Windows also
# '<tool>.exe'; a folder used to be turned into '<folder>/clustalo', which does not exist on Windows).
coal_msa_aligner_for <- function(aa_aln_tool, aa_aln_path = NULL) {
        exe <- if (is.null(aa_aln_path)) NULL else if (dir.exists(aa_aln_path)) tool_command(aa_aln_tool, aa_aln_path) else aa_aln_path
        tryCatch(
                coalr::coal_msa_aligner(aa_aln_tool, path = exe, threads = 1),
                error = function(e) orthologr_abort(class = "tool_missing", conditionMessage(e),
                                                    "\nOr use aa_aln_tool = \"pairwise\" (no external program needed).",
                                                    call. = FALSE)
        )
}

# The status of a cluster in the summary of dNdS_clusters(), from the numbers of its pairs. COAL calls a cluster
# "failed" when none of its pairs has an estimate, also when every pair is only saturated (E. coli K-12 + O157: 37
# clusters "failed" with an empty note and n_failed = 0). Here "failed" means that a pair failed or the cluster
# could not be aligned; a cluster whose pairs are saturated (or have no codons to compare) without any estimate is
# "saturated" ("no_data"), and the note says so. Clusters with an estimate, without pairs, or skipped keep their status.
dnds_cluster_status <- function(summary) {
        need <- c("status", "n_pairs", "n_ok", "n_saturated", "n_no_data", "n_failed")
        if (!is.data.frame(summary) || nrow(summary) == 0 || !all(is.element(need, names(summary))))
                return(summary)
        if (!is.element("note", names(summary)))
                summary$note <- NA_character_
        n0 <- function(x) { x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; x }
        n_pairs <- n0(summary$n_pairs); n_ok <- n0(summary$n_ok); n_sat <- n0(summary$n_saturated)
        n_nd <- n0(summary$n_no_data); n_fail <- n0(summary$n_failed)
        n_nc <- if (is.element("n_not_converged", names(summary))) n0(summary$n_not_converged) else integer(nrow(summary))
        none <- n_pairs > 0 & n_ok == 0 & summary$status %in% "failed"
        empty_note <- is.na(summary$note) | !nzchar(trimws(summary$note))
        words <- function(i) {
                parts <- c(if (n_sat[i] > 0) paste(n_sat[i], "saturated"), if (n_nd[i] > 0) paste(n_nd[i], "without codons to compare"),
                           if (n_nc[i] > 0) paste(n_nc[i], "not converged"), if (n_fail[i] > 0) paste(n_fail[i], "failed"))
                paste0("no pair has an estimate (", n_pairs[i], " pair", if (n_pairs[i] != 1) "s", ": ",
                       paste(parts, collapse = ", "), ")")
        }
        for (i in which(none)) {
                if (n_fail[i] == 0) {
                        summary$status[i] <- if (n_sat[i] > 0) "saturated" else if (n_nd[i] > 0) "no_data" else
                                if (n_nc[i] > 0) "not_converged" else summary$status[i]
                }
                if (empty_note[i])
                        summary$note[i] <- words(i)
        }
        summary
}

# The rule by which dNdS_clusters() and dNdS_orthogroups() summarise a group, in one place:
#   every statistic of dN (dS, dN/dS) is taken over ALL pairs of the group that have a dN (dS, dN/dS) value, whatever
#   their status - a pair that is saturated at synonymous sites still counts with its dN. The columns n_dN, n_dS
#   and n_dNdS give the numbers of pairs behind the statistics.
# Returns a data.table: the columns `by`, n_dN, n_dS, n_dNdS and <value>_<statistic> for the statistics asked for.
dnds_value_stats <- function(pairs, by, stats = c("median", "mean", "q25", "q75")) {
        dN <- dS <- dNdS <- NULL
        pairs <- data.table::as.data.table(pairs)
        fun <- list(
                median = function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE),
                mean   = function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE),
                q1     = function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE)),
                q3     = function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE)),
                q25    = function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE)),
                q75    = function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE)),
                min    = function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE),
                max    = function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
        )[stats]
        pairs[ , {
                out <- list(n_dN = sum(!is.na(dN)), n_dS = sum(!is.na(dS)), n_dNdS = sum(!is.na(dNdS)))
                for (v in c("dN", "dS", "dNdS"))
                        for (f in names(fun))
                                out[[paste0(v, "_", f)]] <- fun[[f]](as.numeric(get(v)))
                out
        }, by = by]
}

# The number of pairs by status for each group of `keys` (a data.frame of the grouping columns, one row per pair, e.g.
# cluster_id and method), `status` the status of every pair (dnds_status()): a data.table with the columns of `keys`
# and n_ok, n_dS_zero, n_saturated, n_no_data, n_not_converged, n_failed - the counts of the summaries of
# dNdS_clusters() and dNdS_orthogroups(), in the words of the column 'status' of their pair tables.
dnds_status_count_table <- function(keys, status) {
        .status <- NULL
        d <- data.table::as.data.table(keys)
        by <- names(d)
        d[[".status"]] <- status
        d[ , list(n_ok = sum(.status == "ok"), n_dS_zero = sum(.status == "dS_zero"),
                  n_saturated = sum(.status == "saturated"), n_no_data = sum(.status == "no_data"),
                  n_not_converged = sum(.status == "not_converged"), n_failed = sum(.status == "failed")), by = by]
}

# The counts by status of the summary of coalr::coal_dnds_clusters() (n_ok, n_not_converged, n_saturated, n_no_data,
# n_failed; those it has) from the status of the pairs as dNdS_clusters() returns them (dnds_status()), so that
# n_ok means what it means in dNdS_orthogroups(): the pairs with status "ok" in the table 'pairs'.
dnds_cluster_status_counts <- function(summary, pairs, status) {
        summary <- as.data.frame(summary, stringsAsFactors = FALSE)
        if (nrow(summary) == 0 || !all(is.element(c("cluster_id", "method"), names(summary))) ||
            !all(is.element(c("cluster_id", "method"), names(pairs))))
                return(summary)
        new <- as.data.frame(dnds_status_count_table(as.data.frame(pairs)[ , c("cluster_id", "method")], status))
        at <- match(paste(summary$cluster_id, summary$method, sep = "\r"), paste(new$cluster_id, new$method, sep = "\r"))
        for (col in intersect(c("n_ok", "n_not_converged", "n_saturated", "n_no_data", "n_failed"), names(summary))) {
                # (a cluster without pairs in the table, e.g. one that was skipped or could not be aligned, keeps
                # coalr's counts)
                v <- suppressWarnings(as.integer(summary[[col]]))
                v[!is.na(at)] <- as.integer(new[[col]][at[!is.na(at)]])
                summary[[col]] <- v
        }
        summary
}

# The summary of coalr::coal_dnds_clusters() with its dN / dS / dN/dS statistics computed by dnds_value_stats() from
# the pairs as dNdS_clusters() returns them (coalr takes them over the pairs with status "ok" only, and before
# orthologr's treatment of Comeron and Li values), and the count columns n_dN, n_dS, n_dNdS added.
dnds_cluster_summary_stats <- function(summary, pairs) {
        summary <- as.data.frame(summary, stringsAsFactors = FALSE)
        if (nrow(summary) == 0 || nrow(pairs) == 0 || !all(is.element(c("cluster_id", "method"), names(summary))) ||
            !is.element("cluster_id", names(pairs)))
                return(summary)
        stat_cols <- grep("^(dN|dS|dNdS)_(median|mean|q1|q3|q25|q75|min|max)$", names(summary), value = TRUE)
        stats <- unique(sub("^(dN|dS|dNdS)_", "", stat_cols))
        if (length(stats) == 0)
                return(summary)
        new <- as.data.frame(dnds_value_stats(pairs, by = c("cluster_id", "method"), stats = stats))
        at <- match(paste(summary$cluster_id, summary$method, sep = "\r"), paste(new$cluster_id, new$method, sep = "\r"))
        for (col in stat_cols)
                summary[[col]] <- new[[col]][at]
        counts <- data.frame(n_dN = new$n_dN[at], n_dS = new$n_dS[at], n_dNdS = new$n_dNdS[at])
        counts[is.na(counts)] <- 0L
        # (newer versions of coalr have a column n_dNdS of their own, counted over status "ok")
        summary <- summary[ , setdiff(names(summary), names(counts)), drop = FALSE]
        first_stat <- match(stat_cols[1], names(summary))
        cbind(summary[ , seq_len(first_stat - 1L), drop = FALSE], counts, summary[ , first_stat:ncol(summary), drop = FALSE])
}

#' @name orthologr_options
#' @rdname orthologr_options
#' @section Options and environment variable of the COAL integration:
#' Read by \code{dNdS(engine = "coal")}, \code{\link{dNdS_coal}}, \code{\link{dNdS_clusters}},
#' \code{dNdS_orthogroups(engine = "coal")}, \code{codon_aln(tool = "coal")} and
#' \code{substitutionrate(kaks_tool = "coal")}.
#' \describe{
#' \item{\code{orthologr.coal_backend}}{\code{NULL} (default: the R package \pkg{coalr} when it is installed, else
#' the program \code{coal-cpp}), \code{"coalr"} or \code{"binary"}. \code{"none"} (internal, for tests) behaves as
#' if COAL were not installed.}
#' \item{\code{orthologr.coal_path}}{the full path of the program \code{coal-cpp} (default \code{""}: the
#' environment variable \code{COAL_CPP}, then the \code{PATH}).}
#' \item{environment variable \code{COAL_CPP}}{the full path of the program \code{coal-cpp}, used when it is not on
#' the \code{PATH} and \code{options(orthologr.coal_path = )} is not set.}
#' }
NULL
