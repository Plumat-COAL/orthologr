# Integration of COAL (https://github.com/Plumat-COAL/COAL): codon alignments and Ka/Ks estimation without Perl
# or KaKs_Calculator, and a batched dN/dS path that hands every gene pair to COAL in one call.
#
# COAL is optional. It is used in process through the R package 'coalr' when that is installed, otherwise through
# the 'coal-cpp' program (part of the Python package 'coal-align').

# The Ka/Ks methods every COAL computes (KaKs_Calculator's names).
coal_kaks_calculator_methods <- function() {
        c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN")
}

# The methods of orthologr itself that newer versions of COAL also compute: "Comeron" (as gestimator) and "Li"
# (as seqinr::kaks()).
coal_native_methods <- function() {
        c("Comeron", "Li")
}

# Does the COAL at hand compute Comeron's and Li's method? coalr says so in coal_info()$methods; the coal-cpp
# program in its help text. FALSE when there is no COAL at all. The answer is remembered for the session.
coal_has_native_methods <- function() {
        backend <- tryCatch(coal_backend(), error = function(e) NULL)
        if (is.null(backend))
                return(FALSE)
        key <- paste0("coal_native\r", backend$type, "\r", if (is.null(backend$path)) "" else backend$path)
        hit <- .orthologr_state[[key]]
        if (!is.null(hit))
                return(hit)
        res <- if (backend$type == "coalr") {
                isTRUE(tryCatch(all(is.element(coal_native_methods(), coalr::coal_info()$methods)), error = function(e) FALSE))
        } else {
                help <- tryCatch(run_tool(basename(backend$path), args = list("dnds", "--help"), path = dirname(backend$path),
                                          error = FALSE), error = function(e) NULL)
                any(grepl("comeron-flavour", c(help$stdout, help$stderr), fixed = TRUE))
        }
        .orthologr_state[[key]] <- res
        res
}

# The methods the COAL at hand computes: KaKs_Calculator's seven and, with a COAL that has them, "Comeron" and "Li".
coal_methods <- function() {
        c(coal_kaks_calculator_methods(), if (coal_has_native_methods()) coal_native_methods())
}

# Stop unless COAL computes all `methods`; `what` names the caller for the message ("engine = \"coal\"").
coal_check_methods <- function(methods, what, other = "engine = \"orthologr\"") {
        if (!is.character(methods) || length(methods) < 1 || anyNA(methods))
                stop("dnds_est.method must be the name of a method.", call. = FALSE)
        known <- c(coal_kaks_calculator_methods(), coal_native_methods())
        unknown <- setdiff(methods, known)
        if (length(unknown) > 0)
                stop(what, " computes dnds_est.method = ", paste0("\"", known, "\"", collapse = ", "), ". \"",
                     paste(unknown, collapse = "\", \""), "\" is only available with ", other, ".", call. = FALSE)
        if (any(is.element(methods, coal_native_methods())) && !coal_has_native_methods()) {
                have <- tryCatch(as.character(utils::packageVersion("coalr")), error = function(e) NULL)
                stop(what, ": dnds_est.method = \"", paste(intersect(methods, coal_native_methods()), collapse = "\", \""),
                     "\" needs a COAL that computes Comeron's and Li's method (coalr >= 0.5.0 whose ",
                     "coalr::coal_info()$methods lists \"Comeron\" and \"Li\"; ",
                     if (is.null(have)) "coalr is not installed" else paste0("the installed coalr ", have, " does not"),
                     "). Update it with remotes::install_github(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\"), ",
                     "choose one of ", paste0("\"", coal_kaks_calculator_methods(), "\"", collapse = ", "),
                     ", or use ", other, ".", call. = FALSE)
        }
        invisible(methods)
}

# The versions of COAL in use, for the record of a run.
coal_versions <- function() {
        v <- list()
        if (requireNamespace("coalr", quietly = TRUE)) {
                v$coalr <- package_version_text("coalr")
                engine <- tryCatch(coalr::coal_info()$version, error = function(e) NULL)
                if (!is.null(engine) && length(engine) == 1)
                        v$coal <- as.character(engine)
        } else {
                exe <- tryCatch(coal_backend()$path, error = function(e) NULL)
                if (!is.null(exe))
                        v$`coal-cpp` <- tool_version_line(exe, "--version")
        }
        v
}

# The arguments that make COAL's Comeron and Li give orthologr's numbers, for a function `fun` of coalr (only the
# arguments that function has - older versions of coalr have none of them):
#   comeron_flavour = "gestimator", li_flavour = "seqinr" : the pair is read as gestimator / seqinr::kaks() read it
#   terminal_stops = "drop"    : orthologr 0.5.0 translates and aligns without the stop codon;
#                    "compare" : with options(orthologr.legacy_pairwise_aln = TRUE), i.e. the alignments of 0.4.3,
#                                in which the two stop codons are an aligned column that gestimator compares
coal_flavour_args <- function(fun) {
        args <- list(
                comeron_flavour = "gestimator",
                li_flavour      = "seqinr",
                terminal_stops  = if (isTRUE(orthologr_option("orthologr.legacy_pairwise_aln", FALSE))) "compare" else "drop"
        )
        args[intersect(names(args), names(formals(fun)))]
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
                "    remotes::install_github(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\")\n",
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
                stop("options(orthologr.coal_backend) must be \"coalr\" or \"binary\".", call. = FALSE)
        if (identical(forced, "none"))
                stop(coal_install_message(), call. = FALSE)

        if (!identical(forced, "binary") && requireNamespace("coalr", quietly = TRUE))
                return(list(type = "coalr", path = NULL))

        if (identical(forced, "coalr") || identical(need, "coalr"))
                stop(
                        "This needs the R package 'coalr' (COAL's C++ engine inside R). Install it with:\n",
                        "    remotes::install_github(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\")",
                        call. = FALSE
                )

        exe <- coal_binary()
        if (is.null(exe))
                stop(coal_install_message(), call. = FALSE)
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
                aln <- seqinr::read.alignment(file_aln, format = "clustal", forceToLower = FALSE)
                rows <- structure(toupper(unlist(aln$seq)), names = sub("\\s.*$", "", aln$nam))
        }
        cds <- read_fasta_strings(file_nuc)
        # PAL2NAL pairs the two files by position; so does orthologr (both files hold the pair in the same order)
        if (length(rows) != length(cds))
                stop("COAL: the protein alignment holds ", length(rows), " sequences and the CDS file ", length(cds),
                     " for ", context, ".", call. = FALSE)
        codon <- tryCatch(
                coalr::coal_codon_align(unname(rows), unname(cds)),
                error = function(e) stop("COAL could not build the codon alignment for ", context, ": ",
                                         conditionMessage(e), call. = FALSE)
        )
        names(codon) <- names(cds)
        write_alignment(codon, tmp.out, format = format)
        publish_file(tmp.out, file.out)
        invisible(file.out)
}

# Ka/Ks of a two-sequence codon alignment file with COAL: the table substitutionrate() returns for the
# KaKs_Calculator methods (query_id, subject_id, dN, dS, dNdS, method).
coal_kaks_file <- function(file, est.method, subst_name = NULL) {
        query_id <- method <- NULL
        coal_check_methods(est.method, "kaks_tool = \"coal\"", other = "kaks_tool = \"kaks_calculator\"")
        backend <- coal_backend()
        context <- pair_context(file, subst_name)
        seqs <- read_fasta_strings(file)
        if (length(seqs) != 2)
                stop("The codon alignment of ", context, " must hold two sequences, but '", file, "' holds ",
                     length(seqs), ".", call. = FALSE)

        if (backend$type == "coalr") {
                k <- do.call(coalr::coal_kaks, c(list(seqs[[1]], seqs[[2]], method = est.method),
                                                 coal_flavour_args(coalr::coal_kaks)))
                if (identical(k$status[1], "error"))
                        stop("COAL could not estimate dN/dS for ", context, ": ", k$note[1], call. = FALSE)
                k <- coal_orthologr_values(k)
                res <- data.table::data.table(query_id = names(seqs)[1], subject_id = names(seqs)[2],
                                              dN = k$dN[1], dS = k$dS[1], dNdS = k$dNdS[1], method = k$method[1])
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
                                c(list("dnds", axt_file, "-o", out_file, "-m", est.method, "-j", "1", "-q",
                                       "--digits", "17"),
                                  if (is.element(est.method, coal_native_methods()))
                                          list("--comeron-flavour", "gestimator", "--li-flavour", "seqinr")),
                                context = context, output = out_file)
                k <- utils::read.delim(out_file, sep = "\t", quote = "", na.strings = "NA", stringsAsFactors = FALSE)
                if (nrow(k) < 1)
                        stop("COAL returned no result for ", context, ".", call. = FALSE)
                if (identical(k$status[1], "error"))
                        stop("COAL could not estimate dN/dS for ", context, ": ", k$note[1], call. = FALSE)
                k <- coal_orthologr_values(k)
                res <- data.table::data.table(query_id = names(seqs)[1], subject_id = names(seqs)[2],
                                              dN = as.numeric(k$dN[1]), dS = as.numeric(k$dS[1]),
                                              dNdS = as.numeric(k$dNdS[1]), method = as.character(k$method[1]))
        }
        # (substitutionrate() returns no 'method' column for its own two methods)
        if (is.element(est.method, coal_native_methods()))
                res[ , method := NULL]
        data.table::setkey(res, query_id)
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
#' installed \pkg{coalr} lists them in \code{coalr::coal_info()$methods} (coalr >= 0.5.0); otherwise an error
#' says so. See Details.
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
#' Li's method. \code{attr(result, "coal")} holds COAL's unrounded values and notes. Differences between the
#' engines that remain come from the protein alignment (co-optimal alignments; see NEWS).
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
#' \code{remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")}.
#' @author Chloe Newton
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
                      sensitivity_mode   = "fast",
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

        best_hit_rule <- orthologr_best_hit_rule(if (missing(best_hit_rule)) NULL else best_hit_rule)

        if (!is.element(aligner, c("diamond", "blast")))
                stop("engine = \"coal\" detects orthologs with aligner = \"diamond\" or \"blast\".", call. = FALSE)
        if (!is.element(ortho_detection, c("BH", "RBH")))
                stop("engine = \"coal\" supports ortho_detection = \"RBH\" or \"BH\".", call. = FALSE)
        if (!identical(seq_type, "cds"))
                stop("dN/dS estimation needs coding sequences: seq_type = \"cds\".", call. = FALSE)
        check_choice(format, "fasta", "format", hint = "dNdS_coal() reads coding sequences from FASTA files.")
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
        end_dropped <- dropped_start()
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
                matrix      = aa_sub_matrix,
                preset      = aa_aln_scoring
        ), coal_flavour_args(coalr::coal_dnds_pairs)))
        # (the attribute "coal" keeps COAL's own values; the result holds them as orthologr's methods define them)
        coal_raw <- coal_tbl
        coal_tbl <- coal_orthologr_values(coal_tbl)

        # (.coal_status / .coal_note: what COAL said about the row; used for the column 'status' and removed)
        dNdS_tbl <- data.table::data.table(
                query_id   = coal_tbl$seq1,
                subject_id = coal_tbl$seq2,
                dN         = coal_tbl$dN,
                dS         = coal_tbl$dS,
                dNdS       = coal_tbl$dNdS,
                method     = coal_tbl$method,
                .coal_status = as.character(coal_tbl$status),
                .coal_note   = as.character(coal_tbl$note)
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
        if (length(dnds_est.method) == 1 && is.element(dnds_est.method, coal_native_methods()))
                dNdS_tbl$method <- NULL

        # joined with the hit table by the pair (query_id and subject_id), as dNdS() does
        res <- dplyr::inner_join(dNdS_tbl, tibble::as_tibble(hit.table), by = c("query_id", "subject_id"))

        # The end of the run, as in dNdS(): the status of every row (last column), the attributes
        # "dropped_sequences", "failed_pairs", "coal" (COAL's own table) and "orthologr_run", one summary message.
        engine_status <- res[[".coal_status"]]
        note <- res[[".coal_note"]]
        res[[".coal_status"]] <- NULL
        res[[".coal_note"]] <- NULL
        res <- dnds_finish(
                res, failed_pairs = failed, dropped = end_dropped(), engine_status = engine_status,
                not_converged = sum(grepl("not converge|did not converge|no convergence", note, ignore.case = TRUE)),
                record = list(fun = "dNdS_coal", args = run_args, engine = "coal",
                              tools = c(run_tool_versions(aligner = aligner, aligner_path = aligner_path), coal_versions()),
                              files = c(query_file, subject_file)),
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
#' \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"} or \code{"MYN"}. \code{"Comeron"} and \code{"Li"} with a \pkg{coalr} that
#' has them (see \code{\link{dNdS_coal}}; values as \code{\link{substitutionrate}} returns them).
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
#' \code{remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")}.
#' A gene, pair or cluster that cannot be used never ends the run: see the \code{status} and \code{note} columns and
#' the \code{failed} table.
#' @author Chloe Newton
#' @return A list of three tibbles:
#' \itemize{
#' \item \code{pairs}: one row per gene pair and method: \code{cluster_id}, \code{query_id}, \code{subject_id},
#' \code{method}, \code{dN}, \code{dS}, \code{dNdS}, standard errors, numbers of sites and differences,
#' \code{codons_compared}, \code{prot_aln_length}, \code{prot_identity}, \code{status} and \code{note};
#' \item \code{summary}: one row per cluster and method: \code{cluster_id}, \code{n_members}, \code{n_used},
#' \code{method}, \code{pair_selection}, \code{aligner}, \code{n_pairs}, \code{n_ok}, \code{n_saturated},
#' \code{n_no_data}, \code{n_failed}, \code{n_dN}, \code{n_dS}, \code{n_dNdS} (the numbers of pairs with a dN, a dS
#' and a dN/dS value), and median, mean, quartiles (\code{_q25}, \code{_q75}), minimum and maximum of dN, dS and
#' dN/dS. \strong{Which pairs enter the statistics:} every statistic of dN (dS, dN/dS) is taken over all pairs of the
#' cluster that have a dN (dS, dN/dS) value, whatever their status - a pair with status \code{"saturated"} that
#' has a dN counts with it. \code{n_dN}, \code{n_dS} and \code{n_dNdS} are the numbers behind the statistics
#' (\code{n_ok} counts the pairs with status \code{"ok"}). \code{\link{dNdS_orthogroups}} uses the same rule.
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
                stop("The following input file(s) were not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
        if (!is.element(aa_aln_tool, c("clustalo", "mafft", "pairwise")))
                stop("aa_aln_tool must be \"clustalo\", \"mafft\" or \"pairwise\".", call. = FALSE)
        if (!is.element(large_clusters, c("representative", "sample", "skip")))
                stop("large_clusters must be \"representative\", \"sample\" or \"skip\".", call. = FALSE)
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
                        stop("The clusters file '", clusters, "' was not found.", call. = FALSE)
        } else if (!is.data.frame(clusters) || ncol(clusters) < 2) {
                stop("'clusters' must be the table returned by deepclust() / deepclust_realign() (representative_id, member_id), ",
                     "a file of that table, or NULL.", call. = FALSE)
        }
        if (is.data.frame(clusters))
                clusters <- as.data.frame(clusters)[ , 1:2]

        geneids <- seqs <- NULL
        cds <- do.call(c, lapply(cds_file, function(f) {
                cds_vector(read.cds(file = f, format = "fasta", delete_corrupt_cds = delete_corrupt_cds))
        }))

        msa_aligner <- NULL
        if (aa_aln_tool != "pairwise") {
                exe <- if (is.null(aa_aln_path)) NULL else {
                        if (dir.exists(aa_aln_path)) file.path(aa_aln_path, aa_aln_tool) else aa_aln_path
                }
                msa_aligner <- tryCatch(
                        coalr::coal_msa_aligner(aa_aln_tool, path = exe, threads = 1),
                        error = function(e) stop(conditionMessage(e),
                                                 "\nOr use aa_aln_tool = \"pairwise\" (no external program needed).",
                                                 call. = FALSE)
                )
        }

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
        ), coal_flavour_args(coalr::coal_dnds_clusters)))

        pairs <- coal_orthologr_values(res$pairs)
        names(pairs)[names(pairs) == "seq1"] <- "query_id"
        names(pairs)[names(pairs) == "seq2"] <- "subject_id"
        # the statistics of the summary by the rule of dNdS_orthogroups(): over all pairs that have the value
        res$summary <- dnds_cluster_summary_stats(res$summary, pairs)

        if (nrow(res$failed) > 0)
                message(nrow(res$failed), " gene(s) or cluster(s) could not be used; see the 'failed' table of the result.")

        # the column 'status' in the words used by all dN/dS functions (see ?dNdS, "Status of a gene pair"; COAL's
        # own words are "ok", "saturated", "no_data", "warning", "error"), one summary message, the record of the run
        if (nrow(pairs) > 0 && all(is.element(c("dN", "dS", "dNdS", "status"), names(pairs)))) {
                not_converged <- if (is.element("note", names(pairs)))
                        sum(grepl("not converge|did not converge|no convergence", pairs$note, ignore.case = TRUE)) else 0
                status <- dnds_status(pairs$dN, pairs$dS, pairs$dNdS, engine_status = pairs$status)
                if (dnds_status_column())
                        pairs$status <- status
                message(dnds_summary_text(status, dS = pairs$dS, not_converged = not_converged, fun = "dNdS_clusters()",
                                          what = paste0("gene pairs in ", length(unique(pairs$cluster_id)), " clusters"),
                                          failed_where = "the column 'note' of the table 'pairs'"))
        } else {
                status <- NULL
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
                files = cds_file, status = status)
        out
}

# The rule by which dNdS_clusters() and dNdS_orthogroups() summarise a group, in one place:
#   every statistic of dN (dS, dN/dS) is taken over ALL pairs of the group that have a dN (dS, dN/dS) value, whatever
#   their status - a pair that is saturated at synonymous sites still counts with its dN. The columns n_dN, n_dS
#   and n_dNdS give the numbers of pairs behind the statistics.
# Returns a data.table: the columns `by`, n_dN, n_dS, n_dNdS and <value>_<statistic> for the statistics asked for.
dnds_value_stats <- function(pairs, by, stats = c("median", "mean", "q1", "q3")) {
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
#' the program \code{coal-cpp}), \code{"coalr"} or \code{"binary"}.}
#' \item{\code{orthologr.coal_path}}{the full path of the program \code{coal-cpp} (default \code{""}: the
#' environment variable \code{COAL_CPP}, then the \code{PATH}).}
#' \item{environment variable \code{COAL_CPP}}{the full path of the program \code{coal-cpp}, used when it is not on
#' the \code{PATH} and \code{options(orthologr.coal_path = )} is not set.}
#' }
NULL
