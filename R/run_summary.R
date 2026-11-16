# What a dN/dS run reports about itself: the outcome of every gene pair (column 'status'), one summary message at
# the end of the run, progress messages, and the record of the run (attribute "orthologr_run").

# ---- status of a gene pair -------------------------------------------------------------------------------------

# The words used in the column 'status', in the order in which they are counted.
dnds_status_levels <- function() {
        c("ok", "dS_zero", "saturated", "no_data", "not_converged", "failed")
}

# The status of gene pairs from their values.
#   dN, dS, dNdS  : the values of the result
#   failed        : logical, TRUE for pairs that could not be processed (a program failed; attr "failed_pairs")
#   engine_status : NULL, or what the engine said about the pair (COAL: "ok", "warning", "not_converged",
#                   "saturated", "no_data", "error"; NA where the engine said nothing)
# Both engines are read by the same rule, so that the same pair has the same status with either:
#   failed        - the pair could not be processed
#   ok            - there is a dN/dS value
#   dS_zero       - dN and dS were estimated and dS is 0: the ratio is not defined
#   no_data       - the engine reports that the pair has no codons to compare
#   not_converged - the engine reports that its iteration (YN, MYN) or optimisation (GY) did not converge, so
#                   there is no estimate (COAL; KaKs_Calculator does not report it)
#   saturated     - dN or dS could not be estimated (the correction for multiple substitutions is not defined)
dnds_status <- function(dN, dS, dNdS, failed = NULL, engine_status = NULL) {
        dN <- suppressWarnings(as.numeric(dN))
        dS <- suppressWarnings(as.numeric(dS))
        dNdS <- suppressWarnings(as.numeric(dNdS))
        n <- length(dNdS)
        status <- rep("saturated", n)
        if (!is.null(engine_status)) {
                engine_status <- as.character(engine_status)
                status[!is.na(engine_status) & engine_status == "no_data"] <- "no_data"
                status[!is.na(engine_status) & engine_status == "not_converged"] <- "not_converged"
        }
        status[!is.na(dN) & !is.na(dS) & dS == 0] <- "dS_zero"
        status[is.finite(dNdS)] <- "ok"
        if (!is.null(engine_status))
                status[!is.na(engine_status) & is.element(engine_status, c("error", "failed"))] <- "failed"
        if (!is.null(failed))
                status[!is.na(failed) & failed] <- "failed"
        status
}

# Is the column 'status' part of the tables? (options(orthologr.status_column = FALSE), or the legacy switch)
dnds_status_column <- function() {
        isTRUE(orthologr_option("orthologr.status_column", TRUE))
}

# Counts by status, every level present.
dnds_status_counts <- function(status) {
        tab <- table(factor(status, levels = dnds_status_levels()))
        stats::setNames(as.integer(tab), names(tab))
}

# Which rows of a pair table are in a table of failed pairs (query_id, subject_id)?
dnds_is_failed <- function(tbl, failed_pairs) {
        if (is.null(failed_pairs) || NROW(failed_pairs) == 0)
                return(rep(FALSE, nrow(tbl)))
        pair <- paste(tbl$query_id, tbl$subject_id, sep = "\r")
        fp_pair <- paste(failed_pairs$query_id, failed_pairs$subject_id, sep = "\r")
        # with several methods per pair (dnds_est.method = "ALL") a failed estimate of one method names the method
        if (is.element("method", names(failed_pairs)) && is.element("method", names(tbl))) {
                whole <- is.na(failed_pairs$method)
                return(is.element(pair, fp_pair[whole]) |
                               is.element(paste(pair, tbl$method, sep = "\r"),
                                          paste(fp_pair[!whole], failed_pairs$method[!whole], sep = "\r")))
        }
        is.element(pair, fp_pair)
}

# The reason of every row of a pair table from the table of failed pairs (NA for rows that did not fail); with a
# column 'method' in both (dnds_est.method = "ALL") a failed estimate of one method is the reason of its row only.
dnds_failed_reason <- function(tbl, failed_pairs) {
        out <- rep(NA_character_, nrow(tbl))
        if (is.null(failed_pairs) || NROW(failed_pairs) == 0)
                return(out)
        pair <- paste(tbl$query_id, tbl$subject_id, sep = "\r")
        fp_pair <- paste(failed_pairs$query_id, failed_pairs$subject_id, sep = "\r")
        if (is.element("method", names(failed_pairs)) && is.element("method", names(tbl))) {
                whole <- is.na(failed_pairs$method)
                out <- failed_pairs$reason[!whole][match(paste(pair, tbl$method, sep = "\r"),
                                                          paste(fp_pair[!whole], failed_pairs$method[!whole], sep = "\r"))]
                w <- match(pair, fp_pair[whole])
                out[!is.na(w)] <- failed_pairs$reason[whole][w[!is.na(w)]]
                return(out)
        }
        failed_pairs$reason[match(pair, fp_pair)]
}

# ---- summary message -------------------------------------------------------------------------------------------

# The text of the summary of a run.
#   status        : the status of every gene pair
#   dS            : their dS values (for the median)
#   dropped       : the table of input sequences left out before the search (file, id, reason), or NULL
#   what          : the unit, e.g. "gene pairs"
#   species_pair  : NULL (two species, dNdS()), or for every gene pair the species it compares (e.g. "A / B"):
#                   the note on saturation then names the species pairs with the highest median dS
dnds_summary_text <- function(status, dS = NULL, dropped = NULL, what = "gene pairs", fun = "dNdS()",
                              failed_where = "attr(<result>, \"failed_pairs\")", trimmed = NULL, ambiguous_rows = 0,
                              ambiguous_where = "attr(<result>, \"ambiguous_sequences\")", stops = NULL,
                              species_pair = NULL) {
        n <- length(status)
        counts <- dnds_status_counts(status)
        parts <- c(
                paste0(counts[["ok"]], " with dN/dS"),
                if (counts[["dS_zero"]] > 0) paste0(counts[["dS_zero"]], " with dS = 0 (no ratio)"),
                if (counts[["saturated"]] > 0) paste0(counts[["saturated"]], " saturated (dN or dS not estimable)"),
                if (counts[["no_data"]] > 0) paste0(counts[["no_data"]], " without codons to compare"),
                if (counts[["not_converged"]] > 0) paste0(counts[["not_converged"]], " not converged (no estimate)"),
                if (counts[["failed"]] > 0) paste0(counts[["failed"]], " failed (reasons in ", failed_where, ")")
        )
        # (one gene pair: "1 gene pair")
        if (n == 1 && identical(what, "gene pairs"))
                what <- "gene pair"
        txt <- paste0(fun, ": ", n, " ", what, " - ", paste(parts, collapse = ", "), ".")
        dS <- suppressWarnings(as.numeric(dS))
        dS_all <- dS
        dS <- dS[!is.na(dS)]
        med <- if (length(dS) > 0) stats::median(dS) else NA_real_
        if (!is.na(med))
                txt <- paste0(txt, " Median dS: ", signif(med, 3), ".")
        if (!is.null(dropped) && NROW(dropped) > 0) {
                tab <- table(dropped$reason)
                txt <- paste0(txt, "\n", nrow(dropped), " input sequence", if (nrow(dropped) > 1) "s were" else " was",
                              " left out before the search (", paste0(tab, " x ", names(tab), collapse = "; "),
                              "). They are listed in attr(<result>, \"dropped_sequences\").")
        }
        if (!is.null(trimmed) && NROW(trimmed) > 0) {
                txt <- paste0(txt, "\n", nrow(trimmed), " coding sequence", if (nrow(trimmed) > 1) "s have" else " has",
                              " a length that is not a multiple of 3 and ", if (nrow(trimmed) > 1) "were" else "was",
                              " trimmed at the 3' end (delete_corrupt_cds = FALSE; a sequence incomplete at its 5' end is ",
                              "then read in the wrong frame). They are listed in attr(<result>, \"trimmed_sequences\"); ",
                              "as read: ", paste0("'", unique(trimmed$corrupted_file), "'", collapse = ", "),
                              if (any(in_tempdir(unique(trimmed$corrupted_file))))
                                      paste0(" (inside tempdir(), which R removes when the session ends; ",
                                             "options(orthologr.corrupted_cds_dir = \"<folder>\") keeps the file)"),
                              ".")
        }
        if (ambiguous_rows > 0)
                txt <- paste0(txt, "\n", ambiguous_rows, if (ambiguous_rows > 1) " rows have" else " row has",
                              " a coding sequence with ambiguous bases (e.g. N; options(orthologr.ambiguous_cds = \"keep\"), ",
                              "the default): codons with such a base were not compared (see ", ambiguous_where, "); ",
                              "options(orthologr.ambiguous_cds = \"drop\") leaves them out, as 0.4.3 did.")
        if (!is.null(stops) && NROW(stops) > 0)
                txt <- paste0(txt, "\n", nrow(stops), " coding sequence", if (nrow(stops) > 1) "s have" else " has",
                              " in-frame stop codons before the last codon (a pseudogene, a wrong frame or another ",
                              "genetic code?); they were used as they are. They are listed in ",
                              "attr(<result>, \"stop_codon_sequences\").")
        if (dnds_saturation_hint(status, med))
                txt <- paste0(txt, "\n", saturation_hint_text(counts[["saturated"]], n, med,
                                                              species_pair = species_pair, status = status, dS = dS_all))
        txt
}

# Which rows of a pair table involve a coding sequence kept with ambiguous bases (table of ambiguous_table())?
dnds_ambiguous_rows <- function(tbl, ambiguous) {
        if (is.null(ambiguous) || NROW(ambiguous) == 0 || nrow(tbl) == 0)
                return(rep(FALSE, nrow(tbl)))
        is.element(as.character(tbl$query_id), ambiguous$id) | is.element(as.character(tbl$subject_id), ambiguous$id)
}

# The end of a dN/dS run of dNdS() / dNdS_coal(): the column 'status' (last column; not with
# options(orthologr.status_column = FALSE)), the attributes, one summary message, the record of the run and the
# citation note.
#   res           : the result table (query_id, subject_id, dN, dS, dNdS, ...)
#   failed_pairs  : data.frame (query_id, subject_id, reason) or NULL
#   dropped       : data.frame (file, id, reason) or NULL
#   engine_status : NULL, or the status the engine gave for every row of `res`
#   record        : list(fun, args, engine, tools, files) for run_record()
#   trimmed       : data.frame (file, id, removed_nt, corrupted_file) of the sequences trimmed to a multiple of 3,
#                   or NULL
dnds_finish <- function(res, failed_pairs, dropped, engine_status = NULL, record,
                        print_citation = TRUE, extra_attr = list(), trimmed = NULL, ambiguous = NULL, stops = NULL) {
        status <- dnds_status(res$dN, res$dS, res$dNdS, failed = dnds_is_failed(res, failed_pairs),
                              engine_status = engine_status)
        if (dnds_status_column())
                res[["status"]] <- status
        if (!is.null(dropped))
                attr(res, "dropped_sequences") <- dropped
        attr(res, "failed_pairs") <- failed_pairs
        if (!is.null(trimmed) && nrow(trimmed) > 0)
                attr(res, "trimmed_sequences") <- trimmed
        # options(orthologr.ambiguous_cds = "keep"): the sequences kept with ambiguous bases, and how many of their
        # codons hold such a base (those codons are not compared)
        if (!is.null(ambiguous) && nrow(ambiguous) > 0)
                attr(res, "ambiguous_sequences") <- ambiguous
        # coding sequences with in-frame stop codons (file, id, stop_codons); their values are unchanged
        if (!is.null(stops) && nrow(stops) > 0)
                attr(res, "stop_codon_sequences") <- stops
        for (n in names(extra_attr))
                attr(res, n) <- extra_attr[[n]]
        attr(res, "orthologr_run") <- run_record(fun = record$fun, args = record$args, engine = record$engine,
                                                 tools = record$tools, files = record$files, status = status,
                                                 settings = record$settings)
        # several methods per gene pair (dnds_est.method = "ALL"): the rows are counted as gene pairs x methods
        what <- "gene pairs"
        if (is.element("method", names(res)) && length(unique(stats::na.omit(res$method))) > 1)
                what <- local({
                        n_pairs <- nrow(unique(as.data.frame(res)[, c("query_id", "subject_id")]))
                        paste0("rows (", n_pairs, if (n_pairs == 1) " gene pair x " else " gene pairs x ",
                               length(unique(stats::na.omit(res$method))), " methods)")
                })
        message(dnds_summary_text(status, dS = res$dS, dropped = dropped,
                                  fun = paste0(record$fun, "()"), trimmed = trimmed, what = what,
                                  ambiguous_rows = sum(dnds_ambiguous_rows(res, ambiguous)), stops = stops))
        citation_message(print_citation)
        res
}

# Which of the paths are inside tempdir()?
in_tempdir <- function(path) {
        tmp <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
        p <- normalizePath(path, winslash = "/", mustWork = FALSE)
        startsWith(p, paste0(tmp, "/"))
}

# Most pairs saturated, or the synonymous distance of the estimated pairs high: the species are too distant.
dnds_saturation_hint <- function(status, median_dS) {
        n <- length(status)
        if (n < 10)
                return(FALSE)
        isTRUE(mean(status == "saturated", na.rm = TRUE) > 0.2) || isTRUE(!is.na(median_dS) && median_dS > 1.5)
}

# The note on saturation.
#   species_pair : NULL (the two species of dNdS()), or the species pair of every gene pair (with `status` and
#                  `dS` of every gene pair): with more than one species pair the note names the (up to three)
#                  species pairs with the highest median dS, and the share of their saturated gene pairs
saturation_hint_text <- function(n_saturated, n, median_dS, species_pair = NULL, status = NULL, dS = NULL) {
        groups <- if (!is.null(species_pair) && length(species_pair) == n && length(status) == n && length(dS) == n)
                unique(stats::na.omit(as.character(species_pair))) else character(0)
        who <- if (length(groups) > 1) {
                sp <- as.character(species_pair)
                med <- vapply(groups, function(g) {
                        d <- dS[!is.na(sp) & sp == g & !is.na(dS)]
                        if (length(d) == 0) NA_real_ else stats::median(d)
                }, 0)
                sat <- vapply(groups, function(g) mean(status[!is.na(sp) & sp == g] == "saturated"), 0)
                # (species pairs with saturated gene pairs only have no median: they come first)
                top <- order(-ifelse(is.na(med), Inf, med), -sat, groups)[seq_len(min(3, length(groups)))]
                paste0("some of the species compared are probably too distant for pairwise dN/dS (highest median dS: ",
                       paste0(groups[top], " (", ifelse(is.na(med[top]), "-", signif(med[top], 3)), ", ",
                              round(100 * sat[top]), " % saturated)", collapse = "; "), "; dS above about 1.5 - 2 is ",
                       "not informative). Compare more closely related species, or use dN only.")
        } else {
                paste0(if (length(groups) == 1) paste0("the two species (", groups, ") are") else "the two species are",
                       " probably too distant for pairwise dN/dS (dS above about 1.5 - 2 is not informative). ",
                       "Compare more closely related species, or use dN only.")
        }
        paste0("Note: ", n_saturated, " of ", n, " pairs (", round(100 * n_saturated / n), " %) are saturated",
               if (!is.na(median_dS)) paste0(" and the median dS of the others is ", signif(median_dS, 3)) else "",
               ". Synonymous sites have changed too often for dS to be estimated reliably: ", who)
}

# ---- messages ------------------------------------------------------------------------------------------------

# A message that is only given when the caller is not quiet.
detail_message <- function(quiet, ...) {
        if (!isTRUE(quiet))
                message(...)
        invisible(NULL)
}

# Run `expr`; when quiet, the progress chatter of the search step (program banners, "creating a database",
# timings) is not shown. Every other message, and every warning, passes.
quiet_search <- function(expr, quiet = TRUE) {
        if (!isTRUE(quiet))
                return(expr)
        chatter <- paste0("^(Running (diamond|blast|DIAMOND|BLAST)|sensitivity mode: |creating a (diamond|blast) database|",
                          "Starting DIAMOND2 search|DIAMOND2 search completed|\\s*$)")
        withCallingHandlers(expr, message = function(m) {
                if (grepl(chatter, conditionMessage(m)))
                        invokeRestart("muffleMessage")
        })
}

# The citation, once per R session.
citation_message <- function(print_citation = TRUE) {
        if (!isTRUE(print_citation) || isTRUE(.orthologr_state$citation_shown))
                return(invisible(FALSE))
        .orthologr_state$citation_shown <- TRUE
        message("Please cite: Drost et al. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns ",
                "in Animal and Plant Embryogenesis. 2015. Mol. Biol. Evol. 32 (5): 1221-1231. ",
                "(This note is shown once per session; print_citation = FALSE switches it off.)")
        invisible(TRUE)
}

# How many gene pairs lie between two progress messages: about a tenth of the run, at least 50 pairs.
progress_step <- function(n) {
        max(50L, as.integer(ceiling(n / 10)))
}

progress_message <- function(done, n, t0, what = "gene pairs") {
        message("dN/dS: ", done, " of ", n, " ", what, " (", format_seconds(proc.time()[["elapsed"]] - t0), ")")
}

format_seconds <- function(s) {
        if (s < 120)
                return(paste0(round(s), " s"))
        if (s < 7200)
                return(paste0(round(s / 60, 1), " min"))
        paste0(round(s / 3600, 1), " h")
}

# ---- the record of a run ---------------------------------------------------------------------------------------

# The arguments of the calling function as a list of plain values (for the record of the run): long or
# structured values are described, not stored.
run_arguments <- function(env = parent.frame(), fun = sys.function(sys.parent())) {
        args <- mget(setdiff(names(formals(fun)), "..."), envir = env, ifnotfound = list(NULL))
        out <- lapply(args, function(v) {
                if (is.null(v))
                        return(NULL)
                if (is.function(v))
                        return("<function>")
                if (is.matrix(v))
                        return(paste0("<", nrow(v), " x ", ncol(v), " matrix>"))
                if (is.data.frame(v))
                        return(paste0("<table with ", nrow(v), " rows>"))
                if (is.atomic(v) && length(v) <= 50)
                        return(if (is.factor(v)) as.character(v) else unname(v))
                if (is.atomic(v))
                        return(paste0("<", typeof(v), " vector of length ", length(v), ">"))
                paste0("<", class(v)[1], ">")
        })
        # where the run starts in the log of search commands (run_record()'s "commands")
        attr(out, "command_mark") <- length(.orthologr_state$command_log)
        out
}

# The search commands (DIAMOND, BLAST+) of the session, as run (run_tool() adds them).
log_command <- function(cmd, args) {
        .orthologr_state$command_log <- c(.orthologr_state$command_log, list(c(cmd, quote_args(args))))
        invisible(NULL)
}

# The user's call: the outermost call of a function of the package that records its runs (dNdS(engine = "coal")
# runs dNdS_coal(), whose record then says how it was called).
RECORDED_FUNCTIONS <- c("dNdS", "dNdS_coal", "dNdS_clusters", "dNdS_orthogroups", "divergence_stratigraphy")
run_user_call <- function() {
        calls <- sys.calls()
        frames <- sys.frames()
        fns <- lapply(seq_along(calls), function(i) sys.function(i))
        ns <- asNamespace("orthologr")
        for (i in seq_along(calls)) {
                nm <- tryCatch(as.character(calls[[i]][[1]]), error = function(e) "")
                nm <- nm[length(nm)]
                if (length(nm) == 1 && is.element(nm, RECORDED_FUNCTIONS) &&
                    identical(environment(fns[[i]]), ns)) {
                        return(list("function" = nm,
                                    call = paste(deparse(calls[[i]], width.cutoff = 500L), collapse = " "),
                                    arguments = run_arguments(env = frames[[i]], fun = fns[[i]])))
                }
        }
        NULL
}

# Every option of the package with the value in force (set values, and the defaults of the ones that are not set).
orthologr_options_in_force <- function() {
        defaults <- list(
                orthologr.legacy                    = FALSE,
                orthologr.best_hit_rule             = "bitscore",
                orthologr.sensitivity_mode          = "default",
                orthologr.legacy_pairwise_aln       = FALSE,
                orthologr.legacy_comeron_saturation = FALSE,
                orthologr.comeron_twosubs           = "comeron",
                orthologr.comeron_min_sites         = 3,
                orthologr.legacy_kaks_zero_dn       = FALSE,
                orthologr.legacy_kaks_saturation    = FALSE,
                orthologr.max_aln_cells             = 5e8,
                orthologr.large_pair_aligner        = "coal",
                orthologr.dS_zero                   = "drop",
                orthologr.ambiguous_cds             = "keep",
                orthologr.status_column             = TRUE,
                orthologr.isoform_warning           = TRUE,
                orthologr.cluster_type              = if (.Platform$OS.type == "windows") "PSOCK" else "FORK",
                orthologr.keep_tmp                  = FALSE
        )
        out <- lapply(stats::setNames(names(defaults), names(defaults)), function(n) orthologr_option(n, defaults[[n]]))
        set <- options()
        set <- set[startsWith(names(set), "orthologr.")]
        # (the working folder of the call is not a setting of the user; whether the note of library(orthologr) is
        # shown, and whether PAL2NAL / KaKs_Calculator handle a chunk of gene pairs per process, change no result)
        set[["orthologr.tmpdir"]] <- NULL
        set[["orthologr.quiet_upgrade"]] <- NULL
        set[["orthologr.batch_tools"]] <- NULL
        for (n in setdiff(names(set), names(out)))
                out[[n]] <- set[[n]]
        out[order(names(out))]
}

# The first line a program prints for `args` ("" when it cannot be started). Asked once per session and program.
tool_version_line <- function(program, args, path = NULL, pattern = NULL) {
        cmd <- tool_command(program, path)
        key <- paste(c("version", cmd, args), collapse = " ")
        if (!is.null(.orthologr_state[[key]]))
                return(.orthologr_state[[key]])
        out <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE, timeout = 20)),
                        error = function(e) character(0))
        out <- trimws(out[nzchar(trimws(out))])
        if (!is.null(pattern) && any(grepl(pattern, out)))
                out <- out[grepl(pattern, out)]
        res <- if (length(out) == 0) "" else substr(out[1], 1, 200)
        .orthologr_state[[key]] <- res
        res
}

# KaKs_Calculator's version: `-h` prints "Program: KaKs_Calculator Toolbox" and, on the next line,
# "Version: 2.0, June. 2009"; the line with the version, else the program line.
kaks_version_line <- function(kaks_calc_path = NULL) {
        line <- tool_version_line("KaKs_Calculator", "-h", kaks_calc_path, pattern = "[Vv]ersion")
        if (grepl("[Vv]ersion", line))
                return(sub("^[Vv]ersion: *", "", line))
        line <- tool_version_line("KaKs_Calculator", "-h", kaks_calc_path, pattern = "KaKs_Calculator")
        if (nzchar(line)) line else "version not reported by the program"
}

package_version_text <- function(pkg) {
        tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

# The versions of the programs and packages a run used.
#   aligner : "diamond" / "blast" / NULL;  aa_aln_tool : "NW", "pairwise", a multiple alignment program or NULL
#   kaks    : TRUE when KaKs_Calculator computed the values
run_tool_versions <- function(aligner = NULL, aligner_path = NULL, aa_aln_tool = NULL, aa_aln_path = NULL,
                              codon_aln_tool = NULL, kaks = FALSE, kaks_calc_path = NULL) {
        v <- list()
        if (identical(aligner, "diamond"))
                v$diamond <- tool_version_line("diamond", "version", aligner_path)
        if (identical(aligner, "blast"))
                v$blast <- tool_version_line("blastp", "-version", aligner_path)
        if (!is.null(aa_aln_tool)) {
                if (is.element(aa_aln_tool, c("NW", "pairwise"))) {
                        v$pwalign <- package_version_text("pwalign")
                } else {
                        v[[aa_aln_tool]] <- switch(
                                aa_aln_tool,
                                clustalo = tool_version_line("clustalo", "--version", aa_aln_path),
                                mafft    = tool_version_line("mafft", "--version", aa_aln_path),
                                muscle   = tool_version_line("muscle", "-version", aa_aln_path),
                                clustalw = tool_version_line("clustalw2", "-help", aa_aln_path, pattern = "CLUSTAL"),
                                t_coffee = tool_version_line("t_coffee", "-version", aa_aln_path, pattern = "[Vv]ersion"),
                                NA_character_)
                }
        }
        if (identical(codon_aln_tool, "pal2nal")) {
                v$pal2nal <- "v14 (shipped with orthologr)"
                v$perl <- tool_version_line("perl", "-v", pattern = "This is perl")
        }
        if (isTRUE(kaks))
                v$KaKs_Calculator <- kaks_version_line(kaks_calc_path)
        v$Biostrings <- package_version_text("Biostrings")
        v$seqinr <- package_version_text("seqinr")
        v$R <- paste(R.version$major, R.version$minor, sep = ".")
        v
}

# Size and MD5 checksum of the input files.
run_input_files <- function(files) {
        files <- as.character(unlist(files))
        files <- files[!is.na(files) & nzchar(files)]
        if (length(files) == 0)
                return(data.frame(file = character(0), size = numeric(0), md5 = character(0), stringsAsFactors = FALSE))
        ok <- file.exists(files) & !dir.exists(files)
        md5 <- rep(NA_character_, length(files))
        md5[ok] <- unname(tools::md5sum(files[ok]))
        data.frame(file = normalizePath(files, mustWork = FALSE), size = as.numeric(file.size(files)), md5 = md5,
                   stringsAsFactors = FALSE)
}

# The record of a run: a list that says which version, arguments, options, programs and input files gave the
# result, and how the gene pairs came out.
#   fun      : the name of the function
#   args     : its arguments (run_arguments())
#   engine   : "orthologr" or "coal"
#   tools    : versions (run_tool_versions())
#   files    : the input files
#   status   : the status of every gene pair, or NULL
#   settings : what the engine was given (dNdS_coal(): the arguments of coalr), or NULL
run_record <- function(fun, args, engine, tools, files, status = NULL, settings = NULL) {
        mark <- attr(args, "command_mark")
        attr(args, "command_mark") <- NULL
        log <- .orthologr_state$command_log
        commands <- if (!is.null(mark) && length(log) > mark)
                vapply(log[(mark + 1):length(log)], paste, "", collapse = " ") else NULL
        user <- run_user_call()
        out <- list(
                orthologr = package_version_text("orthologr"),
                "function" = fun,
                arguments = args,
                call      = user$call,
                called_as = if (!is.null(user) && !identical(user[["function"]], fun))
                        user[c("function", "arguments")] else NULL,
                commands  = commands,
                engine_settings = settings,
                options   = orthologr_options_in_force(),
                engine    = engine,
                tools     = tools,
                date      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                input_files = run_input_files(files),
                status_counts = if (is.null(status)) NULL else as.list(dnds_status_counts(status))
        )
        out[!(vapply(out, is.null, TRUE) & is.element(names(out), c("call", "called_as", "commands", "engine_settings")))]
}

# ---- writing and reading a result ------------------------------------------------------------------------------

#' @title Write and read a dN/dS table together with the record of its run
#' @description \code{write_dNdS()} writes the table returned by \code{\link{dNdS}} (or \code{\link{dNdS_coal}}, the
#' \code{pairs} table of \code{\link{dNdS_clusters}} / \code{\link{dNdS_orthogroups}}, or the table in
#' \code{attr(, "dNdS_tbl")} of a \code{\link{divergence_stratigraphy}} result) as a tab-separated file and, next to
#' it, the file \code{<file>.run.json} with what a plain table file loses: the record of the run
#' (\code{attr(x, "orthologr_run")}), the tables of the run in the attributes \code{"failed_pairs"},
#' \code{"dropped_sequences"}, \code{"trimmed_sequences"}, \code{"stop_codon_sequences"},
#' \code{"ambiguous_sequences"} and \code{"coal"}, and the type of every column. \code{read_dNdS()} reads both
#' files back. \code{dnds_equal()} compares two results without the parts of the record that differ between two
#' runs of the same analysis.
#' @param x a dN/dS table (a \code{data.frame}), or the list returned by \code{dNdS_clusters()} /
#' \code{dNdS_orthogroups()} (its \code{pairs} table is written).
#' @param file the path of the tab-separated file. The record is written to \code{paste0(file, ".run.json")}.
#' @param y a second result, compared with \code{x} by \code{dnds_equal()}.
#' @param ignore the parts of the record of the run (\code{attr(, "orthologr_run")}, see Details) that
#' \code{dnds_equal()} leaves out of the comparison: by default \code{"date"}, \code{"commands"} (the command
#' lines of the searches name files of the working folder of the call, whose names are new for every call) and
#' \code{"call"}. \code{character(0)}: everything is compared, as \code{identical()} does.
#' @details
#' \strong{The record of a run.} Every result of \code{dNdS()}, \code{dNdS_coal()}, \code{dNdS_clusters()},
#' \code{dNdS_orthogroups()} and \code{divergence_stratigraphy()} carries the attribute \code{"orthologr_run"}, a
#' list with
#' \describe{
#' \item{\code{orthologr}}{the version of the package;}
#' \item{\code{function}, \code{arguments}}{the function and the values of its arguments;}
#' \item{\code{call}}{the call as it was written (the outermost of \code{dNdS()}, \code{dNdS_coal()},
#' \code{dNdS_clusters()}, \code{dNdS_orthogroups()}, \code{divergence_stratigraphy()});}
#' \item{\code{called_as}}{when that call is of another function than \code{function} (\code{dNdS(engine = "coal")}
#' runs \code{dNdS_coal()}): its name and the values of all its arguments, also of those the engine does not use;}
#' \item{\code{commands}}{the command lines of the sequence searches the run started (DIAMOND, BLAST+), with the
#' file names of the working folder;}
#' \item{\code{engine_settings}}{with \code{engine = "coal"}: what orthologr gave COAL (the arguments of
#' \pkg{coalr}) and the equivalent \code{coal pipeline} command line (see \code{\link{dNdS_coal}});}
#' \item{\code{options}}{every \code{orthologr.*} option with the value in force (see \code{\link{orthologr_options}});
#' this shows whether a result was computed with the settings of 0.4.3;}
#' \item{\code{engine}}{\code{"orthologr"} or \code{"coal"};}
#' \item{\code{tools}}{the versions of the programs and packages used (DIAMOND or BLAST+, the alignment program or
#' \pkg{pwalign}, PAL2NAL and Perl, KaKs_Calculator, \pkg{Biostrings}, \pkg{seqinr}, R);}
#' \item{\code{date}}{when the run ended;}
#' \item{\code{input_files}}{path, size and MD5 checksum of the input files;}
#' \item{\code{status_counts}}{the number of gene pairs by \code{status} (see \code{\link{dNdS}}, section
#' "Status of a gene pair").}
#' }
#' Attributes are lost by most functions that change a table (\code{dplyr} verbs, \code{merge()}, subsetting a
#' \code{data.frame}) and by \code{write.table()} / \code{readr::write_tsv()}: save the result of a run with
#' \code{write_dNdS()} (or \code{saveRDS()}) before working with it.
#'
#' \strong{Comparing two results.} Two runs of the same analysis give the same table and the same record of the
#' run except for the date, the file names of the working folder in the command lines of the searches, and - when
#' it was written differently - the call; \code{identical()} of the two results is therefore \code{FALSE}.
#' \code{dnds_equal(x, y)} is \code{identical()} without these parts of the record (\code{ignore}); it also
#' compares the results of \code{dNdS_clusters()} / \code{dNdS_orthogroups()} (lists) and of
#' \code{divergence_stratigraphy()}. To see what differs, compare the parts, e.g.
#' \code{all.equal(as.data.frame(x), as.data.frame(y))} and \code{attr(x, "orthologr_run")$tools}.
#'
#' The sidecar file is JSON and needs the package \pkg{jsonlite}. The table file has a header line, tab-separated
#' columns and \code{NA} for missing values; it can be read by any program. Numbers are written with the digits
#' needed to read the same number back (15 or 17 significant digits; a number that R would not read back exactly
#' from 17 digits, as happens on some platforms for very small e-values, in the hexadecimal notation of C99, e.g.
#' \code{0x1.7p-715}, which R and C's \code{strtod()} read exactly).
#' @return \code{write_dNdS()}: the two file names, invisibly. \code{read_dNdS()}: the table (a tibble) with the
#' attributes \code{"orthologr_run"}, \code{"failed_pairs"}, \code{"dropped_sequences"}, \code{"trimmed_sequences"},
#' \code{"stop_codon_sequences"}, \code{"ambiguous_sequences"} and \code{"coal"} as far as they were written (files
#' written by orthologr 0.6.0 and earlier hold the first three only). A table file without a sidecar file is read
#' with a warning. \code{dnds_equal()}: \code{TRUE} or \code{FALSE}.
#' @examples
#' tbl <- data.frame(query_id = c("a", "b"), subject_id = c("x", "y"), dN = c(0.01, 0.02),
#'                   dS = c(0.1, 0.2), dNdS = c(0.1, 0.1), status = "ok")
#' attr(tbl, "orthologr_run") <- list(orthologr = "0.5.0", engine = "orthologr")
#' file <- tempfile(fileext = ".tsv")
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'         write_dNdS(tbl, file)
#'         back <- read_dNdS(file)
#'         attr(back, "orthologr_run")$orthologr
#' }
#' # two runs of the same analysis differ in the date of their records
#' tbl2 <- tbl
#' attr(tbl, "orthologr_run")$date <- "2026-11-16 10:00:00 GMT"
#' attr(tbl2, "orthologr_run")$date <- "2026-11-16 10:05:00 GMT"
#' identical(tbl, tbl2)
#' dnds_equal(tbl, tbl2)
#' @seealso \code{\link{dNdS}}, \code{\link{orthologr_options}}
#' @export
write_dNdS <- function(x, file) {
        need_package("jsonlite", "write_dNdS() (the record of the run is written as JSON)")
        if (!is.character(file) || length(file) != 1 || is.na(file) || !nzchar(file))
                orthologr_abort(class = "bad_argument", "'file' must be the path of the file to write.")
        if (!is.data.frame(x) && is.list(x) && is.data.frame(x$pairs)) {
                tbl <- x$pairs
        } else if (is.data.frame(x)) {
                tbl <- x
        } else {
                orthologr_abort(class = "bad_argument", "'x' must be a dN/dS table (a data.frame) or the list returned by ",
                                "dNdS_clusters() / dNdS_orthogroups().")
        }
        # (an attribute of the list returned by dNdS_clusters() / dNdS_orthogroups(), else of its table)
        attr_of <- function(name) {
                v <- attr(x, name)
                if (is.null(v)) attr(tbl, name) else v
        }
        run <- attr_of("orthologr_run")
        plain <- as.data.frame(tbl, stringsAsFactors = FALSE)
        for (col in names(plain))
                if (is.factor(plain[[col]]))
                        plain[[col]] <- as.character(plain[[col]])
        classes <- lapply(plain, function(v) class(v)[1])
        # numbers as text with the digits needed to read the same number back (15 significant digits where that
        # is enough, else 17)
        for (col in names(plain))
                if (is.double(plain[[col]]))
                        plain[[col]] <- number_text(plain[[col]])
        # (the text is written byte for byte; the encoding marks of non-ASCII text go into the record, see below)
        data.table::fwrite(plain, file, sep = "\t", na = "NA", quote = "auto")
        side <- list(
                format   = "orthologr dN/dS table",
                table    = basename(file),
                rows     = nrow(plain),
                columns  = classes,
                text_encoding = text_marks(plain),
                run      = run
        )
        # every table of the run kept in an attribute: its values (doubles as text that reads back as the same
        # number), the encoding marks of its text, the type of its columns and its class
        tables <- list()
        for (a in dnds_table_attributes()) {
                v <- attr_of(a)
                if (!is.data.frame(v))
                        next
                p <- as.data.frame(v, stringsAsFactors = FALSE)
                for (col in names(p))
                        if (is.factor(p[[col]]))
                                p[[col]] <- as.character(p[[col]])
                cls <- lapply(p, function(col) class(col)[1])
                for (col in names(p))
                        if (is.double(p[[col]]))
                                p[[col]] <- number_text(p[[col]])
                side[[a]] <- p
                side[[paste0(a, "_encoding")]] <- text_marks(p)
                tables[[a]] <- list(columns = cls, class = class(v))
        }
        if (length(tables) > 0)
                side$tables <- tables
        side_file <- paste0(file, ".run.json")
        writeLines(jsonlite::toJSON(json_text(side), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null",
                                    na = "null", dataframe = "columns"), side_file, useBytes = TRUE)
        invisible(c(table = file, run = side_file))
}

#' @rdname write_dNdS
#' @export
read_dNdS <- function(file) {
        if (!is.character(file) || length(file) != 1 || is.na(file) || !file.exists(file))
                orthologr_abort(class = "bad_input", "The file '", paste(file, collapse = "', '"), "' was not found.", call. = FALSE)
        side_file <- paste0(file, ".run.json")
        side <- NULL
        if (file.exists(side_file)) {
                need_package("jsonlite", "read_dNdS() (the record of the run is stored as JSON)")
                side <- jsonlite::fromJSON(side_file, simplifyVector = TRUE, simplifyDataFrame = FALSE)
        } else {
                warning("'", side_file, "' was not found: the table is read without the record of its run.", call. = FALSE)
        }
        # every column is read as text and then given its type: R's own conversion reads a number exactly
        # The text is read byte for byte, as write_dNdS() wrote it, and gets back the encoding marks it had (from
        # the record; without one: unmarked, as dNdS() returns the IDs of the FASTA files). IDs read back thus equal
        # the IDs of the run in every locale (non-ASCII IDs marked as UTF-8 would be unequal to those of dNdS() in a
        # C / POSIX locale).
        tbl <- data.table::fread(file, sep = "\t", header = TRUE, na.strings = "NA", colClasses = "character",
                                 data.table = FALSE, showProgress = FALSE, encoding = "unknown")
        # write_dNdS() quotes a field that holds a quote (or a tab or a line end) and doubles the quotes inside it;
        # fread() removes the outer quotes but keeps the doubled ones ('d"q' is read as 'd""q'). Every '""' in a
        # field therefore is one quote.
        for (col in names(tbl))
                if (is.character(tbl[[col]])) {
                        v <- tbl[[col]]
                        dq <- which(grepl("\"\"", v, fixed = TRUE, useBytes = TRUE))
                        if (length(dq) > 0) {
                                v[dq] <- gsub("\"\"", "\"", v[dq], fixed = TRUE, useBytes = TRUE)
                                Encoding(v[dq]) <- "unknown"
                        }
                        tbl[[col]] <- v
                }
        tbl <- restore_text_marks(tbl, side$text_encoding, from_json = FALSE)
        classes <- if (!is.null(side$columns)) unlist(side$columns) else
                vapply(tbl, function(v) if (all(is.na(v) | !is.na(suppressWarnings(as.numeric(v))))) "numeric" else "character", "")
        for (col in intersect(names(tbl), names(classes)))
                tbl[[col]] <- switch(classes[[col]],
                                     numeric = as.numeric(tbl[[col]]), integer = as.integer(tbl[[col]]),
                                     logical = as.logical(tbl[[col]]), tbl[[col]])
        tbl <- tibble::as_tibble(tbl)
        as_table <- function(l) {
                if (is.null(l))
                        return(NULL)
                l <- lapply(l, function(v) if (length(v) == 0) character(0) else unlist(lapply(v, function(e) if (is.null(e)) NA else e)))
                as.data.frame(l, stringsAsFactors = FALSE)
        }
        if (!is.null(side)) {
                run <- side$run
                if (!is.null(run$input_files)) {
                        run$input_files <- as_table(run$input_files)
                        run$input_files$size <- as.numeric(run$input_files$size)
                }
                attr(tbl, "orthologr_run") <- run
                for (a in dnds_table_attributes()) {
                        if (is.null(side[[a]]))
                                next
                        v <- restore_text_marks(as_table(side[[a]]), side[[paste0(a, "_encoding")]], from_json = TRUE)
                        info <- side$tables[[a]]
                        # (files of 0.6.0 and earlier: failed_pairs and dropped_sequences, without types)
                        if (!is.null(info)) {
                                cls <- unlist(info$columns)
                                for (col in intersect(names(v), names(cls)))
                                        v[[col]] <- switch(cls[[col]],
                                                           numeric = as.numeric(v[[col]]), integer = as.integer(v[[col]]),
                                                           logical = as.logical(v[[col]]),
                                                           character = as.character(v[[col]]), v[[col]])
                                v <- restore_table_class(v, unlist(info$class))
                        }
                        attr(tbl, a) <- v
                }
        }
        tbl
}

# The attributes of a dN/dS result that hold tables of the run (written and read by write_dNdS() / read_dNdS()).
dnds_table_attributes <- function() {
        c("failed_pairs", "dropped_sequences", "trimmed_sequences", "stop_codon_sequences", "ambiguous_sequences", "coal")
}

# A data.frame read back with the class it had: a tibble or a data.table where it was one.
restore_table_class <- function(v, cls) {
        if (is.element("tbl_df", cls))
                return(tibble::as_tibble(v))
        if (is.element("data.table", cls))
                return(data.table::as.data.table(v))
        v
}

#' @rdname write_dNdS
#' @export
dnds_equal <- function(x, y, ignore = c("date", "commands", "call")) {
        if (!is.null(ignore) && (!is.character(ignore) || anyNA(ignore)))
                orthologr_abort(class = "bad_argument", "'ignore' must be the names of parts of the record of the run, ",
                                "e.g. c(\"date\", \"commands\", \"call\"), or character(0).")
        identical(dnds_strip_record(x, ignore), dnds_strip_record(y, ignore))
}

# `x` without the parts `ignore` of every record of a run (attribute "orthologr_run") in it: of `x`, of the
# elements of a list (dNdS_clusters(), dNdS_orthogroups()) and of tables in attributes (divergence_stratigraphy():
# "dNdS_tbl").
dnds_strip_record <- function(x, ignore, depth = 0L) {
        if (length(ignore) == 0 || depth > 3L)
                return(x)
        run <- attr(x, "orthologr_run", exact = TRUE)
        if (is.list(run))
                attr(x, "orthologr_run") <- run[setdiff(names(run), ignore)]
        for (a in setdiff(names(attributes(x)), c("orthologr_run", "names", "row.names", "class", "dim", "dimnames")))
                if (is.list(attr(x, a, exact = TRUE)))
                        attr(x, a) <- dnds_strip_record(attr(x, a, exact = TRUE), ignore, depth + 1L)
        if (is.list(x) && !is.data.frame(x))
                for (i in seq_along(x))
                        if (is.list(x[[i]]))
                                x[[i]] <- dnds_strip_record(x[[i]], ignore, depth + 1L)
        x
}

# The encoding marks of the non-ASCII text of the character columns of a table: list(<column> = list(row = <rows>,
# mark = <mark>)), NULL when all text is ASCII. <mark> is Encoding() ("unknown", "UTF-8", "latin1", "bytes"), and
# "native" for unmarked text that is not valid UTF-8 (text in the encoding of a non-UTF-8 locale). R compares text
# by its marks as well as by its bytes: in a C / POSIX locale an unmarked non-ASCII ID (as dNdS() returns the IDs
# of a FASTA file) differs from the same bytes marked as UTF-8. write_dNdS() stores the marks so that read_dNdS()
# returns the text as it was.
text_marks <- function(tbl) {
        out <- list()
        for (col in names(tbl)) {
                v <- tbl[[col]]
                if (!is.character(v))
                        next
                rows <- which(!is_ascii_text(v))
                if (length(rows) > 0) {
                        mark <- Encoding(v[rows])
                        mark[mark == "unknown" & !validUTF8(v[rows])] <- "native"
                        out[[col]] <- list(row = rows, mark = mark)
                }
        }
        if (length(out) == 0) NULL else out
}

# TRUE for NA and for text of ASCII bytes only (whatever the encoding mark and the locale).
is_ascii_text <- function(v) {
        is.na(v) | !is.na(iconv(v, from = "latin1", to = "ASCII"))
}

# The text of a list (the record written by write_dNdS()) prepared for jsonlite::toJSON(): unmarked non-ASCII
# text that is valid UTF-8 is marked as UTF-8 (the bytes are not changed), so that it is written as it is (in a
# C / POSIX locale jsonlite would write such text as escapes like '<c3><bc>').
json_text <- function(x) {
        if (is.character(x)) {
                u <- which(!is_ascii_text(x) & Encoding(x) == "unknown" & validUTF8(x))
                if (length(u) > 0)
                        Encoding(x[u]) <- "UTF-8"
                return(x)
        }
        if (is.list(x)) {
                a <- attributes(x)
                x <- lapply(x, json_text)
                attributes(x) <- a
        }
        x
}

# The encoding marks recorded by text_marks() given back to the text of a table read by read_dNdS().
#   from_json : the text was read from the JSON record, where jsonlite marks all of it as UTF-8 (text in latin1 or
#               in a native non-UTF-8 encoding was converted to UTF-8 when it was written, and is converted back)
# Non-ASCII text without a record (a file written by an earlier version) is unmarked, as dNdS() returns IDs.
restore_text_marks <- function(tbl, marks, from_json = FALSE) {
        if (is.null(tbl))
                return(tbl)
        for (col in names(tbl)) {
                v <- tbl[[col]]
                if (!is.character(v))
                        next
                rows <- which(!is_ascii_text(v))
                if (length(rows) == 0)
                        next
                mark <- rep("unknown", length(rows))
                rec <- marks[[col]]
                if (!is.null(rec)) {
                        r <- as.integer(unlist(rec$row))
                        m <- as.character(unlist(rec$mark))
                        if (length(r) == length(m)) {
                                hit <- match(rows, r)
                                mark[!is.na(hit)] <- m[hit[!is.na(hit)]]
                        }
                }
                for (m in unique(mark)) {
                        i <- rows[mark == m]
                        if (from_json) {
                                v[i] <- switch(m,
                                               "UTF-8"  = v[i],
                                               latin1   = iconv(v[i], from = "UTF-8", to = "latin1"),
                                               native   = enc2native(v[i]),
                                               { x <- v[i]; Encoding(x) <- "unknown"; x })
                        } else if (is.element(m, c("UTF-8", "latin1", "bytes"))) {
                                Encoding(v[i]) <- m
                        } else {
                                Encoding(v[i]) <- "unknown"
                        }
                }
                tbl[[col]] <- v
        }
        tbl
}

# Numbers as text that reads back as the same number: 15 significant digits where that is enough, else 17, and
# where R does not read even 17 digits back exactly (R's own conversion is not exact for some numbers on platforms
# without an extended-precision long double, e.g. 1.37e-215 on macOS arm64) the hexadecimal notation of C99
# ("%a", e.g. 0x1.5p-3), which R reads exactly.
number_text <- function(x) {
        txt <- as.character(x)
        ok <- !is.na(x) & is.finite(x)
        longer <- ok & suppressWarnings(as.numeric(txt)) != x
        txt[longer] <- sprintf("%.17g", x[longer])
        hex <- which(longer & suppressWarnings(as.numeric(txt)) != x)
        # (only where R reads the hexadecimal text back exactly; not every platform does so for subnormal numbers)
        h <- sprintf("%a", x[hex])
        exact <- !is.na(suppressWarnings(as.numeric(h))) & suppressWarnings(as.numeric(h)) == x[hex]
        txt[hex[exact]] <- h[exact]
        txt[is.na(x)] <- NA_character_
        txt
}

# ---- one switch for the results of 0.4.3 -----------------------------------------------------------------------

# The values of the options that give the results of orthologr 0.4.3.
orthologr_legacy_values <- function() {
        list(
                orthologr.best_hit_rule             = "evalue",
                orthologr.sensitivity_mode          = "fast",
                orthologr.legacy_pairwise_aln       = TRUE,
                orthologr.legacy_comeron_saturation = TRUE,
                orthologr.comeron_twosubs           = "libsequence",
                orthologr.legacy_kaks_zero_dn       = TRUE,
                orthologr.legacy_kaks_saturation    = TRUE,
                orthologr.max_aln_cells             = Inf,
                orthologr.dS_zero                   = "drop",
                orthologr.ambiguous_cds             = "drop",
                orthologr.status_column             = FALSE
        )
}

# The value of an option of the package: what was set with options(); otherwise, with
# options(orthologr.legacy = TRUE), the value that gives the results of 0.4.3; otherwise `default`.
orthologr_option <- function(name, default = NULL) {
        value <- getOption(name, NULL)
        if (!is.null(value))
                return(value)
        if (isTRUE(getOption("orthologr.legacy", FALSE))) {
                legacy <- orthologr_legacy_values()
                if (is.element(name, names(legacy)))
                        return(legacy[[name]])
        }
        default
}

#' @title The settings that reproduce results of orthologr 0.4.3
#' @description Version 0.5.0 changed defaults that change results (how best hits are chosen, how proteins are
#' aligned pairwise, what is returned for saturated gene pairs, the value of dN and dS for pairs without
#' nonsynonymous or without any change) and added a column \code{status} to the dN/dS tables. \code{orthologr_legacy_options()} sets all
#' options that restore the behaviour of 0.4.3 at once, so that an analysis made with 0.4.3 can be repeated.
#' @param set \code{TRUE} (default): the options are set for the session. \code{FALSE}: nothing is set; the list
#' of options is returned.
#' @details The options (each can also be set on its own; see \code{\link{orthologr_options}}):
#' \tabular{ll}{
#' \code{orthologr.best_hit_rule = "evalue"} \tab best hits by e-value, then alignment length \cr
#' \code{orthologr.sensitivity_mode = "fast"} \tab DIAMOND's \code{--fast} mode (the default of the argument \code{sensitivity_mode}) \cr
#' \code{orthologr.legacy_pairwise_aln = TRUE} \tab pairwise protein alignments without substitution matrix, cut to the aligned range, stop codon aligned \cr
#' \code{orthologr.legacy_comeron_saturation = TRUE} \tab Comeron's method returns a finite value for saturated pairs \cr
#' \code{orthologr.comeron_twosubs = "libsequence"} \tab Comeron's method counts codons that differ at two positions along libsequence's branches \cr
#' \code{orthologr.legacy_kaks_zero_dn = TRUE} \tab KaKs_Calculator methods: \code{NA} for a dN or dS that is 0 (pairs without nonsynonymous or without any change) \cr
#' \code{orthologr.legacy_kaks_saturation = TRUE} \tab KaKs_Calculator methods: its capped values (e.g. dS = 99) for saturated estimates, and \code{nan} as \code{NaN} (status \code{"saturated"}, not \code{"failed"}) \cr
#' \code{orthologr.max_aln_cells = Inf} \tab pairwise alignments of any size (default: pairs above 5e8 cells are not aligned) \cr
#' \code{orthologr.dS_zero = "drop"} \tab \code{filter_dNdS()} / \code{divergence_stratigraphy()} leave pairs with dS = 0 out (also the default of 0.5.0) \cr
#' \code{orthologr.ambiguous_cds = "drop"} \tab coding sequences with ambiguous bases (e.g. N) are left out before the search \cr
#' \code{orthologr.status_column = FALSE} \tab the dN/dS tables have no column \code{status}
#' }
#' \code{options(orthologr.legacy = TRUE)} is a shorthand with the same effect: while it is set, each of these
#' options that is not set itself has its legacy value.
#'
#' \strong{Run time:} the pairwise alignments of 0.4.3 (\code{orthologr.legacy_pairwise_aln = TRUE}) take about 0.6 -
#' 0.8 s per gene pair, 14 to 19 times as long as the default (150 gene pairs: about 115 s against 6 s; 1,000 pairs:
#' about 10 minutes per core), as long as with 0.4.3 itself. Reproducing a genome-wide run of 0.4.3 takes hours;
#' use \code{comp_cores}.
#'
#' Not covered by an option, because they are arguments: T-Coffee was run with \code{-mode accurate} in 0.4.3
#' (\code{aa_aln_params = "-mode accurate"} restores it), and \code{aa_aln_params} was ignored for all alignment
#' programs (do not set it). With these settings the 20-gene and the 1,000-gene example of the package give the
#' tables of 0.4.3 (tested against a stored table; see the article "Upgrading from 0.4.3").
#'
#' Only \code{engine = "orthologr"} (the default) gives the numbers of 0.4.3. The COAL engine
#' (\code{dNdS(engine = "coal")}, \code{\link{dNdS_coal}}, \code{dNdS_orthogroups(engine = "coal")},
#' \code{\link{dNdS_clusters}}) makes its own pairwise alignments and cannot make those of 0.4.3; with these
#' settings it warns (class \code{orthologr_warning_legacy_engine}). 20-gene example with these settings: 10 of 20
#' pairs differ from 0.4.3, largest difference in dN/dS 0.147.
#' @return With \code{set = TRUE}: the previous values of the options, invisibly, as \code{options()} returns
#' them - \code{old <- orthologr_legacy_options(); ...; options(old)} restores the settings.
#' With \code{set = FALSE}: the named list of legacy options.
#' @examples
#' orthologr_legacy_options(set = FALSE)
#'
#' old <- orthologr_legacy_options()
#' getOption("orthologr.best_hit_rule")
#' options(old)
#' @seealso \code{\link{orthologr_options}}, \code{\link{dNdS}}
#' @export
orthologr_legacy_options <- function(set = TRUE) {
        check_flag(set, "set")
        legacy <- orthologr_legacy_values()
        if (!set)
                return(legacy)
        old <- do.call(options, legacy)
        message("orthologr: the settings of version 0.4.3 are in force (", paste(names(legacy), collapse = ", "),
                "). options(<the value returned by this call>) restores the previous settings.")
        invisible(old)
}
