# What a dN/dS run reports about itself: the outcome of every gene pair (column 'status'), one summary message at
# the end of the run, progress messages, and the record of the run (attribute "orthologr_run").

# ---- status of a gene pair -------------------------------------------------------------------------------------

# The words used in the column 'status', in the order in which they are counted.
dnds_status_levels <- function() {
        c("ok", "dS_zero", "saturated", "no_data", "failed")
}

# The status of gene pairs from their values.
#   dN, dS, dNdS  : the values of the result
#   failed        : logical, TRUE for pairs that could not be processed (a program failed; attr "failed_pairs")
#   engine_status : NULL, or what the engine said about the pair ("ok", "saturated", "no_data", "error", ...)
# Both engines are read by the same rule, so that the same pair has the same status with either:
#   failed    - the pair could not be processed
#   ok        - there is a dN/dS value
#   dS_zero   - dN and dS were estimated and dS is 0: the ratio is not defined
#   no_data   - the engine reports that the pair has no codons to compare
#   saturated - dN or dS could not be estimated (the correction for multiple substitutions is not defined)
dnds_status <- function(dN, dS, dNdS, failed = NULL, engine_status = NULL) {
        dN <- suppressWarnings(as.numeric(dN))
        dS <- suppressWarnings(as.numeric(dS))
        dNdS <- suppressWarnings(as.numeric(dNdS))
        n <- length(dNdS)
        status <- rep("saturated", n)
        if (!is.null(engine_status)) {
                engine_status <- as.character(engine_status)
                status[!is.na(engine_status) & engine_status == "no_data"] <- "no_data"
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
        is.element(paste(tbl$query_id, tbl$subject_id, sep = "\r"),
                   paste(failed_pairs$query_id, failed_pairs$subject_id, sep = "\r"))
}

# ---- summary message -------------------------------------------------------------------------------------------

# The text of the summary of a run.
#   status        : the status of every gene pair
#   dS            : their dS values (for the median)
#   dropped       : the table of input sequences left out before the search (file, id, reason), or NULL
#   not_converged : number of pairs for which the engine reported that an iteration did not converge
#   what          : the unit, e.g. "gene pairs"
dnds_summary_text <- function(status, dS = NULL, dropped = NULL, not_converged = 0, what = "gene pairs", fun = "dNdS()",
                              failed_where = "attr(<result>, \"failed_pairs\")") {
        n <- length(status)
        counts <- dnds_status_counts(status)
        parts <- c(
                paste0(counts[["ok"]], " with dN/dS"),
                if (counts[["dS_zero"]] > 0) paste0(counts[["dS_zero"]], " with dS = 0 (no ratio)"),
                if (counts[["saturated"]] > 0) paste0(counts[["saturated"]], " saturated (dN or dS not estimable)"),
                if (counts[["no_data"]] > 0) paste0(counts[["no_data"]], " without codons to compare"),
                if (counts[["failed"]] > 0) paste0(counts[["failed"]], " failed (reasons in ", failed_where, ")")
        )
        txt <- paste0(fun, ": ", n, " ", what, " - ", paste(parts, collapse = ", "), ".")
        dS <- suppressWarnings(as.numeric(dS))
        dS <- dS[!is.na(dS)]
        med <- if (length(dS) > 0) stats::median(dS) else NA_real_
        if (!is.na(med))
                txt <- paste0(txt, " Median dS: ", signif(med, 3), ".")
        if (not_converged > 0)
                txt <- paste0(txt, " ", not_converged, " not converged (see the notes of the engine).")
        if (!is.null(dropped) && NROW(dropped) > 0) {
                tab <- table(dropped$reason)
                txt <- paste0(txt, "\n", nrow(dropped), " input sequence", if (nrow(dropped) > 1) "s were" else " was",
                              " left out before the search (", paste0(tab, " x ", names(tab), collapse = "; "),
                              "). They are listed in attr(<result>, \"dropped_sequences\").")
        }
        if (dnds_saturation_hint(status, med))
                txt <- paste0(txt, "\n", saturation_hint_text(counts[["saturated"]], n, med))
        txt
}

# The end of a dN/dS run of dNdS(): the column 'status' (last column; not with
# options(orthologr.status_column = FALSE)), the attributes, one summary message, the record of the run and the
# citation note.
#   res           : the result table (query_id, subject_id, dN, dS, dNdS, ...)
#   failed_pairs  : data.frame (query_id, subject_id, reason) or NULL
#   dropped       : data.frame (file, id, reason) or NULL
#   engine_status : NULL, or the status the engine gave for every row of `res`
#   record        : list(fun, args, engine, tools, files) for run_record()
dnds_finish <- function(res, failed_pairs, dropped, engine_status = NULL, not_converged = 0, record,
                        print_citation = TRUE, extra_attr = list()) {
        status <- dnds_status(res$dN, res$dS, res$dNdS, failed = dnds_is_failed(res, failed_pairs),
                              engine_status = engine_status)
        if (dnds_status_column())
                res[["status"]] <- status
        if (!is.null(dropped))
                attr(res, "dropped_sequences") <- dropped
        attr(res, "failed_pairs") <- failed_pairs
        for (n in names(extra_attr))
                attr(res, n) <- extra_attr[[n]]
        attr(res, "orthologr_run") <- run_record(fun = record$fun, args = record$args, engine = record$engine,
                                                 tools = record$tools, files = record$files, status = status)
        message(dnds_summary_text(status, dS = res$dS, dropped = dropped, not_converged = not_converged,
                                  fun = paste0(record$fun, "()")))
        citation_message(print_citation)
        res
}

# Most pairs saturated, or the synonymous distance of the estimated pairs high: the two species are too distant.
dnds_saturation_hint <- function(status, median_dS) {
        n <- length(status)
        if (n < 10)
                return(FALSE)
        mean(status == "saturated") > 0.2 || (!is.na(median_dS) && median_dS > 1.5)
}

saturation_hint_text <- function(n_saturated, n, median_dS) {
        paste0("Note: ", n_saturated, " of ", n, " pairs (", round(100 * n_saturated / n), " %) are saturated",
               if (!is.na(median_dS)) paste0(" and the median dS of the others is ", signif(median_dS, 3)) else "",
               ". Synonymous sites have changed too often for dS to be estimated reliably: the two species are ",
               "probably too distant for pairwise dN/dS (dS above about 1.5 - 2 is not informative). ",
               "Compare more closely related species, or use dN only.")
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
        lapply(args, function(v) {
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
}

# Every option of the package with the value in force (set values, and the defaults of the ones that are not set).
orthologr_options_in_force <- function() {
        defaults <- list(
                orthologr.legacy                    = FALSE,
                orthologr.best_hit_rule             = "bitscore",
                orthologr.legacy_pairwise_aln       = FALSE,
                orthologr.legacy_comeron_saturation = FALSE,
                orthologr.legacy_kaks_zero_dn       = FALSE,
                orthologr.dS_zero                   = "drop",
                orthologr.status_column             = TRUE,
                orthologr.isoform_warning           = TRUE,
                orthologr.cluster_type              = if (.Platform$OS.type == "windows") "PSOCK" else "FORK",
                orthologr.keep_tmp                  = FALSE
        )
        out <- lapply(stats::setNames(names(defaults), names(defaults)), function(n) orthologr_option(n, defaults[[n]]))
        set <- options()
        set <- set[startsWith(names(set), "orthologr.")]
        # (the working folder of the call is not a setting of the user)
        set[["orthologr.tmpdir"]] <- NULL
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
                v$KaKs_Calculator <- {
                        line <- tool_version_line("KaKs_Calculator", "-h", kaks_calc_path, pattern = "[Vv]ersion|KaKs_Calculator")
                        if (nzchar(line)) line else "version not reported by the program"
                }
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
#   engine   : "orthologr"
#   tools    : versions (run_tool_versions())
#   files    : the input files
#   status   : the status of every gene pair, or NULL
run_record <- function(fun, args, engine, tools, files, status = NULL) {
        list(
                orthologr = package_version_text("orthologr"),
                "function" = fun,
                arguments = args,
                options   = orthologr_options_in_force(),
                engine    = engine,
                tools     = tools,
                date      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                input_files = run_input_files(files),
                status_counts = if (is.null(status)) NULL else as.list(dnds_status_counts(status))
        )
}

# ---- writing and reading a result ------------------------------------------------------------------------------

#' @title Write and read a dN/dS table together with the record of its run
#' @description \code{write_dNdS()} writes the table returned by \code{\link{dNdS}} (or the
#' \code{pairs} table of \code{\link{dNdS_orthogroups}}, or the table in
#' \code{attr(, "dNdS_tbl")} of a \code{\link{divergence_stratigraphy}} result) as a tab-separated file and, next to
#' it, the file \code{<file>.run.json} with what a plain table file loses: the record of the run
#' (\code{attr(x, "orthologr_run")}), the tables \code{attr(x, "failed_pairs")} and
#' \code{attr(x, "dropped_sequences")}, and the type of every column. \code{read_dNdS()} reads both files back.
#' @param x a dN/dS table (a \code{data.frame}), or the list returned by \code{dNdS_orthogroups()} (its
#' \code{pairs} table is written).
#' @param file the path of the tab-separated file. The record is written to \code{paste0(file, ".run.json")}.
#' @details
#' \strong{The record of a run.} Every result of \code{dNdS()},
#' \code{dNdS_orthogroups()} and \code{divergence_stratigraphy()} carries the attribute \code{"orthologr_run"}, a
#' list with
#' \describe{
#' \item{\code{orthologr}}{the version of the package;}
#' \item{\code{function}, \code{arguments}}{the function and the values of its arguments;}
#' \item{\code{options}}{every \code{orthologr.*} option with the value in force (see \code{\link{orthologr_options}});
#' this shows whether a result was computed with the settings of 0.4.3;}
#' \item{\code{engine}}{\code{"orthologr"};}
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
#' The sidecar file is JSON and needs the package \pkg{jsonlite}. The table file has a header line, tab-separated
#' columns and \code{NA} for missing values; it can be read by any program. Numbers are written with the digits
#' needed to read the same number back.
#' @return \code{write_dNdS()}: the two file names, invisibly. \code{read_dNdS()}: the table (a tibble) with the
#' attributes \code{"orthologr_run"}, \code{"failed_pairs"} and \code{"dropped_sequences"} as far as they were
#' written. A table file without a sidecar file is read with a warning.
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
#' @seealso \code{\link{dNdS}}, \code{\link{orthologr_options}}
#' @export
write_dNdS <- function(x, file) {
        need_package("jsonlite", "write_dNdS() (the record of the run is written as JSON)")
        if (!is.character(file) || length(file) != 1 || is.na(file) || !nzchar(file))
                stop("'file' must be the path of the file to write.", call. = FALSE)
        run <- attr(x, "orthologr_run")
        if (!is.data.frame(x) && is.list(x) && is.data.frame(x$pairs)) {
                tbl <- x$pairs
        } else if (is.data.frame(x)) {
                tbl <- x
        } else {
                stop("'x' must be a dN/dS table (a data.frame) or the list returned by dNdS_orthogroups().",
                     call. = FALSE)
        }
        failed <- attr(x, "failed_pairs")
        dropped <- attr(x, "dropped_sequences")
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
        data.table::fwrite(plain, file, sep = "\t", na = "NA", quote = "auto")
        side <- list(
                format   = "orthologr dN/dS table",
                table    = basename(file),
                rows     = nrow(plain),
                columns  = classes,
                run      = run,
                failed_pairs      = if (is.data.frame(failed)) failed else NULL,
                dropped_sequences = if (is.data.frame(dropped)) dropped else NULL
        )
        side_file <- paste0(file, ".run.json")
        writeLines(jsonlite::toJSON(side, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null",
                                    dataframe = "columns"), side_file)
        invisible(c(table = file, run = side_file))
}

#' @rdname write_dNdS
#' @export
read_dNdS <- function(file) {
        if (!is.character(file) || length(file) != 1 || is.na(file) || !file.exists(file))
                stop("The file '", paste(file, collapse = "', '"), "' was not found.", call. = FALSE)
        side_file <- paste0(file, ".run.json")
        side <- NULL
        if (file.exists(side_file)) {
                need_package("jsonlite", "read_dNdS() (the record of the run is stored as JSON)")
                side <- jsonlite::fromJSON(side_file, simplifyVector = TRUE, simplifyDataFrame = FALSE)
        } else {
                warning("'", side_file, "' was not found: the table is read without the record of its run.", call. = FALSE)
        }
        # every column is read as text and then given its type: R's own conversion reads a number exactly
        tbl <- data.table::fread(file, sep = "\t", header = TRUE, na.strings = "NA", colClasses = "character",
                                 data.table = FALSE, showProgress = FALSE)
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
                if (!is.null(side$failed_pairs))
                        attr(tbl, "failed_pairs") <- as_table(side$failed_pairs)
                if (!is.null(side$dropped_sequences))
                        attr(tbl, "dropped_sequences") <- as_table(side$dropped_sequences)
        }
        tbl
}

# Numbers as text that reads back as the same number.
number_text <- function(x) {
        txt <- as.character(x)
        ok <- !is.na(x) & is.finite(x)
        longer <- ok & suppressWarnings(as.numeric(txt)) != x
        txt[longer] <- sprintf("%.17g", x[longer])
        txt[is.na(x)] <- NA_character_
        txt
}

# ---- one switch for the results of 0.4.3 -----------------------------------------------------------------------

# The values of the options that give the results of orthologr 0.4.3.
orthologr_legacy_values <- function() {
        list(
                orthologr.best_hit_rule             = "evalue",
                orthologr.legacy_pairwise_aln       = TRUE,
                orthologr.legacy_comeron_saturation = TRUE,
                orthologr.legacy_kaks_zero_dn       = TRUE,
                orthologr.dS_zero                   = "drop",
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
#' aligned pairwise, what is returned for saturated gene pairs, the value of dN for pairs without nonsynonymous
#' change) and added a column \code{status} to the dN/dS tables. \code{orthologr_legacy_options()} sets all
#' options that restore the behaviour of 0.4.3 at once, so that an analysis made with 0.4.3 can be repeated.
#' @param set \code{TRUE} (default): the options are set for the session. \code{FALSE}: nothing is set; the list
#' of options is returned.
#' @details The options (each can also be set on its own; see \code{\link{orthologr_options}}):
#' \tabular{ll}{
#' \code{orthologr.best_hit_rule = "evalue"} \tab best hits by e-value, then alignment length \cr
#' \code{orthologr.legacy_pairwise_aln = TRUE} \tab pairwise protein alignments without substitution matrix, cut to the aligned range, stop codon aligned \cr
#' \code{orthologr.legacy_comeron_saturation = TRUE} \tab Comeron's method returns a finite value for saturated pairs \cr
#' \code{orthologr.legacy_kaks_zero_dn = TRUE} \tab KaKs_Calculator methods: \code{dN = NA} for a pair without nonsynonymous change \cr
#' \code{orthologr.dS_zero = "drop"} \tab \code{filter_dNdS()} / \code{divergence_stratigraphy()} leave pairs with dS = 0 out (also the default of 0.5.0) \cr
#' \code{orthologr.status_column = FALSE} \tab the dN/dS tables have no column \code{status}
#' }
#' \code{options(orthologr.legacy = TRUE)} is a shorthand with the same effect: while it is set, each of these
#' options that is not set itself has its legacy value.
#'
#' Not covered by an option, because they are arguments: T-Coffee was run with \code{-mode accurate} in 0.4.3
#' (\code{aa_aln_params = "-mode accurate"} restores it), and \code{aa_aln_params} was ignored for all alignment
#' programs (do not set it). With these settings the 20-gene and the 1,000-gene example of the package give the
#' tables of 0.4.3 (tested against a stored table; see the article "Upgrading from 0.4.3").
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
