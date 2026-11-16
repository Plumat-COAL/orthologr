# Internal helpers for running external programs (BLAST+, DIAMOND, PAL2NAL, KaKs_Calculator, ClustalW, ...).
#
# All external programs are started through run_tool(): arguments are quoted (file names with spaces work), the
# exit status is checked, the expected output files must exist afterwards, and the program's own messages are
# quoted in the error. Nothing here changes the working directory or the PATH of the R session.

.orthologr_state <- new.env(parent = emptyenv())

# The command to start `program`: '<path>/<program>' when a folder is given, else the bare name (looked up on PATH).
tool_command <- function(program, path = NULL) {
        if (is.null(path) || !nzchar(path))
                return(program)
        cmd <- file.path(path, program)
        if (.Platform$OS.type == "windows" && !file.exists(cmd) && file.exists(paste0(cmd, ".exe")))
                cmd <- paste0(cmd, ".exe")
        cmd
}

# Split a free-text parameter string (blast_params = "-seg no -word_size 3") into its arguments, the way a shell
# would split it into words - but without a shell: white space separates arguments; single quotes keep everything
# between them as it is; double quotes do the same, with \" and \\ for a quote and a backslash inside; outside
# quotes a backslash takes the next character literally. Nothing else has a meaning: ';', '|', '&', '$(...)', '`',
# '>', '*' are ordinary characters of an argument and are handed to the program as they are.
split_params <- function(x) {
        if (is.null(x) || length(x) == 0)
                return(character(0))
        x <- paste(as.character(x[!is.na(x)]), collapse = " ")
        chars <- strsplit(x, "", fixed = TRUE)[[1]]
        out <- character(0)
        cur <- character(0)
        open_tok <- FALSE
        quote <- ""
        i <- 1L
        n <- length(chars)
        while (i <= n) {
                ch <- chars[i]
                if (quote == "'") {
                        if (ch == "'") quote <- "" else cur <- c(cur, ch)
                } else if (quote == "\"") {
                        if (ch == "\"") {
                                quote <- ""
                        } else if (ch == "\\" && i < n && is.element(chars[i + 1L], c("\"", "\\"))) {
                                i <- i + 1L
                                cur <- c(cur, chars[i])
                        } else {
                                cur <- c(cur, ch)
                        }
                } else if (ch == "'" || ch == "\"") {
                        quote <- ch
                        open_tok <- TRUE
                } else if (ch == "\\" && i < n) {
                        i <- i + 1L
                        cur <- c(cur, chars[i])
                        open_tok <- TRUE
                } else if (grepl("[[:space:]]", ch)) {
                        if (open_tok || length(cur) > 0)
                                out <- c(out, paste(cur, collapse = ""))
                        cur <- character(0)
                        open_tok <- FALSE
                } else {
                        cur <- c(cur, ch)
                }
                i <- i + 1L
        }
        if (nzchar(quote))
                orthologr_abort(class = "bad_argument", "Unbalanced quote (", quote, ") in the parameter string: ", x)
        if (open_tok || length(cur) > 0)
                out <- c(out, paste(cur, collapse = ""))
        out
}

# Quote arguments for the command line. An argument wrapped in I() is a free-text parameter string users hand over
# (e.g. blast_params = "-seg no -word_size 3"): it is split into its arguments (split_params()) and every one of
# them is quoted like all other arguments, so that the command line the shell sees holds no character with a
# meaning for the shell.
quote_args <- function(args) {
        if (length(args) == 0)
                return(character(0))
        if (!is.list(args))
                args <- as.list(args)
        type <- if (.Platform$OS.type == "windows") "cmd" else "sh"
        unlist(lapply(args, function(a) {
                if (is.null(a) || length(a) == 0)
                        return(character(0))
                if (inherits(a, "AsIs")) {
                        a <- split_params(unclass(a))
                        if (length(a) == 0)
                                return(character(0))
                }
                shQuote(as.character(a), type = type)
        }), use.names = FALSE)
}

# Run an external program.
#   program : executable name, e.g. "blastp"
#   args    : list / character vector of arguments; each one is quoted unless wrapped in I()
#   path    : optional folder that holds the executable. It is also put first on the PATH while the program
#             runs (so that helper programs next to it are found) and the PATH is restored afterwards.
#   stdout  : file that receives the program's standard output (NULL: a log file that is removed again)
#   tool    : name used in messages; context: what the program was run for
#   output  : files that must exist and be non-empty when the program has finished
#   error   : stop on failure (TRUE) or return the result for the caller to inspect (FALSE)
#   echo    : show the program's messages (stderr) with message() also when it succeeded
#   record  : FALSE: not listed in the commands of the record of a run (version checks of tool_available())
#   limit   : the default time limit of this run in seconds, in place of that of its program in
#             tool_timeout_defaults() (options(orthologr.tool_timeout) still decide when they name it); NULL: the
#             default of the program
# Returns (invisibly) list(status, stdout = <lines or NULL when sent to a file>, stderr = <lines>).
run_tool <- function(program, args = list(), path = NULL, stdout = NULL, tool = program, context = NULL,
                     output = NULL, error = TRUE, echo = FALSE, record = TRUE, limit = NULL) {
        cmd <- tool_command(program, path)
        if (!is.null(path) && nzchar(path)) {
                old_path <- Sys.getenv("PATH")
                on.exit(Sys.setenv(PATH = old_path), add = TRUE)
                Sys.setenv(PATH = paste(path, old_path, sep = .Platform$path.sep))
        }
        if (record && is.element(program, c("diamond", "blastp", "makeblastdb")))
                log_command(cmd, args)
        err_file <- tempfile("orthologr_stderr_")
        out_file <- if (is.null(stdout)) tempfile("orthologr_stdout_") else stdout
        on.exit(unlink(c(err_file, if (is.null(stdout)) out_file)), add = TRUE)
        # (a program that does not finish within its time limit is stopped: tool_timeout(), status 124)
        limit <- tool_timeout(tool, program, default = limit)
        timed_out <- FALSE
        t0 <- proc.time()[["elapsed"]]
        status <- tryCatch(
                withCallingHandlers(
                        system2(cmd, args = quote_args(args), stdout = out_file, stderr = err_file,
                                timeout = if (is.finite(limit)) limit else 0),
                        warning = function(w) {
                                if (grepl("timed out", conditionMessage(w), fixed = TRUE))
                                        timed_out <<- TRUE
                                invokeRestart("muffleWarning")
                        }),
                error = function(e) {
                        writeLines(conditionMessage(e), err_file)
                        127L
                }
        )
        # (the warning of system2() may be translated: status 124 after the time limit is a timeout as well)
        if (is.finite(limit) && identical(as.integer(status), 124L) &&
            (timed_out || proc.time()[["elapsed"]] - t0 >= limit - 0.5)) {
                # (counted, so that a caller that catches the error knows that a program timed out:
                # compute_dnds() then retries the pairs of a chunk with a shorter limit, timeout_cap())
                .orthologr_state$timeouts <- tool_timeouts() + 1L
                orthologr_abort(class = c("tool_timeout", "tool_failed"),
                                tool, " did not finish within ", format(limit), " seconds",
                                if (!is.null(context) && nzchar(context)) paste0(" for ", context), " and was stopped. ",
                                "The time limit is options(orthologr.tool_timeout) (see ?orthologr_options), e.g. ",
                                "options(orthologr.tool_timeout = c(\"", tool_timeout_key(tool, program), "\" = ",
                                format(limit * 4), ")); Inf: no limit.", call. = FALSE)
        }
        res <- list(
                status = status,
                stdout = if (is.null(stdout) && file.exists(out_file)) readLines(out_file, warn = FALSE) else NULL,
                stderr = if (file.exists(err_file)) readLines(err_file, warn = FALSE) else character(0)
        )
        if (error) {
                # the program's own words: stderr, or stdout for programs that report there
                log <- err_file
                if (!length(res$stderr[nzchar(trimws(res$stderr))]) && is.null(stdout))
                        log <- out_file
                check_external(status, tool = tool, context = context, log = log, output = output)
        }
        if (echo && length(res$stderr) > 0)
                message(paste(res$stderr, collapse = "\n"))
        invisible(res)
}

# The time limit in seconds for one run of an external program (run_tool()): options(orthologr.tool_timeout), either
# one number for every program or numbers named by program (names as tool_timeout_defaults(), or the name of the
# program in messages, e.g. "ClustalW"; case is ignored) that replace the defaults for those programs. Inf (or 0,
# or NA): no limit. Programs without a limit by default: the alignment programs, OrthoFinder, version checks.
# A limit set by timeout_cap() for the program (or its key) applies when it is shorter. `default`: the default limit
# of this run in place of that of its program (diamond_makedb(): scaled with the size of the input).
tool_timeout <- function(tool, program = tool, default = NULL) {
        limits <- tool_timeout_defaults()
        if (!is.null(default)) {
                key <- tool_timeout_key(tool, program)
                hit <- match(tolower(key), tolower(names(limits)))
                limits[[if (is.na(hit)) key else names(limits)[hit]]] <- default
        }
        opt <- getOption("orthologr.tool_timeout", NULL)
        if (!is.null(opt)) {
                if (!is.numeric(opt) || length(opt) == 0 || any(!is.na(opt) & opt < 0))
                        orthologr_abort(class = "bad_argument", "options(orthologr.tool_timeout = ) must be a number of seconds ",
                                        "(Inf: no limit), or numbers named by program, e.g. c(KaKs_Calculator = 600).",
                                        call. = FALSE)
                if (is.null(names(opt))) {
                        if (length(opt) != 1)
                                orthologr_abort(class = "bad_argument", "options(orthologr.tool_timeout = ) must be one number of seconds ",
                                                "or numbers named by program, e.g. c(KaKs_Calculator = 600).", call. = FALSE)
                        return(tool_timeout_capped(tool_timeout_value(opt), tool, program))
                }
                nm <- names(opt)
                if (any(is.na(nm) | !nzchar(nm)))
                        orthologr_abort(class = "bad_argument", "options(orthologr.tool_timeout = ): every number of a named limit needs ",
                                        "the name of its program, e.g. c(KaKs_Calculator = 600, DIAMOND = Inf); an element without ",
                                        "a name is not used for any program.", call. = FALSE)
                # (a name that is one of the defaults but for case replaces it)
                known <- match(tolower(nm), tolower(names(limits)))
                nm[!is.na(known)] <- names(limits)[known[!is.na(known)]]
                limits[nm] <- opt
                # a limit for DIAMOND is also the limit of the steps of DIAMOND that have defaults of their own,
                # unless they are named as well
                steps <- c("diamond makedb", "diamond deepclust", "diamond realign")
                if (is.element("DIAMOND", nm))
                        limits[setdiff(steps, nm)] <- limits[["DIAMOND"]]
        }
        keys <- tolower(names(limits))
        hit <- match(tolower(c(tool, tool_timeout_key(tool, program))), keys)
        hit <- hit[!is.na(hit)]
        value <- if (length(hit) == 0) Inf else tool_timeout_value(limits[[hit[1]]])
        tool_timeout_capped(value, tool, program)
}

# A limit in whole seconds as system2(timeout = ) takes it: fractions are rounded up (0.5 s would be 0 s, which is
# no limit at all), and Inf, 0, NA and numbers too large for an integer mean no limit.
tool_timeout_value <- function(x) {
        x <- as.numeric(x)
        if (is.na(x) || x <= 0 || x > .Machine$integer.max) Inf else ceiling(x)
}

# The default time limits (seconds) of one run: KaKs_Calculator and PAL2NAL (one gene pair, or one chunk of up to
# 100 pairs) 1 hour; building a DIAMOND database (diamond makedb) 1 hour, or longer for a large input
# (makedb_time_limit(), which diamond_makedb() gives run_tool()); DIAMOND, BLAST+ and coal-cpp (a whole
# search or pipeline, which may take many hours on large proteomes) and DIAMOND's clustering of a database
# (deepclust(), realign) 7 days. (DIAMOND 2.2 can hang instead of exiting on some input errors, without using the
# processor: a limit ends such a run.)
tool_timeout_defaults <- function() {
        c(KaKs_Calculator = 3600, PAL2NAL = 3600, "diamond makedb" = 3600, DIAMOND = 604800, BLAST = 604800,
          COAL = 604800, "diamond deepclust" = 604800, "diamond realign" = 604800)
}

# The name of the program in tool_timeout_defaults() for a run of 'program' called 'tool' in messages (the name in
# the messages when none fits). The steps of DIAMOND with limits of their own are matched before DIAMOND.
tool_timeout_key <- function(tool, program = tool) {
        x <- tolower(paste(tool, basename(as.character(program))))
        if (grepl("kaks", x, fixed = TRUE)) return("KaKs_Calculator")
        if (grepl("pal2nal", x, fixed = TRUE)) return("PAL2NAL")
        for (step in c("makedb", "deepclust", "realign"))
                if (grepl("diamond", x, fixed = TRUE) && grepl(step, x, fixed = TRUE)) return(paste("diamond", step))
        if (grepl("diamond", x, fixed = TRUE)) return("DIAMOND")
        if (grepl("blast", x, fixed = TRUE)) return("BLAST")
        if (grepl("coal", x, fixed = TRUE)) return("COAL")
        tool
}

# How many runs of external programs were stopped at their time limit in this R process (run_tool()).
tool_timeouts <- function() {
        n <- .orthologr_state$timeouts
        if (is.null(n)) 0L else n
}

# Evaluate `expr` with the time limit of the program `key` (a name of tool_timeout_defaults()) at most `seconds`
# (or at most a shorter limit that an enclosing timeout_cap() set).
timeout_cap <- function(key, seconds, expr) {
        old <- .orthologr_state$timeout_cap
        cap <- old
        k <- tolower(key)
        cap[[k]] <- if (is.null(old[[k]])) seconds else min(old[[k]], seconds)
        .orthologr_state$timeout_cap <- cap
        on.exit(.orthologr_state$timeout_cap <- old, add = TRUE)
        expr
}

tool_timeout_capped <- function(value, tool, program = tool) {
        cap <- .orthologr_state$timeout_cap
        if (length(cap) == 0)
                return(value)
        c1 <- cap[[tolower(tool_timeout_key(tool, program))]]
        if (!is.null(c1) && is.finite(c1))
                value <- min(value, ceiling(c1))
        value
}

# Is `program` available? It is started once per R session (and per `path`, `version_args`, `pattern` and
# `ok_status`) with `version_args`; when it is there, the first line that matches `pattern` (the version line) is
# remembered, so that functions called for every gene pair do not start the program again just to see that it is
# there. Stops with `not_found` when it is not (that answer is not remembered: a program installed later in the
# session is found). `ok_status`: the exit statuses of a program that runs (NULL: any; the version line decides);
# a program that prints the pattern but exits with another status does not run (e.g. a Python program that stops
# with a traceback naming its own folder). The check is not part of the commands of the record of a run.
tool_available <- function(program, version_args, pattern, path = NULL, not_found = NULL, ok_status = NULL) {
        key <- paste(c(program, if (is.null(path)) "" else path, version_args, pattern, ok_status), collapse = "\r")
        hit <- .orthologr_state[[key]]
        if (!is.null(hit))
                return(hit)
        res <- run_tool(program, args = version_args, path = path, error = FALSE, record = FALSE)
        lines <- c(res$stdout, res$stderr)
        found <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
        status <- as.integer(res$status)
        if (!is.null(ok_status) && length(found) > 0 && !identical(status, 127L) &&
            !(length(status) == 1 && is.element(status, ok_status)))
                found <- character(0)
        if (length(found) == 0 && length(status) == 1 && !is.na(status) && !is.element(status, c(0L, 127L))) {
                # the program is there but does not run (e.g. it crashes, exit status 139): that is not a missing program
                said <- utils::tail(lines[nzchar(trimws(lines))], 3)
                orthologr_abort("The program '", program, "' was found but could not be started ('", program, " ",
                                paste(version_args, collapse = " "), "' exited with status ", status, ")",
                                if (length(said) > 0) paste0(": ", paste(said, collapse = " / ")) else ".",
                                class = "tool_failed")
        }
        if (identical(status, 127L) || length(found) == 0) {
                if (is.null(not_found))
                        not_found <- paste0("The program '", program, "' was not found.")
                orthologr_abort(not_found, class = "tool_missing")
        }
        .orthologr_state[[key]] <- found[1]
        found[1]
}

# Forget which programs were found (for tests, or after installing a program in a running session).
tool_cache_reset <- function() {
        rm(list = setdiff(ls(.orthologr_state, all.names = TRUE), c("seq_cache", "dropped", "current_file")), envir = .orthologr_state)
        invisible(TRUE)
}

# A folder of its own for one call, inside the pipeline folder (or tempdir()). Returns the path; use
# on.exit(call_dir_end(dir)) to remove it (kept with options(orthologr.keep_tmp = TRUE)).
call_dir_start <- function(prefix) {
        d <- tempfile(prefix, tmpdir = orthologr_tmpdir())
        dir.create(d, recursive = TRUE)
        d
}

call_dir_end <- function(dir) {
        if (!keep_tmp_setting()$keep)
                unlink(dir, recursive = TRUE, force = TRUE)
        invisible(dir)
}

# Move a finished result file to its final place. The file is written under a name of its own first, so a file
# left by an earlier or failed call is never taken for the result of this one.
publish_file <- function(from, to) {
        if (file.exists(to))
                unlink(to)
        if (!suppressWarnings(file.rename(from, to))) {
                if (!file.copy(from, to, overwrite = TRUE))
                        orthologr_abort(class = "bad_input", "The result file '", to, "' could not be written.")
                unlink(from)
        }
        invisible(to)
}

# A file name without characters that the database builders cannot handle (spaces, quotes, ...).
safe_file_name <- function(x) {
        gsub("[^A-Za-z0-9._+-]", "_", x)
}

# Number of threads / processes to use: a positive whole number; when more cores are asked for than the machine
# has, the number of available cores is used and a warning says so.
check_comp_cores <- function(comp_cores) {
        if (!is.numeric(comp_cores) || length(comp_cores) != 1 || is.na(comp_cores) || !is.finite(comp_cores) ||
            comp_cores < 1 || comp_cores != round(comp_cores))
                orthologr_abort("comp_cores must be a positive whole number (1, 2, ...).", class = "bad_argument")
        cores <- parallel::detectCores()
        if (!is.na(cores) && comp_cores > cores) {
                orthologr_warn(class = "comp_cores", "comp_cores = ", comp_cores, " is more than the ", cores, " core", if (cores != 1) "s",
                               " available on this machine: using comp_cores = ", cores, ".", call. = FALSE)
                comp_cores <- cores
        }
        comp_cores
}

# Copy a result file into the folder (or to the file) given as save.output.
save_output_copy <- function(from, save.output, name) {
        if (is.null(save.output))
                return(invisible(NULL))
        to <- if (dir.exists(save.output)) file.path(save.output, name) else save.output
        file.copy(from, to, overwrite = TRUE)
        invisible(to)
}

# Stop with an installation hint when a package that is only needed by some functions (DESCRIPTION: Suggests) is
# not installed.
need_package <- function(pkg, what, bioc = FALSE) {
        if (requireNamespace(pkg, quietly = TRUE))
                return(invisible(TRUE))
        orthologr_abort(
                class = "package_missing",
                "The package '", pkg, "' is needed for ", what, " but is not installed. Please install it with ",
                if (bioc) paste0("BiocManager::install(\"", pkg, "\")") else paste0("install.packages(\"", pkg, "\")"),
                ".",
                call. = FALSE
        )
}

# Evaluate every argument of the calling function now. A function that hands a closure to worker processes
# (parallel::parLapply() on a PSOCK cluster, the default on Windows) must call this before the cluster is used:
# an argument that was not evaluated yet travels to the workers as a promise - the expression the user typed plus
# a reference to the environment it was typed in - and is evaluated there, where a variable of the user's
# workspace does not exist ("object 'kp' not found") or means something else (`quiet = qt` found stats::qt).
force_args <- function() {
        env <- parent.frame()
        fn <- sys.function(sys.parent())
        for (arg in setdiff(names(formals(fn)), "..."))
                get(arg, envir = env, inherits = FALSE)
        invisible(TRUE)
}
