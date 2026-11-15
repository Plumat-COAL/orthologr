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
                stop("Unbalanced quote (", quote, ") in the parameter string: ", x, call. = FALSE)
        if (open_tok || length(cur) > 0)
                out <- c(out, paste(cur, collapse = ""))
        out
}

# Quote arguments for the command line. An argument wrapped in I() is a free-text parameter string users hand over
# (e.g. blast_params = "-seg no -word_size 3"): it is split into its arguments (split_params()) and every one of
# them is quoted like all other arguments, so that the command line the shell sees holds no character with a
# meaning for the shell. (Until 0.5.0 such strings were pasted into the command line as they were:
# params = "--force; touch X" ran 'touch X'.)
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
# Returns (invisibly) list(status, stdout = <lines or NULL when sent to a file>, stderr = <lines>).
run_tool <- function(program, args = list(), path = NULL, stdout = NULL, tool = program, context = NULL,
                     output = NULL, error = TRUE, echo = FALSE) {
        cmd <- tool_command(program, path)
        if (!is.null(path) && nzchar(path)) {
                old_path <- Sys.getenv("PATH")
                on.exit(Sys.setenv(PATH = old_path), add = TRUE)
                Sys.setenv(PATH = paste(path, old_path, sep = .Platform$path.sep))
        }
        err_file <- tempfile("orthologr_stderr_")
        out_file <- if (is.null(stdout)) tempfile("orthologr_stdout_") else stdout
        on.exit(unlink(c(err_file, if (is.null(stdout)) out_file)), add = TRUE)
        status <- tryCatch(
                suppressWarnings(system2(cmd, args = quote_args(args), stdout = out_file, stderr = err_file)),
                error = function(e) {
                        writeLines(conditionMessage(e), err_file)
                        127L
                }
        )
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

# Is `program` available? It is started once per R session (and per `path`) with `version_args`; the answer and the
# first line that matches `pattern` (the version line) are remembered, so that functions called for every gene
# pair do not start the program again just to see that it is there. Stops with `not_found` when it is not.
tool_available <- function(program, version_args, pattern, path = NULL, not_found = NULL) {
        key <- paste0(program, "\r", if (is.null(path)) "" else path)
        hit <- .orthologr_state[[key]]
        if (!is.null(hit))
                return(hit)
        res <- run_tool(program, args = version_args, path = path, error = FALSE)
        lines <- c(res$stdout, res$stderr)
        found <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
        if (identical(as.integer(res$status), 127L) || length(found) == 0) {
                if (is.null(not_found))
                        not_found <- paste0("The program '", program, "' was not found.")
                stop(not_found, call. = FALSE)
        }
        .orthologr_state[[key]] <- found[1]
        found[1]
}

# Forget which programs were found (for tests, or after installing a program in a running session).
tool_cache_reset <- function() {
        rm(list = setdiff(ls(.orthologr_state, all.names = TRUE), "seq_cache"), envir = .orthologr_state)
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
        if (!isTRUE(getOption("orthologr.keep_tmp", FALSE)))
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
                        stop("The result file '", to, "' could not be written.", call. = FALSE)
                unlink(from)
        }
        invisible(to)
}

# A file name without characters that the database builders cannot handle (spaces, quotes, ...).
safe_file_name <- function(x) {
        gsub("[^A-Za-z0-9._+-]", "_", x)
}

# Number of threads / processes to use: when more cores are asked for than the machine has, the number of
# available cores is used and a warning says so (this used to be an error).
check_comp_cores <- function(comp_cores) {
        if (!is.numeric(comp_cores) || length(comp_cores) != 1 || is.na(comp_cores) || comp_cores < 1)
                stop("comp_cores must be a positive number.", call. = FALSE)
        cores <- parallel::detectCores()
        if (!is.na(cores) && comp_cores > cores) {
                warning("comp_cores = ", comp_cores, " is more than the ", cores, " core", if (cores != 1) "s",
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
        stop(
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
