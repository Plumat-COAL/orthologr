#' @title Install the DIAMOND aligner into the user cache
#' @description Downloads a pre-compiled DIAMOND2 binary from the official
#' GitHub releases (\url{https://github.com/bbuchfink/diamond/releases}) that
#' matches the operating system of the user, checks the SHA-256 checksum of the download, and stores the
#' program in a user-level cache directory (see \code{\link[tools]{R_user_dir}}). No administrator rights are
#' required. Once installed, \code{\link{diamond}}, \code{\link{diamond_best}},
#' \code{\link{diamond_rec}}, \code{\link{set_diamond}}, \code{\link{deepclust}} and
#' \code{\link{deepclust_realign}} find the cached binary automatically whenever
#' \code{path = NULL} and no \code{diamond} is available on the \code{PATH}.
#' @param version a DIAMOND release tag, e.g. \code{"v2.2.8"}. The default \code{"latest"} is the newest release
#' whose checksum is pinned in this version of orthologr and that has a binary for this operating system
#' (\code{"v2.2.8"}; see Details) - not necessarily the newest release of DIAMOND.
#' @param force logical; re-download even if DIAMOND is already cached.
#' @param quiet logical; suppress progress messages.
#' @param sha256 \code{NULL} (default): the download must match the checksum pinned in orthologr for
#' \code{version}. Or the SHA-256 checksum (64 hexadecimal digits) the download must have, for a release that is
#' not pinned; take it from a source you trust (e.g. the \code{sha256:} digest that GitHub shows next to the file
#' on the release page).
#' @param timeout the most seconds the download may take (R's option \code{timeout} during the download; at
#' least 600 by default). The check that the program runs is given 60 seconds.
#' @details Pre-compiled binaries are provided upstream for Linux (x86_64), macOS (one universal binary for
#' x86_64 and arm64) and Windows (x86_64).
#' On other platforms (e.g. Linux on ARM, Windows on ARM; error of class \code{orthologr_error_tool_missing})
#' please install DIAMOND via conda or build it from source:
#' \code{conda install -c bioconda diamond}.
#'
#' \strong{Checksums.} The DIAMOND releases do not publish checksum files, so orthologr carries the SHA-256
#' checksums of the release files it installs (internal table \code{diamond_release_checksums()}): DIAMOND
#' v2.2.8 (19 September 2026) for Linux, macOS and Windows. They were computed from the files of the GitHub
#' release, downloaded twice, and are compared with the digests GitHub records for the release files by the
#' continuous integration of orthologr. A download that does not match is deleted and nothing is installed
#' (error of class \code{orthologr_error_checksum_mismatch}). A release without a pinned checksum is not
#' downloaded unless \code{sha256} is given (error of class \code{orthologr_error_bad_argument}): this fails
#' closed rather than installing a file that was not verified. DIAMOND v2.2.9 (14 November 2026) is not pinned:
#' when this version of orthologr was made its release had no macOS binary yet and its Linux file had been
#' rebuilt a day after the release, so a checksum taken then might not stay valid.
#'
#' The checksum is computed with \code{tools::sha256sum()} (R >= 4.5.0), else with the program
#' \code{sha256sum}, \code{shasum} or \code{certutil} (Windows); when none of them is available the downloaded
#' file is deleted without being unpacked (error of class \code{orthologr_error_checksum_unavailable}).
#'
#' \strong{Errors.} A download that fails or takes longer than \code{timeout} seconds is an error of class
#' \code{orthologr_error_download_failed} (and \code{orthologr_error_tool_failed}); the partial file is deleted.
#' The archive is unpacked into a temporary folder; the program is copied into the cache folder under a temporary
#' name and run there (not in the temporary folder, which may be mounted \code{noexec}): \code{diamond version}
#' must exit with status 0 within 60 seconds and print a line \code{diamond version <number>} (otherwise an error
#' of class \code{orthologr_error_tool_failed} with the program's last lines of output, and the temporary copy is
#' removed). Only then is it renamed over the cached copy, so a failed installation leaves an existing copy in
#' place. Temporary copies left by an interrupted installation are removed by the next one.
#'
#' With a cached copy and \code{force = FALSE} nothing is downloaded; a \code{version} other than
#' \code{"latest"} or a \code{sha256} is then not used, with a warning.
#' @return Invisibly, the directory containing the \code{diamond} executable. This can be passed
#' to the \code{path} argument of the DIAMOND functions.
#' @note Pre-built DIAMOND2 releases exist for Linux x86_64, macOS and Windows; on other platforms (e.g. Linux
#' aarch64) the function stops with a hint to install DIAMOND2 with a package manager. The archive is fetched over
#' HTTPS from the GitHub releases of \url{https://github.com/bbuchfink/diamond} and verified against its pinned
#' SHA-256 checksum (or the one given as \code{sha256}) before it is unpacked.
#' @examples
#' # NULL when DIAMOND is on the PATH, the folder of a copy installed by diamond_install(),
#' # an error otherwise
#' tryCatch(diamond_path(), error = conditionMessage)
#' \dontrun{
#' # download DIAMOND into that folder (needs an internet connection)
#' diamond_install()
#' # a release that orthologr has no checksum for: give the checksum yourself
#' diamond_install(version = "v2.1.25", sha256 = "<the 64 hexadecimal digits>")
#' }
#' @seealso \code{\link{diamond_path}}, \code{\link{diamond}}
#' @export
diamond_install <- function(version = "latest", force = FALSE, quiet = FALSE, sha256 = NULL,
                            timeout = max(600, getOption("timeout"))) {

        if (!is.character(version) || length(version) != 1 || is.na(version) ||
            !(identical(version, "latest") || grepl("^v?[0-9]+(\\.[0-9]+)*$", version)))
                orthologr_abort(class = "bad_argument", "'version' must be one DIAMOND release tag such as \"v2.2.8\", or \"latest\".")
        if (!is.numeric(timeout) || length(timeout) != 1 || is.na(timeout) || timeout <= 0)
                orthologr_abort(class = "bad_argument", "'timeout' must be a positive number of seconds.")
        if (!is.null(sha256)) {
                if (!is.character(sha256) || length(sha256) != 1 || is.na(sha256) || !grepl("^[0-9A-Fa-f]{64}$", sha256))
                        orthologr_abort(class = "bad_argument", "'sha256' must be NULL or a SHA-256 checksum: 64 hexadecimal digits.")
                sha256 <- tolower(sha256)
        }

        platform <- diamond_platform()
        asset <- diamond_platform_asset(platform[["sysname"]], platform[["machine"]])
        if (is.na(asset))
                orthologr_abort(class = "tool_missing",
                                "No pre-compiled DIAMOND binary is available for ", platform[["sysname"]], " (", platform[["machine"]], "). ",
                                "Please install it via conda: 'conda install -c bioconda diamond' (or build it from source).")

        exe_name <- diamond_exe_name()
        dest     <- diamond_cache_dir()
        exe      <- file.path(dest, exe_name)

        if (file.exists(exe) && !force) {
                if (!identical(version, "latest") || !is.null(sha256))
                        orthologr_warn(class = "tool_version",
                                       "DIAMOND is already installed in ", dest, "; 'version' and 'sha256' were not used. ",
                                       "Use force = TRUE to install the release asked for.")
                else if (!quiet)
                        message("DIAMOND is already installed in ", dest, " (use force = TRUE to reinstall).")
                return(invisible(dest))
        }

        pins <- diamond_release_checksums()
        pins <- pins[pins$asset == asset, , drop = FALSE]
        if (identical(version, "latest")) {
                if (nrow(pins) == 0)
                        orthologr_abort(class = "bad_argument", "orthologr has no pinned checksum for the DIAMOND file '", asset,
                                        "'. Give a release and its checksum: diamond_install(version = \"v...\", sha256 = \"...\").")
                version <- pins$version[order(numeric_version(sub("^v", "", pins$version)), decreasing = TRUE)][1]
        } else {
                version <- sub("^v?", "v", version)
        }
        pin <- pins[pins$version == version, , drop = FALSE]
        expected <- if (nrow(pin) == 1) pin$sha256 else NA_character_
        if (!is.null(sha256)) {
                if (!is.na(expected) && !identical(sha256, expected))
                        orthologr_abort(class = "bad_argument", "'sha256' (", sha256, ") is not the checksum orthologr has pinned for DIAMOND ",
                                        version, " (", asset, "): ", expected, ".")
                expected <- sha256
        }
        if (is.na(expected))
                orthologr_abort(class = "bad_argument",
                                "orthologr has no pinned SHA-256 checksum for DIAMOND ", version, " (", asset, "), so it is not ",
                                "downloaded unverified. Pinned: ", paste(unique(pins$version), collapse = ", "),
                                " (version = \"latest\"). For another release give its checksum from a source you trust: ",
                                "diamond_install(version = \"", version, "\", sha256 = \"<64 hexadecimal digits>\"), ",
                                "or install DIAMOND via conda: 'conda install -c bioconda diamond'.")

        url <- paste0("https://github.com/bbuchfink/diamond/releases/download/", version, "/", asset)

        work <- tempfile("diamond_install_")
        dir.create(work)
        on.exit(unlink(work, recursive = TRUE), add = TRUE)
        archive <- file.path(work, asset)

        if (!quiet) message("Downloading ", url, " ...")
        result <- tryCatch(
                {
                        diamond_download(url, archive, timeout = timeout, quiet = quiet)
                        NULL
                },
                error = function(e) conditionMessage(e),
                warning = function(w) conditionMessage(w)
        )
        if (!is.null(result) || !file.exists(archive)) {
                unlink(archive)
                timed_out <- !is.null(result) && grepl("timeout|timed out", result, ignore.case = TRUE)
                orthologr_abort(class = c("download_failed", "tool_failed"),
                                "Could not download DIAMOND from ", url, ": ",
                                if (timed_out) paste0("the download did not finish within ", timeout, " seconds (argument 'timeout'). ") else "",
                                if (is.null(result)) "no file was written." else result)
        }

        actual <- diamond_sha256(archive)
        if (!identical(actual, expected)) {
                unlink(archive)
                orthologr_abort(class = c("checksum_mismatch", "tool_failed"),
                                "The download of DIAMOND ", version, " (", url, ") does not have the expected SHA-256 checksum ",
                                "and was deleted; nothing was installed. Expected ", expected, ", got ", actual, ". ",
                                "The file may have been changed on the server or in transit. Install DIAMOND via conda instead: ",
                                "'conda install -c bioconda diamond'.")
        }
        if (!quiet) message("SHA-256 checksum verified: ", actual)

        unpacked <- file.path(work, "unpacked")
        dir.create(unpacked)
        if (grepl("\\.zip$", asset)) {
                utils::unzip(archive, exdir = unpacked)
        } else {
                utils::untar(archive, exdir = unpacked)
        }
        new_exe <- file.path(unpacked, exe_name)
        if (!file.exists(new_exe))
                orthologr_abort(class = c("download_failed", "tool_failed"),
                                "The downloaded archive did not contain '", exe_name, "'; nothing was installed.")

        # Into the cache folder under a temporary name (ending in the program's name, so that Windows runs it),
        # checked there - not in the temporary folder, which may be mounted noexec - and then renamed over the old
        # copy in one step. Temporary copies that an interrupted installation left behind are removed first.
        dir.create(dest, recursive = TRUE, showWarnings = FALSE)
        diamond_remove_stale_staged(dest, exe_name)
        staged <- file.path(dest, paste0(".new-", Sys.getpid(), "-", exe_name))
        on.exit(unlink(staged), add = TRUE)
        if (!file.copy(new_exe, staged, overwrite = TRUE))
                orthologr_abort(class = "tool_failed", "Could not copy DIAMOND into ", dest, "; nothing was installed.")
        Sys.chmod(staged, mode = "0755")

        # the program must run: exit status 0 and its version line
        out <- diamond_binary_version(staged)
        status <- attr(out, "status")
        version_line <- grep("^diamond version [0-9]", out, value = TRUE)
        if (is.null(out) || !(is.null(status) || identical(as.integer(status), 0L)) || length(version_line) == 0) {
                last <- utils::tail(out[nzchar(trimws(out))], 5)
                how <- if (is.null(out)) "could not be started"
                       else if (identical(as.integer(status), 124L)) "did not finish within 60 seconds"
                       else paste0("exited with status ", if (is.null(status)) 0L else status)
                orthologr_abort(class = "tool_failed",
                                "DIAMOND ", version, " was downloaded and verified but does not run here ('diamond version' ", how,
                                if (length(last)) paste0("; its last output: ", paste(last, collapse = " | ")) else "; no output",
                                "; it was run from ", dest, "). Nothing was installed. ",
                                if (identical(as.integer(status), 126L)) "Status 126 means it may not be executed there (e.g. a folder on a filesystem mounted noexec; the cache folder is tools::R_user_dir(\"orthologr\", \"cache\"), set R_USER_CACHE_DIR to move it). " else "",
                                "Otherwise please install DIAMOND via conda: 'conda install -c bioconda diamond'.")
        }

        if (!file.rename(staged, exe))
                orthologr_abort(class = "tool_failed", "Could not move DIAMOND into place in ", dest, "; the copy there is unchanged.")

        if (!quiet) message("Installed ", version_line[1], " in ", dest)
        invisible(dest)
}

# The SHA-256 checksums of the DIAMOND release files that diamond_install() downloads (one row per release and
# file). The DIAMOND releases publish no checksum files: these were computed (sha256sum) from the files of the
# GitHub release bbuchfink/diamond, downloaded twice on 15 November 2026 over separate connections (identical);
# the CI job 'diamond-checksums' (.github/workflows/R-CMD-check.yaml) compares them with the digests GitHub
# records for the release files and with a fresh download. Add a release here only with checksums obtained the
# same way.
diamond_release_checksums <- function() {
        data.frame(
                version = c("v2.2.8", "v2.2.8", "v2.2.8"),
                asset   = c("diamond-linux64.tar.gz", "diamond-macos.tar.gz", "diamond-windows.zip"),
                size    = c(10221459, 7469225, 3366672),
                sha256  = c("33a765d192ae4f5b388f15fe2e7c23d3aaea0f9e3f6872efe8c324aec261fbc7",
                            "7d17c9c72b927d74ae1e72b4b7ca44f6100af02e9f753ad1ecded9abb30a3b22",
                            "8035ece14562e13576da6e9583c74ab1aaff9d0c96086d93ab337e95cd288ce1"),
                stringsAsFactors = FALSE
        )
}

diamond_platform <- function() {
        info <- Sys.info()
        c(sysname = info[["sysname"]], machine = info[["machine"]])
}

# the name of the release file for an operating system, NA when DIAMOND has none
# (macOS: one universal binary for x86_64 and arm64; Windows and Linux: x86_64 only)
diamond_platform_asset <- function(sysname, machine) {
        x64 <- tolower(machine) %in% c("x86_64", "amd64", "x86-64", "x64")
        switch(
                sysname,
                Linux   = if (x64) "diamond-linux64.tar.gz" else NA_character_,
                Darwin  = "diamond-macos.tar.gz",
                Windows = if (x64) "diamond-windows.zip" else NA_character_,
                NA_character_
        )
}

# Remove copies '.new-<pid>-diamond[.exe]' that an installation interrupted between the copy and the rename left in
# the cache folder (older than an hour, so that an installation running in another R process is not disturbed).
diamond_remove_stale_staged <- function(dest, exe_name) {
        f <- list.files(dest, pattern = paste0("^\\.new-[0-9]+-", gsub(".", "\\.", exe_name, fixed = TRUE), "$"),
                        all.files = TRUE, full.names = TRUE)
        old <- f[difftime(Sys.time(), file.mtime(f), units = "hours") > 1]
        unlink(old)
        invisible(old)
}

# Download 'url' to 'file' within 'timeout' seconds (R's option 'timeout' applies to the whole transfer).
diamond_download <- function(url, file, timeout, quiet = FALSE, download = utils::download.file) {
        old <- options(timeout = timeout)
        on.exit(options(old), add = TRUE)
        download(url, file, mode = "wb", quiet = quiet)
}

# The SHA-256 checksum of a file as 64 lower-case hexadecimal digits: tools::sha256sum() (R >= 4.5.0) or a
# program of the operating system. An error when none is available: the download is then not used.
diamond_sha256 <- function(file, use_r = TRUE) {
        if (use_r && exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
                h <- unname(get("sha256sum", envir = asNamespace("tools"))(file))
                if (length(h) == 1 && !is.na(h))
                        return(tolower(h))
        }
        path <- normalizePath(file, mustWork = TRUE)
        run <- function(cmd, args) {
                if (!nzchar(Sys.which(cmd)))
                        return(NA_character_)
                out <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE, timeout = 120)),
                                error = function(e) character())
                # sha256sum / shasum: "<digest>  <file>"; certutil: the digest on a line of its own (in pairs
                # separated by spaces on older Windows versions)
                out <- trimws(out)
                cand <- c(sub("[[:space:]].*$", "", out), gsub(" ", "", out[grepl("^[0-9A-Fa-f ]+$", out)]))
                hex <- cand[grepl("^[0-9A-Fa-f]{64}$", cand)]
                if (length(hex) >= 1) tolower(hex[1]) else NA_character_
        }
        h <- run("sha256sum", shQuote(path))
        if (is.na(h)) h <- run("shasum", c("-a", "256", shQuote(path)))
        if (is.na(h) && .Platform$OS.type == "windows") h <- run("certutil", c("-hashfile", shQuote(path), "SHA256"))
        if (is.na(h))
                orthologr_abort(class = c("checksum_unavailable", "tool_missing"),
                                "Cannot compute the SHA-256 checksum of the DIAMOND download (needs R >= 4.5.0, or one of the ",
                                "programs sha256sum, shasum, certutil), so it is not installed. Install DIAMOND via conda: ",
                                "'conda install -c bioconda diamond'.")
        h
}

# the output of 'diamond version' (stdout and stderr, through a file, so that it is kept whatever the exit status)
# with the attribute "status" (the exit status; 124 when it took more than 60 seconds), NULL when the program
# cannot be started at all
diamond_binary_version <- function(exe) {
        log <- tempfile("diamond_version_")
        on.exit(unlink(log), add = TRUE)
        status <- tryCatch(suppressWarnings(system2(exe, "version", stdout = log, stderr = log, timeout = 60)),
                           error = function(e) NULL)
        if (is.null(status))
                return(NULL)
        out <- if (file.exists(log)) suppressWarnings(readLines(log)) else character()
        attr(out, "status") <- as.integer(status)
        out
}

#' @title Locate the DIAMOND executable used by orthologr
#' @description Returns the directory holding the \code{diamond} executable that
#' \code{orthologr} will use. A \code{diamond} on the \code{PATH} takes precedence
#' (returned as \code{NULL}, i.e. no \code{path} needed); otherwise the copy installed with
#' \code{\link{diamond_install}} is returned.
#' @return \code{NULL} if \code{diamond} is on the \code{PATH}, the cache directory if only the
#' cached copy exists, and an error otherwise.
#' @seealso \code{\link{diamond_install}}
#' @export
diamond_path <- function() {
        if (nzchar(Sys.which("diamond"))) return(NULL)
        dir <- diamond_cache_dir()
        if (file.exists(file.path(dir, diamond_exe_name()))) return(dir)
        orthologr_abort(class = "tool_missing", "DIAMOND was not found. Run diamond_install() or install it via conda: ",
             "'conda install -c bioconda diamond'.", call. = FALSE)
}

diamond_exe_name <- function() {
        if (.Platform$OS.type == "windows") "diamond.exe" else "diamond"
}

diamond_cache_dir <- function() {
        file.path(tools::R_user_dir("orthologr", which = "cache"), "diamond")
}

# internal: if the user gave no path and diamond is not on the PATH, fall back
# to the copy installed by diamond_install()
diamond_resolve_path <- function(path) {
        if (!is.null(path) || nzchar(Sys.which("diamond"))) return(path)
        dir <- diamond_cache_dir()
        if (file.exists(file.path(dir, diamond_exe_name()))) return(dir)
        path
}
