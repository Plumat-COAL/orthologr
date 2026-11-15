# ---------------------------------------------------------------------------
# Offline tests: use a fake cache directory (R_USER_CACHE_DIR) and a fake
# 'diamond' file, so no download happens. Safe for GitHub Actions / CRAN.
# ---------------------------------------------------------------------------

# run `code` with a fresh, empty user cache and (optionally) an empty PATH
with_fake_cache <- function(code, empty_path = FALSE) {
        old <- Sys.getenv(c("R_USER_CACHE_DIR", "PATH"), unset = NA)
        on.exit({
                for (v in names(old)) {
                        if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
                }
        }, add = TRUE)
        fake_cache <- tempfile("fake_cache_")
        dir.create(fake_cache)
        on.exit(unlink(fake_cache, recursive = TRUE), add = TRUE)
        Sys.setenv(R_USER_CACHE_DIR = fake_cache)
        if (empty_path) Sys.setenv(PATH = with_coreutils(fake_cache))
        force(code)
}

place_fake_diamond <- function() {
        dir <- orthologr:::diamond_cache_dir()
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        file.create(file.path(dir, orthologr:::diamond_exe_name()))
        dir
}

test_that("diamond_cache_dir() lives in the R user cache", {
        with_fake_cache({
                dir <- orthologr:::diamond_cache_dir()
                dir.create(dir, recursive = TRUE)
                expect_true(startsWith(normalizePath(dir),
                                       normalizePath(Sys.getenv("R_USER_CACHE_DIR"))))
                expect_match(dir, "orthologr")
        })
})

test_that("diamond_install() does nothing if DIAMOND is already cached", {
        with_fake_cache({
                dir <- place_fake_diamond()
                expect_message(res <- diamond_install(), "already installed")
                expect_equal(res, dir)
                expect_silent(diamond_install(quiet = TRUE))
        })
})

test_that("diamond_resolve_path() respects a user supplied path", {
        with_fake_cache({
                place_fake_diamond()
                expect_equal(orthologr:::diamond_resolve_path("/some/path"), "/some/path")
        }, empty_path = TRUE)
})

test_that("diamond_resolve_path() falls back to the cached copy", {
        with_fake_cache({
                expect_null(orthologr:::diamond_resolve_path(NULL))
                dir <- place_fake_diamond()
                expect_equal(orthologr:::diamond_resolve_path(NULL), dir)
        }, empty_path = TRUE)
})

test_that("diamond_path() finds the cached copy or errors", {
        with_fake_cache({
                expect_error(diamond_path(), "diamond_install")
                dir <- place_fake_diamond()
                expect_equal(diamond_path(), dir)
        }, empty_path = TRUE)
})

# ---------------------------------------------------------------------------
# Download, checksum and timeout with a mocked download (no network): the archive "downloaded" is made here, and
# the platform, the pinned checksums and the run of the program are replaced where a test needs it.
# ---------------------------------------------------------------------------

zip_program <- function() {
        z <- Sys.getenv("R_ZIPCMD")
        if (nzchar(z)) z else "zip"
}

# a .tar.gz with a file 'diamond' (or a .zip with 'diamond.exe'), as the DIAMOND releases have
fake_diamond_archive <- function(zip = FALSE) {
        src <- tempfile("fake_diamond_src_")
        dir.create(src)
        exe <- if (zip) "diamond.exe" else "diamond"
        writeLines(c("#!/bin/sh", "echo 'diamond version 2.2.8'"), file.path(src, exe))
        Sys.chmod(file.path(src, exe), "755")
        archive <- tempfile("fake_diamond_", fileext = if (zip) ".zip" else ".tar.gz")
        old <- setwd(src)
        on.exit(setwd(old), add = TRUE)
        if (zip) utils::zip(archive, exe, flags = "-q", zip = zip_program()) else utils::tar(archive, exe, compression = "gzip", tar = "internal")
        archive
}

# mock the download (copy 'archive' or call 'fail'), the platform (Linux x86_64 unless zip), the checksum table
# (the archive's checksum for 'version', or 'pinned_sha') and the run of the program
with_mocked_download <- function(code, archive, version = "v2.2.8", pinned_sha = NULL, fail = NULL, zip = FALSE) {
        calls <- new.env()
        calls$n <- 0L
        asset <- if (zip) "diamond-windows.zip" else "diamond-linux64.tar.gz"
        sha <- if (is.null(pinned_sha)) orthologr:::diamond_sha256(archive) else pinned_sha
        local_mocked_bindings(
                diamond_platform = function() if (zip) c(sysname = "Windows", machine = "x86-64") else c(sysname = "Linux", machine = "x86_64"),
                diamond_exe_name = function() if (zip) "diamond.exe" else "diamond",
                diamond_release_checksums = function() data.frame(version = version, asset = asset, size = file.size(archive),
                                                                  sha256 = sha, stringsAsFactors = FALSE),
                diamond_download = function(url, file, timeout, quiet = FALSE) {
                        calls$n <- calls$n + 1L
                        calls$url <- url
                        calls$timeout <- timeout
                        if (!is.null(fail)) return(fail(url, file, timeout))
                        file.copy(archive, file)
                        invisible(0L)
                },
                diamond_binary_version = function(exe) "diamond version 2.2.8",
                .package = "orthologr"
        )
        force(code)
        calls
}

test_that("the pinned checksums are well-formed and name the files diamond_install() downloads", {
        pins <- orthologr:::diamond_release_checksums()
        expect_true(nrow(pins) >= 3)
        expect_true(all(grepl("^[0-9a-f]{64}$", pins$sha256)))
        expect_true(all(grepl("^v[0-9]+\\.[0-9]+\\.[0-9]+$", pins$version)))
        expect_false(anyDuplicated(pins[, c("version", "asset")]) > 0)
        assets <- c(orthologr:::diamond_platform_asset("Linux", "x86_64"), orthologr:::diamond_platform_asset("Darwin", "arm64"),
                    orthologr:::diamond_platform_asset("Windows", "x86-64"))
        expect_true(all(assets %in% pins$asset))
        expect_true(all(pins$size > 1e6))
        expect_true(is.na(orthologr:::diamond_platform_asset("Linux", "aarch64")))
})

test_that("diamond_sha256() gives the SHA-256 digest (R and the programs of the operating system)", {
        f <- tempfile()
        writeBin(charToRaw("abc"), f)
        abc <- "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        expect_identical(orthologr:::diamond_sha256(f), abc)
        if (nzchar(Sys.which("sha256sum")) || nzchar(Sys.which("shasum")) || .Platform$OS.type == "windows")
                expect_identical(orthologr:::diamond_sha256(f, use_r = FALSE), abc)
})

test_that("diamond_install() installs a download whose checksum matches", {
        archive <- fake_diamond_archive()
        with_fake_cache({
                calls <- with_mocked_download(
                        msgs <- capture_messages(res <- diamond_install()),
                        archive)
                expect_match(msgs, "^Downloading https://github.com/bbuchfink/diamond/releases/download/v2.2.8/", all = FALSE)
                expect_match(msgs, "SHA-256 checksum verified", all = FALSE)
                expect_match(msgs, "^Installed diamond version 2.2.8 in ", all = FALSE)
                expect_equal(calls$n, 1L)
                expect_identical(calls$url, "https://github.com/bbuchfink/diamond/releases/download/v2.2.8/diamond-linux64.tar.gz")
                expect_equal(calls$timeout, max(600, getOption("timeout")))
                expect_equal(res, orthologr:::diamond_cache_dir())
                expect_true(file.exists(file.path(res, "diamond")))
                # the temporary archive and the unpacked copy are gone
                expect_length(list.files(tempdir(), pattern = "^diamond_install_"), 0)
        })
})

test_that("diamond_install() installs a Windows zip whose checksum matches", {
        skip_if_not(nzchar(Sys.which(zip_program())), "no zip program")
        archive <- fake_diamond_archive(zip = TRUE)
        with_fake_cache({
                with_mocked_download(res <- diamond_install(quiet = TRUE), archive, zip = TRUE)
                expect_true(file.exists(file.path(res, "diamond.exe")))
        })
})

test_that("a download with another checksum is deleted and nothing is installed", {
        archive <- fake_diamond_archive()
        wrong <- strrep("0", 64)
        with_fake_cache({
                e <- NULL
                with_mocked_download(e <- tryCatch(diamond_install(quiet = TRUE), error = function(e) e), archive, pinned_sha = wrong)
                expect_s3_class(e, "orthologr_error_checksum_mismatch")
                expect_s3_class(e, "orthologr_error_tool_failed")
                expect_match(conditionMessage(e), wrong)
                expect_match(conditionMessage(e), orthologr:::diamond_sha256(archive))
                expect_false(dir.exists(orthologr:::diamond_cache_dir()) &&
                             file.exists(file.path(orthologr:::diamond_cache_dir(), "diamond")))
                expect_length(list.files(tempdir(), pattern = "^diamond_install_"), 0)
        })
})

test_that("a failed reinstallation (bad checksum) leaves the cached copy in place", {
        archive <- fake_diamond_archive()
        with_fake_cache({
                dir <- place_fake_diamond()
                writeLines("old copy", file.path(dir, "diamond"))
                with_mocked_download(expect_error(diamond_install(force = TRUE, quiet = TRUE), class = "orthologr_error_checksum_mismatch"),
                                     archive, pinned_sha = strrep("a", 64))
                expect_identical(readLines(file.path(dir, "diamond")), "old copy")
        })
})

test_that("a download that times out is an error that names the timeout, and the partial file is removed", {
        archive <- fake_diamond_archive()
        partial <- NULL
        with_fake_cache({
                e <- NULL
                calls <- with_mocked_download(e <- tryCatch(diamond_install(timeout = 7, quiet = TRUE), error = function(e) e),
                                              archive, fail = function(url, file, timeout) {
                                                      writeBin(as.raw(1:10), file)          # a partial file
                                                      partial <<- file
                                                      stop("Timeout of ", timeout, " seconds was reached")
                                              })
                expect_equal(calls$timeout, 7)
                expect_s3_class(e, "orthologr_error_download_failed")
                expect_match(conditionMessage(e), "did not finish within 7 seconds")
                expect_false(file.exists(partial))
                expect_false(file.exists(file.path(orthologr:::diamond_cache_dir(), "diamond")))
                # a warning from download.file (e.g. a truncated transfer) is a failed download too
                e2 <- NULL
                with_mocked_download(e2 <- tryCatch(diamond_install(timeout = 7, quiet = TRUE), error = function(e) e),
                                     archive, fail = function(url, file, timeout) {
                                             warning("downloaded length 10 != reported length 10221459")
                                     })
                expect_s3_class(e2, "orthologr_error_download_failed")
                expect_match(conditionMessage(e2), "downloaded length")
        })
})

test_that("the timeout of the download is R's option 'timeout' during the transfer only", {
        seen <- NULL
        old <- getOption("timeout")
        orthologr:::diamond_download("https://example.invalid/x", tempfile(), timeout = 12, quiet = TRUE,
                                     download = function(url, destfile, ...) { seen <<- getOption("timeout"); 0L })
        expect_equal(seen, 12)
        expect_equal(getOption("timeout"), old)
        # a real transfer (a local file) with the default download function
        src <- tempfile()
        writeLines("x", src)
        dst <- tempfile()
        orthologr:::diamond_download(paste0("file:///", sub("^/", "", normalizePath(src, winslash = "/"))), dst, timeout = 5, quiet = TRUE)
        expect_identical(readLines(dst), "x")
        expect_equal(getOption("timeout"), old)
})

test_that("a release without a pinned checksum is not downloaded unless 'sha256' is given", {
        archive <- fake_diamond_archive()
        with_fake_cache({
                calls <- with_mocked_download(
                        expect_error(diamond_install(version = "v2.1.25", quiet = TRUE), class = "orthologr_error_bad_argument"),
                        archive)
                expect_equal(calls$n, 0L)
                e <- tryCatch(with_mocked_download(diamond_install(version = "2.1.25", quiet = TRUE), archive), error = function(e) e)
                expect_match(conditionMessage(e), "no pinned SHA-256 checksum for DIAMOND v2.1.25")
                expect_match(conditionMessage(e), "sha256 = ")
                # with the checksum given, the release is installed
                sha <- orthologr:::diamond_sha256(archive)
                calls <- with_mocked_download(res <- diamond_install(version = "v2.1.25", sha256 = toupper(sha), quiet = TRUE),
                                              archive, pinned_sha = strrep("b", 64))
                expect_identical(calls$url, "https://github.com/bbuchfink/diamond/releases/download/v2.1.25/diamond-linux64.tar.gz")
                expect_true(file.exists(file.path(res, "diamond")))
        })
})

test_that("'sha256' must agree with a pinned checksum and be 64 hexadecimal digits; bad arguments are refused before a download", {
        archive <- fake_diamond_archive()
        with_fake_cache({
                calls <- with_mocked_download({
                        expect_error(diamond_install(sha256 = strrep("c", 64), quiet = TRUE), "is not the checksum orthologr has pinned",
                                     class = "orthologr_error_bad_argument")
                        expect_error(diamond_install(sha256 = "abc", quiet = TRUE), class = "orthologr_error_bad_argument")
                        expect_error(diamond_install(sha256 = c(strrep("c", 64), strrep("d", 64)), quiet = TRUE), class = "orthologr_error_bad_argument")
                        expect_error(diamond_install(timeout = 0, quiet = TRUE), class = "orthologr_error_bad_argument")
                        expect_error(diamond_install(timeout = NA, quiet = TRUE), class = "orthologr_error_bad_argument")
                        expect_error(diamond_install(version = NA_character_, quiet = TRUE), class = "orthologr_error_bad_argument")
                }, archive)
                expect_equal(calls$n, 0L)
        })
})

test_that("version = 'latest' is the newest pinned release that has a file for the platform", {
        pins <- data.frame(version = c("v2.2.8", "v2.10.1", "v2.9.0"), asset = c("diamond-linux64.tar.gz", "diamond-linux64.tar.gz", "diamond-macos.tar.gz"),
                           size = 1, sha256 = c(strrep("1", 64), strrep("2", 64), strrep("3", 64)), stringsAsFactors = FALSE)
        seen <- NULL
        local_mocked_bindings(diamond_platform = function() c(sysname = "Linux", machine = "x86_64"),
                              diamond_release_checksums = function() pins,
                              diamond_download = function(url, file, timeout, quiet = FALSE) { seen <<- url; stop("offline") },
                              .package = "orthologr")
        with_fake_cache({
                expect_error(diamond_install(quiet = TRUE), class = "orthologr_error_download_failed")
                expect_match(seen, "/v2.10.1/diamond-linux64.tar.gz$")
        })
})

test_that("a binary that does not run is not copied into the cache", {
        archive <- fake_diamond_archive()
        with_fake_cache({
                local_mocked_bindings(diamond_binary_version = function(exe) NULL, .package = "orthologr")
                local_mocked_bindings(diamond_platform = function() c(sysname = "Linux", machine = "x86_64"),
                                      diamond_exe_name = function() "diamond",
                                      diamond_release_checksums = function() data.frame(version = "v2.2.8", asset = "diamond-linux64.tar.gz",
                                                                                        size = 1, sha256 = orthologr:::diamond_sha256(archive),
                                                                                        stringsAsFactors = FALSE),
                                      diamond_download = function(url, file, timeout, quiet = FALSE) file.copy(archive, file),
                                      .package = "orthologr")
                expect_error(diamond_install(quiet = TRUE), "could not be executed")
                expect_false(file.exists(file.path(orthologr:::diamond_cache_dir(), "diamond")))
        })
})

test_that("the real 'diamond version' check is given a time limit", {
        skip_on_os("windows")
        bin <- tempfile("slow_diamond_")
        dir.create(bin)
        writeLines(c("#!/bin/sh", "echo 'diamond version 2.2.8'"), file.path(bin, "diamond"))
        Sys.chmod(file.path(bin, "diamond"), "755")
        expect_match(orthologr:::diamond_binary_version(file.path(bin, "diamond"))[1], "diamond version 2.2.8")
        expect_match(paste(deparse(body(orthologr:::diamond_binary_version)), collapse = " "), "timeout = 60")
})

# ---------------------------------------------------------------------------
# Online tests: download the real binary from GitHub. They only run
# interactively / on a local machine (devtools::test(), devtools::check()),
# never on CRAN or CI, because they need network access and a supported OS.
# ---------------------------------------------------------------------------

# On CI only where ORTHOLOGR_TEST_DIAMOND_DOWNLOAD=true (set in .github/workflows/R-CMD-check.yaml): the real,
# pinned release is downloaded, verified and run there.
skip_if_no_download <- function() {
        skip_on_cran()
        if (!identical(Sys.getenv("ORTHOLOGR_TEST_DIAMOND_DOWNLOAD"), "true"))
                skip_on_ci()
        skip_if_offline("github.com")
        info <- Sys.info()
        supported <- info[["sysname"]] %in% c("Darwin", "Windows") ||
                (info[["sysname"]] == "Linux" && info[["machine"]] %in% c("x86_64", "amd64"))
        skip_if_not(supported, "No pre-compiled DIAMOND binary for this platform")
}

test_that("diamond_install() downloads a working DIAMOND binary", {
        skip_if_no_download()
        with_fake_cache({
                dir <- suppressMessages(diamond_install(quiet = TRUE))
                exe <- file.path(dir, orthologr:::diamond_exe_name())
                expect_true(file.exists(exe))
                out <- system2(exe, "version", stdout = TRUE, stderr = TRUE)
                expect_match(out[1], "diamond version")
                # the DIAMOND wrappers now find it even without diamond on the PATH
                expect_true(is_installed_diamond(diamond_exec_path = dir))
        })
})

test_that("diamond_install(force = TRUE) re-downloads", {
        skip_if_no_download()
        with_fake_cache({
                dir <- place_fake_diamond()
                exe <- file.path(dir, orthologr:::diamond_exe_name())
                expect_equal(file.size(exe), 0)
                suppressMessages(diamond_install(force = TRUE, quiet = TRUE))
                expect_gt(file.size(exe), 0)
        })
})

test_that("diamond_install() fails clearly for a non-existing version", {
        skip_if_no_download()
        with_fake_cache({
                # (refused before any download: no pinned checksum; with one given, the download fails)
                expect_error(diamond_install(version = "v0.0.0", quiet = TRUE), class = "orthologr_error_bad_argument")
                expect_error(suppressWarnings(diamond_install(version = "v0.0.0", sha256 = strrep("0", 64), quiet = TRUE)),
                             "Could not download DIAMOND", class = "orthologr_error_download_failed")
        })
})
