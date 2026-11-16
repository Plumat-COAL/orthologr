# The fixes of REVIEW-OL-3 in the COAL, method and status code: one table of the dN/dS methods (dnds_method_registry());
# "GY" with the COAL engine, kaks_tool = "coal", dNdS_clusters() and dNdS_orthogroups(engine = "coal") when the COAL
# at hand has it (M1); the status "not_converged" (M3); the messages about input sequences of orthogroups() and
# dNdS_orthogroups() (M4); the summary of dNdS_orthogroups() on the statistics helper of dNdS_clusters(); the
# multiple alignment program in a folder; the codon aligner of dNdS_orthogroups(engine = "orthologr").

ol3b_coalr <- requireNamespace("coalr", quietly = TRUE)
ol3b_gy <- ol3b_coalr && is.element("gy94_ml", tryCatch(coalr::coal_info()$features, error = function(e) character(0)))

# A pair whose YN iteration does not converge in COAL (found by a random search; KaKs_Calculator returns the
# values of its last round for it).
ol3b_nc <- c(q = "AACGTTTGCGGTATTCACCAC", s = "AACCAATTCGCTATTCGAGAC")
ol3b_nc_file <- function() {
        f <- tempfile("ol3b_nc_", fileext = ".fasta")
        writeLines(c(">q", ol3b_nc[["q"]], ">s", ol3b_nc[["s"]]), f)
        f
}
ol3b_skip_unless_nc <- function() {
        skip_if_not(ol3b_coalr, "coalr is not installed")
        st <- tryCatch(coalr::coal_kaks(ol3b_nc[["q"]], ol3b_nc[["s"]], method = "YN")$status[1], error = function(e) NA)
        skip_if_not(identical(st, "not_converged"), "this COAL converges for the test pair")
}

# The orthologr options and the environment variables a test changes, restored by the function returned.
ol3b_state <- function() {
        ours <- function() grep("^orthologr\\.", names(options()), value = TRUE)
        opt <- options()[ours()]
        env <- Sys.getenv(c("COAL_CPP", "FAKE_COAL_MODE"), unset = NA, names = TRUE)
        function() {
                added <- setdiff(ours(), names(opt))
                if (length(added) > 0)
                        options(stats::setNames(rep(list(NULL), length(added)), added))
                options(opt)
                for (n in names(env)) {
                        if (is.na(env[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(env[[n]]), n))
                }
        }
}

ol3b_fake_coal <- function(mode = "ok") {
        testthat::skip_on_os("windows")
        testthat::skip_if_not(nzchar(Sys.which("sh")), "no POSIX shell")
        dir <- new_test_dir("ol3b_fake_coal_")
        exe <- file.path(dir, "coal-cpp")
        file.copy(testthat::test_path("fixtures", "coverage", "fake-coal-cpp.sh"), exe)
        Sys.chmod(exe, "0755")
        options(orthologr.coal_backend = "binary", orthologr.coal_path = exe)
        Sys.setenv(FAKE_COAL_MODE = mode)
        normalizePath(exe)
}

ol3b_messages <- function(expr) {
        msgs <- character(0)
        value <- withCallingHandlers(expr, message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
        })
        list(value = value, messages = msgs)
}

# ---- one table of the methods (quality 1) ---------------------------------------------------------------------

test_that("every list of dN/dS methods is a view of dnds_method_registry()", {
        r <- dnds_method_registry()
        expect_identical(lengths(r), stats::setNames(rep(length(r$method), length(r)), names(r)))
        expect_false(anyDuplicated(r$method) > 0)
        expect_identical(dnds_methods(), r$method)
        # the lists that were written out by hand before, as sets
        expect_setequal(kaks_calculator_methods(),
                        c("MA", "MS", "NG", "LWL", "LPB", "MLWL", "MLPB", "YN", "MYN", "GY", "GMYN", "kaks_calc", "ALL"))
        expect_setequal(dnds_methods(coal = "always"), c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"))
        expect_identical(coal_native_methods(), c("Comeron", "Li"))
        expect_setequal(coal_kaks_methods(), c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN", "GY"))
        expect_setequal(kaks_lwl_methods(), c("LWL", "MLWL", "LPB", "MLPB"))
        expect_setequal(kaks_yn_methods(), c("YN", "MYN"))
        expect_setequal(kaks_counted_methods(), c("NG", "YN", "MYN", "GMYN"))
        # is.dnds_est_method(): the names of the table, and the list without an argument
        expect_true(all(vapply(r$method, is.dnds_est_method, TRUE)))
        expect_false(is.dnds_est_method("yn"))
        expect_message(expect_false(is.dnds_est_method()), "The following methods are provided: Comeron, Li")
        expect_identical(dnds_method_hint("gy"), " Did you mean \"GY\"?")
        # the message of engine = "orthologr" names every method once
        txt <- dnds_method_list(kaks_calc = TRUE)
        for (m in r$method)
                expect_true(grepl(paste0("\"", m, "\""), txt, fixed = TRUE), info = m)
})

# ---- GY with COAL (M1) ----------------------------------------------------------------------------------------

test_that("coal_has_gy() asks the COAL at hand: coalr's features, or the help text of coal-cpp", {
        restore <- ol3b_state()
        on.exit(restore(), add = TRUE)
        if (ol3b_coalr)
                expect_identical(coal_has_gy(), is.element("gy94_ml", coalr::coal_info()$features))
        ol3b_fake_coal("ok")
        expect_true(coal_has_gy())
        expect_true(is.element("GY", coal_methods()))
        expect_invisible(coal_check_methods("GY", "kaks_tool = \"coal\""))
        ol3b_fake_coal("nogy")
        expect_false(coal_has_gy())
        expect_true(coal_has_native_methods())
        expect_false(is.element("GY", coal_methods()))
        err <- expect_error(coal_check_methods(c("YN", "GY"), "kaks_tool = \"coal\"", other = "kaks_tool = \"kaks_calculator\""),
                            class = "orthologr_error_package_missing")
        expect_match(conditionMessage(err), "dnds_est.method = \"GY\" needs a COAL that computes Goldman & Yang's", fixed = TRUE)
        expect_match(conditionMessage(err), "coal-cpp program in use does not", fixed = TRUE)
        expect_match(conditionMessage(err), coalr_install_hint(), fixed = TRUE)
})

test_that("without the feature gy94_ml, GY is refused by every COAL entry point with a message that says why", {
        skip_if_not(ol3b_coalr, "coalr is not installed")
        testthat::local_mocked_bindings(coal_has_gy = function() FALSE, .package = "orthologr")
        f <- ol3b_nc_file()
        for (call in list(
                quote(dNdS_coal(thal_cds(), lyra_cds(), dnds_est.method = "GY")),
                quote(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "GY", engine = "coal")),
                quote(dNdS_clusters(thal_cds(), clusters = data.frame(a = "x", b = "y"), dnds_est.method = "GY")),
                quote(dnds_groups_coal_check("GY", "pairwise")),
                quote(substitutionrate(f, est.method = "GY", kaks_tool = "coal")),
                quote(coal_kaks_file(f, "GY")))) {
                err <- expect_error(eval(call), class = "orthologr_error_package_missing")
                expect_match(conditionMessage(err), "gy94_ml", fixed = TRUE)
        }
        # (the version is never asked: a development coalr reports 0.6.0.9000 without GY)
        expect_false(grepl("packageVersion", paste(deparse(coal_has), collapse = "\n"), fixed = TRUE))
})

test_that("GY with kaks_tool = \"coal\": COAL's values, in process and with the coal-cpp program", {
        skip_if_not(ol3b_gy, "coalr without GY (gy94_ml)")
        f <- ol3b_nc_file()
        r <- substitutionrate(f, est.method = "GY", kaks_tool = "coal", quiet = TRUE)
        k <- coalr::coal_kaks(ol3b_nc[["q"]], ol3b_nc[["s"]], method = "GY")
        expect_identical(r$method, "GY")
        expect_identical(c(r$dN, r$dS, r$dNdS), c(k$dN, k$dS, k$dNdS))
        expect_identical(attr(r, "engine_status"), k$status)
})

test_that("GY with engine = \"coal\" and kaks_tool = \"coal\" give the same values (20 genes)", {
        skip_if_not(ol3b_gy && diamond_available, "coalr with GY or DIAMOND is not installed")
        a <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "GY", engine = "coal", print_citation = FALSE))
        b <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "GY", kaks_tool = "coal", print_citation = FALSE))
        expect_identical(a$method[1], "GY")
        expect_true(all(is.element(a$status, dnds_status_levels())))
        a <- a[order(a$query_id), ]
        b <- b[order(b$query_id), ]
        expect_identical(a$query_id, b$query_id)
        expect_equal(a$dNdS, b$dNdS, tolerance = 1e-10)
        expect_identical(a$status, b$status)
        expect_false(is.element(".engine_status", names(b)))
        expect_identical(attr(b, "orthologr_run")$arguments$kaks_tool, "coal")
})

test_that("GY in dNdS_clusters() and dNdS_orthogroups(engine = \"coal\")", {
        skip_if_not(ol3b_gy, "coalr without GY (gy94_ml)")
        cds <- c(a = ol3b_nc[["q"]], b = ol3b_nc[["s"]], c = "ATGGCTAAAGGTCTGGAACGT", d = "ATGGCGAAAGGCCTGGATCGT")
        f <- tempfile("ol3b_cl_", fileext = ".fasta")
        writeLines(as.vector(rbind(paste0(">", names(cds)), cds)), f)
        res <- suppressMessages(dNdS_clusters(f, clusters = data.frame(representative_id = c("c", "c"), member_id = c("c", "d")),
                                              dnds_est.method = c("YN", "GY"), aa_aln_tool = "pairwise"))
        expect_setequal(res$pairs$method, c("YN", "GY"))
        expect_true(all(is.element(res$pairs$status, dnds_status_levels())))
        got <- dnds_groups_coal(pairs = data.table::data.table(orthogroup = "OG1", query_id = "c", subject_id = "d"),
                                cds = cds, dnds_est.method = "GY", aa_aln_tool = "pairwise")
        expect_identical(got$method, "GY")
        expect_identical(got$engine_status, "ok")
})

# ---- not_converged (M3) ---------------------------------------------------------------------------------------

test_that("dnds_status(): COAL's not_converged is a status of its own, counted and summarised", {
        expect_true(is.element("not_converged", dnds_status_levels()))
        st <- dnds_status(c(NA, NA, 0.1, NA, NA), c(NA, NA, 0.2, NA, NA), c(NA, NA, 0.5, NA, NA),
                          engine_status = c("not_converged", "saturated", "ok", "no_data", "error"))
        expect_identical(st, c("not_converged", "saturated", "ok", "no_data", "failed"))
        # a failed pair stays failed; without the engine's word the values decide (KaKs_Calculator)
        expect_identical(dnds_status(NA, NA, NA, failed = TRUE, engine_status = "not_converged"), "failed")
        expect_identical(dnds_status(NA, NA, NA), "saturated")
        expect_identical(dnds_status_counts(st)[["not_converged"]], 1L)
        txt <- dnds_summary_text(st, dS = c(NA, NA, 0.2, NA, NA))
        expect_match(txt, "1 not converged (no estimate)", fixed = TRUE)
        expect_false(grepl("see the notes of the engine", txt, fixed = TRUE))
        # the count of a note grep is gone
        expect_false(is.element("not_converged", names(formals(dnds_summary_text))))
        expect_false(is.element("not_converged", names(formals(dnds_finish))))
})

test_that("kaks_tool = \"coal\": COAL's status leaves coal_kaks_file() and reaches the column 'status' of dNdS()", {
        ol3b_skip_unless_nc()
        f <- ol3b_nc_file()
        r <- substitutionrate(f, est.method = "YN", kaks_tool = "coal", quiet = TRUE)
        expect_identical(attr(r, "engine_status"), "not_converged")
        expect_true(all(is.na(c(r$dN, r$dS, r$dNdS))))
        # compute_dnds() carries it as '.engine_status'
        tbl <- data.table::data.table(query_id = "q", subject_id = "s", query_cds = ol3b_nc[["q"]], subject_cds = ol3b_nc[["s"]],
                                      query_aa = transl(ol3b_nc[["q"]]), subject_aa = transl(ol3b_nc[["s"]]))
        old <- options(orthologr.tmpdir = new_test_dir("ol3b_cd_"))
        on.exit(options(old), add = TRUE)
        res <- suppressMessages(compute_dnds(tbl, aa_aln_type = "pairwise", aa_aln_tool = "NW", codon_aln_tool = "coal",
                                             dnds_est.method = "YN", kaks_tool = "coal", quiet = TRUE))
        expect_identical(res$.engine_status, "not_converged")
})

test_that("dNdS(kaks_tool = \"coal\") reports a pair that did not converge as not_converged, without a helper column", {
        skip_if_not(ol3b_coalr && diamond_available, "coalr or DIAMOND is not installed")
        # every pair "does not converge"
        testthat::local_mocked_bindings(coal_kaks_file = function(file, est.method, subst_name = NULL) {
                res <- data.table::data.table(query_id = "q", subject_id = "s", dN = NA_real_, dS = NA_real_, dNdS = NA_real_,
                                              method = est.method)
                data.table::setattr(res, "engine_status", "not_converged")
                res
        }, .package = "orthologr")
        run <- ol3b_messages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "YN", kaks_tool = "coal", print_citation = FALSE))
        res <- run$value
        expect_true(all(res$status == "not_converged"))
        expect_false(is.element(".engine_status", names(res)))
        expect_identical(names(res)[ncol(res)], "status")
        expect_identical(attr(res, "orthologr_run")$status_counts$not_converged, nrow(res))
        expect_match(run$messages[length(run$messages)], paste0(nrow(res), " not converged (no estimate)"), fixed = TRUE)
})

test_that("kaks_tool = \"coal\" with KaKs_Calculator's methods only: the other routes have no engine status", {
        skip_if_not(ol3b_coalr && diamond_available, "coalr or DIAMOND is not installed")
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_false(is.element(".engine_status", names(res)))
        expect_false(is.element("not_converged", res$status))
})

test_that("the binary backend carries COAL's status too", {
        restore <- ol3b_state()
        on.exit(restore(), add = TRUE)
        ol3b_fake_coal("notconv")
        r <- coal_kaks_file(ol3b_nc_file(), "YN")
        expect_identical(attr(r, "engine_status"), "not_converged")
        expect_true(is.na(r$dNdS))
        ol3b_fake_coal("ok")
        expect_identical(attr(coal_kaks_file(ol3b_nc_file(), "GY"), "engine_status"), "ok")
})

test_that("dNdS_clusters(): the pairs, the summary counts and the message agree on not_converged", {
        ol3b_skip_unless_nc()
        cds <- c(a = ol3b_nc[["q"]], b = ol3b_nc[["s"]])
        f <- tempfile("ol3b_cl_", fileext = ".fasta")
        writeLines(as.vector(rbind(paste0(">", names(cds)), cds)), f)
        run <- ol3b_messages(dNdS_clusters(f, clusters = data.frame(representative_id = c("a", "a"), member_id = c("a", "b")),
                                           dnds_est.method = "YN", aa_aln_tool = "pairwise"))
        res <- run$value
        st <- coalr::coal_dnds_pairs(query_cds = cds["a"], subject_cds = cds["b"], pairs = data.frame(query_id = "a", subject_id = "b"),
                                     methods = "YN", preset = "pwalign")$status
        skip_if_not(identical(st, "not_converged"), "COAL aligns the test pair so that YN converges")
        expect_identical(res$pairs$status, "not_converged")
        expect_identical(res$summary$n_not_converged, 1L)
        expect_identical(res$summary$n_saturated, 0L)
        expect_identical(res$summary$status, "not_converged")
        expect_match(run$messages[length(run$messages)], "1 not converged (no estimate)", fixed = TRUE)
})

test_that("dnds_cluster_status_counts(): the counts of the summary are those of the column 'status' of the pairs", {
        summary <- data.frame(cluster_id = c("c1", "c1", "c2", "c3"), method = c("YN", "NG", "YN", "YN"),
                              n_ok = c(9L, 9L, 9L, 4L), n_not_converged = 0L, n_saturated = 0L, n_no_data = 0L,
                              n_failed = c(0L, 0L, 0L, 2L), stringsAsFactors = FALSE)
        pairs <- data.frame(cluster_id = c("c1", "c1", "c1", "c2"), method = c("YN", "YN", "NG", "YN"))
        out <- dnds_cluster_status_counts(summary, pairs, c("ok", "not_converged", "saturated", "failed"))
        expect_identical(out$n_ok, c(1L, 0L, 0L, 4L))
        expect_identical(out$n_not_converged, c(1L, 0L, 0L, 0L))
        expect_identical(out$n_saturated, c(0L, 1L, 0L, 0L))
        expect_identical(out$n_failed, c(0L, 0L, 1L, 2L))       # c3 has no pairs in the table: coalr's counts stay
})

test_that("dNdS_orthogroups(engine = \"coal\"): COAL's status comes from its column, not from the note", {
        ol3b_skip_unless_nc()
        got <- dnds_groups_coal(pairs = data.table::data.table(orthogroup = "OG1", query_id = "a", subject_id = "b"),
                                cds = c(a = ol3b_nc[["q"]], b = ol3b_nc[["s"]]), dnds_est.method = "YN", aa_aln_tool = "pairwise")
        expect_true(is.element(got$engine_status, c("ok", "warning", "not_converged", "saturated", "no_data")))
        # a pair COAL returns no row for is an error of COAL
        testthat::local_mocked_bindings(coal_orthologr_values = function(tbl) tbl[0, , drop = FALSE], .package = "orthologr")
        lost <- dnds_groups_coal(pairs = data.table::data.table(orthogroup = "OG1", query_id = "a", subject_id = "b"),
                                 cds = c(a = ol3b_nc[["q"]], b = ol3b_nc[["s"]]), dnds_est.method = "YN", aa_aln_tool = "pairwise")
        expect_identical(lost$engine_status, "error")
        expect_identical(dnds_status(lost$dN, lost$dS, lost$dNdS, engine_status = lost$engine_status), "failed")
})

# ---- the summary of dNdS_orthogroups() (quality 5) ------------------------------------------------------------

test_that("dnds_groups_summary(): the statistics of dnds_value_stats(), the counts by status, the columns of 0.6.0", {
        set.seed(11)
        n <- 60
        pairs <- data.table::data.table(orthogroup = sample(sprintf("OG%02d", 1:7), n, TRUE), method = "YN",
                                        dN = round(stats::runif(n, 0, 0.5), 6), dS = round(stats::runif(n, 0, 3), 6))
        pairs$dN[sample(n, 8)] <- NA
        pairs$dS[sample(n, 8)] <- NA
        pairs$dS[sample(n, 4)] <- 0
        pairs$dNdS <- ifelse(is.finite(pairs$dN / pairs$dS), pairs$dN / pairs$dS, NA_real_)
        members <- data.table::data.table(orthogroup = rep(sprintf("OG%02d", 1:7), each = 3),
                                          species = rep(c("A", "B", "C"), 7), gene = paste0("g", 1:21))
        status <- dnds_status(pairs$dN, pairs$dS, pairs$dNdS)
        status[c(1, 2)] <- "not_converged"
        s <- dnds_groups_summary(pairs, members, dS_saturation = 2, status = status)
        expect_identical(names(s)[1:15], c("orthogroup", "method", "n_genes", "n_species", "n_pairs", "n_ok", "n_na", "n_dS_zero",
                                           "n_saturated", "n_no_data", "n_not_converged", "n_failed", "n_dS_above", "n_dN", "n_dS"))
        # the statistics as written out in 0.6.0 (quantile() for the median): the same numbers
        q <- function(v, p) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, p, na.rm = TRUE, names = FALSE))
        for (g in unique(s$orthogroup)) {
                p <- pairs[pairs$orthogroup == g]
                row <- s[s$orthogroup == g]
                for (v in c("dN", "dS", "dNdS")) {
                        expect_identical(row[[paste0(v, "_median")]], q(p[[v]], 0.5))
                        expect_identical(row[[paste0(v, "_q25")]], q(p[[v]], 0.25))
                        expect_identical(row[[paste0(v, "_q75")]], q(p[[v]], 0.75))
                        expect_identical(row[[paste0(v, "_mean")]], if (all(is.na(p[[v]]))) NA_real_ else mean(p[[v]], na.rm = TRUE))
                }
                expect_identical(row$n_ok, sum(status[pairs$orthogroup == g] == "ok"))
                expect_identical(row$n_na, sum(is.na(p$dNdS)))
        }
        expect_identical(sum(s$n_not_converged), 2L)
        expect_identical(sum(s$n_pairs), as.integer(n))
})

# ---- messages about input sequences (M4) ----------------------------------------------------------------------

test_that("orthogroups() and dNdS_orthogroups(): one line about the input sequences, the tables in the result", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        files <- c(thal = system.file("seqs/ortho_thal_cds_1000.fasta", package = "orthologr"),
                   lyra = system.file("seqs/ortho_lyra_cds_1000.fasta", package = "orthologr"))
        run <- ol3b_messages(orthogroups(files, comp_cores = 2))
        og <- run$value
        tr <- attr(og, "trimmed_sequences")
        expect_true(is.data.frame(tr) && nrow(tr) > 0)
        expect_false(any(grepl("There seem to be", run$messages, fixed = TRUE)))
        expect_identical(sum(grepl("^orthogroups\\(\\): input coding sequences - ", run$messages)), 1L)
        expect_lte(length(run$messages), 6)
        # quiet = FALSE: the messages of every file as before
        loud <- ol3b_messages(orthogroups(files, comp_cores = 2, quiet = FALSE))
        expect_true(any(grepl("There seem to be", loud$messages, fixed = TRUE)))
        # dNdS_orthogroups() from the files: no message of read.cds(), the trimmed sequences counted in the summary
        skip_if_not(ol3b_coalr, "coalr is not installed")
        run <- ol3b_messages(dNdS_orthogroups(files, engine = "coal", aa_aln_tool = "pairwise", comp_cores = 2))
        expect_false(any(grepl("There seem to be|All corrupted CDS", run$messages)))
        expect_lte(length(run$messages), 6)
        expect_true(is.data.frame(attr(run$value, "trimmed_sequences")))
        expect_match(run$messages[length(run$messages)], "have a length that is not a multiple of 3", fixed = TRUE)
        expect_identical(names(run$value), c("pairs", "summary", "failed"))
})

test_that("input_report_text() says nothing when there is nothing to say", {
        expect_null(input_report_text(list(owner = TRUE), "orthogroups()"))
        rep <- list(owner = TRUE, trimmed = data.frame(file = c("a", "a"), id = c("x", "y")),
                    stops = data.frame(file = "b", id = "z", stop_codons = 1L))
        txt <- input_report_text(rep, "orthogroups()")
        expect_match(txt, "^orthogroups\\(\\): input coding sequences - 2 with a length that is not a multiple of 3")
        expect_match(txt, "1 with in-frame stop codons", fixed = TRUE)
        res <- input_report_attributes(list(a = 1), rep)
        expect_identical(attr(res, "trimmed_sequences"), rep$trimmed)
        expect_null(attr(res, "dropped_sequences"))
})

# ---- lows ------------------------------------------------------------------------------------------------------

test_that("coal_msa_aligner_for(): a folder is searched as for every other program (tool_command())", {
        skip_if_not(ol3b_coalr, "coalr is not installed")
        seen <- NULL
        testthat::local_mocked_bindings(tool_command = function(program, path = NULL) {
                seen <<- c(program, path)
                file.path(path, program)
        }, .package = "orthologr")
        dir <- new_test_dir("ol3b_msa_")
        tryCatch(coal_msa_aligner_for("clustalo", dir), error = function(e) NULL)
        expect_identical(seen, c("clustalo", dir))
        err <- expect_error(coal_msa_aligner_for("clustalo", file.path(dir, "not-there")), class = "orthologr_error_tool_missing")
        expect_match(conditionMessage(err), "aa_aln_tool = \"pairwise\"", fixed = TRUE)
})

test_that("dNdS_orthogroups(engine = \"orthologr\") builds the codon alignments with the codon aligner of dNdS()", {
        restore <- ol3b_state()
        on.exit(restore(), add = TRUE)
        options(orthologr.codon_aln_tool = "pal2nal")
        expect_identical(dnds_groups_codon_aln_tool(), "pal2nal")
        options(orthologr.codon_aln_tool = "coal", orthologr.coal_backend = "none")
        expect_identical(dnds_groups_codon_aln_tool(), "pal2nal")
        options(orthologr.coal_backend = NULL, orthologr.codon_aln_tool = NULL)
        expect_identical(dnds_groups_codon_aln_tool(), if (ol3b_coalr) "coal" else "pal2nal")
        expect_false(any(grepl("\"pal2nal\"", deparse(dNdS_orthogroups), fixed = TRUE)))
})

test_that("user-facing errors of the method and status code have orthologr's classes", {
        f <- ol3b_nc_file()
        expect_error(substitutionrate(f, est.method = "Comeron", format = "clustal"), class = "orthologr_error_bad_argument")
        expect_error(kaks_calc_args("kaks_calc"), class = "orthologr_error_bad_argument")
        expect_error(kaks_result_values(data.frame()), class = "orthologr_error_tool_failed")
        expect_error(write_orthogroups(list(), tempdir()), class = "orthologr_error_bad_argument")
})

test_that("the coalr of DESCRIPTION's Remotes and of the install hints is COAL's release v0.6.0", {
        expect_identical(coalr_install_hint(), "remotes::install_github(\"Plumat-COAL/COAL/packaging/r/coalr@v0.6.0\")")
        remotes <- utils::packageDescription("orthologr")$Remotes
        skip_if(is.null(remotes), "the installed DESCRIPTION has no Remotes field")
        expect_match(remotes, "coalr=Plumat-COAL/COAL/packaging/r/coalr@v0.6.0", fixed = TRUE)
        expect_match(coal_install_message(), coalr_install_hint(), fixed = TRUE)
})
