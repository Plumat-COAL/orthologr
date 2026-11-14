# Coverage of the COAL glue (R/coal.R): the backend choice and its options, the "binary" backend (the coal-cpp
# program, played here by a POSIX shell script in fixtures/coverage), the small helpers, and the argument checks
# and error paths of dNdS_coal() / dNdS_clusters(). The tests of the in-process backend need the R package coalr;
# the tests of the binary backend need a POSIX shell (skipped on Windows). No other program is needed.

cov_coalr <- requireNamespace("coalr", quietly = TRUE)

# The orthologr options and the environment variables these tests change, restored by the function returned:
#     restore <- cov_state(); on.exit(restore(), add = TRUE)
cov_state <- function(vars = c("COAL_CPP", "FAKE_COAL_MODE")) {
        ours <- function() grep("^orthologr\\.", names(options()), value = TRUE)
        opt <- options()[ours()]
        env <- Sys.getenv(vars, unset = NA, names = TRUE)
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

# A copy of the fake coal-cpp in a folder of its own (each test gets a fresh path, so the answers orthologr
# remembers per program path never leak between tests), chosen as the COAL backend. Use after cov_state().
local_fake_coal <- function(mode = "ok") {
        testthat::skip_on_os("windows")
        testthat::skip_if_not(nzchar(Sys.which("sh")), "no POSIX shell")
        dir <- new_test_dir("fake_coal_")
        exe <- file.path(dir, "coal-cpp")
        file.copy(testthat::test_path("fixtures", "coverage", "fake-coal-cpp.sh"), exe)
        Sys.chmod(exe, "0755")
        options(orthologr.coal_backend = "binary", orthologr.coal_path = exe)
        Sys.setenv(FAKE_COAL_MODE = mode)
        normalizePath(exe)
}

cov_write_fasta <- function(seqs, ids, file = tempfile("cov_", fileext = ".fasta")) {
        writeLines(as.vector(rbind(paste0(">", ids), seqs)), file)
        file
}

cov_cds_a <- "ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACC"
cov_cds_b <- "ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTGAAGGCGTGGTTGACC"

# ---------------------------------------------------------------------------------------------------------------
# choosing the backend

test_that("coal_backend() checks options(orthologr.coal_backend) and says how to install COAL", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        options(orthologr.coal_backend = "fortran")
        expect_error(coal_backend(), class = "orthologr_error_bad_argument")
        expect_false(is_installed_coal())
        expect_false(coal_has_native_methods())
        expect_identical(coal_versions(), list())

        options(orthologr.coal_backend = "none")
        err <- expect_error(coal_backend(), class = "orthologr_error_package_missing")
        expect_s3_class(err, "orthologr_error_tool_missing")
        expect_match(conditionMessage(err), "pip install coal-align", fixed = TRUE)
        expect_identical(conditionMessage(err), coal_install_message())
        # methods without COAL: only KaKs_Calculator's seven
        expect_identical(coal_methods(), coal_kaks_calculator_methods())
})

test_that("coal_backend(): the binary backend is used when forced, and 'coalr' is asked for by name", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        options(orthologr.coal_backend = "binary", orthologr.coal_path = "")
        Sys.setenv(COAL_CPP = "")
        # (a coal-cpp on the PATH would be found: only then is there a binary)
        if (!nzchar(Sys.which("coal-cpp"))) {
                expect_null(coal_binary())
                expect_error(coal_backend(), "COAL is not available", class = "orthologr_error_tool_missing")
        }
        # need = "coalr" with the binary forced: the package is asked for by name
        expect_error(coal_backend(need = "coalr"), "needs the R package 'coalr'", class = "orthologr_error_package_missing")
})

test_that("coal_binary(): options(orthologr.coal_path), then COAL_CPP; folders and missing files are not taken", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        dir <- new_test_dir("coal_bin_")
        exe <- file.path(dir, "coal-cpp")
        writeLines("#!/bin/sh", exe)
        Sys.setenv(COAL_CPP = "")
        options(orthologr.coal_path = exe)
        expect_identical(coal_binary(), normalizePath(exe))
        # a folder is not a program
        options(orthologr.coal_path = dir)
        if (!nzchar(Sys.which("coal-cpp")))
                expect_null(coal_binary())
        # the environment variable when the option names nothing
        options(orthologr.coal_path = file.path(dir, "not-there"))
        Sys.setenv(COAL_CPP = exe)
        expect_identical(coal_binary(), normalizePath(exe))
        # a malformed option is passed over
        options(orthologr.coal_path = c("a", "b"))
        expect_identical(coal_binary(), normalizePath(exe))
})

test_that("coal_backend(forced = 'coalr') without coalr asks for the package", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        skip_if(cov_coalr, "coalr is installed (the error is only given without it)")
        options(orthologr.coal_backend = "coalr")
        expect_error(coal_backend(), "needs the R package 'coalr'", class = "orthologr_error_package_missing")
})

# ---------------------------------------------------------------------------------------------------------------
# the binary backend (fake coal-cpp)

test_that("the binary backend: version, native methods (from the help text) and the record of a run", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        exe <- local_fake_coal("ok")
        expect_true(is_installed_coal())
        expect_identical(coal_backend(), list(type = "binary", path = exe))
        expect_true(coal_has_native_methods())
        # (the answer is remembered: asked again it comes from the cache)
        expect_true(coal_has_native_methods())
        expect_identical(coal_methods(), c(coal_kaks_calculator_methods(), coal_native_methods()))
        v <- coal_versions()
        expect_identical(v$coal_backend, "binary")
        expect_identical(v$`coal-cpp`, "coal-cpp 0.0.0-fake")
        expect_invisible(coal_check_methods(c("YN", "Comeron"), "kaks_tool = \"coal\""))
})

test_that("the binary backend: a coal-cpp without Comeron and Li says what is needed", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        local_fake_coal("nonative")
        expect_false(coal_has_native_methods())
        expect_identical(coal_methods(), coal_kaks_calculator_methods())
        expect_error(coal_check_methods("Li", "engine = \"coal\""), "needs a COAL that computes Comeron's and Li's method",
                     class = "orthologr_error_package_missing")
})

test_that("coal_check_methods() checks the methods it is given", {
        expect_error(coal_check_methods(1, "x"), class = "orthologr_error_bad_argument")
        expect_error(coal_check_methods(character(0), "x"), class = "orthologr_error_bad_argument")
        expect_error(coal_check_methods(NA_character_, "x"), class = "orthologr_error_bad_argument")
        expect_error(coal_check_methods(c("YN", "GY", "ML"), "engine = \"coal\"", other = "the other engine"),
                     "\"GY\", \"ML\" is only available with the other engine", class = "orthologr_error_bad_argument")
        expect_invisible(coal_check_methods(coal_kaks_calculator_methods(), "x"))
})

test_that("the binary backend: substitutionrate's Ka/Ks of one codon alignment (coal_kaks_file)", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        local_fake_coal("ok")
        f <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("q-1", "s-1"))
        res <- coal_kaks_file(f, "YN", subst_name = "p1")
        expect_s3_class(res, "data.table")
        expect_identical(names(res), c("query_id", "subject_id", "dN", "dS", "dNdS", "method"))
        expect_identical(res$query_id, "q-1")
        expect_identical(res$subject_id, "s-1")
        expect_equal(c(res$dN, res$dS, res$dNdS), c(0.1, 0.5, 0.2))
        expect_identical(res$method, "YN")
        # Comeron and Li: no 'method' column, as substitutionrate() returns them
        res <- coal_kaks_file(f, "Comeron")
        expect_identical(names(res), c("query_id", "subject_id", "dN", "dS", "dNdS"))
        # the input and output files of the call are removed
        expect_length(list.files(file.path(orthologr_tmpdir(), "_calculation"), pattern = "coal\\.tsv$"), 0)
})

test_that("the binary backend: Comeron's undefined correction is NA unless the legacy option is set", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        local_fake_coal("undefined")
        f <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("q", "s"))
        res <- coal_kaks_file(f, "Comeron")
        expect_equal(res$dN, 0.1235)
        expect_true(is.na(res$dS) && is.na(res$dNdS))
        options(orthologr.legacy_comeron_saturation = TRUE)
        res <- coal_kaks_file(f, "Comeron")
        expect_equal(c(res$dN, res$dS, res$dNdS), c(0.1235, 1.5, 0.0823))
})

test_that("the binary backend: errors of coal-cpp name the gene pair", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        f <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("qq", "ss"))
        local_fake_coal("error")
        expect_error(coal_kaks_file(f, "YN"), "COAL could not estimate dN/dS for gene pair 'qq' / 'ss': sequences of unequal length")
        Sys.setenv(FAKE_COAL_MODE = "empty")
        expect_error(coal_kaks_file(f, "YN"), "COAL returned no result for gene pair 'qq' / 'ss'")
        Sys.setenv(FAKE_COAL_MODE = "fail")
        expect_error(coal_kaks_file(f, "YN"), "failed on purpose", class = "orthologr_error_tool_failed")
        Sys.setenv(FAKE_COAL_MODE = "nooutput")
        expect_error(coal_kaks_file(f, "YN"), "wrote no output", class = "orthologr_error_alignment_failed")
        # a file that does not hold two sequences
        three <- cov_write_fasta(c(cov_cds_a, cov_cds_b, cov_cds_a), c("a", "b", "c"))
        Sys.setenv(FAKE_COAL_MODE = "ok")
        expect_error(coal_kaks_file(three, "YN"), "must hold two sequences")
        # an unknown method is refused before COAL is started
        expect_error(coal_kaks_file(f, "GY"), "kaks_tool = \"kaks_calculator\"", class = "orthologr_error_bad_argument")
})

test_that("the binary backend: codon alignments with PAL2NAL's command line (coal_codon_aln_file)", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        local_fake_coal("ok")
        aln <- cov_write_fasta(c("MKTAYIAKQRLEGDLKAWLT", "MKTAYIAKQRLDGDLKAWLT"), c("a", "b"))
        nuc <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("a", "b"))
        out <- file.path(new_test_dir("cov_aln_"), "pair.aln")
        writeLines("an alignment left by an earlier call", out)
        expect_identical(coal_codon_aln_file(aln, nuc, out), out)
        expect_identical(readLines(out), readLines(nuc))
        # the temporary file is gone
        expect_identical(list.files(dirname(out)), "pair.aln")
        Sys.setenv(FAKE_COAL_MODE = "fail")
        expect_error(coal_codon_aln_file(aln, nuc, out, context = "the test pair"), "for the test pair",
                     class = "orthologr_error_alignment_failed")
        # the old result was removed, not left behind as if it were the new one
        expect_false(file.exists(out))
})

test_that("codon_aln(tool = \"coal\") runs the binary backend", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        local_fake_coal("ok")
        aln <- cov_write_fasta(c("MKTAYIAKQRLEGDLKAWLT", "MKTAYIAKQRLDGDLKAWLT"), c("a", "b"))
        nuc <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("a", "b"))
        res <- suppressMessages(codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", tool = "coal",
                                          get_aln = TRUE, quiet = TRUE))
        expect_false(is.null(res))
})

# ---------------------------------------------------------------------------------------------------------------
# small helpers

test_that("read_fasta_strings() and write_alignment() (FASTA and CLUSTAL)", {
        f <- tempfile(fileext = ".fasta")
        writeLines(c(">id1 a description", "ACGT", "ACG", ">id2", "aaaa"), f)
        x <- read_fasta_strings(f)
        expect_identical(x, c(id1 = "ACGTACG", id2 = "aaaa"))

        out <- tempfile(fileext = ".fasta")
        expect_identical(write_alignment(x, out), out)
        expect_identical(readLines(out), c(">id1", "ACGTACG", ">id2", "aaaa"))

        long <- c(seq_one = strrep("ACG", 30), s2 = strrep("ACT", 30))
        expect_identical(write_alignment(long, out, format = "clustal"), out)
        lines <- readLines(out)
        expect_identical(lines[1], "CLUSTAL W multiple sequence alignment")
        # 90 columns: two blocks of 60 and 30, each name padded to the longest + 8 characters
        expect_identical(lines[3], paste0(formatC("seq_one", width = -15), substr(long[[1]], 1, 60)))
        expect_identical(lines[7], paste0(formatC("s2", width = -15), substr(long[[2]], 61, 90)))
        expect_length(lines, 8)
})

test_that("cds_vector() turns a CDS table into the named vector COAL takes", {
        dt <- data.table::data.table(geneids = c("g1", "g2"), seqs = c("ATG", "TTT"))
        expect_identical(cds_vector(dt), c(g1 = "ATG", g2 = "TTT"))
})

test_that("check_aa_aln_scoring() and the matrices of the COAL engine", {
        expect_silent(check_aa_aln_scoring("pwalign"))
        expect_silent(check_aa_aln_scoring("needle"))
        expect_error(check_aa_aln_scoring("smith-waterman"), "aa_aln_scoring")
        expect_identical(coal_sub_matrix_names(), c("BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250"))
        # only pwalign scoring with BLOSUM62 may become pwalign's copy of the matrix
        expect_identical(coal_matrix_arg("BLOSUM45", "pwalign"), "BLOSUM45")
        expect_identical(coal_matrix_arg("BLOSUM62", "needle"), "BLOSUM62")
        have_pw <- is.element("BLOSUM62-PWALIGN", tryCatch(coalr::coal_info()$matrices, error = function(e) character(0)))
        expect_identical(coal_matrix_arg("BLOSUM62", "pwalign"), if (have_pw) "BLOSUM62-PWALIGN" else "BLOSUM62")
})

test_that("coal_legacy_warning() warns only with the 0.4.3 alignments", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        options(orthologr.legacy_pairwise_aln = FALSE)
        expect_silent(coal_legacy_warning("x()"))
        options(orthologr.legacy_pairwise_aln = TRUE)
        expect_warning(coal_legacy_warning("x()"), "x\\(\\): the COAL engine cannot reproduce",
                       class = "orthologr_warning_legacy_engine")
})

test_that("coal_flavour_args() warns when a coalr without the Comeron options is asked for Comeron", {
        f_mid <- function(a, b, comeron_flavour, li_flavour, terminal_stops) NULL
        expect_warning(a <- coal_flavour_args(f_mid, methods = "Comeron"), class = "orthologr_warning_coal_version")
        expect_identical(names(a), c("comeron_flavour", "li_flavour", "terminal_stops"))
        expect_silent(coal_flavour_args(f_mid, methods = "YN"))
        f_new <- function(a, comeron_flavour, comeron_twosubs, comeron_min_sites, li_flavour, terminal_stops) NULL
        expect_silent(a <- coal_flavour_args(f_new, methods = "Comeron"))
        expect_identical(sort(names(a)), sort(c("comeron_flavour", "comeron_twosubs", "comeron_min_sites", "li_flavour",
                                               "terminal_stops")))
})

test_that("coal_engine_settings() records the coalr call and the matching COAL command line", {
        f_tie <- function(query_cds, methods, preset, matrix, tie_break, comeron_flavour, li_flavour, terminal_stops) NULL
        f_old <- function(query_cds, methods, preset) NULL

        s <- coal_engine_settings("coal_dnds_pairs", f_tie, list(methods = c("YN", "NG"), preset = "pwalign"))
        expect_identical(s$coalr_function, "coal_dnds_pairs")
        expect_identical(s$tie_break, "pwalign")
        expect_identical(s$matrix, "BLOSUM62-PWALIGN")
        expect_identical(s$comeron_flavour, "gestimator")
        expect_identical(s$terminal_stops, "drop")
        expect_match(s$coal_command, "^coal pipeline pairs --query QUERY_CDS")
        expect_match(s$coal_command, "-m YN -m NG --preset pwalign --matrix BLOSUM62-PWALIGN --tie-break pwalign",
                     fixed = TRUE)
        expect_match(s$coal_command, "--comeron-flavour gestimator --li-flavour seqinr --terminal-stops drop -o OUT_DIR",
                     fixed = TRUE)

        s <- coal_engine_settings("coal_dnds_clusters", f_tie, list(methods = "YN", preset = "needle", matrix = "PAM250"),
                                  mode = "clusters")
        expect_identical(s$tie_break, "leftmost")
        expect_identical(s$matrix, "PAM250")
        expect_match(s$coal_command, "^coal pipeline clusters --cds CDS --clusters CLUSTERS.tsv -m YN --preset needle --matrix PAM250")

        s <- coal_engine_settings("coal_dnds_pairs", f_old, list(methods = "YN", preset = "pwalign"))
        expect_identical(s$tie_break, "leftmost (coalr without tie_break)")
        expect_identical(s$matrix, "BLOSUM62")
        expect_null(s$comeron_flavour)
        expect_false(grepl("--tie-break|--comeron-flavour", s$coal_command))
        expect_error(coal_engine_settings("f", f_old, list(methods = "YN", preset = "x"), mode = "trios"))
})

test_that("coal_orthologr_values(): tables without rows or methods, and Li's saturation marker", {
        expect_identical(nrow(coal_orthologr_values(data.frame(method = character(0)))), 0L)
        no_method <- data.frame(dN = 1, dS = 2)
        expect_identical(coal_orthologr_values(no_method), no_method)
        tbl <- data.frame(method = c("Li", "Li", "Li", "YN", NA), dN = c(0.1, -0.2, 0.3, 5, 1),
                          dS = c(0.2, 0.4, 10, -1, 1), dNdS = c(0.5, -0.5, 0.03, -5, 1), stringsAsFactors = FALSE)
        r <- coal_orthologr_values(tbl)
        expect_equal(r$dN, c(0.1, NA, 0.3, 5, 1))
        expect_equal(r$dS, c(0.2, 0.4, NA, -1, 1))
        expect_equal(r$dNdS, c(0.5, NA, NA, -5, 1))
        # Comeron without a 'note' column: rounded to 4 significant digits only
        r <- coal_orthologr_values(data.frame(method = "Comeron", dN = 0.123456789, dS = 0.987654321, dNdS = 0.125))
        expect_equal(c(r$dN, r$dS, r$dNdS), c(0.1235, 0.9877, 0.125))
        # Comeron with a missing note and a correction undefined at nonsynonymous sites
        r <- coal_orthologr_values(data.frame(method = c("Comeron", "Comeron"), dN = c(2, 0.2), dS = c(1, 1),
                                              dNdS = c(2, 0.2), note = c("Ka: correction undefined", NA),
                                              stringsAsFactors = FALSE))
        expect_equal(r$dN, c(NA, 0.2))
        expect_equal(r$dNdS, c(NA, 0.2))
})

test_that("dnds_cluster_status() words clusters without an estimate by what their pairs are", {
        expect_identical(dnds_cluster_status(data.frame()), data.frame())
        partial <- data.frame(status = "failed", n_pairs = 1)
        expect_identical(dnds_cluster_status(partial), partial)
        sm <- data.frame(cluster_id = c("c1", "c2", "c3", "c4", "c5", "c6"),
                         status = c("ok", "failed", "failed", "failed", "failed", "skipped"),
                         n_pairs = c(3, 2, 1, 2, 3, 0), n_ok = c(3, 0, 0, 0, 0, 0),
                         n_saturated = c(0, 2, 0, 0, 1, 0), n_no_data = c(0, 0, 1, 0, 0, 0),
                         n_failed = c(0, 0, 0, 0, 1, 0), n_not_converged = c(0, 0, 0, 2, 1, 0),
                         stringsAsFactors = FALSE)
        out <- dnds_cluster_status(sm)
        expect_identical(out$status, c("ok", "saturated", "no_data", "not_converged", "failed", "skipped"))
        expect_identical(out$note[2], "no pair has an estimate (2 pairs: 2 saturated)")
        expect_identical(out$note[3], "no pair has an estimate (1 pair: 1 without codons to compare)")
        expect_identical(out$note[5], "no pair has an estimate (3 pairs: 1 saturated, 1 not converged, 1 failed)")
        expect_true(is.na(out$note[1]))
        # a note that is there is kept
        sm$note <- c(NA, "COAL's words", NA, NA, NA, NA)
        expect_identical(dnds_cluster_status(sm)$note[2], "COAL's words")
})

test_that("dnds_cluster_summary_stats() leaves summaries it cannot use as they are", {
        pairs <- data.frame(cluster_id = c("c1", "c1"), method = "YN", dN = c(0.1, NA), dS = c(0.2, 0.4),
                            dNdS = c(0.5, NA))
        sm <- data.frame(cluster_id = "c1", method = "YN", n_pairs = 2, dN_median = 9, dS_median = 9,
                         dNdS_median = 9, dN_min = 9, dS_max = 9)
        expect_identical(dnds_cluster_summary_stats(sm[0, ], pairs), sm[0, ])
        expect_identical(dnds_cluster_summary_stats(sm, pairs[0, ]), sm)
        expect_identical(dnds_cluster_summary_stats(sm[ , c("cluster_id", "n_pairs")], pairs), sm[ , c("cluster_id", "n_pairs")])
        expect_identical(dnds_cluster_summary_stats(sm[ , 1:3], pairs), sm[ , 1:3])
        out <- dnds_cluster_summary_stats(sm, pairs)
        expect_identical(names(out), c("cluster_id", "method", "n_pairs", "n_dN", "n_dS", "n_dNdS", "dN_median",
                                       "dS_median", "dNdS_median", "dN_min", "dS_max"))
        expect_equal(unlist(out[1, c("n_dN", "n_dS", "n_dNdS")]), c(n_dN = 1, n_dS = 2, n_dNdS = 1))
        expect_equal(c(out$dN_median, out$dS_median, out$dN_min, out$dS_max), c(0.1, 0.3, 0.1, 0.4))
        # a cluster without pairs gets counts of 0
        sm2 <- rbind(sm, transform(sm, cluster_id = "c2"))
        out <- dnds_cluster_summary_stats(sm2, pairs)
        expect_equal(out$n_dN, c(1, 0))
        expect_true(is.na(out$dN_median[2]))
})

test_that("dnds_value_stats() computes every statistic it knows, NA for groups without values", {
        pairs <- data.frame(g = c("a", "a", "a", "b"), dN = c(0.1, 0.2, 0.3, NA), dS = c(1, 2, 3, NA),
                            dNdS = c(0.1, 0.1, 0.1, NA))
        st <- as.data.frame(dnds_value_stats(pairs, by = "g", stats = c("median", "mean", "q1", "q3", "q25", "q75",
                                                                         "min", "max")))
        expect_equal(st$dS_q1, c(1.5, NA))
        expect_equal(st$dS_q3, c(2.5, NA))
        expect_equal(st$dN_min, c(0.1, NA))
        expect_equal(st$dN_max, c(0.3, NA))
        expect_equal(st$n_dN, c(3L, 0L))
})

# ---------------------------------------------------------------------------------------------------------------
# argument checks and error paths of the entry points (before any search)

test_that("dNdS_coal() checks its arguments before the search", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        q <- thal_cds()
        s <- lyra_cds()
        run <- function(...) dNdS_coal(query_file = q, subject_file = s, print_citation = FALSE, ...)
        expect_error(run(aligner = "mmseqs2"), "aligner = \"diamond\" or \"blast\"", class = "orthologr_error_bad_argument")
        expect_error(run(ortho_detection = "OrthoFinder"), "ortho_detection", class = "orthologr_error_bad_argument")
        expect_error(run(seq_type = "protein"), "seq_type = \"cds\"", class = "orthologr_error_bad_argument")
        expect_error(run(format = "gbk"), "format")
        expect_error(run(aa_aln_scoring = "x"), "aa_aln_scoring")
        expect_error(run(aa_sub_matrix = "BLOSUM50"), "aa_sub_matrix")
        options(orthologr.coal_backend = "none")
        expect_error(run(), "COAL is not available")
})

test_that("dNdS_clusters() checks its arguments", {
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        cds <- thal_cds()
        expect_error(dNdS_clusters(c(cds, file.path(tempdir(), "no-such.fasta"))), "no-such.fasta",
                     class = "orthologr_error_bad_input")
        expect_error(dNdS_clusters(cds, aa_aln_tool = "muscle"), "aa_aln_tool", class = "orthologr_error_bad_argument")
        expect_error(dNdS_clusters(cds, large_clusters = "all"), "large_clusters", class = "orthologr_error_bad_argument")
        expect_error(dNdS_clusters(cds, aa_aln_scoring = "x"), "aa_aln_scoring")
        options(orthologr.coal_backend = "none")
        expect_error(dNdS_clusters(cds, clusters = data.frame(a = 1, b = 2)), "COAL is not available")
})

test_that("dNdS_clusters() checks the clusters it is given (coalr)", {
        skip_if_not(cov_coalr, "the R package coalr is not installed")
        cds <- thal_cds()
        expect_error(dNdS_clusters(cds, clusters = file.path(tempdir(), "no-such-clusters.tsv")), "was not found",
                     class = "orthologr_error_bad_input")
        expect_error(dNdS_clusters(cds, clusters = 1:3), "'clusters' must be", class = "orthologr_error_bad_argument")
        expect_error(dNdS_clusters(cds, clusters = data.frame(a = "x")), "'clusters' must be",
                     class = "orthologr_error_bad_argument")
        expect_error(dNdS_clusters(cds, clusters = data.frame(r = "x", m = "x"), dnds_est.method = "GY"),
                     "only available with dNdS\\(\\)", class = "orthologr_error_bad_argument")
})

# ---------------------------------------------------------------------------------------------------------------
# the in-process backend (coalr)

test_that("coal_codon_aln_file() with coalr: CLUSTAL input, and errors with the class of PAL2NAL's", {
        skip_if_not(cov_coalr, "the R package coalr is not installed")
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        options(orthologr.coal_backend = "coalr")
        prot <- c(a = "MKTAYIAKQRLEGDLKAWLT", b = "MKTAYIAKQRLDGDLKAWLT")
        aln_fa <- cov_write_fasta(prot, names(prot))
        # (a CLUSTAL file as ClustalW writes it, with a conservation line)
        aln_cl <- tempfile("cov_", fileext = ".aln")
        writeLines(c("CLUSTAL W (1.83) multiple sequence alignment", "", "",
                     "a        MKTAYIAKQRLEGDLKAWLT", "b        MKTAYIAKQRLDGDLKAWLT",
                     "         *********** ********", ""), aln_cl)
        nuc <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("a", "b"))
        out_fa <- tempfile("cov_out_", fileext = ".fasta")
        out_cl <- tempfile("cov_out_", fileext = ".fasta")
        coal_codon_aln_file(aln_fa, nuc, out_fa)
        coal_codon_aln_file(aln_cl, nuc, out_cl)
        expect_identical(readLines(out_cl), readLines(out_fa))
        expect_identical(read_fasta_strings(out_fa), c(a = cov_cds_a, b = cov_cds_b))
        # the two files must hold the same number of sequences
        nuc3 <- cov_write_fasta(c(cov_cds_a, cov_cds_b, cov_cds_a), c("a", "b", "c"))
        expect_error(coal_codon_aln_file(aln_fa, nuc3, out_fa), "holds 2 sequences and the CDS file 3",
                     class = "orthologr_error_alignment_failed")
        # a protein alignment that is not the translation of the CDS
        bad <- cov_write_fasta(c("WWWWWWWWWWWWWWWWWWWW", "WWWWWWWWWWWWWWWWWWWW"), c("a", "b"))
        res <- tryCatch(coal_codon_aln_file(bad, nuc, out_fa), error = function(e) e)
        if (inherits(res, "error")) {
                expect_s3_class(res, "orthologr_error_alignment_failed")
                expect_match(conditionMessage(res), "COAL could not build the codon alignment for gene pair 'a' / 'b'")
        } else {
                expect_true(file.exists(out_fa))
        }
})

cov_toy_clusters <- function(toy, families = 1:4) {
        tr <- toy$truth[toy$truth$family %in% families, ]
        rep_id <- tapply(tr$gene, tr$family, function(g) g[1])
        data.frame(representative_id = unname(rep_id[as.character(tr$family)]), member_id = tr$gene,
                   stringsAsFactors = FALSE)
}

test_that("dNdS_clusters() with clusters given as a table or a file, aligned pairwise by COAL (coalr)", {
        skip_if_not(cov_coalr, "the R package coalr is not installed")
        toy <- make_toy_species()
        # (unnamed: a named vector of files is a known bug, see the test below)
        toy$files <- unname(toy$files)
        cl <- cov_toy_clusters(toy)
        res <- suppressMessages(dNdS_clusters(toy$files, clusters = tibble::as_tibble(cl), aa_aln_tool = "pairwise",
                                              dnds_est.method = c("YN", "NG")))
        expect_identical(names(res), c("pairs", "summary", "failed"))
        # four families in four species: 6 pairs each, two methods
        expect_identical(nrow(res$pairs), 48L)
        expect_true(all(is.element(c("query_id", "subject_id", "dN", "dS", "dNdS", "status"), names(res$pairs))))
        expect_identical(nrow(res$summary), 8L)
        expect_true(all(is.element(c("n_dN", "n_dS", "n_dNdS"), names(res$summary))))
        expect_identical(nrow(res$failed), 0L)
        rec <- attr(res, "orthologr_run")
        expect_false(is.null(rec))

        # the same clusters from a tab-separated file
        f <- tempfile("cov_clusters_", fileext = ".tsv")
        utils::write.table(cl, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        res_f <- suppressMessages(dNdS_clusters(toy$files, clusters = f, aa_aln_tool = "pairwise",
                                                dnds_est.method = c("YN", "NG")))
        expect_equal(res_f$pairs$dN, res$pairs$dN)

        # a member without a coding sequence is reported in 'failed', with a message
        cl2 <- rbind(cl, data.frame(representative_id = cl$representative_id[1], member_id = "no-such-gene"))
        msgs <- testthat::capture_messages(res2 <- dNdS_clusters(toy$files, clusters = cl2, aa_aln_tool = "pairwise"))
        expect_true(any(grepl("could not be used", msgs)))
        expect_true(nrow(res2$failed) >= 1)
})

test_that("dNdS_clusters(): a multiple alignment program that cannot be found suggests aa_aln_tool = \"pairwise\" (coalr)", {
        skip_if_not(cov_coalr, "the R package coalr is not installed")
        toy <- make_toy_species()
        toy$files <- unname(toy$files)
        cl <- cov_toy_clusters(toy, 1)
        missing_tool <- file.path(new_test_dir("cov_no_mafft_"), "mafft")
        expect_error(suppressMessages(dNdS_clusters(toy$files, clusters = cl, aa_aln_tool = "mafft",
                                                    aa_aln_path = missing_tool)),
                     "aa_aln_tool = \"pairwise\"")
        # a folder for aa_aln_path is the folder of the program
        expect_error(suppressMessages(dNdS_clusters(toy$files, clusters = cl, aa_aln_tool = "clustalo",
                                                    aa_aln_path = dirname(missing_tool))),
                     "aa_aln_tool = \"pairwise\"")
})

# ---------------------------------------------------------------------------------------------------------------
# bugs found by the coverage pass (B1, B2; fixed in R/coal.R)

test_that("B1: dNdS_clusters() with a named vector of CDS files finds the coding sequences (coalr)", {
        skip_if_not(cov_coalr, "the R package coalr is not installed")
        toy <- make_toy_species()
        cl <- cov_toy_clusters(toy, 1:2)
        named <- suppressMessages(dNdS_clusters(toy$files, clusters = cl, aa_aln_tool = "pairwise"))
        unnamed <- suppressMessages(dNdS_clusters(unname(toy$files), clusters = cl, aa_aln_tool = "pairwise"))
        expect_false(is.null(names(toy$files)))
        expect_identical(NROW(named$failed), 0L)
        expect_gt(nrow(named$pairs), 0)
        expect_identical(named$pairs, unnamed$pairs)
})

test_that("B2: coal_codon_aln_file() reads a CLUSTAL protein alignment without a conservation line (coalr)", {
        skip_if_not(cov_coalr, "the R package coalr is not installed")
        restore <- cov_state()
        on.exit(restore(), add = TRUE)
        options(orthologr.coal_backend = "coalr")
        prot <- c(a = "MKTAYIAKQRLEGDLKAWLT", b = "MKTAYIAKQRLDGDLKAWLT")
        # the CLUSTAL format orthologr itself writes (as PAL2NAL; T-Coffee and MUSCLE write no conservation line either)
        aln <- write_alignment(prot, tempfile("cov_", fileext = ".aln"), format = "clustal")
        nuc <- cov_write_fasta(c(cov_cds_a, cov_cds_b), c("a", "b"))
        out <- tempfile("cov_out_", fileext = ".fasta")
        coal_codon_aln_file(aln, nuc, out)
        expect_identical(read_fasta_strings(out), c(a = cov_cds_a, b = cov_cds_b))
})
