# The COAL parts of the third review: Comeron and Li through COAL, one alignment scoring for every entry point
# (R15), the rule of the summaries (R16), attributes of dNdS(engine = "coal") (R21), dN = 0 (R14).

coalr_available <- requireNamespace("coalr", quietly = TRUE)
native_in_coal <- coalr_available && coal_has_native_methods()
frame <- function(x, cols = c("query_id", "subject_id", "dN", "dS", "dNdS")) {
        x <- as.data.frame(x)[ , cols]
        x <- x[order(x$query_id, x$subject_id), ]
        rownames(x) <- NULL
        x
}

test_that("coal_orthologr_values() applies orthologr's rules to COAL's Comeron and Li values", {
        tbl <- data.frame(
                method = c("Comeron", "Comeron", "Comeron", "Comeron", "Li", "Li", "Li", "YN"),
                dN   = c(0.123456, 0.2, 0.3, NA, 0.05, -0.001, 0.1, 0.123456),
                dS   = c(0.654321, 0.5, 0.6, 0.4, 0.5, 0.2, 0, 0.654321),
                dNdS = c(0.123456 / 0.654321, 0.4, 0.5, NA, 0.1, -0.005, Inf, 0.188678),
                status = c("ok", "warning", "warning", "saturated", "ok", "warning", "ok", "ok"),
                note = c("", "Ka: correction undefined, gestimator's value (undefined logarithms set to 0)",
                         "Ks: correction undefined, gestimator's value (undefined logarithms set to 0)", "Ka saturated", "", "", "", ""),
                stringsAsFactors = FALSE)
        r <- coal_orthologr_values(tbl)
        expect_equal(r$dN[1:4], c(0.1235, NA, 0.3, NA))            # 4 significant digits; NA where Ka is undefined
        expect_equal(r$dS[1:4], c(0.6543, 0.5, NA, 0.4))
        expect_equal(r$dNdS[1:4], c(signif(0.123456 / 0.654321, 4), NA, NA, NA))
        expect_equal(r$dN[5:7], c(0.05, NA, 0.1))                  # Li: a negative value is NA ...
        expect_equal(r$dNdS[5:7], c(0.1, NA, NA))                  # ... and so is a ratio that is not finite
        expect_equal(unlist(r[8, c("dN", "dS", "dNdS")], use.names = FALSE), c(0.123456, 0.654321, 0.188678))  # other methods untouched
        old <- options(orthologr.legacy_comeron_saturation = TRUE)
        on.exit(options(old), add = TRUE)
        r <- coal_orthologr_values(tbl)
        expect_equal(r$dN[1:3], c(0.1235, 0.2, 0.3))               # the finite values of 0.4.3
        expect_equal(r$dS[1:3], c(0.6543, 0.5, 0.6))
        expect_equal(r$dNdS[2:3], c(0.4, 0.5))
})

test_that("a COAL without Comeron and Li gives a message that names what is needed", {
        testthat::local_mocked_bindings(coal_has_native_methods = function() FALSE, .package = "orthologr")
        expect_identical(coal_methods(), coal_kaks_calculator_methods())
        expect_error(coal_check_methods("Comeron", "engine = \"coal\""),
                     "engine = \"coal\": dnds_est.method = \"Comeron\" needs a COAL that computes Comeron's and Li's method \\(coalr >= 0.5.0")
        expect_error(coal_check_methods(c("YN", "Li"), "engine = \"coal\""), "\"Li\" needs a COAL")
        expect_error(coal_check_methods("MA", "engine = \"coal\""), "\"MA\" is only available with engine = \"orthologr\"")
        expect_silent(coal_check_methods(c("YN", "NG"), "engine = \"coal\""))
})

test_that("coal_flavour_args() asks for orthologr's numbers and only passes what the function takes", {
        f_new <- function(a, comeron_flavour, li_flavour, terminal_stops) NULL
        f_old <- function(a, yn_flavour) NULL
        expect_identical(coal_flavour_args(f_new), list(comeron_flavour = "gestimator", li_flavour = "seqinr", terminal_stops = "drop"))
        expect_identical(coal_flavour_args(f_old), structure(list(), names = character(0)))
        old <- options(orthologr.legacy_pairwise_aln = TRUE)
        on.exit(options(old), add = TRUE)
        expect_identical(coal_flavour_args(f_new)$terminal_stops, "compare")
})

test_that("substitutionrate(kaks_tool = \"coal\") returns the table of the built-in Comeron and Li code", {
        skip_if_not(native_in_coal, "coalr with Comeron and Li is not installed")
        d <- new_test_dir("coal_native_")
        files <- c(
                pal2nal = system.file("seqs/pal2nal.aln", package = "orthologr"),
                # one synonymous substitution: dS is saturated in the smallest possible alignment
                tiny = file.path(d, "tiny.fa"),
                # 30 codons with a saturated nonsynonymous class (R17)
                short = file.path(d, "short.fa"),
                # a class without sites and an undefined logarithm: not a saturation (R17)
                nosites = file.path(d, "nosites.fa"))
        writeLines(c(">a", "ATGGCTAAA", ">b", "ATGGCCAAA"), files[["tiny"]])
        writeLines(c(">a", "ACGCCAAAACTCGTTCTGTGGAAAAGACTAACACATGTCAATTGGCAGACCAGCACCTCCTTTCAGTTACAGGCGGCTTTAGTTCAAAGG",
                     ">b", "ACGCCAAAACTCGTTCTGCGGAAAACACAAACACATGTCAAGTGGCAGACCAGCACCTCCTTTCCGTTACAGGCGGCATATGTTCAAAGG"), files[["short"]])
        writeLines(c(">a", "GGTGCCCTAAACGTCGAAAAAGAGTGCGGA", ">b", "AAAGCCCTAAATGGGGTCCTAGTGTGCTCA"), files[["nosites"]])
        for (legacy in c(FALSE, TRUE)) {
                old <- options(orthologr.legacy_comeron_saturation = legacy)
                for (n in names(files)) for (m in c("Comeron", "Li")) {
                        native <- tryCatch(substitutionrate(files[[n]], est.method = m, quiet = TRUE), error = function(e) NULL)
                        if (is.null(native)) next
                        coal <- substitutionrate(files[[n]], est.method = m, quiet = TRUE, kaks_tool = "coal")
                        expect_identical(names(coal), names(native), info = paste(n, m, legacy))
                        expect_equal(as.data.frame(coal), as.data.frame(native), tolerance = 1e-12, info = paste(n, m, legacy))
                }
                options(old)
        }
        # the NA pattern that is expected with the default rule
        expect_true(is.na(substitutionrate(files[["tiny"]], est.method = "Comeron", quiet = TRUE, kaks_tool = "coal")$dS))
        expect_true(is.na(substitutionrate(files[["short"]], est.method = "Comeron", quiet = TRUE, kaks_tool = "coal")$dN))
        expect_equal(substitutionrate(files[["nosites"]], est.method = "Comeron", quiet = TRUE, kaks_tool = "coal")$dN, 0.6615)
})

test_that("dNdS() with COAL's tools on its own alignments equals the built-in Comeron and Li for every pair", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        skip_if_not(native_in_coal, "coalr with Comeron and Li is not installed")
        run <- function(...) suppressMessages(dNdS(thal_cds(), lyra_cds(), quiet = TRUE, print_citation = FALSE, ...))
        for (m in c("Comeron", "Li")) {
                native <- run(dnds_est.method = m)
                tools <- run(dnds_est.method = m, codon_aln_tool = "coal", kaks_tool = "coal")
                expect_identical(names(tools), names(native), info = m)
                expect_equal(as.data.frame(tools), as.data.frame(native), tolerance = 1e-12, ignore_attr = TRUE, info = m)
                # the batched engine: the same table (columns, types, attributes); values equal wherever COAL's
                # protein alignment is the one pairwise_aln() makes (co-optimal alignments can differ)
                engine <- run(dnds_est.method = m, engine = "coal")
                expect_identical(names(engine), names(native), info = m)
                expect_identical(vapply(engine, function(v) class(v)[1], ""), vapply(native, function(v) class(v)[1], ""), info = m)
                expect_identical(frame(engine)[ , 1:2], frame(native)[ , 1:2], info = m)
                same <- abs(frame(engine)$dN - frame(native)$dN) < 1e-12 & abs(frame(engine)$dS - frame(native)$dS) < 1e-12
                expect_true(sum(same, na.rm = TRUE) >= 10, info = m)
                expect_true(all(c("failed_pairs", "dropped_sequences", "coal") %in% names(attributes(engine))), info = m)
        }
        # the default method of dNdS() works with the COAL engine
        def <- run(engine = "coal")
        expect_identical(frame(def), frame(run(dnds_est.method = "Comeron", engine = "coal")))
        expect_false(is.element("method", names(def)))
})

test_that("R21: dNdS(engine = \"coal\") carries the attributes of dNdS() and joins by both IDs", {
        skip_on_cran()
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        d <- new_test_dir("r21_")
        # one coding sequence whose length is not a multiple of 3: it is left out and reported, as by dNdS()
        q <- readLines(thal_cds())
        first_seq <- which(!startsWith(q, ">"))[1]
        q[first_seq] <- substr(q[first_seq], 1, nchar(q[first_seq]) - 1)
        qf <- file.path(d, "q.fa"); writeLines(q, qf)
        res <- suppressWarnings(suppressMessages(dNdS(qf, lyra_cds(), dnds_est.method = "YN", engine = "coal", delete_corrupt_cds = TRUE,
                                                      print_citation = FALSE)))
        expect_true(all(c("failed_pairs", "dropped_sequences", "coal") %in% names(attributes(res))))
        dropped <- attr(res, "dropped_sequences")
        expect_true(is.data.frame(dropped))
        if (perl_available) {
                ref <- suppressWarnings(suppressMessages(dNdS(qf, lyra_cds(), dnds_est.method = "Comeron", delete_corrupt_cds = TRUE,
                                                              print_citation = FALSE)))
                expect_identical(dropped, attr(ref, "dropped_sequences"))
                expect_identical(frame(res)[ , 1:2], frame(ref)[ , 1:2])
                hit_cols <- setdiff(names(ref), c("dN", "dS", "dNdS", "method"))
                expect_identical(frame(res, hit_cols), frame(ref, hit_cols))
        }
        expect_false(any(duplicated(paste(res$query_id, res$subject_id))))
})

test_that("R15: the same gene pair gets the same values from every COAL entry point", {
        skip_on_cran()
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        methods <- c("YN", if (native_in_coal) c("Comeron", "Li"))
        og <- suppressMessages(orthogroups(c(thal = thal_cds(), lyra = lyra_cds())))
        for (m in methods) {
                a <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = m, engine = "coal", print_citation = FALSE))
                b <- suppressMessages(dNdS_coal(thal_cds(), lyra_cds(), dnds_est.method = m, print_citation = FALSE))
                expect_identical(frame(a), frame(b), info = m)
                g <- suppressMessages(dNdS_orthogroups(og, engine = "coal", aa_aln_tool = "pairwise", dnds_est.method = m))$pairs
                k <- match(paste(a$query_id, a$subject_id), paste(g$query_id, g$subject_id))
                expect_false(anyNA(k), info = m)
                expect_equal(g$dN[k], a$dN, tolerance = 1e-12, info = m)
                expect_equal(g$dS[k], a$dS, tolerance = 1e-12, info = m)
                expect_equal(g$dNdS[k], a$dNdS, tolerance = 1e-12, info = m)
                # the pairs as clusters of two
                cl <- data.frame(representative_id = rep(a$query_id, 2), member_id = c(a$query_id, a$subject_id))
                cc <- suppressMessages(dNdS_clusters(c(thal_cds(), lyra_cds()), clusters = cl, aa_aln_tool = "pairwise", dnds_est.method = m))$pairs
                key <- function(x, y) ifelse(x < y, paste(x, y), paste(y, x))
                k <- match(key(a$query_id, a$subject_id), key(cc$query_id, cc$subject_id))
                expect_false(anyNA(k), info = m)
                expect_equal(cc$dN[k], a$dN, tolerance = 1e-12, info = m)
                expect_equal(cc$dS[k], a$dS, tolerance = 1e-12, info = m)
        }
        # "needle" is another scoring: some pairs get another alignment (the scoring both functions used before)
        a <- suppressMessages(dNdS_coal(thal_cds(), lyra_cds(), dnds_est.method = "YN", print_citation = FALSE))
        n <- suppressMessages(dNdS_coal(thal_cds(), lyra_cds(), dnds_est.method = "YN", aa_aln_scoring = "needle", print_citation = FALSE))
        g <- suppressMessages(dNdS_orthogroups(og, engine = "coal", aa_aln_tool = "pairwise", dnds_est.method = "YN", aa_aln_scoring = "needle"))$pairs
        k <- match(paste(n$query_id, n$subject_id), paste(g$query_id, g$subject_id))
        expect_equal(g$dS[k], n$dS, tolerance = 1e-12)
        expect_identical(eval(formals(dNdS)$aa_aln_scoring), "pwalign")
        expect_identical(eval(formals(dNdS_coal)$aa_aln_scoring), "pwalign")
        expect_identical(eval(formals(dNdS_clusters)$aa_aln_scoring), "pwalign")
        expect_identical(eval(formals(dNdS_orthogroups)$aa_aln_scoring), "pwalign")
})

test_that("R16: both summaries take their statistics over all pairs with a value and say how many", {
        pairs <- data.frame(cluster_id = c("c1", "c1", "c1", "c2"), method = "YN",
                            dN = c(0.1, 0.3, 0.5, NA), dS = c(0.2, NA, 0.6, NA), dNdS = c(0.5, NA, 0.8333, NA),
                            status = c("ok", "saturated", "ok", "no_data"), stringsAsFactors = FALSE)
        st <- as.data.frame(dnds_value_stats(pairs, by = c("cluster_id", "method"), stats = c("median", "mean", "min", "max")))
        expect_equal(st$n_dN, c(3, 0)); expect_equal(st$n_dS, c(2, 0)); expect_equal(st$n_dNdS, c(2, 0))
        expect_equal(st$dN_median, c(0.3, NA))          # the saturated pair counts with its dN
        expect_equal(st$dS_mean, c(0.4, NA))
        expect_equal(st$dN_max, c(0.5, NA))
        # coalr's summary (statistics over status "ok" only) is replaced, its other columns kept
        summary <- data.frame(cluster_id = c("c2", "c1"), n_members = c(2, 3), method = "YN", n_pairs = c(1, 3), n_ok = c(0, 2),
                              dN_median = c(NA, 0.3), dN_mean = c(NA, 0.3), dS_median = c(NA, 0.4), dNdS_max = c(NA, 0.8333))
        new <- dnds_cluster_summary_stats(summary, pairs)
        expect_identical(names(new), c("cluster_id", "n_members", "method", "n_pairs", "n_ok", "n_dN", "n_dS", "n_dNdS",
                                       "dN_median", "dN_mean", "dS_median", "dNdS_max"))
        expect_equal(new$n_dN, c(0, 3))
        expect_equal(new$dN_mean, c(NA, 0.3))
        expect_equal(new$dS_median, c(NA, 0.4))
        # dNdS_orthogroups(): the same rule
        m <- data.table::data.table(orthogroup = c("OG1", "OG1"), gene = c("a", "b"), species = c("x", "y"))
        p <- data.table::data.table(orthogroup = "OG1", method = "Comeron", dN = c(0.1, 0.3), dS = c(0.2, NA), dNdS = c(0.5, NA))
        s <- dnds_groups_summary(p, m)
        expect_equal(c(s$n_dN, s$n_dS, s$n_dNdS), c(2, 1, 1))
        expect_equal(s$dN_median, 0.2)
        expect_equal(s$dS_median, 0.2)
        
        skip_on_cran()
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        cl <- suppressMessages(deepclust(c(thal_cds(), lyra_cds()), seq_type = "cds"))
        res <- suppressMessages(dNdS_clusters(c(thal_cds(), lyra_cds()), clusters = cl, aa_aln_tool = "pairwise"))
        expect_true(all(c("n_dN", "n_dS", "n_dNdS") %in% names(res$summary)))
        again <- as.data.frame(dnds_value_stats(as.data.frame(res$pairs), by = c("cluster_id", "method"), stats = c("median", "mean")))
        k <- match(again$cluster_id, res$summary$cluster_id)
        expect_equal(res$summary$dN_median[k], again$dN_median)
        expect_equal(res$summary$dS_mean[k], again$dS_mean)
        expect_equal(res$summary$n_dNdS[k], again$n_dNdS)
})

test_that("R14: KaKs_Calculator and COAL agree on a gene pair without nonsynonymous change", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        skip_if_not(is_installed_coal(), "COAL is not installed")
        a   <- "ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACCGGTCTGGAAGCGAAACTGGGTATTCCGGAAGCG"
        syn <- "ATGAAGACCGCCTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGTTGACCGGTCTGGAAGCGAAGCTGGGTATTCCGGAAGCG"
        f <- file.path(new_test_dir("r14c_"), "syn.fasta"); writeLines(c(">q", a, ">s", syn), f)
        for (m in c("NG", "YN", "LWL", "MLWL", "LPB", "MLPB")) {
                k <- substitutionrate(f, est.method = m, quiet = TRUE)
                cc <- substitutionrate(f, est.method = m, quiet = TRUE, kaks_tool = "coal")
                expect_equal(k$dN, 0, info = m)
                expect_equal(cc$dN, 0, info = m)
                expect_equal(k$dNdS, cc$dNdS, info = m)
                expect_equal(k$dS, cc$dS, tolerance = 1e-5, info = m)
                expect_equal(nrow(suppressMessages(filter_dNdS(as.data.frame(k)))), nrow(suppressMessages(filter_dNdS(as.data.frame(cc)))), info = m)
        }
})
