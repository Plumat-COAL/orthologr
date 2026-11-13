# Regression tests for the findings of the third review (R01 - R23); one block per finding.

# Evaluate `call` (quoted) in the global environment with the variables `vars` defined there, as a user at the
# prompt would: this is the situation in which an unevaluated argument reaches a PSOCK worker as a promise.
with_globals <- function(vars, call) {
        for (n in names(vars))
                assign(n, vars[[n]], envir = globalenv())
        on.exit(rm(list = names(vars), envir = globalenv()), add = TRUE)
        eval(call, globalenv())
}

psock_possible <- function() {
        # PSOCK workers load the installed package
        !isTRUE(tryCatch(pkgload::is_dev_package("orthologr"), error = function(e) FALSE)) && parallel::detectCores() >= 2
}

# ---------------------------------------------------------------- R08: PSOCK workers and arguments given as variables

test_that("R08: force_args() evaluates every argument of the calling function", {
        f <- function(a, b = 2, ...) {
                force_args()
                environment()
        }
        n <- 0
        e <- f({n <- n + 1; 5})
        expect_equal(n, 1)
        g <- function(a, b) { force_args(); TRUE }
        expect_error(g(stop("evaluated now")), "evaluated now")
        expect_error(g(1, stop("second too")), "second too")
})

test_that("R08: dNdS() on PSOCK workers with arguments passed as global variables", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        skip_if_not(psock_possible(), "PSOCK workers need the installed package and two cores")
        old <- options(orthologr.cluster_type = "PSOCK")
        on.exit(options(old), add = TRUE)
        ref <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        # quiet = qt: the worker used to find stats::qt
        res <- suppressMessages(with_globals(
                list(orthologr_q0 = thal_cds(), orthologr_s0 = lyra_cds(), qt = TRUE, cf = FALSE, sl = FALSE),
                quote(orthologr::dNdS(orthologr_q0, orthologr_s0, comp_cores = 2, quiet = qt, clean_folders = cf,
                                      store_locally = sl, print_citation = FALSE))))
        expect_equal(nrow(attr(res, "failed_pairs")), 0)
        expect_equal(as.data.frame(res), as.data.frame(ref), ignore_attr = TRUE)
        # a variable local to a calling function
        f <- function() {
                quiet_here <- TRUE
                suppressMessages(dNdS(thal_cds(), lyra_cds(), comp_cores = 2, quiet = quiet_here, print_citation = FALSE))
        }
        expect_equal(as.data.frame(f()), as.data.frame(ref), ignore_attr = TRUE)
})

test_that("R08: kaks_calc_path and aa_aln_path as global variables on PSOCK workers", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        skip_if_not(psock_possible(), "PSOCK workers need the installed package and two cores")
        old <- options(orthologr.cluster_type = "PSOCK")
        on.exit(options(old), add = TRUE)
        if (kaks_available) {
                res <- suppressMessages(with_globals(
                        list(orthologr_q0 = thal_cds(), orthologr_s0 = lyra_cds(), kp = dirname(Sys.which("KaKs_Calculator"))),
                        quote(orthologr::dNdS(orthologr_q0, orthologr_s0, comp_cores = 2, dnds_est.method = "YN",
                                              kaks_calc_path = kp, print_citation = FALSE))))
                expect_equal(nrow(attr(res, "failed_pairs")), 0)
                expect_true(sum(!is.na(res$dNdS)) >= 18)
        }
        if (clustalo_available) {
                res <- suppressMessages(with_globals(
                        list(orthologr_q0 = thal_cds(), orthologr_s0 = lyra_cds(), ap = dirname(Sys.which("clustalo")),
                             apar = NULL, at = "clustalo"),
                        quote(orthologr::dNdS(orthologr_q0, orthologr_s0, comp_cores = 2, aa_aln_type = "multiple",
                                              aa_aln_tool = at, aa_aln_path = ap, aa_aln_params = apar,
                                              print_citation = FALSE))))
                expect_equal(nrow(attr(res, "failed_pairs")), 0)
                expect_true(sum(!is.na(res$dNdS)) >= 18)
        }
        skip_if_not(kaks_available || clustalo_available, "neither KaKs_Calculator nor Clustal Omega is installed")
})

test_that("R08: divergence_stratigraphy() and dNdS_orthogroups() on PSOCK workers with global variables", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        skip_if_not(psock_possible(), "PSOCK workers need the installed package and two cores")
        old <- options(orthologr.cluster_type = "PSOCK")
        on.exit(options(old), add = TRUE)
        ds <- suppressWarnings(suppressMessages(with_globals(
                list(orthologr_q0 = thal_cds(), orthologr_s0 = lyra_cds(), qt = TRUE),
                quote(orthologr::divergence_stratigraphy(orthologr_q0, orthologr_s0, comp_cores = 2, quiet = qt)))))
        tbl <- attr(ds, "dNdS_tbl")
        expect_equal(nrow(attr(tbl, "failed_pairs")), 0)
        expect_true(sum(!is.na(tbl$dNdS)) >= 18)
        
        og <- suppressMessages(orthogroups(c(a = thal_cds(), b = lyra_cds())))
        # pairwise: the pipeline of dNdS() on the workers
        res <- suppressMessages(with_globals(
                list(orthologr_og = og, tool = "pairwise", qt = TRUE, kp = NULL),
                quote(orthologr::dNdS_orthogroups(orthologr_og, aa_aln_tool = tool, comp_cores = 2, quiet = qt,
                                                  kaks_calc_path = kp))))
        expect_true(all(is.na(res$pairs$note)))
        expect_true(sum(!is.na(res$pairs$dNdS)) >= 18)
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        res <- suppressMessages(with_globals(
                list(orthologr_og = og, ap = dirname(Sys.which("clustalo")), apar = NULL),
                quote(orthologr::dNdS_orthogroups(orthologr_og, aa_aln_tool = "clustalo", aa_aln_path = ap,
                                                  aa_aln_params = apar, comp_cores = 2))))
        expect_true(all(is.na(res$pairs$note)))
        expect_true(sum(!is.na(res$pairs$dNdS)) >= 18)
})

test_that("R13: the workers of an unfinished parallel run are ended", {
        skip_on_cran()
        skip_on_os("windows")
        skip_if(parallel::detectCores() < 2, "one core only")
        cl <- dnds_cluster(2)
        pids <- attr(cl, "orthologr_pids")
        expect_length(pids, 2)
        alive <- function(p) vapply(p, function(x) file.exists(file.path("/proc", x)) &&
                                            !grepl(") Z ", readLines(file.path("/proc", x, "stat"), warn = FALSE)[1], fixed = TRUE), logical(1))
        skip_if_not(dir.exists("/proc"), "no /proc")
        expect_true(all(alive(pids)))
        dnds_cluster_stop(cl, finished = FALSE)
        Sys.sleep(0.5)
        expect_false(any(alive(pids)))
})

# ---------------------------------------------------------------- R03: annotation lines of more than 65,000 characters

test_that("R03: retrieve_longest_isoforms() keeps features annotated in lines of more than 65,000 characters", {
        d <- new_test_dir("r03_")
        a <- suppressWarnings(make_test_annotation(thal_aa(), "Athaliana", function(x) sub("\\..*", "", x), d))
        long_note <- strrep("A", 70000)
        for (fmt in c("gff", "gtf")) {
                f <- a[[fmt]]
                x <- readLines(f)
                # the CDS and transcript rows of the long isoform of the first gene get a very long attribute
                hit <- which(grepl(a$ids[1], x, fixed = TRUE) & !grepl("_iso2", x, fixed = TRUE) & !grepl("\tgene\t", x, fixed = TRUE))
                expect_true(length(hit) >= 2)
                x[hit] <- if (fmt == "gff") paste0(x[hit], ";Note=", long_note) else paste0(x[hit], " note \"", long_note, "\";")
                long_file <- file.path(d, paste0("long.", fmt))
                writeLines(x, long_file)
                expect_true(all(nchar(readLines(long_file)[hit]) > 65000))
                out <- file.path(d, paste0("longest_", fmt, ".fa"))
                msgs <- character(0)
                res <- withCallingHandlers(
                        retrieve_longest_isoforms(a$proteome, long_file, out, annotation_format = fmt),
                        message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
                got <- sub("\\s.*", "", names(Biostrings::readAAStringSet(out)))
                expect_setequal(got, a$ids)                    # the long isoform of every gene, also of the first
                expect_false(any(grepl("not in the annotation", msgs)))
                expect_false(any(grepl("outlier", msgs)))
                # the same result as for the file without the long attribute
                out2 <- file.path(d, paste0("ref_", fmt, ".fa"))
                suppressMessages(retrieve_longest_isoforms(a$proteome, f, out2, annotation_format = fmt))
                expect_identical(readLines(out), readLines(out2))
        }
})

test_that("R03: check_annotation() reports over-long lines and removes them only on request, with a warning", {
        d <- new_test_dir("r03b_")
        f <- file.path(d, "a.gff")
        lines <- c("##gff-version 3", "chr1\tx\tgene\t1\t20\t.\t+\t.\tID=g1",
                   paste0("chr1\tx\tCDS\t1\t20\t.\t+\t0\tID=c1;Parent=g1;Note=", strrep("A", 70000)))
        writeLines(lines, f)
        expect_warning(res <- suppressMessages(check_annotation(f)), "orthologr read such lines")
        expect_identical(res, f)
        expect_warning(copy <- suppressMessages(check_annotation(f, remove_annotation_outliers = TRUE)),
                       "features annotated there .* are missing")
        expect_identical(readLines(copy), lines[1:2])
        expect_identical(readLines(f), lines)
        expect_error(check_annotation(f, remove_annotation_outliers = "yes"), "remove_annotation_outliers")
})

test_that("R03: NCBI TAIR10.1 - no protein is lost and AT5G38410 gets its longest isoform (needs the files)", {
        dl <- Sys.getenv("ORTHOLOGR_TEST_TAIR_DIR", "")
        faa <- file.path(dl, "GCF_000001735.4_TAIR10.1_protein.faa.gz")
        skip_if_not(nzchar(dl) && file.exists(faa), "set ORTHOLOGR_TEST_TAIR_DIR to the folder with the NCBI TAIR10.1 files")
        for (fmt in c("gff", "gtf")) {
                out <- tempfile(fileext = ".fa")
                msgs <- character(0)
                withCallingHandlers(
                        retrieve_longest_isoforms(faa, file.path(dl, paste0("GCF_000001735.4_TAIR10.1_genomic.", fmt, ".gz")), out, fmt),
                        message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
                got <- sub("\\s.*", "", names(Biostrings::readAAStringSet(out)))
                expect_true(is.element("NP_001119331.1", got))
                expect_false(is.element("NP_198657.1", got))
                # 27,562 loci from the GFF3 file, 27,561 from the GTF file (one protein is not in the GTF file), as an
                # independent parser finds; 27,560 / 27,551 when the long lines were removed
                expect_length(got, if (fmt == "gff") 27562 else 27561)
                expect_false(any(grepl("NP_001119331.1", msgs, fixed = TRUE)))
        }
})

# ---------------------------------------------------------------- R14: dN = 0 with the KaKs_Calculator methods

test_that("R14: kaks_zero_dn() sets dN = 0 only where KaKs_Calculator printed NA for a known zero", {
        #                 zero dN   identical  failed  dS < 1e-6  normal  -0     tiny Ka           dS NA
        dN   <- c(NA,      NA,       NaN,    0.0147,    0.0143, NA,    NA,               NA)
        dS   <- c(0.2046,  NA,       NaN,    NA,        0.2144, 0.189, 0.5,              NA)
        dNdS <- c(0,       NA,       NaN,    NA,        0.0665, -0,    1e-7,             0)
        r <- kaks_zero_dn(dN, dS, dNdS)
        expect_equal(r$changed, c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))
        expect_equal(r$dN[c(1, 6)], c(0, 0))
        expect_equal(r$dN[7], 5e-8)
        expect_true(all(is.na(r$dN[c(2, 3, 8)])))
        expect_equal(r$dN[4:5], dN[4:5])
        expect_identical(sign(1 / r$dNdS[6]), 1)        # -0 became 0
        # a table as KaKs_Calculator writes it
        tbl <- data.frame(Sequence = c("a", "b", "c"), Method = "YN", Ka = c("NA", "NA", "0.0143"), Ks = c("0.204638", "NA", "0.214"),
                          `Ka/Ks` = c("-0", "NA", "0.0665"), check.names = FALSE)
        res <- kaks_result_table(tbl, c("q", "s"))
        expect_equal(res$dN, c(0, NA, 0.0143))
        expect_equal(res$dNdS, c(0, NA, 0.0665))
        old <- options(orthologr.legacy_kaks_zero_dn = TRUE)
        on.exit(options(old), add = TRUE)
        res <- kaks_result_table(tbl, c("q", "s"))
        expect_equal(res$dN, c(NA, NA, 0.0143))
        expect_equal(res$dNdS, c(0, NA, 0.0665))
})

test_that("R14: a gene pair with synonymous changes only has dN = 0 with KaKs_Calculator and survives filter_dNdS()", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        a   <- "ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACCGGTCTGGAAGCGAAACTGGGTATTCCGGAAGCG"
        syn <- "ATGAAGACCGCCTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGTTGACCGGTCTGGAAGCGAAGCTGGGTATTCCGGAAGCG"
        d <- new_test_dir("r14_")
        f <- file.path(d, "syn.fasta"); writeLines(c(">q", a, ">s", syn), f)
        g <- file.path(d, "ident.fasta"); writeLines(c(">q", a, ">s", a), g)
        for (m in c("NG", "YN", "LWL")) {
                r <- substitutionrate(f, est.method = m, quiet = TRUE)
                expect_equal(r$dN, 0, info = m)
                expect_true(r$dS > 0.1, info = m)
                expect_equal(r$dNdS, 0, info = m)
                kept <- suppressMessages(filter_dNdS(as.data.frame(r), dnds.threshold = 2))
                expect_equal(nrow(kept), 1, info = m)
                # identical sequences: KaKs_Calculator gives no value at all; nothing is invented
                i <- substitutionrate(g, est.method = m, quiet = TRUE)
                expect_true(is.na(i$dN) && is.na(i$dS) && is.na(i$dNdS), info = m)
        }
        old <- options(orthologr.legacy_kaks_zero_dn = TRUE)
        on.exit(options(old), add = TRUE)
        r <- substitutionrate(f, est.method = "YN", quiet = TRUE)
        expect_true(is.na(r$dN))
        expect_equal(r$dNdS, 0)
})

# ---------------------------------------------------------------- R23: the output contract of divergence_stratigraphy()

test_that("R23: divergence_stratigraphy() returns a divergence map also when the binning fails", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        testthat::local_mocked_bindings(divergence_map = function(...) stop("no strata today"), .package = "orthologr")
        for (sid in c(FALSE, TRUE)) {
                expect_warning(
                        ds <- suppressMessages(divergence_stratigraphy(thal_cds(), lyra_cds(), quiet = TRUE, subject.id = sid)),
                        "could not be sorted into divergence strata \\(no strata today\\).*without rows")
                expect_s3_class(ds, "data.frame")
                expect_identical(names(ds), c("DS", "query_id", if (sid) "subject_id"))
                expect_equal(nrow(ds), 0)
                expect_equal(nrow(attr(ds, "dNdS_tbl")), 20)          # the run is not lost
                expect_true(is.data.frame(attr(ds, "left_out")))
        }
})

test_that("R23: the defaults of dS_zero are \"drop\" (0.4.3 behaviour) in filter_dNdS() and divergence_stratigraphy()", {
        old <- options(orthologr.dS_zero = NULL)
        on.exit(options(old), add = TRUE)
        expect_identical(eval(formals(filter_dNdS)$dS_zero), "drop")
        expect_identical(eval(formals(divergence_stratigraphy)$dS_zero), "drop")
        x <- data.frame(query_id = c("a", "b"), subject_id = c("x", "y"), dN = c(0, 0.1), dS = c(0, 0.5), dNdS = c(NA, 0.2))
        res <- suppressMessages(filter_dNdS(x))
        expect_identical(res$query_id, "b")
        expect_identical(names(res), names(x))                        # no extra column
})
