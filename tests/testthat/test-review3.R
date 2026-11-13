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
