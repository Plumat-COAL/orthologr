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
        !loaded_with_load_all() && parallel::detectCores() >= 2
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

# ---------------------------------------------------------------- R07: parameter strings are argument lists, not shell text

test_that("R07: split_params() splits like a shell and gives no character a meaning", {
        expect_identical(split_params(NULL), character(0))
        expect_identical(split_params(""), character(0))
        expect_identical(split_params("  -seg no   -word_size 3 "), c("-seg", "no", "-word_size", "3"))
        expect_identical(split_params("--force; touch INJ"), c("--force;", "touch", "INJ"))
        expect_identical(split_params("-a 'b c' \"d e\" f\\ g"), c("-a", "b c", "d e", "f g"))
        expect_identical(split_params("-a \"x\\\"y\" ''"), c("-a", "x\"y", ""))
        expect_identical(split_params("$(touch X) `touch Y` | cat > f && ls *"),
                         c("$(touch", "X)", "`touch", "Y`", "|", "cat", ">", "f", "&&", "ls", "*"))
        expect_identical(split_params("-a\t-b\n-c"), c("-a", "-b", "-c"))
        expect_identical(split_params(c("-a 1", "-b 2")), c("-a", "1", "-b", "2"))
        expect_error(split_params("-a 'open"), "Unbalanced quote")
        # every token is quoted for the shell
        skip_on_os("windows")
        expect_identical(quote_args(list("-i", I("--id 10; touch X"))), c("'-i'", "'--id'", "'10;'", "'touch'", "'X'"))
        expect_identical(quote_args(list(I(""))), character(0))
        res <- run_tool("echo", args = list(I("a; echo INJECTED $(echo sub) `echo bq` > /dev/null | cat")), error = FALSE)
        expect_identical(res$stdout, "a; echo INJECTED $(echo sub) `echo bq` > /dev/null | cat")
})

test_that("R07: an injection attempt in any parameter string runs no command", {
        skip_on_cran()
        skip_on_os("windows")
        d <- new_test_dir("r07_")
        n <- 0
        attempt <- function(expr_fun) {
                n <<- n + 1
                marker <- file.path(d, paste0("INJ", n))
                # (the injected command makes a folder: a program that takes the path for the name of an output
                # file of its own - T-Coffee writes its log there - is not mistaken for an executed command)
                for (inj in c(paste0("; mkdir ", marker), paste0("$(mkdir ", marker, ")"), paste0("`mkdir ", marker, "`"),
                              paste0("| mkdir ", marker), paste0("&& mkdir ", marker), paste0("\n mkdir ", marker))) {
                        try(suppressWarnings(suppressMessages(utils::capture.output(expr_fun(inj)))), silent = TRUE)
                        expect_false(dir.exists(marker), info = paste(n, inj))
                        unlink(marker)
                }
        }
        aln <- system.file("seqs/aa_seqs.fasta", package = "orthologr")
        if (clustalo_available) attempt(function(p) multi_aln(aln, tool = "clustalo", params = paste("--force", p), quiet = TRUE))
        if (mafft_available) attempt(function(p) multi_aln(aln, tool = "mafft", params = paste("--quiet", p), quiet = TRUE))
        if (muscle_available) attempt(function(p) multi_aln(aln, tool = "muscle", params = p, quiet = TRUE))
        if (clustalw_available) attempt(function(p) multi_aln(aln, tool = "clustalw", params = paste("-quiet", p), quiet = TRUE))
        if (tcoffee_available) attempt(function(p) multi_aln(aln, tool = "t_coffee", params = paste("-quiet", p), quiet = TRUE))
        if (perl_available) attempt(function(p) codon_aln(
                file_aln = system.file("seqs/aa_seqs.aln", package = "orthologr"),
                file_nuc = system.file("seqs/dna_seqs.fasta", package = "orthologr"),
                format = "clustal", codon_aln_name = "r07", params = p, quiet = TRUE))
        if (kaks_available) attempt(function(p) substitutionrate(
                system.file("seqs/pal2nal.aln", package = "orthologr"), est.method = "YN", kaks_calc.params = paste("-c 1", p), quiet = TRUE))
        if (diamond_available) {
                attempt(function(p) diamond(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = paste("--id 10", p)))
                attempt(function(p) diamond_best(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = p))
                attempt(function(p) diamond_rec(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = p))
                attempt(function(p) deepclust(thal_aa(), direct_input = TRUE, deepclust_params = p))
                if (perl_available) {
                        attempt(function(p) dNdS(thal_cds(), lyra_cds(), aligner_params = p, print_citation = FALSE))
                        if (clustalo_available) attempt(function(p) dNdS(
                                thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "clustalo",
                                aa_aln_params = paste("--force", p), print_citation = FALSE))
                }
        }
        if (blast_available) {
                attempt(function(p) blast(thal_aa(), lyra_aa(), seq_type = "protein", blast_params = paste("-seg no", p)))
                attempt(function(p) blast_best(thal_aa(), lyra_aa(), seq_type = "protein", blast_params = p))
                attempt(function(p) blast_rec(thal_aa(), lyra_aa(), seq_type = "protein", blast_params = p))
        }
        expect_length(list.dirs(d, recursive = FALSE), 0)
        expect_true(n > 0 || !perl_available)
})

test_that("R07: ordinary parameter strings still reach the program", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        a <- suppressMessages(diamond(thal_aa(), lyra_aa(), seq_type = "protein"))
        b <- suppressMessages(diamond(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = "--id 99.5  --max-hsps 1"))
        expect_true(nrow(b) < nrow(a))
        expect_true(all(b$perc_identity >= 99.5))
})

# ---------------------------------------------------------------- R09: a missing program is an error before the search

test_that("R09: a missing program stops dNdS() and dNdS_orthogroups() before the search", {
        skip_on_cran()
        skip_on_os("windows")
        skip_if_not(perl_available, "perl is not installed")
        # a PATH with perl only
        bin <- new_test_dir("r09_bin_")
        file.symlink(Sys.which("perl"), file.path(bin, "perl"))
        old <- Sys.getenv("PATH")
        on.exit({ Sys.setenv(PATH = old); tool_cache_reset() }, add = TRUE)
        Sys.setenv(PATH = bin)
        tool_cache_reset()
        searched <- FALSE
        testthat::local_mocked_bindings(
                diamond_rec = function(...) { searched <<- TRUE; stop("the search was started") },
                .package = "orthologr")
        expect_error(suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "YN", print_citation = FALSE)),
                     "KaKs_Calculator.*needed for dnds_est.method = \"YN\"")
        expect_error(suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "clustalo", print_citation = FALSE)),
                     "Clustal.Omega")
        expect_error(suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "mafft", print_citation = FALSE)),
                     "MAFFT")
        expect_error(suppressMessages(divergence_stratigraphy(thal_cds(), lyra_cds(), dnds_est.method = "YN")), "KaKs_Calculator")
        expect_false(searched)
        # no perl: PAL2NAL cannot run
        Sys.setenv(PATH = new_test_dir("r09_empty_"))
        tool_cache_reset()
        expect_error(suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE)), "'perl' was not found")
        expect_false(searched)
        # what is needed and there passes
        Sys.setenv(PATH = old); tool_cache_reset()
        expect_true(dnds_require_programs(aa_aln_tool = "NW", dnds_est.method = "Comeron"))
        expect_true(dnds_require_programs(aa_aln_tool = "pairwise", dnds_est.method = "Li"))
})

# ---------------------------------------------------------------- R10: longest_per_group_do() and a column 'pos'

test_that("R10: longest_per_group_do() is not misled by a column named 'pos'", {
        set.seed(10)
        for (k in 1:20) {
                n <- sample(5:60, 1)
                x <- tibble::tibble(gene_id = sample(paste0("g", 1:8), n, replace = TRUE), id = paste0("t", seq_len(n)),
                                    width = sample(50:60, n, replace = TRUE), pos = sample(n))
                a <- longest_per_group(x, "gene_id")
                b <- longest_per_group_do(x, "gene_id")
                expect_identical(as.data.frame(a), as.data.frame(b))
                expect_equal(nrow(b), length(unique(x$gene_id)))
        }
        # the path on which the reference code is used: missing widths
        x$width[1] <- NA
        res <- longest_per_group(x, "gene_id")
        expect_equal(nrow(res), length(unique(x$gene_id)))
        best <- vapply(split(x, x$gene_id), function(g) g$id[which.max(g$width)[1]], character(1))
        expect_setequal(res$id, best)
})

# ---------------------------------------------------------------- R11: the pkgdown reference index is complete

test_that("R11: every exported function and every help topic is in the reference index of _pkgdown.yml", {
        yml <- test_path("..", "..", "_pkgdown.yml")
        man <- test_path("..", "..", "man")
        skip_if_not(file.exists(yml) && dir.exists(man), "the package sources (_pkgdown.yml, man/) are not available")
        lines <- readLines(yml, warn = FALSE)
        from <- grep("^reference:", lines)
        to <- c(grep("^[a-z]+:", lines), length(lines) + 1L)
        lines <- lines[from:(min(to[to > from]) - 1L)]
        listed <- trimws(sub("^\\s*-\\s+", "", grep("^\\s+-\\s+[A-Za-z._0-9]+\\s*$", lines, value = TRUE)))
        expect_false(any(duplicated(listed)))
        topics <- list()
        for (f in list.files(man, pattern = "\\.Rd$", full.names = TRUE)) {
                rd <- readLines(f, warn = FALSE)
                if (any(grepl("^\\\\keyword\\{internal\\}", rd)))
                        next
                name <- sub("^\\\\name\\{(.*)\\}\\s*$", "\\1", grep("^\\\\name\\{", rd, value = TRUE)[1])
                topics[[name]] <- sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", grep("^\\\\alias\\{", rd, value = TRUE))
        }
        # pkgdown stops with "topics missing from index" unless every topic is listed (by name or alias)
        is_listed <- vapply(names(topics), function(n) any(is.element(c(n, topics[[n]]), listed)), logical(1))
        expect_identical(names(topics)[!is_listed], character(0))
        # every exported function is documented in a listed topic
        exported <- getNamespaceExports("orthologr")
        aliases <- unlist(topics, use.names = FALSE)
        internal_ok <- character(0)
        for (f in list.files(man, pattern = "\\.Rd$", full.names = TRUE)) {
                rd <- readLines(f, warn = FALSE)
                if (any(grepl("^\\\\keyword\\{internal\\}", rd)))
                        internal_ok <- c(internal_ok, sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", grep("^\\\\alias\\{", rd, value = TRUE)))
        }
        expect_identical(sort(setdiff(exported, c(aliases, internal_ok))), character(0))
        # and nothing is listed that does not exist
        expect_identical(setdiff(listed, c(names(topics), aliases)), character(0))
})
