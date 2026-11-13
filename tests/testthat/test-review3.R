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

# ---------------------------------------------------------------- R06: the same sequence ID in query and subject

test_that("R06: dNdS() works when query and subject sequences have the same IDs, on every tool path", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        d <- new_test_dir("r06_")
        q <- Biostrings::readDNAStringSet(thal_cds()); s <- Biostrings::readDNAStringSet(lyra_cds())
        ref <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        ref_ids <- sub("\\s.*", "", names(q))
        for (scheme in c("same", "numeric", "minus")) {
                ids <- switch(scheme, same = sprintf("gene%02d", 1:20), numeric = as.character(1:20), minus = as.character(-(1:20)))
                x <- q; names(x) <- ids; y <- s; names(y) <- ids
                qf <- file.path(d, paste0(scheme, "_q.fa")); sf <- file.path(d, paste0(scheme, "_s.fa"))
                Biostrings::writeXStringSet(x, qf); Biostrings::writeXStringSet(y, sf)
                run <- function(...) suppressMessages(dNdS(qf, sf, print_citation = FALSE, ...))
                res <- run()
                expect_equal(nrow(attr(res, "failed_pairs")), 0, info = scheme)
                expect_identical(res$query_id, res$subject_id, info = scheme)       # gene i pairs with gene i
                # the values of the run with the original IDs
                m <- match(ref$query_id, ref_ids[match(res$query_id, ids)])
                expect_equal(res$dN[m], ref$dN, info = scheme)
                expect_equal(res$dS[m], ref$dS, info = scheme)
                expect_equal(run(dnds_est.method = "Li")$dN[m], suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "Li", print_citation = FALSE))$dN,
                             info = scheme)
                if (kaks_available) {
                        r <- run(dnds_est.method = "YN")
                        expect_equal(nrow(attr(r, "failed_pairs")), 0, info = scheme)
                        expect_true(sum(!is.na(r$dS)) >= 18, info = scheme)
                }
                for (tool in c("clustalo", "mafft", "muscle", "clustalw")) {
                        if (!get(paste0(tool, "_available"))) next
                        r <- run(aa_aln_type = "multiple", aa_aln_tool = tool)
                        expect_equal(nrow(attr(r, "failed_pairs")), 0, info = paste(scheme, tool))
                        expect_true(all(is.element(r$query_id, ids)) && sum(!is.na(r$dN)) >= 18, info = paste(scheme, tool))
                }
        }
})

test_that("R06: the reason of a failed pair names the IDs of the pair, and stored alignments come with an index", {
        skip_if_not(perl_available, "perl is not installed")
        expect_identical(dnds_pair_reason("PAL2NAL failed for gene pair 'query' / 'subject' (q3): x", "A|1", "A|1"),
                         "PAL2NAL failed for gene pair 'A|1' / 'A|1' (q3): x")
        tbl <- data.table::data.table(
                query_id = c("same", "bad"), subject_id = c("same", "bad"),
                query_cds   = rep("ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACC", 2),
                # the coding sequence of the subject of 'bad' lacks its last nucleotide
                subject_cds = c("ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTGAAGGCGTGGTTGACC",
                                "ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTGAAGGCGTGGTTGAC"),
                query_aa    = c("MKTAYIAKQRLEGDLKAWLT", "MKTAYIAKQRLEGDLKAWLT"),
                subject_aa  = c("MKTAYIAKQRLDGDLKAWLT", "MKTAYIAKQRLDGDLKAWLT"))
        expect_warning(
                res <- suppressMessages(compute_dnds(tbl, aa_aln_type = "pairwise", aa_aln_tool = "NW", dnds_est.method = "Comeron", quiet = TRUE)),
                "1 of 2 gene pairs could not be processed")
        expect_identical(res$query_id, c("bad", "same"))
        expect_false(is.na(res$dN[2]))
        fp <- attr(res, "failed_pairs")
        expect_match(fp$reason, "gene pair 'bad' / 'bad'")
        expect_false(grepl("'query' / 'subject'", fp$reason, fixed = TRUE))
        
        wd <- new_test_dir("r06_store_")
        old <- setwd(wd); on.exit(setwd(old), add = TRUE)
        dir.create("orthologr_alignment_files")
        suppressMessages(compute_dnds(tbl[1], aa_aln_type = "pairwise", aa_aln_tool = "NW", dnds_est.method = "Comeron",
                                      quiet = TRUE, store_locally = TRUE))
        idx <- utils::read.delim(file.path("orthologr_alignment_files", "_alignment", "pair_ids.tsv"), colClasses = "character")
        expect_identical(idx, data.frame(file_prefix = "q1", query = "same", subject = "same"))
        expect_identical(fasta_ids(file.path("orthologr_alignment_files", "_alignment", "pairwise_aln", "q1_NW_AA.aln")), c("query", "subject"))
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

# ---------------------------------------------------------------- R17: Comeron saturation and classes without sites

test_that("R17: an undefined logarithm in a class without sites does not make Comeron's estimate NA", {
        gest <- function(a, b, legacy = FALSE) {
                f <- tempfile(fileext = ".fa"); o <- tempfile()
                on.exit(unlink(c(f, o)), add = TRUE)
                writeLines(c(">a", a, ">b", b), f)
                gestimator(f, o, 3L, FALSE, FALSE, legacy)
                v <- as.numeric(scan(o, what = "", sep = "\t", quiet = TRUE)[3:5])
                v[v == 999] <- NA
                v
        }
        # three pairs of the review's simulation that were NA in dN although the class concerned has no sites:
        # the values are those of the uncorrected (legacy) computation
        cases <- list(
                c("GGTGCCCTAAACGTCGAAAAAGAGTGCGGA", "AAAGCCCTAAATGGGGTCCTAGTGTGCTCA", 0.6615, 0.6216),
                c("TACTTCCCGACTACCGCAGCTTACTATTTG", "CAGATGCCGTATTCCGAAGCTTTCCTAATG", 0.8430, 0.2085),
                c("TTGTTGTTGGACGTTGTAGTGTTGCTTTCT", "TTGTTGCTGGACGATATGGTGTTACGTTCC", 0.1814, 0.4013))
        for (x in cases) {
                v <- gest(x[1], x[2])
                expect_equal(v[1:2], as.numeric(x[3:4]), info = x[1])
                expect_equal(v, gest(x[1], x[2], legacy = TRUE), info = x[1])
        }
        # a class WITH sites that is saturated still gives NA (30 codons, 8 differences), and the legacy value on request
        a <- "ACGCCAAAACTCGTTCTGTGGAAAAGACTAACACATGTCAATTGGCAGACCAGCACCTCCTTTCAGTTACAGGCGGCTTTAGTTCAAAGG"
        b <- "ACGCCAAAACTCGTTCTGCGGAAAACACAAACACATGTCAAGTGGCAGACCAGCACCTCCTTTCCGTTACAGGCGGCATATGTTCAAAGG"
        v <- gest(a, b)
        expect_true(is.na(v[1]) && is.na(v[3]))
        expect_equal(v[2], 0.06677)
        expect_equal(gest(a, b, legacy = TRUE), c(0.09579, 0.06677, 1.435))
        # the documented example of saturation at synonymous sites is unchanged
        expect_true(is.na(gest("ATGGCTAAA", "ATGGCCAAA")[2]))
})
