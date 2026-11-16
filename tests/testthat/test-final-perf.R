# orthologr 0.6.0, last round: the result tables of a chunk of gene pairs built at once (KaKs_Calculator and
# Comeron routes), gestimator with the codon alignments as strings, and the fixes GMYN without nonsynonymous change,
# GMYN chunk with a pair without rows, a pair named "Sequence", stop codons written with ambiguity codes, and time
# limits of external programs.

final_kaks_header <- paste("Sequence", "Method", "Ka", "Ks", "Ka/Ks", "P-Value(Fisher)", "Length", "S-Sites", "N-Sites",
                           "Fold-Sites(0:2:4)", "Substitutions", "S-Substitutions", "N-Substitutions",
                           "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)", "Divergence-Time",
                           "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
                           "Akaike-Weight", "Model", sep = "\t")
# KaKs_Calculator 2.0 -m GMYN for ATGGCTAAA / ATGGCCAAA (one synonymous difference) and ATGGCTAAAGATGAA /
# ATGGATAAAGATGAA (one nonsynonymous difference)
final_gmyn_rows <- c(
        "short\tGMYN\tNA\t0.873157\t1.2715e-15\tNA\t9\t1.93846\t7.06154\tNA\t1\t1\tNA\tNA\tNA\t0.188065\t2:2:1:1:1:1\t0.388889(0.333333:0.333333:0.5)\tNA\tNA\tNA\tNA",
        "nosyn\tGMYN\t0.0826775\tNA\tNA\tNA\t15\t1.83744\t13.1626\tNA\t1\tNA\t1\tNA\tNA\t0.0725499\t2:2:1:1:1:1\t0.3(0.6:0.1:0.2)\tNA\tNA\tNA\tNA")

# A stand-in for a program: a shell script 'name' in a new folder; returns the folder.
final_fake_program <- function(name, body) {
        d <- new_test_dir("fake_tool_")
        writeLines(c("#!/bin/sh", body), file.path(d, name))
        Sys.chmod(file.path(d, name), "755")
        d
}

# A FASTA codon alignment of two sequences in a file
final_aln_file <- function(a, b, names = c("query", "subject")) {
        f <- tempfile("aln_", fileext = ".fasta")
        writeLines(c(paste0(">", names[1]), a, paste0(">", names[2]), b), f)
        f
}

test_that("GMYN without nonsynonymous difference: dN and dN/dS are exactly 0, dS = 0 without synonymous difference", {
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy_kaks_zero_dn = NULL, orthologr.legacy = NULL)
        tbl <- utils::read.csv(text = c(final_kaks_header, final_gmyn_rows), sep = "\t", header = TRUE)
        res <- kaks_result_table(tbl[1, ], c("q", "s"))
        expect_identical(res$dN, 0)
        expect_identical(res$dNdS, 0)
        expect_identical(res$dS, 0.873157)
        res2 <- kaks_result_table(tbl[2, ], c("q", "s"))
        expect_identical(res2$dS, 0)
        expect_true(is.na(res2$dNdS))
        expect_identical(res2$dN, 0.0826775)
        expect_identical(dnds_status(res$dN, res$dS, res$dNdS), "ok")
        expect_identical(dnds_status(res2$dN, res2$dS, res2$dNdS), "dS_zero")
        # other values are left as they are: a GMYN row with nonsynonymous differences; the ML methods (GY, MS) print
        # expected numbers of differences and are not touched
        x <- tbl[1, ]
        x$N.Substitutions <- 2
        expect_identical(kaks_result_table(x, c("q", "s"))$dN, 1.2715e-15 * 0.873157)
        y <- tbl[1, ]
        y$Method <- "GY-HKY"
        expect_identical(kaks_result_table(y, c("q", "s"))$dN, 1.2715e-15 * 0.873157)
        # options(orthologr.legacy_kaks_zero_dn = TRUE): as printed
        rlang::local_options(orthologr.legacy_kaks_zero_dn = TRUE)
        expect_true(is.na(kaks_result_table(tbl[1, ], c("q", "s"))$dN))
})

test_that("the values of a table of several gene pairs are those of each pair's table (kaks_result_values())", {
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy_kaks_zero_dn = NULL, orthologr.legacy = NULL)
        lines <- readLines(test_path("fixtures", "codon_aln", "lwl_exact_zero_scer_cgla.fa"))
        cgla <- lines[!startsWith(lines, ">")]
        rows <- c(final_gmyn_rows,
                  "q8\tYN\t0.081132\t99\t0.000819515\t0\t384\t70.7133\t313.287\tNA\t109\t84.9412\t24.0588\tNA\tNA\t18.297\t0.662209:0.662209:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
                  "q39\tYN\tNA\tNA\tNA\tNA\t171\t40.0171\t130.983\tNA\t0\tNA\tNA\tNA\tNA\tNA\t2:2:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
                  "q7\tMYN\tnan\t0.3\tnan\tNA\t171\t40.0171\t130.983\tNA\t9\t3\t6\tNA\tNA\tNA\t2:2:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
                  "q13\tLWL\t0.143407\t0.741265\t0.193462\t8.94205e-93\t1452\t324.333\t1127.67\t925:304:223\t415\t240.417\t161.715\t40.5833:106.25:64.6667\t68.8333:37.1667:97.5\t0.27695\t2.07862:2.07862:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA")
        tbl <- utils::read.csv(text = c(final_kaks_header, rows), sep = "\t", header = TRUE)
        aln <- rep(list(NULL, cgla), length.out = nrow(tbl))
        ctx <- paste0("gene pair ", seq_len(nrow(tbl)), " (method x)")
        v <- kaks_result_values(tbl, context = ctx, aln = aln)
        for (i in seq_len(nrow(tbl))) {
                one <- kaks_result_table(tbl[i, ], c("q", "s"), context = ctx[i], aln = aln[[i]])
                expect_identical(v$dN[i], one$dN, info = i)
                expect_identical(v$dS[i], one$dS, info = i)
                expect_identical(v$dNdS[i], one$dNdS, info = i)
                f <- attr(one, "kaks_failed")
                expect_identical(v$failed[i], !is.null(f), info = i)
                if (!is.null(f))
                        expect_identical(v$reason[i], f$reason, info = i)
        }
        expect_identical(kaks_failed_reason(c("A", "B", "A"), c(FALSE, TRUE, FALSE), c(NA, "r", NA)), c(NA, "r", NA))
        expect_identical(kaks_failed_reason("A", FALSE, NA_character_), NA_character_)
})

test_that("a gene pair named 'Sequence' with header-less KaKs_Calculator output is read (not taken for the header)", {
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy_kaks_zero_dn = NULL, orthologr.legacy = NULL)
        row <- sub("^short", "Sequence", final_gmyn_rows[1])
        expect_identical(kaks_output_lines(row), c(final_kaks_header, row))
        expect_identical(kaks_output_lines(c(final_kaks_header, row)), c(final_kaks_header, row))
        expect_identical(kaks_output_lines(c(row, final_kaks_header)), c(final_kaks_header, row))
        skip_on_os("windows")
        # a stand-in for KaKs_Calculator that writes the row without the header line ('-o' is the 4th argument)
        out <- tempfile("gmyn_row_")
        writeLines(row, out)
        bin <- final_fake_program("KaKs_Calculator", paste0('cp ', shQuote(out), ' "$4"'))
        f <- final_aln_file("ATGGCTAAA", "ATGGCCAAA")
        res <- suppressMessages(substitutionrate(f, est.method = "GMYN", kaks_calc_path = bin, subst_name = "Sequence", quiet = TRUE))
        expect_identical(nrow(res), 1L)
        expect_identical(res$dN, 0)
        expect_identical(res$dS, 0.873157)
        expect_identical(res$method, "GMYN")
})

test_that("kaks_calc_batch(): a record without rows (GMYN of one codon) is left to its own call, the other pairs are read", {
        skip_on_cran()
        skip_on_os("windows")
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        rlang::local_options(orthologr.legacy = NULL)
        files <- c(final_aln_file("ATGGCTAAAGATGAA", "ATGGATAAAGCTGAA"), final_aln_file("ATG", "ATG"),
                   final_aln_file("ATGGCTAAAGATGAACCC", "ATGGCCAAAGATGAGCCA"))
        old <- suppressMessages(lapply(files, function(f) tryCatch(substitutionrate(f, est.method = "GMYN", quiet = TRUE),
                                                                    error = function(e) e)))
        expect_s3_class(old[[2]], "error")
        res <- kaks_calc_batch(files, subst_name = c("q1", "q2", "q3"), est.method = "GMYN")
        expect_null(res[[2]])
        for (k in c(1, 3)) {
                expect_s3_class(res[[k]], "dnds_values")
                expect_identical(res[[k]]$dN, old[[k]]$dN, info = k)
                expect_identical(res[[k]]$dS, old[[k]]$dS, info = k)
                expect_identical(res[[k]]$dNdS, old[[k]]$dNdS, info = k)
                expect_identical(res[[k]]$method, old[[k]]$method, info = k)
        }
        # YN gives a row for every record (the single codon: 'nan'): every pair is read from the batch
        expect_false(any(vapply(kaks_calc_batch(files, subst_name = c("q3", "q2", "q1"), est.method = "YN"), is.null, TRUE)))
})

test_that("stop codons written with an ambiguity code (TAR, TRA) are counted as in-frame stop codons", {
        expect_identical(ambiguous_stop_codons(c("ATGTARTRAAAA", "ATGtarAAATAR", "ATGTAA", "ATGTRG", "TAR", NA, "ATGNNNTAR")),
                         c(2L, 1L, 0L, 0L, 0L, 0L, 0L))
        end <- dropped_start(quiet = TRUE)
        on.exit(end(), add = TRUE)
        cds <- c("ATGTARAAATAA", "ATGAAATAA", "ATGTAAATGTRATAG")
        stop_codon_record(c("a", "b", "c"), translate_cds(cds), "f.fa", cds = cds)
        expect_identical(stop_codon_table(), data.frame(file = "f.fa", id = c("a", "c"), stop_codons = c(1L, 2L),
                                                        stringsAsFactors = FALSE))
})

test_that("an external program that does not finish in time is stopped (options(orthologr.tool_timeout))", {
        # the limits
        rlang::local_options(orthologr.tool_timeout = NULL)
        expect_identical(tool_timeout("KaKs_Calculator"), 3600)
        expect_identical(tool_timeout("PAL2NAL", "perl"), 3600)
        expect_identical(tool_timeout("diamond blastp", "diamond"), 604800)
        expect_identical(tool_timeout("BLAST", "blastp"), 604800)
        expect_identical(tool_timeout("makeblastdb"), 604800)
        expect_identical(tool_timeout("COAL (coal-cpp)", "/x/coal-cpp"), 604800)
        expect_identical(tool_timeout("diamond deepclust", "diamond"), Inf)
        expect_identical(tool_timeout("ClustalW", "clustalw2"), Inf)
        rlang::local_options(orthologr.tool_timeout = 5)
        expect_identical(tool_timeout("ClustalW", "clustalw2"), 5)
        expect_identical(tool_timeout("DIAMOND", "diamond"), 5)
        rlang::local_options(orthologr.tool_timeout = c(kaks_calculator = 7, ClustalW = 9, DIAMOND = Inf))
        expect_identical(tool_timeout("KaKs_Calculator"), 7)
        expect_identical(tool_timeout("ClustalW", "clustalw2"), 9)
        expect_identical(tool_timeout("diamond makedb", "diamond"), Inf)
        expect_identical(tool_timeout("PAL2NAL"), 3600)
        rlang::local_options(orthologr.tool_timeout = 0)
        expect_identical(tool_timeout("KaKs_Calculator"), Inf)
        for (bad in list("10", -1, c(1, 2)))
                local({
                        rlang::local_options(orthologr.tool_timeout = bad)
                        expect_error(tool_timeout("KaKs_Calculator"), class = "orthologr_error_bad_argument")
                })

        skip_on_cran()
        skip_on_os("windows")
        # stand-ins that hang: every program run through run_tool() is stopped at its limit
        rlang::local_options(orthologr.tool_timeout = 1)
        for (p in list(c("KaKs_Calculator", "KaKs_Calculator"), c("perl", "PAL2NAL"), c("diamond", "DIAMOND"),
                       c("blastp", "BLAST"), c("coal-cpp", "COAL (coal-cpp)"))) {
                bin <- final_fake_program(p[1], "exec sleep 30")
                t0 <- proc.time()[["elapsed"]]
                e <- tryCatch(run_tool(p[1], args = list("x"), path = bin, tool = p[2], context = "a test"), error = function(e) e)
                expect_lt(proc.time()[["elapsed"]] - t0, 15)
                expect_s3_class(e, "orthologr_error_tool_timeout")
                expect_s3_class(e, "orthologr_error_tool_failed")
                expect_match(conditionMessage(e), paste0(p[2], " did not finish within 1 seconds for a test"), fixed = TRUE)
        }
        # a program that ends in time with status 124 is not a timeout
        bin <- final_fake_program("KaKs_Calculator", "exit 124")
        e <- tryCatch(run_tool("KaKs_Calculator", path = bin), error = function(e) e)
        expect_false(inherits(e, "orthologr_error_tool_timeout"))
        expect_s3_class(e, "orthologr_error_tool_failed")
        # substitutionrate() with a KaKs_Calculator that hangs
        bin <- final_fake_program("KaKs_Calculator", "exec sleep 30")
        f <- final_aln_file("ATGGCTAAAGATGAA", "ATGGATAAAGCTGAA")
        expect_error(substitutionrate(f, est.method = "YN", kaks_calc_path = bin, quiet = TRUE),
                     class = "orthologr_error_tool_timeout")
})

test_that("dNdS(): gene pairs whose KaKs_Calculator hangs fail with the time limit as reason, the run ends", {
        skip_on_cran()
        skip_on_os("windows")
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        rlang::local_options(orthologr.tool_timeout = c(KaKs_Calculator = 1))
        q <- c("ATGGCTAAAGATGAACCCGGG", "ATGAAACCCGGGTTTAAA")
        s <- c("ATGGCCAAAGATGAGCCCGGA", "ATGAAGCCCGGGTTCAAA")
        tbl <- data.table::data.table(query_id = c("a", "b"), subject_id = c("x", "y"), query_cds = q, subject_cds = s,
                                      query_aa = translate_cds(q), subject_aa = translate_cds(s))
        bin <- final_fake_program("KaKs_Calculator", "exec sleep 30")
        t0 <- proc.time()[["elapsed"]]
        res <- suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                             codon_aln_tool = "coal", dnds_est.method = "YN", kaks_calc_path = bin,
                                             quiet = TRUE))
        # one chunk (1 s), then each pair on its own (1 s each)
        expect_lt(proc.time()[["elapsed"]] - t0, 30)
        fp <- attr(res, "failed_pairs")
        expect_identical(fp$query_id, c("a", "b"))
        expect_match(fp$reason, "KaKs_Calculator did not finish within 1 seconds for gene pair 'a' / 'x'", fixed = TRUE, all = FALSE)
        expect_true(all(is.na(res$dN)))
})

test_that("gestimator_pairs() gives the numbers gestimator() writes for the same alignment", {
        set.seed(7)
        cod <- c("TTT", "TTC", "TTA", "TTG", "CTT", "ATG", "GCT", "GCC", "GCA", "AAA", "AAG", "GAT", "GAC", "TGG",
                 "CGA", "AGA", "TAA", "---", "NNN", "ACN")
        a <- b <- character(0)
        for (i in 1:120) {
                n <- sample(c(1:3, 10, 50, 200), 1)
                x <- sample(cod, n, TRUE, prob = c(rep(5, 16), 0.3, 1, 0.3, 0.3))
                y <- x
                ch <- sample(n, max(1, round(n * runif(1, 0, 0.9))))
                y[ch] <- sample(cod, length(ch), TRUE, prob = c(rep(5, 16), 0.3, 1, 0.3, 0.3))
                a <- c(a, paste(x, collapse = ""))
                b <- c(b, paste(y, collapse = ""))
        }
        for (cfg in list(list(FALSE, FALSE, 3), list(TRUE, TRUE, 0), list(FALSE, TRUE, 3))) {
                got <- gestimator_pairs(a, b, legacy_saturation = cfg[[1]], twosubs_libsequence = cfg[[2]], min_sites = cfg[[3]])
                for (i in seq_along(a)) {
                        f <- final_aln_file(a[i], b[i])
                        out <- tempfile()
                        gestimator(f, file_out = out, legacy_saturation = cfg[[1]], twosubs_libsequence = cfg[[2]],
                                   min_sites = cfg[[3]])
                        expect_identical(got[i], sub("^query\tsubject\t", "", readLines(out)), info = i)
                        unlink(c(f, out))
                }
        }
        expect_identical(gestimator_pairs(c("ATGA", "AT", "ATG"), c("ATGA", "AT", "ATGAAA")), rep(NA_character_, 3))
        expect_error(gestimator_pairs("ATG", "ATG", maxHits = 4), "maxHits")
})

test_that("compute_dnds(), Comeron: the gene pairs computed with strings give the table of the route through files", {
        skip_on_cran()
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        q <- suppressMessages(read.cds(thal_cds(), format = "fasta", delete_corrupt_cds = FALSE))
        s <- suppressMessages(read.cds(lyra_cds(), format = "fasta", delete_corrupt_cds = FALSE))
        tbl <- data.table::data.table(query_id = q$geneids, subject_id = s$geneids, query_cds = q$seqs, subject_cds = s$seqs)
        # a pair of identical sequences, a pair with N, a protein that does not fit its CDS, a single codon
        tbl$subject_cds[2] <- tbl$query_cds[2]
        substr(tbl$query_cds[3], 40, 40) <- "N"
        tbl$query_cds[5] <- "ATG"
        tbl$subject_cds[5] <- "ATG"
        tbl$query_aa <- translate_cds(tbl$query_cds)
        tbl$subject_aa <- translate_cds(tbl$subject_cds)
        tbl$query_aa[4] <- tbl$query_aa[6]
        run <- function(keep, cores = 1) {
                old <- options(orthologr.keep_tmp = keep)
                on.exit(options(old), add = TRUE)
                n <- 0L
                f <- check_codon_alignment_file
                testthat::local_mocked_bindings(check_codon_alignment_file = function(...) { n <<- n + 1L; f(...) },
                                                .package = "orthologr")
                r <- suppressMessages(suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise",
                                                                    aa_aln_tool = "NW", codon_aln_tool = "coal",
                                                                    dnds_est.method = "Comeron", quiet = TRUE,
                                                                    comp_cores = cores)))
                attr(r, "files_read") <- n
                r
        }
        strings <- run(FALSE)
        # (keep_tmp keeps the files of every pair: the route through the files)
        files <- run(TRUE)
        # (the reasons name the temporary files of a call)
        clean <- function(r) { attr(r, "files_read") <- NULL; f <- attr(r, "failed_pairs")
                f$reason <- gsub("'[^']*(_alignment|_calculation)[^']*'", "<file>", f$reason); attr(r, "failed_pairs") <- f; r }
        expect_identical(clean(strings), clean(files))
        expect_identical(attr(strings, "files_read"), 0L)
        expect_gt(attr(files, "files_read"), 10L)
        skip_if(parallel::detectCores() < 2 || .Platform$OS.type == "windows", "fewer than 2 cores")
        expect_identical(clean(run(FALSE, cores = 2)), clean(files))
})

test_that("dnds_bind(): one table for consecutive records, in order, with the columns of the tables per pair", {
        rec <- function(q, n = 1, method = NULL, reason = rep(NA_character_, n)) {
                r <- dnds_values(dN = seq_len(n) / 10, dS = rep(0.5, n), dNdS = seq_len(n) / 5,
                                 method = if (!is.null(method)) rep(method, n), reason = reason)
                r$query_id <- q
                r$subject_id <- paste0(q, "_s")
                r
        }
        tab <- function(q, method = NULL) {
                t <- data.table::data.table(query_id = q, subject_id = paste0(q, "_s"), dN = NA_real_, dS = NA_real_, dNdS = NA_real_)
                if (!is.null(method)) t$method <- method
                t$.reason <- "failed"
                t
        }
        per_pair <- function(r) if (inherits(r, "dnds_values")) {
                t <- data.table::data.table(query_id = rep(r$query_id, length(r$dN)), subject_id = rep(r$subject_id, length(r$dN)),
                                            dN = r$dN, dS = r$dS, dNdS = r$dNdS)
                if (!is.null(r$method)) t$method <- r$method
                t$.reason <- r$reason
                t
        } else r
        for (m in list(NULL, "YN")) {
                res <- list(rec("a", method = m), rec("b", 2, m, c(NA, "x")), tab("c", m), rec("d", method = m), tab("e", m), tab("f", m))
                b <- dnds_bind(res)
                expect_identical(length(b), 5L)
                expect_identical(data.table::rbindlist(b, use.names = TRUE, fill = TRUE),
                                 data.table::rbindlist(lapply(res, per_pair), use.names = TRUE, fill = TRUE))
        }
        expect_identical(dnds_bind(list()), list())
        only <- list(tab("a"), tab("b"))
        expect_identical(dnds_bind(only), only)
})
