# What is read from KaKs_Calculator's table (default engine): values it prints as NA although they are 0 (O03), and
# capped values it prints for saturated estimates (O02). The rows are real output of KaKs_Calculator 2.0 (S.
# cerevisiae / C. glabrata and E. coli K-12 / O157 gene pairs); the COAL engine gives the same status for each.

kaks_lines <- c(
        paste("Sequence", "Method", "Ka", "Ks", "Ka/Ks", "P-Value(Fisher)", "Length", "S-Sites", "N-Sites",
              "Fold-Sites(0:2:4)", "Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)",
              "Fold-N-Substitutions(0:2:4)", "Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)",
              "GC(1:2:3)", "ML-Score", "AICc", "Akaike-Weight", "Model", sep = "\t"),
        # 1: YN, Ks capped at 99 (P + Q > 1)
        "q8\tYN\t0.081132\t99\t0.000819515\t0\t384\t70.7133\t313.287\tNA\t109\t84.9412\t24.0588\tNA\tNA\t18.297\t0.662209:0.662209:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
        # 2: YN, Ks at JC69's cap 0.75 log(S-Sites)
        "q2\tYN\t0.281572\t4.19758\t0.0670796\t3.2505e-74\t1233\t269.555\t963.445\tNA\t451\t225.146\t225.854\tNA\tNA\t1.13768\t1.47965:1.47965:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
        # 3: YN, an estimate (dS 2.49)
        "q1\tYN\t0.0829333\t2.49034\t0.0333021\t3.57675e-25\t255\t72.926\t182.074\tNA\t67\t52.718\t14.282\tNA\tNA\t0.771412\t2.02383:2.02383:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
        # 4: LWL, 4-fold class saturated (1 - 2P - Q < 0): Ks underestimated by KaKs_Calculator
        "q13\tLWL\t0.143407\t0.741265\t0.193462\t8.94205e-93\t1452\t324.333\t1127.67\t925:304:223\t415\t240.417\t161.715\t40.5833:106.25:64.6667\t68.8333:37.1667:97.5\t0.27695\t2.07862:2.07862:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
        # 5: YN, identical sequences
        "q39\tYN\tNA\tNA\tNA\tNA\t171\t40.0171\t130.983\tNA\t0\tNA\tNA\tNA\tNA\tNA\t2:2:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
        # 6: YN, two nonsynonymous differences only
        "q115\tYN\t0.00514085\tNA\tNA\tNA\t528\t137.561\t390.439\tNA\t2\tNA\t2\tNA\tNA\t0.00380149\t2.16641:2.16641:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA",
        # 7: LPB, one nonsynonymous difference (a transversion at a 2-fold site)
        "q91\tLPB\t0.00554893\tNA\tNA\tNA\t216\tnan\t185\t146:39:31\t1\tNA\t1.02655\tNA\t0:1:0\t0.00475256\t0:0:1:1:1:1\t0(0:0:0)\tNA\tNA\tNA\tNA"
)
kaks_example_table <- function() utils::read.csv(text = kaks_lines, sep = "\t", header = TRUE)

test_that("capped values of KaKs_Calculator are NA (saturated) and zeros it prints as NA are 0 (dS_zero)", {
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy_kaks_zero_dn = NULL, orthologr.legacy = NULL)
        tbl <- kaks_example_table()
        res <- do.call(rbind, lapply(seq_len(nrow(tbl)), function(i) kaks_result_table(tbl[i, ], c("q", "s"))))
        st <- dnds_status(res$dN, res$dS, res$dNdS)
        expect_identical(st, c("saturated", "saturated", "ok", "saturated", "dS_zero", "dS_zero", "dS_zero"))
        expect_true(all(is.na(res$dS[c(1, 2, 4)])))
        expect_true(all(is.na(res$dNdS[c(1, 2, 4, 5, 6, 7)])))
        expect_equal(res$dN[1:4], c(0.081132, 0.281572, 0.0829333, 0.143407))
        expect_equal(res$dS[3], 2.49034)
        expect_equal(res$dN[5:7], c(0, 0.00514085, 0.00554893))
        expect_equal(res$dS[5:7], c(0, 0, 0))
        # the helpers on the whole table at once
        sat <- kaks_saturated(tbl)
        expect_identical(sat$Ks, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
        expect_identical(sat$Ka, rep(FALSE, 7))
        z <- kaks_zero_ds(tbl)
        expect_identical(z$identical, c(rep(FALSE, 4), TRUE, FALSE, FALSE))
        expect_identical(z$no_synonymous, c(rep(FALSE, 5), TRUE, TRUE))
        # a table without the columns of KaKs_Calculator (only the first five) is read as before
        expect_null(kaks_saturated(tbl[, 1:5]))
        expect_null(kaks_zero_ds(tbl[, 1:5]))
})

# Two LWL-family pairs whose 4-fold class has 1 - 2P - Q = 0 exactly (2 Si + Vi = L in twelfths): KaKs_Calculator's
# floating-point result decides, and its printed counts (6 digits) cannot tell which way it went (V02).
#   q1966 (S. cerevisiae / C. glabrata, LWL): KaKs_Calculator computed the class (Ks 13.4); the printed counts give
#     an argument of -1e-6 (taken for saturated before)
#   q2376 (S. cerevisiae / K. lactis, MLPB): KaKs_Calculator found the class undefined and counted 0 substitutions
#     (Ks 0.967, an underestimate); the printed counts give +7e-7 (taken for an estimate before)
lwl_zero_lines <- c(
        kaks_lines[1],
        "q1966\tLWL\t0.310604\t13.4228\t0.02314\t0\t423\t88.5\t334.5\t268.5:99:55.5\t160\t1187.92\t103.897\t28.6667:34:14.1667\t36.0833:19.9167:27.1667\t3.05394\t1.8477:1.8477:1:1:1:1\t0.399351(0.480519:0.418831:0.298701)\tNA\tNA\tNA\tNA",
        "q2376\tMLPB\t0.670748\t0.967005\t0.693635\t5.15685e-40\t1179\t413.5\t950.788\t765.5:273.5:140\t614\t399.856\t637.739\t110.833:89.75:36.5833\t225.042:84.9583:66.8333\t0.880064\t1.25874:1.25874:1:1:1:1\t0.367492(0.394144:0.34009:0.368243)\tNA\tNA\tNA\tNA"
)
lwl_zero_aln <- function(name) {
        f <- test_path("fixtures", "codon_aln", paste0("lwl_exact_zero_", name, ".fa"))
        x <- readLines(f)
        x[!startsWith(x, ">")]
}

test_that("the saturation test of the LWL family follows KaKs_Calculator's own numbers where the printed counts cannot tell (V02)", {
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy = NULL)
        tbl <- utils::read.csv(text = lwl_zero_lines, sep = "\t", header = TRUE)
        cgla <- lwl_zero_aln("scer_cgla")
        klac <- lwl_zero_aln("scer_klac")
        # the counts computed from the codon alignments are those KaKs_Calculator printed
        k <- kaks_columns(tbl)
        e1 <- kaks_lwl_counts(cgla[1], cgla[2], "LWL")
        e2 <- kaks_lwl_counts(klac[1], klac[2], "MLPB")
        expect_true(kaks_counts_agree(e1, k$L[1, ], k$Si[1, ], k$Vi[1, ]))
        expect_true(kaks_counts_agree(e2, k$L[2, ], k$Si[2, ], k$Vi[2, ]))
        expect_equal(e1$L, c(268.5, 99, 55.5))
        expect_equal(2 * e2$Si[3] + e2$Vi[3], 140)
        # KaKs_Calculator's test on its own numbers: q1966 defined (> 0), q2376 not (<= 0)
        expect_gt(1 - 2 * e1$P[3] - e1$Q[3], 0)
        expect_lte(1 - 2 * e2$P[3] - e2$Q[3], 0)
        expect_identical(kaks_saturated(tbl[1, ], aln = cgla)$Ks, FALSE)
        expect_identical(kaks_saturated(tbl[2, ], aln = klac)$Ks, TRUE)
        expect_identical(kaks_saturated(tbl[2, ], aln = klac)$Ka, FALSE)
        res <- kaks_result_table(tbl[2, ], c("q", "s"), aln = klac)
        expect_true(is.na(res$dS))
        expect_equal(res$dN, 0.670748)
        # without the alignment (or with another one) the printed counts decide, as before
        expect_identical(kaks_saturated(tbl)$Ks, c(TRUE, FALSE))
        expect_identical(kaks_saturated(tbl[1, ], aln = klac)$Ks, TRUE)
        # another genetic code: no table of it here, so the printed counts decide; '-c' is read from kaks_calc.params
        expect_identical(kaks_genetic_code(NULL), 1L)
        expect_identical(kaks_genetic_code("-c 11"), 11L)
        expect_identical(kaks_genetic_code(c("-m", "LWL", "-c", "2")), 2L)
        expect_true(is.na(kaks_genetic_code("-c x")))
        expect_identical(kaks_saturated(tbl[2, ], aln = klac, genetic_code = NA)$Ks, FALSE)
        # the tables of the standard code: the classes of a few codons
        tab <- kaks_lwl_tables(1, modified = FALSE)
        expect_equal(unname(tab$L[c(1, 2, 3, 64), ]), matrix(c(2, 1, 0, 2, 1, 0, 1, 2, 0, 2, 0, 1), ncol = 3, byrow = TRUE))
})

test_that("KaKs_Calculator and orthologr agree on the LWL-family pairs at 1 - 2P - Q = 0 (Linux x86-64)", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        # (the sign of an exact 0 is a matter of KaKs_Calculator's rounding; builds for other processors may round
        # differently, e.g. with fused multiply-add on arm64)
        skip_if_not(Sys.info()[["sysname"]] == "Linux" && R.version$arch == "x86_64", "the expected rounding is that of x86-64")
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy = NULL)
        klac <- suppressMessages(substitutionrate(test_path("fixtures", "codon_aln", "lwl_exact_zero_scer_klac.fa"), est.method = "MLPB", quiet = TRUE))
        expect_true(is.na(klac$dS))
        expect_true(is.finite(klac$dN))
        cgla <- suppressMessages(substitutionrate(test_path("fixtures", "codon_aln", "lwl_exact_zero_scer_cgla.fa"), est.method = "LWL", quiet = TRUE))
        expect_equal(cgla$dS, 13.4228, tolerance = 1e-4)
})

test_that("a KaKs_Calculator estimate printed as 'nan' is a failed estimate, not a saturated one (V07)", {
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy = NULL)
        # MYN rows of S. cerevisiae / C. glabrata (KaKs_Calculator 2.0): Ka and Ks 'nan', and Ks 'nan' with a Ka
        lines <- c(kaks_lines[1],
                   "q2\tMYN\t-nan\t-nan\t-nan\t0\t216\t-nan\t-nan\tNA\t87\t25\t5\tNA\tNA\t-nan\t2.70104:inf:1:1:1:1\t0.484018(0.424658:0.424658:0.60274)\tNA\tNA\tNA\tNA",
                   "q7\tMYN\t0.0546929\t-nan\t-nan\t0.267293\t144\t29.2971\t114.703\tNA\t9\t3\t6\tNA\tNA\t-nan\t2:2:1:1:1:1\t0.177083(0.197917:0.21875:0.114583)\tNA\tNA\tNA\tNA",
                   kaks_lines[4])
        tbl <- utils::read.csv(text = lines, sep = "\t", header = TRUE)
        res <- kaks_result_table(tbl, c("a", "b"), "gene pair 'a' / 'b' (method MYN)")
        expect_true(all(is.na(res$dN[1:2])) && all(is.na(res$dS[1:2])) && all(is.na(res$dNdS[1:2])))
        expect_false(any(is.nan(res$dN)))
        expect_equal(res$dS[3], 2.49034)
        f <- attr(res, "kaks_failed")
        expect_identical(f$method, c("MYN", "MYN"))
        expect_match(f$reason[1], "printed 'nan' for Ka and Ks (method MYN) of gene pair 'a' / 'b': its estimate failed", fixed = TRUE)
        expect_match(f$reason[2], "printed 'nan' for Ks (method MYN)", fixed = TRUE)
        # status: failed (with the table of failed pairs), not saturated
        fp <- data.frame(query_id = "a", subject_id = "b", reason = f$reason[1], stringsAsFactors = FALSE)
        st <- dnds_status(res$dN[1], res$dS[1], res$dNdS[1], failed = dnds_is_failed(data.frame(query_id = "a", subject_id = "b"), fp))
        expect_identical(st, "failed")
        # several methods per pair: only the row of the method that failed
        rows <- data.frame(query_id = "a", subject_id = "b", method = c("YN", "MYN"))
        fpm <- data.frame(query_id = "a", subject_id = "b", reason = "r", method = "MYN", stringsAsFactors = FALSE)
        expect_identical(dnds_is_failed(rows, fpm), c(FALSE, TRUE))
        expect_identical(dnds_failed_reason(rows, fpm), c(NA, "r"))
        fpm$method <- NA_character_
        expect_identical(dnds_is_failed(rows, fpm), c(TRUE, TRUE))
        # options(orthologr.legacy_kaks_saturation = TRUE): as printed (NaN), without the note
        rlang::local_options(orthologr.legacy_kaks_saturation = TRUE)
        old <- kaks_result_table(tbl, c("a", "b"))
        expect_true(is.nan(old$dS[1]) && is.nan(old$dN[1]))
        expect_equal(old$dN[2], 0.0546929)
        expect_null(attr(old, "kaks_failed"))
})

test_that("options(orthologr.legacy_kaks_saturation / legacy_kaks_zero_dn = TRUE) return the values as printed (0.4.3)", {
        rlang::local_options(orthologr.legacy_kaks_saturation = TRUE, orthologr.legacy_kaks_zero_dn = TRUE)
        tbl <- kaks_example_table()
        res <- do.call(rbind, lapply(seq_len(nrow(tbl)), function(i) kaks_result_table(tbl[i, ], c("q", "s"))))
        num <- function(x) suppressWarnings(as.numeric(as.character(x)))
        expect_equal(res$dN, num(tbl$Ka))
        expect_equal(res$dS, num(tbl$Ks))
        expect_equal(res$dNdS, num(tbl$Ka.Ks))
        # the shorthand of the legacy switch sets both
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy_kaks_zero_dn = NULL, orthologr.legacy = TRUE)
        res2 <- do.call(rbind, lapply(seq_len(nrow(tbl)), function(i) kaks_result_table(tbl[i, ], c("q", "s"))))
        expect_equal(res2$dS, res$dS)
        expect_equal(res2$dN, res$dN)
})

test_that("identical sequences and nonsynonymous-only pairs with KaKs_Calculator: dS = 0 as with the COAL engine", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        rlang::local_options(orthologr.legacy_kaks_saturation = NULL, orthologr.legacy_kaks_zero_dn = NULL, orthologr.legacy = NULL)
        a   <- "ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACCGGTCTGGAAGCGAAACTGGGTATTCCGGAAGCG"
        # GCG -> ACG (Ala -> Thr) and AAA -> CAA (Lys -> Gln): nonsynonymous only
        non <- "ATGAAAACCACGTATATTGCGCAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACCGGTCTGGAAGCGAAACTGGGTATTCCGGAAGCG"
        d <- new_test_dir("kaks_zero_")
        g <- file.path(d, "ident.fasta"); writeLines(c(">q", a, ">s", a), g)
        h <- file.path(d, "non.fasta"); writeLines(c(">q", a, ">s", non), h)
        for (m in c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN")) {
                i <- substitutionrate(g, est.method = m, quiet = TRUE)
                expect_equal(c(i$dN, i$dS), c(0, 0), info = m)
                expect_true(is.na(i$dNdS), info = m)
                expect_identical(dnds_status(i$dN, i$dS, i$dNdS), "dS_zero")
                n <- substitutionrate(h, est.method = m, quiet = TRUE)
                expect_true(n$dN > 0 && n$dS == 0 && is.na(n$dNdS), info = m)
                if (requireNamespace("coalr", quietly = TRUE) && is.element("kaks_tool", names(formals(substitutionrate)))) {
                        cn <- substitutionrate(h, est.method = m, quiet = TRUE, kaks_tool = "coal")
                        expect_equal(cn$dS, 0, info = m)
                        expect_equal(n$dN, cn$dN, tolerance = 1e-5, info = m)
                }
        }
        old <- options(orthologr.legacy_kaks_zero_dn = TRUE)
        on.exit(options(old), add = TRUE)
        i <- substitutionrate(g, est.method = "YN", quiet = TRUE)
        expect_true(is.na(i$dN) && is.na(i$dS) && is.na(i$dNdS))
})
