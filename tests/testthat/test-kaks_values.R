# What is read from KaKs_Calculator's table (default engine): values it prints as NA although they are 0 (O03), and
# capped values it prints for saturated estimates (O02). The rows are real output of KaKs_Calculator 2.0 (S.
# cerevisiae / C. glabrata and E. coli K-12 / O157 gene pairs).

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

test_that("identical sequences and nonsynonymous-only pairs with KaKs_Calculator: dS = 0", {
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
        }
        old <- options(orthologr.legacy_kaks_zero_dn = TRUE)
        on.exit(options(old), add = TRUE)
        i <- substitutionrate(g, est.method = "YN", quiet = TRUE)
        expect_true(is.na(i$dN) && is.na(i$dS) && is.na(i$dNdS))
})
