# The counts of the LWL family recomputed by orthologr (R/lwl_counts.R). Expected values are those KaKs_Calculator
# 2.0 printed for the same input ('Fold-Sites', 'Fold-S-Substitutions', 'Fold-N-Substitutions'; 6 significant
# digits) - see CLEANROOM.md of branch chore/gpl2-cleanroom.

lwl_printed <- function(e) {
        g <- function(v) sub("^\\s+", "", formatC(v, digits = 6, format = "g"))
        fold <- function(v, na_when_zero) if (na_when_zero && all(v == 0)) "NA" else paste(g(v), collapse = ":")
        c(fold(e$L, TRUE), fold(e$Si, TRUE), fold(e$Vi, TRUE))
}

lwl_cases <- list(
        # the arginine and isoleucine codon pairs of Tzeng et al. (2004), with a background
        arg_ile = c("CGAAGACGGAGGATAATAATTATGGCTTTTAAACATGCTGGCCCAGTGTTTAAAGATCATGCTGGCCCAGTGTTTAAAGATCATGCTGGCCCAGTGTTTAAAGATCAT",
                    "AGACGAAGGCGGATTATGATAATACCTTTCCAGCACGCCGGCCCAGTGTTTAAAGATCATGCCGGCCCAGTGTTTAAAGATCATGCCGGCCCAGTGTTTAAAGATCAT"),
        # stop codon, lower case, gap, N, ambiguity code
        odd = c("ATGTAAgct---NNNGCTTTTAGRCATGGGTTTCCCAAAGGGTTTCCC",
                "ATGCAAGCCGCTGCTGCATTCAGACAAGGATTCCCGAAGGAATTTCCA"),
        # codons that differ at two and three positions (pathways, some through stop codons)
        multi = c("ATGGCTTGGCCCAAAGATGAGCGTTTATCACATCAGGGGTGGATAACTGACATG",
                  "ATGGAGTGTCAAAGGTAAGAAAGCTCAGTCCGTCATGGCTGCATGACCTGGATG"))

lwl_expected <- utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
code case method L S V
1 arg_ile LWL 69:24:15 1:6:3 2:4:0
1 arg_ile MLWL 69:24:15 0:9:3 3:1:0
1 arg_ile LPB 69:24:15 1:6:3 2:4:0
1 arg_ile MLPB 69:24:15 0:9:3 3:1:0
1 odd LWL 25:5.5:5.5 1:3.5:2.5 0:1:3
1 odd MLPB 25:5.5:5.5 1:3.5:2.5 0:1:3
1 multi LWL 37:8.5:5.5 5.5:3:1.5 6.16667:3.75:3.08333
1 multi MLWL 37:8.5:5.5 5:2.5:1.5 6.66667:4.25:3.08333
1 multi LPB 37:8.5:5.5 5.5:3:1.5 6.16667:3.75:3.08333
1 multi MLPB 37:8.5:5.5 5:2.5:1.5 6.66667:4.25:3.08333
2 arg_ile LWL 64:19:13 0:5:3 2:2:0
2 arg_ile MLWL 64:19:13 0:5:3 2:2:0
2 odd LWL 24:6.5:5.5 1:3.5:2.5 0:1:3
2 multi LWL 31.5:11:5.5 4:2.5:1.5 4.91667:5:3.08333
2 multi MLWL 31.5:11:5.5 4:1.5:1.5 4.91667:6:3.08333
13 arg_ile LWL 70:23:15 1:6:3 3:3:0
13 arg_ile MLWL 70:23:15 2:7:3 2:2:0
13 odd MLWL 22.5:8:5.5 1:3.5:2.5 0:1:3
13 multi LPB 32.5:13:5.5 5:3.5:1.5 4.66667:5.25:3.08333
13 multi MLPB 32.5:13:5.5 5:2.5:1.5 4.66667:6.25:3.08333
")

test_that("kaks_lwl_counts() gives the counts KaKs_Calculator prints (fixed cases, codes 1, 2, 13)", {
        for (r in seq_len(nrow(lwl_expected))) {
                s <- lwl_cases[[lwl_expected$case[r]]]
                e <- kaks_lwl_counts(s[1], s[2], lwl_expected$method[r], lwl_expected$code[r])
                expect_identical(lwl_printed(e), unlist(lwl_expected[r, c("L", "S", "V")], use.names = FALSE),
                                 info = paste(lwl_expected[r, 1:3], collapse = " "))
                expect_identical(e$P, e$Si / e$L)
                expect_identical(e$Q, e$Vi / e$L)
        }
})

test_that("the genetic codes and site classes are NCBI's, for the codes KaKs_Calculator has", {
        t1 <- kaks_transl_table(1)
        expect_identical(paste(t1, collapse = ""), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
        expect_identical(names(t1)[c(1, 2, 3, 64)], c("TTT", "TTC", "TTA", "GGG"))
        expect_identical(unname(kaks_transl_table(2)[c("AGA", "ATA", "TGA")]), c("*", "M", "W"))
        for (code in c(7, 8, 17, 24, 33, NA))
                expect_null(kaks_transl_table(code))
        expect_length(names(kaks_genetic_codes()), 17)
        # the threefold third position of ATA counts as twofold; stop codons have no sites
        L <- kaks_lwl_tables(1, modified = FALSE)$L
        expect_equal(unname(L["ATA", ]), c(2, 1, 0))
        expect_equal(unname(L["TAA", ]), c(0, 0, 0))
        expect_equal(unname(L["CGA", ]), c(1, 1, 1))
        expect_equal(unname(kaks_lwl_tables(2, modified = FALSE)$L["ATA", ]), c(2, 1, 0))
        expect_null(kaks_lwl_tables(NA, modified = FALSE))
})

test_that("kaks_lwl_counts() returns NULL where it cannot compute", {
        expect_null(kaks_lwl_counts("ATG", "ATG", "YN"))
        expect_null(kaks_lwl_counts("ATG", "ATGA", "LWL"))
        expect_null(kaks_lwl_counts("ATGA", "ATGA", "LWL"))
        expect_null(kaks_lwl_counts("ATG", "ATG", "LWL", genetic_code = 7))
        expect_null(kaks_lwl_counts("ATG", "ATG", "LWL", genetic_code = 24))
        expect_null(kaks_lwl_counts("", "", "LWL"))
        expect_null(kaks_lwl_counts(c("ATG", "ATG"), "ATG", "LWL"))
        # RNA letters: left to the counts KaKs_Calculator printed
        expect_null(kaks_lwl_counts("AUGTTT", "ATGTTC", "LWL"))
        expect_null(kaks_lwl_counts("ATGTTT", "augttc", "MLPB"))
        # a stop codon pair and an identical pair
        e <- kaks_lwl_counts("TAAATG", "TAAATG", "MLWL")
        expect_identical(e$L, c(3, 0, 0))
        expect_identical(e$Si + e$Vi, c(0, 0, 0))
})

test_that("kaks_lwl_counts() adds the counts of the codon pairs in alignment order, in double precision", {
        # Exact-zero cases: 1 - 2P - Q (or 1 - 2Q) of one class is 0 in exact arithmetic, so its sign comes from the
        # rounding of the running totals. 'kc' is the decision KaKs_Calculator 2.0 took (from its printed Ka/Ks).
        # Sequential double-precision addition in alignment order reproduces it in every case; sum() (extended
        # precision) gives the other decision in each of these cases.
        x <- utils::read.delim(test_path("fixtures", "lwl_boundary.tsv"), colClasses = "character",
                               na.strings = character())
        expect_gte(nrow(x), 10)
        for (r in seq_len(nrow(x))) {
                e <- kaks_lwl_counts(x$s1[r], x$s2[r], x$method[r], as.integer(x$code[r]))
                k <- as.integer(x$target[r]) + 1
                defined <- (1 - 2 * e$P[k] - e$Q[k]) > 0 && (1 - 2 * e$Q[k]) > 0
                expect_identical(if (defined) "nonsat" else "sat", x$kc[r],
                                 info = paste("row", r, x$method[r], "class", x$target[r]))
        }
        # a codon pair that differs at three positions has counts in twelfths, which are not exact in binary
        e <- kaks_lwl_counts("GCTGCTGCT", "ATCGCCGCC", "LWL")
        expect_identical(e$Si, Reduce(`+`, list(kaks_lwl_counts("GCT", "ATC", "LWL")$Si,
                                                 kaks_lwl_counts("GCT", "GCC", "LWL")$Si,
                                                 kaks_lwl_counts("GCT", "GCC", "LWL")$Si), c(0, 0, 0)))
})

test_that("kaks_lwl_counts() agrees with KaKs_Calculator on random pairs of every genetic code", {
        skip_on_cran()
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        dir <- tempfile("lwl")
        dir.create(dir)
        on.exit(unlink(dir, recursive = TRUE), add = TRUE)
        set.seed(20261116)
        codons <- kaks_codons()
        for (code in as.integer(names(kaks_genetic_codes()))) {
                n <- 8
                s1 <- s2 <- character(n)
                for (i in seq_len(n)) {
                        a <- sample(codons, sample(5:120, 1), replace = TRUE)
                        b <- vapply(a, function(x) {
                                y <- strsplit(x, "")[[1]]
                                k <- runif(3) < runif(1, 0, 0.8)
                                y[k] <- sample(c("T", "C", "A", "G"), sum(k), replace = TRUE)
                                paste(y, collapse = "")
                        }, "")
                        s1[i] <- paste(a, collapse = "")
                        s2[i] <- paste(b, collapse = "")
                }
                axt <- file.path(dir, paste0("c", code, ".axt"))
                out <- file.path(dir, paste0("c", code, ".kaks"))
                writeLines(as.vector(rbind(paste0("p", seq_len(n)), s1, s2, "")), axt)
                system2(Sys.which("KaKs_Calculator"), c("-i", axt, "-o", out, "-c", code, "-m", "LWL", "-m", "MLWL",
                                                       "-m", "LPB", "-m", "MLPB"), stdout = FALSE, stderr = FALSE)
                tbl <- utils::read.delim(out, colClasses = "character", na.strings = character(), check.names = FALSE)
                for (r in seq_len(nrow(tbl))) {
                        i <- as.integer(sub("^p", "", tbl[r, 1]))
                        e <- kaks_lwl_counts(s1[i], s2[i], tbl[r, 2], code)
                        expect_identical(lwl_printed(e), unname(unlist(tbl[r, c(10, 14, 15)])),
                                         info = paste("code", code, tbl[r, 1], tbl[r, 2]))
                }
        }
})
