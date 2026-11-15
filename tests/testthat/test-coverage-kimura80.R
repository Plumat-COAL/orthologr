# Coverage of src/Kimura80.h through substitutionrate(est.method = "Comeron") (gestimator, src/Comeron95.h).
# Comeron's method asks for Kimura's two-parameter distance (K80) of the whole pair when its nonsynonymous
# transversion term Ba is not finite (1 - 2 Qa <= 0: at least half of the nonsynonymous sites differ by a
# transversion). A K80 distance below 1 sets Ba to 0 (libsequence's ad-hoc rule); orthologr then calls the pair
# saturated at nonsynonymous sites (dN NA) unless options(orthologr.legacy_comeron_saturation = TRUE).
# The three pairs reach the three outcomes of Kimura80::Compute() / K(): 1 - 2P - Q = 0 (distance 0), a finite
# distance below 1 (with a transition, so both Ts and Tv are counted), and a distance that is not finite (999).
# (Kimura80's message for sequences of unequal length cannot be reached through gestimator, and a negative
# distance cannot arise.) No external program is needed.

cov_k80 <- function(q, s, legacy = FALSE) {
        old <- options(orthologr.legacy_comeron_saturation = legacy)
        on.exit(options(old), add = TRUE)
        f <- tempfile("cov_k80_", fileext = ".fasta")
        writeLines(c(">query", q, ">subject", s), f)
        as.data.frame(substitutionrate(f, est.method = "Comeron", quiet = TRUE))
}

test_that("Comeron / K80: every site a transversion (1 - 2P - Q = 0, distance 0 < 1, Ba set to 0)", {
        # GGG (Gly) -> CCC (Pro): transversions at all three positions
        q <- strrep("GGG", 12)
        s <- strrep("CCC", 12)
        r <- cov_k80(q, s)
        expect_identical(names(r), c("query_id", "subject_id", "dN", "dS", "dNdS"))
        expect_true(is.na(r$dN))
        expect_true(is.na(r$dNdS))
        # libsequence's value: Ba = 0 because K80 < 1
        r_legacy <- cov_k80(q, s, legacy = TRUE)
        expect_equal(r_legacy$dN, 1.54)
        expect_true(is.na(r_legacy$dS))
})

test_that("Comeron / K80: transversions at the nondegenerate sites of most codons and a transition (finite K80 < 1)", {
        # 8 codons GGG -> CCG (transversions at the two nondegenerate positions), 5 identical codons, and one
        # synonymous transition GGG -> GGA
        q <- strrep("GGG", 14)
        s <- paste0(strrep("CCG", 8), strrep("GGG", 5), "GGA")
        r <- cov_k80(q, s)
        r_legacy <- cov_k80(q, s, legacy = TRUE)
        expect_true(is.na(r$dN))
        expect_true(is.na(r$dNdS))
        expect_true(is.finite(r$dS))
        expect_true(is.finite(r_legacy$dN) && r_legacy$dN > 0)
        expect_equal(r_legacy$dS, r$dS)
        expect_equal(r_legacy$dNdS, signif(r_legacy$dN / r_legacy$dS, 4), tolerance = 1e-3)
})

test_that("Comeron / K80: more than half of all sites transversions (K80 not finite: 999, Ba stays undefined)", {
        # 10 of 12 codons GGG -> CCC: 1 - 2Q < 0
        q <- strrep("GGG", 12)
        s <- paste0(strrep("CCC", 10), strrep("GGG", 2))
        for (legacy in c(FALSE, TRUE)) {
                r <- cov_k80(q, s, legacy = legacy)
                expect_true(is.na(r$dN))
                expect_true(is.na(r$dNdS))
        }
})
