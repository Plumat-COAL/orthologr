# Coverage of src/Kimura80.h through substitutionrate(est.method = "Comeron") (gestimator, src/Comeron95.h).
# Comeron's method asks for Kimura's two-parameter distance of the whole pair when its nonsynonymous transversion
# term is not finite (1 - 2 Qa <= 0, i.e. at least half of the nonsynonymous sites differ by a transversion):
# a K80 distance below 1 sets that term to 0 (libsequence's ad-hoc rule); with options(orthologr.legacy_comeron_saturation
# = FALSE), the default, orthologr then calls the pair saturated at nonsynonymous sites (dN NA).
# The pairs below reach the three outcomes of Kimura80::Compute(): 1 - 2P - Q == 0 (distance 0), a finite
# distance < 1, and a distance that is not finite (K() returns 999). No external program is needed.

cov_k80_pair <- function(q, s) {
        f <- tempfile("cov_k80_", fileext = ".fasta")
        writeLines(c(">query", q, ">subject", s), f)
        f
}
cov_k80 <- function(q, s, legacy = FALSE) {
        old <- options(orthologr.legacy_comeron_saturation = legacy)
        on.exit(options(old), add = TRUE)
        as.data.frame(substitutionrate(cov_k80_pair(q, s), est.method = "Comeron", quiet = TRUE))
}

test_that("Comeron: every nonsynonymous site a transversion (K80: 1 - 2P - Q = 0, distance 0)", {
        # GGG (Gly) -> CCC (Pro): transversions at all positions; position 3 is fourfold degenerate in both
        q <- strrep("GGG", 12)
        s <- strrep("CCC", 12)
        r <- cov_k80(q, s)
        expect_identical(names(r), c("query_id", "subject_id", "dN", "dS", "dNdS"))
        expect_true(is.na(r$dN))
        expect_true(is.na(r$dNdS))
        # the legacy rule keeps libsequence's value (the term set to 0 because K80 < 1)
        r_legacy <- cov_k80(q, s, legacy = TRUE)
        expect_true(is.finite(r_legacy$dN))
        expect_equal(r_legacy$dN, r_legacy$dN)
})

test_that("Comeron: half of the nonsynonymous sites transversions (K80 finite and below 1)", {
        # 6 codons GGG -> CCG (two transversions at nondegenerate positions, third position kept) and 6 identical
        # codons; one synonymous transition (GGG -> GGA) and one gapped codon for the other branches of Kimura80
        q <- paste0(strrep("GGG", 6), strrep("GGG", 6), "GGG", "GGG")
        s <- paste0(strrep("CCG", 6), strrep("GGG", 6), "GGA", "---")
        r <- cov_k80(q, s)
        r_legacy <- cov_k80(q, s, legacy = TRUE)
        expect_true(is.finite(r_legacy$dN))
        expect_true(is.finite(r_legacy$dS))
        expect_true(is.na(r$dN) || isTRUE(all.equal(r$dN, r_legacy$dN)))
        expect_equal(r$dS, r_legacy$dS)
})

test_that("Comeron: most sites transversions (K80 not finite: 999, the term stays undefined)", {
        # 10 of 12 codons GGG -> CCC: more than half of all sites differ by a transversion, 1 - 2Q < 0
        q <- strrep("GGG", 12)
        s <- paste0(strrep("CCC", 10), strrep("GGG", 2))
        r <- cov_k80(q, s)
        expect_true(is.na(r$dN))
        expect_true(is.na(r$dNdS))
        r_legacy <- cov_k80(q, s, legacy = TRUE)
        expect_identical(names(r_legacy), names(r))
})
