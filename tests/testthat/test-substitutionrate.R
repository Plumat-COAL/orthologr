codon_aln_example <- system.file("seqs/pal2nal.aln", package = "orthologr")

test_that("substitutionrate() computes Comeron's dN, dS and dN/dS for the example codon alignment", {
        res <- substitutionrate(codon_aln_example, est.method = "Comeron", quiet = TRUE)
        expect_s3_class(res, "data.table")
        expect_identical(names(res), c("query_id", "subject_id", "dN", "dS", "dNdS"))
        expect_identical(res$query_id, "AT1G01030.1")
        expect_identical(res$subject_id, "470180|PACid_16054974")
        expect_equal(res$dN, 0.01499, tolerance = 1e-4)
        expect_equal(res$dS, 0.1265, tolerance = 1e-4)
        expect_equal(res$dNdS, 0.1185, tolerance = 1e-3)
})

test_that("substitutionrate() computes Li's dN and dS with seqinr::kaks()", {
        res <- substitutionrate(codon_aln_example, est.method = "Li", quiet = TRUE)
        k <- seqinr::kaks(seqinr::read.alignment(codon_aln_example, format = "fasta"))
        expect_equal(res$dN, as.numeric(k$ka))
        expect_equal(res$dS, as.numeric(k$ks))
        expect_equal(res$dNdS, as.numeric(k$ka) / as.numeric(k$ks))
})

test_that("substitutionrate() checks method and format", {
        expect_error(substitutionrate(codon_aln_example, est.method = "XX"), "dNdS estimation method")
        expect_error(substitutionrate(codon_aln_example, est.method = "Comeron", format = "nexus"), "format")
        expect_error(substitutionrate(codon_aln_example, est.method = "Comeron", format = "clustal"), "fasta")
})

test_that("substitutionrate() runs the KaKs_Calculator models", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        # every documented KaKs_Calculator method (MLPB returned NULL up to this fix)
        for (m in c("NG", "YN", "LWL", "MLWL", "LPB", "MLPB", "MYN", "GY", "GMYN", "MA", "MS")) {
                res <- substitutionrate(codon_aln_example, est.method = m, quiet = TRUE)
                expect_s3_class(res, "data.table")
                # (KaKs_Calculator names some methods with their model, e.g. "GY-HKY")
                expect_true(startsWith(res$method, m), info = m)
                expect_identical(res$query_id, "AT1G01030.1")
                expect_true(res$dN > 0 && res$dN < 0.05, info = m)
                expect_true(res$dS > 0.05 && res$dS < 0.3, info = m)
                expect_equal(res$dNdS, res$dN / res$dS, tolerance = 1e-3, info = m)
        }
        both <- substitutionrate(codon_aln_example, est.method = "kaks_calc", kaks_calc.params = "-m NG -m YN", quiet = TRUE)
        expect_setequal(both$method, c("NG", "YN"))
        all_m <- substitutionrate(codon_aln_example, est.method = "ALL", quiet = TRUE)
        expect_gt(nrow(all_m), 10)
        expect_true(all(c("NG", "YN", "MLPB", "MYN") %in% all_m$method))
})

test_that("every method that is.dnds_est_method() accepts gives a table in substitutionrate() (never NULL)", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        methods <- setdiff(c("Comeron", "Li", kaks_calculator_methods()), "kaks_calc")
        for (m in methods) {
                expect_true(is.dnds_est_method(m), info = m)
                res <- substitutionrate(codon_aln_example, est.method = m, quiet = TRUE)
                expect_s3_class(res, "data.table")
                expect_true(nrow(res) >= 1 && all(c("dN", "dS", "dNdS") %in% names(res)), info = m)
        }
})
