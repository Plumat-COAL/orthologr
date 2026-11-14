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
        # every documented KaKs_Calculator method (MLPB returned NULL up to this fix). On GitHub's Ubuntu runners
        # KaKs_Calculator now and then writes no result line for one of the model-based methods (GY, GMYN, MA, MS,
        # and so ALL; runs 43557520713 and 43563309087): that must be an error naming the method, never NULL.
        model_based <- c("GY", "GMYN", "MA", "MS")
        for (m in c("NG", "YN", "LWL", "MLWL", "LPB", "MLPB", "MYN", model_based)) {
                res <- tryCatch(substitutionrate(codon_aln_example, est.method = m, quiet = TRUE),
                                error = function(e) if (is.element(m, model_based)) e else stop(e))
                if (inherits(res, "error")) {
                        expect_match(conditionMessage(res), paste0("KaKs_Calculator.*method ", m), info = m)
                        next
                }
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
        all_m <- tryCatch(substitutionrate(codon_aln_example, est.method = "ALL", quiet = TRUE), error = function(e) e)
        if (inherits(all_m, "error")) {
                expect_match(conditionMessage(all_m), "KaKs_Calculator.*method ALL")
        } else {
                expect_gt(nrow(all_m), 10)
                expect_true(all(c("NG", "YN", "MLPB", "MYN") %in% all_m$method))
        }
})

test_that("every method that is.dnds_est_method() accepts gives a table in substitutionrate() (never NULL)", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        methods <- setdiff(c("Comeron", "Li", kaks_calculator_methods()), "kaks_calc")
        for (m in methods) {
                expect_true(is.dnds_est_method(m), info = m)
                # (see above: a run of KaKs_Calculator without result line is an error naming the method, never NULL)
                res <- tryCatch(substitutionrate(codon_aln_example, est.method = m, quiet = TRUE), error = function(e) e)
                if (inherits(res, "error")) {
                        expect_match(conditionMessage(res), paste0("KaKs_Calculator.*method ", m), info = m)
                        next
                }
                expect_s3_class(res, "data.table")
                expect_true(nrow(res) >= 1 && all(c("dN", "dS", "dNdS") %in% names(res)), info = m)
        }
})


# Comeron, decisions D8 / D9 of 2026-11-14: codons that differ at two positions are counted along Comeron's two
# pathways by default (libsequence's branches made the estimate depend on which sequence is first;
# options(orthologr.comeron_twosubs = "libsequence") keeps them), and a site class with fewer than
# options(orthologr.comeron_min_sites) (3) sites cannot make dN or dS NA. The COAL engine gets the same settings.
comeron_random_pairs <- function(n, codons, seed) {
        set.seed(seed)
        bases <- c("T", "C", "A", "G")
        stops <- c("TAA", "TAG", "TGA")
        sense <- setdiff(apply(expand.grid(bases, bases, bases)[, 3:1], 1, paste0, collapse = ""), stops)
        mut2 <- function(c) repeat {
                d <- strsplit(c, "")[[1]]
                for (p in sample(3, 2)) d[p] <- sample(setdiff(bases, d[p]), 1)
                d <- paste0(d, collapse = "")
                if (!is.element(d, stops)) return(d)
        }
        replicate(n, {
                x <- sample(sense, codons, TRUE)
                y <- vapply(x, function(c) if (stats::runif(1) < 0.3) mut2(c) else c, "")
                c(paste0(x, collapse = ""), paste0(y, collapse = ""))
        })
}

comeron_native <- function(a, b) {
        f <- tempfile(fileext = ".fasta")
        on.exit(unlink(f))
        writeLines(c(">a", a, ">b", b), f)
        r <- substitutionrate(f, est.method = "Comeron", quiet = TRUE, format = "fasta")
        c(dN = r$dN, dS = r$dS)
}

test_that("Comeron: the estimate does not depend on which sequence is first; libsequence's branches as an option (D8)", {
        p <- comeron_random_pairs(40, 150, 20261114)
        ab <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
        ba <- t(apply(p, 2, function(x) comeron_native(x[2], x[1])))
        expect_identical(ab, ba)
        old <- options(orthologr.comeron_twosubs = "libsequence")
        on.exit(options(old), add = TRUE)
        lab <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
        lba <- t(apply(p, 2, function(x) comeron_native(x[2], x[1])))
        expect_gt(sum(lab[, "dN"] != lba[, "dN"]), 20)
        expect_false(identical(lab, ab))
        options(orthologr.comeron_twosubs = "both")
        expect_error(comeron_native(p[1, 1], p[2, 1]), class = "orthologr_error_bad_argument")
})

test_that("Comeron: the native engine and COAL agree to the printed digits with both pathway orders (D8)", {
        skip_if_not_installed("coalr")
        skip_if_not(is.element("comeron_twosubs", names(formals(coalr::coal_kaks))), "coalr predates comeron_twosubs")
        p <- comeron_random_pairs(40, 150, 7)
        for (ts in c("comeron", "libsequence")) {
                old <- options(orthologr.comeron_twosubs = ts)
                nat <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
                co <- do.call(coalr::coal_kaks, c(list(p[1, ], p[2, ], method = "Comeron"),
                                                  orthologr:::coal_flavour_args(coalr::coal_kaks, "Comeron")))
                co <- orthologr:::coal_orthologr_values(co)
                expect_identical(unname(nat[, "dN"]), co$dN, info = ts)
                expect_identical(unname(nat[, "dS"]), co$dS, info = ts)
                options(old)
        }
})

test_that("Comeron: a site class with fewer than 3 sites cannot make dN NA; options(orthologr.comeron_min_sites = 0) restores it (D9)", {
        one <- comeron_native("GCT", "GAT")                 # one codon, nonsynonymous: the 2V class has < 1 site
        expect_false(is.na(one[["dN"]]))
        old <- options(orthologr.comeron_min_sites = 0)
        on.exit(options(old), add = TRUE)
        expect_true(is.na(comeron_native("GCT", "GAT")[["dN"]]))
        # never changes a value that is not NA otherwise
        p <- comeron_random_pairs(30, 30, 11)
        v0 <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
        options(orthologr.comeron_min_sites = 3)
        v3 <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
        expect_identical(v3[!is.na(v0)], v0[!is.na(v0)])
        options(orthologr.comeron_min_sites = -1)
        expect_error(comeron_native("GCT", "GAT"), class = "orthologr_error_bad_argument")
        skip_if_not_installed("coalr")
        skip_if_not(is.element("comeron_min_sites", names(formals(coalr::coal_kaks))), "coalr predates comeron_min_sites")
        # the COAL engine: the same rule (its note for such a class is not "correction undefined")
        for (ms in c(0, 3)) {
                options(orthologr.comeron_min_sites = ms)
                co <- do.call(coalr::coal_kaks, c(list("GCT", "GAT", method = "Comeron"),
                                                  orthologr:::coal_flavour_args(coalr::coal_kaks, "Comeron")))
                co <- orthologr:::coal_orthologr_values(co)
                expect_identical(co$dN, comeron_native("GCT", "GAT")[["dN"]], info = ms)
        }
})

# D9 refinement (COAL #76): the rule of small classes only applies where it gives a distance > 0. A small class's
# undefined logarithm set to 0 could leave dS at exactly 0 although there was a synonymous change (ATGGCTAAA /
# ATGGCCAAA), which reads as "no synonymous change" (dS_zero); such a value stays NA, as with min_sites = 0.
comeron_short_pairs <- function(per_length, seed) {
        set.seed(seed)
        bases <- c("T", "C", "A", "G")
        stops <- c("TAA", "TAG", "TGA")
        sense <- setdiff(apply(expand.grid(bases, bases, bases)[, 3:1], 1, paste0, collapse = ""), stops)
        ts <- c("AG", "GA", "CT", "TC")
        out <- list()
        for (n in c(5, 10, 15, 20)) for (k in seq_len(per_length)) {
                anc <- sample(sense, n, TRUE)
                s <- anc
                need <- max(1, round(stats::runif(1, 0.02, 0.35) * n))
                while (sum(s != anc) < need) {
                        i <- sample(n, 1); pos <- sample(3, 1)
                        d <- strsplit(s[i], "")[[1]]
                        opts <- setdiff(bases, d[pos])
                        d[pos] <- sample(opts, 1, prob = ifelse(paste0(d[pos], opts) %in% ts, 2, 1))
                        m <- paste0(d, collapse = "")
                        if (!is.element(m, stops)) s[i] <- m
                }
                out[[length(out) + 1]] <- c(paste0(anc, collapse = ""), paste0(s, collapse = ""))
        }
        do.call(cbind, out)
}

test_that("Comeron: the min-sites rule only gives values > 0; a 0 from a small class stays NA (D9 refinement)", {
        old <- options(orthologr.comeron_min_sites = 3)
        on.exit(options(old), add = TRUE)
        # one synonymous change; the fourfold class has 1 site: dS was 0 (= "no synonymous change"), now NA
        tiny <- comeron_native("ATGGCTAAA", "ATGGCCAAA")
        expect_true(is.na(tiny[["dS"]]))
        expect_equal(tiny[["dN"]], 0)
        p <- comeron_short_pairs(40, 20261114)
        v3 <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
        options(orthologr.comeron_min_sites = 0)
        v0 <- t(apply(p, 2, function(x) comeron_native(x[1], x[2])))
        # a value that is not NA without the rule is unchanged; a value the rule gives in place of NA is > 0
        expect_identical(v3[!is.na(v0)], v0[!is.na(v0)])
        recovered <- is.na(v0) & !is.na(v3)
        expect_gt(sum(recovered), 0)
        expect_true(all(v3[recovered] > 0))
        # the legacy switch is not affected: 0.4.3's finite values, 0 included
        options(orthologr.comeron_min_sites = 3, orthologr.legacy_comeron_saturation = TRUE)
        expect_equal(comeron_native("ATGGCTAAA", "ATGGCCAAA")[["dS"]], 0)
        options(orthologr.legacy_comeron_saturation = FALSE)
        # the COAL engine (coalr with COAL #76) gives the same NA pattern
        skip_if_not_installed("coalr")
        skip_if_not(is.element("comeron_min_sites", names(formals(coalr::coal_kaks))), "coalr predates comeron_min_sites")
        co <- do.call(coalr::coal_kaks, c(list("ATGGCTAAA", "ATGGCCAAA", method = "Comeron"),
                                          orthologr:::coal_flavour_args(coalr::coal_kaks, "Comeron")))
        skip_if_not(is.na(orthologr:::coal_orthologr_values(co)$dS), "coalr predates COAL's D9 refinement (#76)")
        co <- do.call(coalr::coal_kaks, c(list(p[1, ], p[2, ], method = "Comeron"),
                                          orthologr:::coal_flavour_args(coalr::coal_kaks, "Comeron")))
        co <- orthologr:::coal_orthologr_values(co)
        expect_identical(is.na(co$dS), is.na(unname(v3[, "dS"])))
        expect_identical(is.na(co$dN), is.na(unname(v3[, "dN"])))
})
