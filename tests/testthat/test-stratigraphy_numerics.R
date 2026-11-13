# Divergence strata with ties, n_quantile edge cases, a finished run that is never lost, the dS = 0 policy,
# promoter pairing, input checks and saturation of Comeron's method. (Second audit: N07, N09, N10, N18, N19, N26,
# N27, N35, N37.)

# the binning of orthologr <= 0.4.3 (with the corrected probabilities)
old_strata <- function(x, n) {
        q <- stats::quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)
        if (any(duplicated(q)))
                stop("duplicate threshold")
        base::cut(x, breaks = q, include.lowest = TRUE, labels = FALSE)
}

gest <- function(a, b, ..., ids = c("a", "b")) {
        f <- tempfile(fileext = ".fa"); o <- tempfile()
        writeLines(c(paste0(">", ids[1]), a, paste0(">", ids[2]), b), f)
        gestimator(f, o, ...)
        x <- read_gestimator_table(o)
        for (j in c("dN", "dS", "dNdS")) x[[j]][x[[j]] == 999] <- NA
        as.data.frame(x)
}

random_codons <- function(n, seed) {
        set.seed(seed)
        cod <- apply(expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T")), 1, paste, collapse = "")
        sample(setdiff(cod, c("TAA", "TAG", "TGA")), n, TRUE)
}

# ---------------------------------------------------------------- divergence strata (N07, N26, N27)

test_that("without ties the strata are those of the quantile classes of earlier versions, for every n_quantile ...", {
        set.seed(11)
        x <- round(stats::rexp(1000, 3), 6)
        x <- x[!duplicated(x)]
        for (n in c(2:30, 49, 50, 98, 100)) {
                expect_identical(divergence_strata_of(x, n), old_strata(x, n), info = n)
                # the rank rule alone gives the same classes
                r <- rank(x, ties.method = "min")
                k <- pmax(1L, as.integer(((r - 1) * n + (length(x) - 2)) %/% (length(x) - 1)))
                expect_identical(k, old_strata(x, n), info = n)
        }
})

test_that("n_quantile = 49 and 98 give no NA stratum (N26) ...", {
        tb <- data.frame(query_id = paste0("q", 1:1000), subject_id = paste0("s", 1:1000), dNdS = seq(0.001, 1, length.out = 1000))
        for (n in 2:100) {
                r <- divergence_map(tb, n_quantile = n)
                expect_false(anyNA(r$divergence_strata), info = n)
                expect_equal(length(unique(r$divergence_strata)), n, info = n)
        }
})

test_that("many equal values: equal values share a stratum, strata stay ordered, a warning names the empty strata (N07) ...", {
        tb <- tibble::tibble(query_id = paste0("q", 30:1), subject_id = paste0("s", 30:1), dN = 0.1, dS = 0.2,
                             dNdS = c(rep(0, 6), seq(0.01, 0.24, length.out = 24)))
        expect_warning(r <- divergence_map(tb), "the most frequent value, 0, occurs 6 times among 30 values.*1 of the 10 strata is empty \\(2\\)")
        r <- as.data.frame(r)
        m <- merge(r, as.data.frame(tb), by = "query_id")
        expect_true(all(m$divergence_strata[m$dNdS == 0] == 1L))
        expect_equal(nrow(r), 30)
        expect_false(anyNA(r$divergence_strata))
        # monotone: a larger dNdS never has a smaller stratum
        o <- order(m$dNdS)
        expect_false(is.unsorted(m$divergence_strata[o]))
        expect_setequal(unique(m$divergence_strata), c(1L, 3:10))
        # the documented example: 280 of 1,000 genes with dNdS = 0
        x <- c(rep(0, 280), seq(0.001, 0.72, length.out = 720))
        k <- suppressWarnings(divergence_strata_of(x, 10))
        expect_identical(unique(k[1:280]), 1L)
        expect_identical(k[281], 3L)
        expect_false(2L %in% k)
        expect_identical(max(k), 10L)
})

test_that("divergence_map(): NA values, edge cases of n_quantile, and the table of the caller is left alone (N27) ...", {
        tb <- tibble::tibble(query_id = paste0("q", 30:1), subject_id = paste0("s", 30:1), dNdS = c(NA, NA, seq(0.01, 0.28, length.out = 28)))
        before <- tb
        r <- as.data.frame(divergence_map(tb, n_quantile = 4, subject.id = TRUE))
        expect_identical(tb, before)
        expect_s3_class(tb, "tbl_df")
        expect_identical(tb$query_id[1], "q30")
        expect_equal(sum(is.na(r$divergence_strata)), 2)
        expect_setequal(r$query_id[is.na(r$divergence_strata)], c("q30", "q29"))
        expect_equal(as.vector(table(r$divergence_strata)), c(7, 7, 7, 7))
        
        expect_identical(as.data.frame(divergence_map(tb, n_quantile = 1))$divergence_strata[1:3], rep(1L, 3))
        # more strata than values: every value its own stratum, the others empty
        small <- tb[3:6, ]
        r2 <- divergence_map(small, n_quantile = 10)
        expect_false(anyDuplicated(r2$divergence_strata) > 0)
        expect_identical(range(r2$divergence_strata), c(1L, 10L))
        for (bad in list(0, 2.5, NA, "10", c(5, 10), -1))
                expect_error(divergence_map(tb, n_quantile = bad), "n_quantile must be one whole number >= 1")
        expect_error(divergence_map(tb[ , c("query_id", "dNdS")]), "columns 'query_id', 'subject_id' and 'dNdS'")
        one <- tb[3, ]
        expect_identical(divergence_map(one)$divergence_strata, 1L)
})

test_that("divergence_stratigraphy() keeps the dN/dS table of the run and survives many identical genes (N07, N09) ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        # two 'genomes' in which 8 of 20 genes are identical (dN = dS = 0) and one differs by one nonsynonymous change
        s <- Biostrings::readDNAStringSet(lyra_cds())
        names(s) <- paste0("s", seq_along(s))
        q <- s
        names(q) <- paste0("q", seq_along(q))
        other <- Biostrings::readDNAStringSet(thal_cds())
        ref <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        # genes 1-11: the A. thaliana orthologs (diverged); genes 12-19: copies; gene 20: one changed codon
        pairs <- ref[seq_len(min(11, nrow(ref))), ]
        q_div <- other[match(pairs$query_id, sub(" .*", "", names(other)))]
        keep_s <- match(pairs$subject_id, sub(" .*", "", names(Biostrings::readDNAStringSet(lyra_cds()))))
        rest <- setdiff(seq_along(s), keep_s)
        q_same <- s[rest]
        x <- as.character(q_same[[length(q_same)]])
        substr(x, 4, 6) <- if (substr(x, 4, 6) == "GAT") "AAA" else "GAT"
        q_seqs <- c(as.character(q_div), as.character(q_same[-length(q_same)]), x)
        qf <- tempfile(fileext = ".fa"); sf <- tempfile(fileext = ".fa")
        writeLines(as.vector(rbind(paste0(">q", seq_along(q_seqs)), q_seqs)), qf)
        Biostrings::writeXStringSet(s, sf)
        
        w <- NULL
        ds <- withCallingHandlers(
                suppressMessages(divergence_stratigraphy(qf, sf, quiet = TRUE, subject.id = TRUE)),
                warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }
        )
        expect_s3_class(ds, "data.frame")
        expect_identical(names(ds), c("DS", "query_id", "subject_id"))
        full <- attr(ds, "dNdS_tbl")
        expect_true(all(c("query_id", "subject_id", "dN", "dS", "dNdS", "evalue") %in% names(full)))
        expect_equal(nrow(full), 20)
        n_zero <- sum(full$dS == 0, na.rm = TRUE)
        expect_gte(n_zero, 8)
        # kept by default: the genes with dS = 0 are rows with DS = NA
        expect_equal(nrow(ds), sum(!is.na(full$dN) & !is.na(full$dS) & (full$dS == 0 | full$dNdS <= 2), na.rm = TRUE))
        expect_equal(sum(is.na(ds$DS)), n_zero)
        # the behaviour of orthologr <= 0.4.3
        old <- suppressWarnings(suppressMessages(divergence_stratigraphy(qf, sf, quiet = TRUE, dS_zero = "drop")))
        expect_equal(nrow(old), nrow(ds) - n_zero)
        expect_false(anyNA(old$DS))
        # the run can be binned again from the attribute
        again <- suppressWarnings(suppressMessages(divergence_map(filter_dNdS(full), n_quantile = 3)))
        expect_equal(nrow(again), nrow(ds))
        
        # checked before the run
        expect_error(divergence_stratigraphy(qf, sf, n_quantile = 0), "n_quantile must be one whole number")
        expect_error(divergence_stratigraphy(qf, sf, dS_zero = "x"), "\"keep\" or \"drop\"")
})

# ---------------------------------------------------------------- gestimator: input checks (N18, N35, N37)

test_that("gestimator() stops with an R error on input that is not a codon alignment (N18) ...", {
        f <- tempfile(fileext = ".fa"); o <- tempfile()
        wr <- function(a, b) writeLines(c(">a", a, ">b", b), f)
        wr("AT", "AT")
        expect_error(gestimator(f, o), "shorter than one codon")
        wr("", "")
        expect_error(gestimator(f, o), "shorter than one codon|at least two sequences")
        wr(paste(rep("ATGGCTAAA", 4000), collapse = ""), "ATGGCT")
        expect_error(gestimator(f, o), "unequal lengths \\(36000 and 6\\)")
        wr("ATGGCTAA", "ATGGCTAA")
        expect_error(gestimator(f, o), "not a multiple of 3")
        writeLines(c(">a", "ATGGCTAAA"), f)
        expect_error(gestimator(f, o), "at least two sequences")
        expect_error(gestimator(tempfile(), o), "not found")
        wr("ATGGCTAAA", "ATGGCTAAA")
        expect_error(gestimator(f, o, maxHits = 7), "maxHits")
        
        # the same through substitutionrate(): the message names the pair
        wr(paste(rep("ATGGCTAAA", 10), collapse = ""), "ATGGCT")
        expect_error(substitutionrate(f, est.method = "Comeron", quiet = TRUE),
                     "Comeron's method \\(gestimator\\) failed for gene pair 'a' / 'b': the aligned sequences have unequal lengths")
        wr("ATGGCTAA", "ATGGCTAA")
        expect_error(substitutionrate(f, est.method = "Comeron", quiet = TRUE), "not a multiple of 3")
        expect_error(check_codon_alignment_file(f), "not a multiple of 3")
        wr("AT", "AT")
        expect_error(check_codon_alignment_file(f), "shorter than one codon")
})

test_that("gestimator() reads CRLF files, wrapped sequences and headers with descriptions (N35) ...", {
        a <- paste(random_codons(60, 1), collapse = "")
        b <- a; substr(b, 10, 12) <- "GCC"; substr(b, 100, 102) <- "AAA"
        ref <- gest(a, b)
        f <- tempfile(fileext = ".fa"); o <- tempfile()
        wrap <- function(x) substring(x, seq(1, nchar(x), 60), pmin(nchar(x), seq(60, nchar(x) + 59, 60)))
        con <- file(f, "wb")
        writeLines(c(">a desc>x more", wrap(a), ">b\tmore", wrap(b)), con, sep = "\r\n")
        close(con)
        gestimator(f, o)
        x <- read_gestimator_table(o)
        expect_identical(c(x$query_id, x$subject_id), c("a", "b"))
        expect_equal(x$dN, ref$dN)
        expect_equal(x$dS, ref$dS)
        y <- substitutionrate(f, est.method = "Comeron", quiet = TRUE)
        expect_equal(y$dN, ref$dN)
})

test_that("gestimator() frees what it allocates (N37) ...", {
        skip_on_cran()
        skip_on_os("windows")
        a <- paste(rep("ATGGCTAAACTGGGT", 400), collapse = "")
        f <- tempfile(fileext = ".fa"); o <- tempfile()
        writeLines(c(">a", a, ">b", a), f)
        rss <- function() as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)) / 1024
        for (i in 1:500) gestimator(f, o)
        r0 <- rss()
        for (i in 1:4000) gestimator(f, o)
        # it grew by about 6 MB per 4,000 calls before
        expect_lt(rss() - r0, 3)
})

# ---------------------------------------------------------------- Comeron: saturation (N19)

test_that("Comeron: an undefined correction gives NA, not a finite number; the option restores the old values (N19) ...", {
        # one synonymous transition at the only fourfold site: 1 - 2P - Q < 0
        new <- gest("ATGGCTAAA", "ATGGCCAAA")
        old <- gest("ATGGCTAAA", "ATGGCCAAA", legacy_saturation = TRUE)
        expect_true(is.na(new$dS))
        expect_true(is.na(new$dNdS))
        expect_equal(new$dN, 0)
        expect_equal(old$dS, 0)          # 'no synonymous divergence' for a pair with a synonymous difference
        
        # unrelated sequences
        a <- paste(random_codons(300, 1), collapse = "")
        b <- paste(random_codons(300, 2), collapse = "")
        new <- gest(a, b); old <- gest(a, b, legacy_saturation = TRUE)
        expect_true(is.na(new$dS) || is.na(new$dN))
        
        # nothing changes for pairs without saturation
        b2 <- a; substr(b2, 10, 12) <- "GCC"; substr(b2, 100, 102) <- "AAA"; substr(b2, 301, 303) <- "TTT"
        expect_identical(gest(a, b2), gest(a, b2, legacy_saturation = TRUE))
        expect_identical(gest(a, a), gest(a, a, legacy_saturation = TRUE))
        expect_equal(gest(a, a)$dN, 0)
        expect_equal(gest(a, a)$dS, 0)
        
        # the option is what substitutionrate() and dNdS() use
        f <- tempfile(fileext = ".fa")
        writeLines(c(">a", "ATGGCTAAA", ">b", "ATGGCCAAA"), f)
        expect_true(is.na(substitutionrate(f, est.method = "Comeron", quiet = TRUE)$dS))
        opt <- options(orthologr.legacy_comeron_saturation = TRUE)
        expect_equal(substitutionrate(f, est.method = "Comeron", quiet = TRUE)$dS, 0)
        options(opt)
})

test_that("Comeron on the example data: the option reproduces the values of earlier versions, the default flags saturation only ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        new <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        opt <- options(orthologr.legacy_comeron_saturation = TRUE)
        old <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        options(opt)
        expect_identical(new$query_id, old$query_id)
        same <- !is.na(new$dS) & !is.na(new$dN)
        expect_identical(new$dN[same], old$dN[same])
        expect_identical(new$dS[same], old$dS[same])
        # the 20 example pairs are far from saturation
        expect_identical(as.data.frame(new), as.data.frame(old))
})

# ---------------------------------------------------------------- promoters (N10)

test_that("promotor_index() finds promoters by gene ID, also for IDs with '_' ...", {
        nm <- c("G1_prom", "NM_001_1000bp", "NM_002_1000bp", "LOC_Os01g01010", "AT1G01010")
        expect_identical(promotor_index(c("G1", "NM_002", "NM_001", "LOC_Os01g01010", "AT1G01010", "NM", "G9"), nm),
                         c(1L, 3L, 2L, 4L, 5L, NA, NA))
})

test_that("promotor_divergence_of_orthologous_genes() pairs promoters by ortholog pair, not by position (N10) ...", {
        skip_if_not_installed("ape")
        d <- new_test_dir("prom_")
        set.seed(3)
        rs <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
        mut <- function(s, k) { x <- strsplit(s, "")[[1]]; i <- sample(length(x), k); x[i] <- sample(c("A", "C", "G", "T"), k, TRUE); paste(x, collapse = "") }
        genes <- paste0("G", 1:6)
        qp <- vapply(genes, function(g) rs(sample(180:220, 1)), "")
        sp <- vapply(qp, mut, "", k = 10)
        names(sp) <- paste0("S", 1:6)
        dir.create(file.path(d, "prom")); dir.create(file.path(d, "orth"))
        # the query lacks the promoter of G2, the subject that of S5: both files hold 5 sequences
        Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(qp[-2], paste0(genes[-2], "_prom"))), file.path(d, "prom", "Qsp_promotors.fa"))
        Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(sp[-5], paste0(names(sp)[-5], "_prom"))), file.path(d, "prom", "Ssp_promotors.fa"))
        ot <- data.frame(query_species = "Qsp", subject_species = "Ssp", query_id = paste0(genes, ".1"), query_gene_locus_id = genes,
                         subject_id = paste0(names(sp), ".1"), subject_gene_locus_id = names(sp), dN = 0.01, dS = 0.1, dNdS = 0.1, evalue = 0,
                         bit_score = 100, perc_identity = 90, num_ident_matches = 90, alig_length = 100L, mismatches = 1L, gap_openings = 0L,
                         n_gaps = 0, pos_match = 90, ppos = 90, q_start = 1L, q_end = 100L, q_len = 100, qcov = 100, qcovhsp = 100,
                         s_start = 1L, s_end = 100L, s_len = 100, score_raw = 100)
        utils::write.table(ot, file.path(d, "orth", "Qsp_vs_Ssp.csv"), sep = ";", row.names = FALSE, quote = FALSE)
        res <- NULL
        utils::capture.output(res <- suppressMessages(suppressWarnings(
                promotor_divergence_of_orthologous_genes(file.path(d, "prom"), file.path(d, "orth"),
                                                         ortholog_promotor_seq_output = file.path(d, "pout")))))
        expect_equal(nrow(res), 4)
        expect_identical(res$query_gene_locus_id, c("G1", "G3", "G4", "G6"))
        expect_identical(res$subject_gene_locus_id, c("S1", "S3", "S4", "S6"))
        # promoters of orthologs differ at about 10 of 200 sites; shifted pairs would be unrelated sequences
        expect_true(all(res$promotor_dist < 0.1))
})
