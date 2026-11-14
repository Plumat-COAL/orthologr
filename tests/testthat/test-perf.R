# orthologr 0.6.0 performance: the batched code paths of compute_dnds() give exactly what the code that handled
# one gene pair at a time gave. The per-pair code is still there (pairwise_aln(), codon_aln()) and is the reference:
# it is used for every pair when pairwise_aln_rows() returns nothing.

# The gene pairs of the 20-gene example in the form dNdS() hands to compute_dnds() (pairs by position; the IDs are
# made unique per pair).
example_pair_table <- function(n = 20) {
        q <- suppressMessages(read.cds(thal_cds(), format = "fasta", delete_corrupt_cds = FALSE))
        s <- suppressMessages(read.cds(lyra_cds(), format = "fasta", delete_corrupt_cds = FALSE))
        n <- min(n, nrow(q), nrow(s))
        data.table::data.table(query_id = q$geneids[seq_len(n)], subject_id = s$geneids[seq_len(n)],
                               query_cds = q$seqs[seq_len(n)], subject_cds = s$seqs[seq_len(n)],
                               query_aa = toupper(translate_cds(q$seqs[seq_len(n)])),
                               subject_aa = toupper(translate_cds(s$seqs[seq_len(n)])))
}

# compute_dnds() of the pair table with the batched paths, and with every pair through pairwise_aln() / codon_aln()
# (pairwise_aln_rows() returning NULL for each pair, as it does when it cannot handle a pair)
both_ways <- function(tbl, ...) {
        new <- suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                             quiet = TRUE, ...))
        local({
                testthat::local_mocked_bindings(pairwise_aln_rows = function(aa1, ...) vector("list", length(aa1)),
                                                .package = "orthologr")
                old <- suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                                     quiet = TRUE, ...))
                # (a failure of PAL2NAL names the temporary file it was to write, a new name in every call)
                tmp_names <- function(x) {
                        f <- attr(x, "failed_pairs")
                        f$reason <- gsub("aln_[0-9a-f]+\\.tmp", "aln_<tmp>.tmp", f$reason)
                        attr(x, "failed_pairs") <- f
                        x
                }
                list(new = tmp_names(new), old = tmp_names(old))
        })
}

test_that("first_field() gives what unlist(stringr::str_split(x, sep))[1] gave", {
        expect_identical(first_field("Ath.cds.fa", "."), "Ath")
        expect_identical(first_field("Ath_thal_cds", "_"), "Ath")
        expect_identical(first_field("AT1G01010.1 gene=x", " "), "AT1G01010.1")
        expect_identical(first_field("nosep", "."), "nosep")
        expect_identical(first_field("", "."), "")
        expect_identical(first_field(".hidden", "."), "")
        expect_identical(first_field(NA_character_, "."), NA_character_)
        expect_identical(first_field(c("a.b", "c.d.e"), "."), c("a", "c"))
        # the separator is taken literally ('.' is not "any character")
        expect_identical(first_field("ab", "."), "ab")
})

test_that("write_pair_fasta() writes the bytes seqinr::write.fasta() wrote for a gene pair", {
        f_old <- tempfile(); f_new <- tempfile()
        seqs <- list(paste(rep("ATGAAACCC", 50), collapse = ""), "ATGTAA")
        seqinr::write.fasta(sequences = seqs, names = as.list(c("query", "subject")), file.out = f_old)
        write_pair_fasta(seqs, c("query", "subject"), f_new)
        expect_identical(unname(tools::md5sum(f_new)), unname(tools::md5sum(f_old)))
        # an empty sequence and NA
        seqinr::write.fasta(sequences = list("", NA_character_), names = as.list(c("query", "subject")), file.out = f_old)
        write_pair_fasta(list("", NA_character_), c("query", "subject"), f_new)
        expect_identical(unname(tools::md5sum(f_new)), unname(tools::md5sum(f_old)))
})

test_that("aligned_rows_list() gives aligned_rows() of each alignment of a vectorised pwalign call", {
        tbl <- example_pair_table()
        a <- sub("\\*$", "", tbl$query_aa)
        b <- sub("\\*$", "", rev(tbl$subject_aa))     # also unrelated pairs
        aln <- pwalign::pairwiseAlignment(Biostrings::AAStringSet(a), Biostrings::AAStringSet(b),
                                          type = "global", substitutionMatrix = "BLOSUM62")
        rows <- aligned_rows_list(aln, a, b)
        for (i in seq_along(a)) {
                one <- pwalign::pairwiseAlignment(Biostrings::AAString(a[i]), Biostrings::AAString(b[i]),
                                                  type = "global", substitutionMatrix = "BLOSUM62")
                expect_identical(c(rows[[1]][i], rows[[2]][i]), aligned_rows(one, a[i], b[i]), info = i)
                expect_identical(pwalign::score(aln)[i], pwalign::score(one), info = i)
        }
})

test_that("pairwise_aln_rows() gives the rows pairwise_aln() writes for each pair", {
        tbl <- example_pair_table()
        pairs <- list(q = tbl$query_aa, s = tbl$subject_aa)
        pairs$s[3] <- rev(tbl$subject_aa)[3]          # an unrelated pair
        for (m in c("BLOSUM62", "PAM250", "BLOSUM45")) {
                rows <- pairwise_aln_rows(pairs$q, pairs$s, m)
                for (i in seq_along(rows)) {
                        f <- tempfile(fileext = ".fasta")
                        seqinr::write.fasta(list(pairs$q[i], pairs$s[i]), names = list("query", "subject"), file.out = f)
                        pairwise_aln(f, tool = "NW", seq_type = "protein", pairwise_aln_name = "perf_ref", quiet = TRUE,
                                     aa_sub_matrix = m)
                        written <- read_fasta_strings(file.path(orthologr_tmpdir(), "_alignment", "pairwise_aln",
                                                                "perf_ref_NW_AA.aln"))
                        expect_identical(rows[[i]], unname(written), info = paste(m, i))
                }
        }
})

test_that("pairwise_aln_rows() leaves to pairwise_aln() what it cannot handle", {
        good <- "MKTAYIAKQRLEGDLKAWLT"
        # no substitution matrix (the scoring of 0.4.3): nothing
        expect_identical(pairwise_aln_rows(c(good, good), c(good, good), NULL), list(NULL, NULL))
        # characters other than A - Z and '*', NA, an empty sequence: that pair only
        r <- pairwise_aln_rows(c(good, "mkta", NA, "", good), c(good, good, good, good, good), "BLOSUM62")
        expect_identical(vapply(r, is.null, TRUE), c(FALSE, TRUE, TRUE, TRUE, FALSE))
        # a stop codon symbol at the end is removed, one alone is not
        expect_identical(pairwise_aln_rows(paste0(good, "*"), good, "BLOSUM62")[[1]], c(good, good))
        # above options(orthologr.max_aln_cells): that pair goes to pairwise_aln() (COAL or an error there)
        withr_opts <- options(orthologr.max_aln_cells = 1000)
        on.exit(options(withr_opts), add = TRUE)
        long <- paste(rep(good, 2), collapse = "")
        r <- pairwise_aln_rows(c(good, long), c(good, long), "BLOSUM62")
        expect_identical(vapply(r, is.null, TRUE), c(FALSE, TRUE))
        # an option that is not valid: nothing (pairwise_aln() reports it for every pair)
        options(orthologr.max_aln_cells = "many")
        expect_identical(pairwise_aln_rows(good, good, "BLOSUM62"), list(NULL))
        options(withr_opts)
        # a matrix pwalign does not know: nothing
        expect_identical(pairwise_aln_rows(good, good, "NOMATRIX"), list(NULL))
        expect_identical(pairwise_aln_rows(character(0), character(0), "BLOSUM62"), list())
})

test_that("compute_dnds(): the batched alignments give the table of the per-pair code (Comeron, Li; COAL and PAL2NAL)", {
        skip_on_cran()
        tbl <- example_pair_table()
        # with an unrelated pair, a pair whose protein does not fit its CDS (fails) and a short pair
        tbl$subject_aa[5] <- tbl$subject_aa[6]
        tbl$subject_cds[5] <- tbl$subject_cds[6]
        tbl$query_aa[7] <- paste0(substr(tbl$query_aa[7], 1, 10), "WWWWWWWWWWWW", substr(tbl$query_aa[7], 11, 400))
        tools <- c(if (requireNamespace("coalr", quietly = TRUE)) "coal", if (perl_available) "pal2nal")
        skip_if(length(tools) == 0, "neither coalr nor perl is installed")
        for (tool in tools)
                for (m in c("Comeron", "Li")) {
                        r <- both_ways(tbl, codon_aln_tool = tool, dnds_est.method = m)
                        expect_identical(r$new, r$old, info = paste(tool, m))
                        expect_identical(nrow(attr(r$new, "failed_pairs")), 1L, info = paste(tool, m))
                }
})

test_that("compute_dnds(): KaKs_Calculator and COAL's KaKs methods give the table of the per-pair code", {
        skip_on_cran()
        tbl <- example_pair_table(8)
        if (kaks_available) {
                r <- both_ways(tbl, dnds_est.method = "YN")
                expect_identical(r$new, r$old)
        }
        if (requireNamespace("coalr", quietly = TRUE)) {
                r <- both_ways(tbl, codon_aln_tool = "coal", dnds_est.method = "YN", kaks_tool = "coal")
                expect_identical(r$new, r$old)
        }
})

test_that("compute_dnds(): the table does not depend on the chunk size or the number of processes", {
        skip_on_cran()
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        tbl <- example_pair_table()
        ref <- compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                            codon_aln_tool = "coal", dnds_est.method = "Comeron", quiet = TRUE)
        for (size in c(1L, 3L, 7L)) {
                local({
                        testthat::local_mocked_bindings(dnds_chunk_size = function() size, .package = "orthologr")
                        expect_identical(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                                      codon_aln_tool = "coal", dnds_est.method = "Comeron", quiet = TRUE),
                                         ref, info = size)
                })
        }
        expect_identical(dnds_chunks(1:10, 1), list(1:10))
        expect_identical(dnds_chunks(1:10, 3), list(1:4, 5:8, 9:10))
        expect_identical(dnds_chunks(integer(0), 2), list())
        expect_identical(lengths(dnds_chunks(1:250, 1)), c(100L, 100L, 50L))
        skip_if(parallel::detectCores() < 2 || .Platform$OS.type == "windows", "no forked processes")
        expect_identical(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                      codon_aln_tool = "coal", dnds_est.method = "Comeron", quiet = TRUE, comp_cores = 2),
                         ref)
})

test_that("compute_dnds(): progress messages after the same numbers of pairs as before", {
        skip_on_cran()
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        tbl <- example_pair_table()
        testthat::local_mocked_bindings(progress_step = function(n) 6L, .package = "orthologr")
        msgs <- character(0)
        withCallingHandlers(
                compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                             codon_aln_tool = "coal", dnds_est.method = "Comeron", quiet = FALSE),
                message = function(m) {
                        msgs <<- c(msgs, conditionMessage(m))
                        invokeRestart("muffleMessage")
                })
        done <- as.integer(sub("^\\D*(\\d+) of 20.*$", "\\1", grep(" of 20 ", msgs, value = TRUE)))
        expect_identical(done, c(6L, 12L, 18L))
})

test_that("the sequence checks with PCRE on bytes answer as the default regular expressions did", {
        set.seed(1)
        al <- c("A", "C", "G", "T", "a", "c", "g", "t", "N", "n", "R", "y", "X", "-", "*", "\n", " ", "ü", "(", "|")
        x <- c(vapply(1:5000, function(i) paste(sample(al, sample(0:12, 1), TRUE, prob = c(rep(10, 8), rep(1, 12))),
                                                collapse = ""), ""), NA, "", "ACGT\n", "\nACGT")
        expect_identical(is.dnaSequence(x), !grepl("[^ACGTacgt]", x))
        expect_identical(vapply(x, is.dnaSequence, TRUE, USE.NAMES = FALSE), is.dnaSequence(x))
        old_fuzzy <- grepl("^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]*$", x)
        expect_identical(grepl("^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]*\\z", x, perl = TRUE, useBytes = TRUE), old_fuzzy)
        dna <- x[!is.na(x) & nchar(x) %% 3 == 0 & nchar(x) > 0 & !grepl("[\n ()|ü]", x)]
        expect_identical(translate_cds(dna), vapply(dna, function(s) translate_cds(s), "", USE.NAMES = FALSE))
})
