# The vectorised code paths return what the per-gene / per-group code returned (audit 2026-11-11, problem 8).
# The old implementations are kept as *_do() / transl() and serve as the reference.

thal_1000 <- function() system.file("seqs/ortho_thal_cds_1000.fasta", package = "orthologr")
lyra_1000 <- function() system.file("seqs/ortho_lyra_cds_1000.fasta", package = "orthologr")

# ---------------------------------------------------------------- translation

test_that("translate_cds() returns what transl() (seqinr) returns for every sequence of the example files", {
        for (f in c(thal_cds(), lyra_cds(), thal_1000(), lyra_1000())) {
                cds <- suppressMessages(read.cds(f, format = "fasta", delete_corrupt_cds = FALSE))
                expect_identical(translate_cds(cds$seqs), vapply(cds$seqs, transl, "", USE.NAMES = FALSE), info = basename(f))
        }
})

test_that("translate_cds() agrees with transl() on alternative start codons, ambiguity codes and partial codons", {
        set.seed(2)
        r <- vapply(1:500, function(i) paste(sample(c("a", "c", "g", "t", "A", "C", "G", "T"),
                                                    sample(c(3, 6, 30, 31, 32, 300), 1), TRUE), collapse = ""), "")
        r <- c(r, "ttgaaa", "ctgaaa", "attaaa", "gtgaaa", "atgnnnaaa", "atga", "atgaartaa", "atgaa", "taatagtga")
        expect_identical(translate_cds(r), vapply(r, transl, "", USE.NAMES = FALSE))
        # TTG, CTG, ... at the start are not turned into M
        expect_identical(translate_cds(c("ttgaaa", "ctgaaa")), c("LK", "LK"))
        expect_identical(translate_cds(character(0)), character(0))
})

test_that("cds2aa() gives the per-gene result, also for IDs that occur twice", {
        dt <- suppressMessages(read.cds(thal_cds(), format = "fasta"))
        new <- cds2aa(data.table::copy(dt))
        old <- data.table::copy(dt)[ , aa := transl(seqs), by = geneids]
        expect_identical(new$aa, old$aa)
        # an ID that occurs twice was an error of the per-gene translation and still is
        dup <- data.table::data.table(geneids = c("a", "a", "b"), seqs = c("atgaaa", "atgccc", "atgggg"))
        expect_error(data.table::copy(dup)[ , aa := transl(seqs), by = geneids])
        expect_error(cds2aa(data.table::copy(dup)), "could not be translated")
})

test_that("write_protein_fasta() writes the file seqinr::write.fasta() wrote", {
        dt <- protein_table(thal_cds())
        f_old <- tempfile(); f_new <- tempfile()
        seqinr::write.fasta(sequences = as.list(dt$aa), names = dt$geneids, nbchar = 80, open = "w", file.out = f_old)
        write_protein_fasta(dt, f_new)
        expect_identical(unname(tools::md5sum(f_new)), unname(tools::md5sum(f_old)))
})

# ---------------------------------------------------------------- reading each file once

test_that("inside dNdS() / diamond_rec() a file is read once; the tables are the ones read before", {
        f <- thal_cds()
        direct_cds <- read.cds(f, format = "fasta", delete_corrupt_cds = FALSE)
        direct_pt  <- protein_table_read(f, seq_type = "cds", delete_corrupt_cds = FALSE)
        expect_false(seq_cache_active())
        end <- seq_cache_start()
        expect_true(seq_cache_active())
        inner_end <- seq_cache_start()          # a nested call does not end the outer one
        a <- read_cds_once(f, delete_corrupt_cds = FALSE)
        inner_end()
        expect_true(seq_cache_active())
        b <- read_cds_once(f, delete_corrupt_cds = FALSE)
        p1 <- protein_table(f, seq_type = "cds", delete_corrupt_cds = FALSE)
        p2 <- protein_table(f, seq_type = "cds", delete_corrupt_cds = FALSE)
        expect_equal(a, direct_cds)
        expect_equal(b, direct_cds)
        expect_equal(p1, direct_pt)
        # copies: changing one does not change the next one handed out
        data.table::setnames(p1, "aa", "changed")
        expect_true("aa" %in% names(p2))
        expect_true("aa" %in% names(protein_table(f, seq_type = "cds", delete_corrupt_cds = FALSE)))
        expect_equal(cds_protein_table_file(f, delete_corrupt_cds = FALSE),
                     cds_protein_table(direct_cds, delete_corrupt_cds = FALSE))
        end()
        expect_false(seq_cache_active())
})

test_that("a file that changes between two calls is read again", {
        d <- new_test_dir()
        f <- file.path(d, "cds.fasta")
        writeLines(c(">g1", "ATGAAATTTTAA"), f)
        end <- seq_cache_start()
        on.exit(end(), add = TRUE)
        expect_identical(read_cds_once(f)$geneids, "g1")
        writeLines(c(">g1", "ATGAAATTTTAA", ">g2", "ATGCCCGGGTAA"), f)
        expect_identical(read_cds_once(f)$geneids, c("g1", "g2"))
})

# ---------------------------------------------------------------- best hits

hit_table_example <- function() {
        set.seed(11)
        n <- 400
        tibble::tibble(
                query_id    = sample(paste0("q", 1:60), n, TRUE),
                subject_id  = sample(paste0("s", 1:80), n, TRUE),
                alig_length = sample(c(50L, 100L, 150L), n, TRUE),
                evalue      = sample(c(0, 1e-50, 1e-20, 1e-5), n, TRUE),
                locus       = sample(c("L3", "L1", "L2", NA), n, TRUE),
                bit_score   = stats::runif(n)
        )
}

test_that("best_hits() is identical to the per-query dplyr::do() result, with ties in e-value and length", {
        h <- hit_table_example()
        expect_identical(best_hits(h, "query_id"), best_hits_do(h, "query_id"))
        # a grouping column that is not the first one, with missing values
        expect_identical(best_hits(h, "locus"), best_hits_do(h, "locus"))
        # the form used until 0.4.3
        query_id <- NULL
        expect_identical(best_hits(h), dplyr::do(dplyr::group_by(h, query_id), filter_best_hits(.)))
        # all hits tied: the first row of each query wins
        h$evalue <- 0; h$alig_length <- 10L
        expect_identical(best_hits(h), best_hits_do(h))
        # missing e-values take the old path
        h$evalue[3] <- NA
        expect_identical(best_hits(h), best_hits_do(h))
})

test_that("diamond_best() and diamond_rec() on the 1,000-gene example are identical to the dplyr::do() result", {
        skip_on_cran()
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        fwd <- suppressMessages(diamond(thal_1000(), lyra_1000()))
        rev <- suppressMessages(diamond(lyra_1000(), thal_1000()))
        expect_gt(nrow(fwd), 1000)
        expect_identical(best_hits(fwd), best_hits_do(fwd))
        expect_identical(best_hits(rev), best_hits_do(rev))
        expect_identical(suppressMessages(diamond_best(thal_1000(), lyra_1000())), best_hits_do(fwd))
        # reciprocal best hits as they were computed before
        a <- best_hits_do(fwd); b <- best_hits_do(rev)
        colnames(b)[1:2] <- c("subject_id", "query_id")
        expect_identical(suppressMessages(diamond_rec(thal_1000(), lyra_1000())),
                         dplyr::semi_join(a, b, by = c("query_id", "subject_id")))
})

test_that("blast_best() is identical to the dplyr::do() result", {
        skip_if_not(blast_available, "BLAST+ is not installed")
        pars <- "-best_hit_score_edge 0.05 -best_hit_overhang 0.25"
        h <- suppressMessages(blast(thal_cds(), lyra_cds(), blast_params = pars))
        expect_identical(suppressMessages(blast_best(thal_cds(), lyra_cds())), best_hits_do(h))
})

# ---------------------------------------------------------------- annotation helpers

test_that("extract_features_all() is identical to one extract_features() call per gene locus (gtf and gff)", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        for (fmt in c("gtf", "gff")) {
                ann <- tibble::as_tibble(rtracklayer::import(x$query[[fmt]]))
                col <- if (fmt == "gtf") "gene_id" else "Name"
                expect_identical(extract_features_all(ann, fmt, col), suppressWarnings(extract_features_do(ann, fmt, col)), info = fmt)
                set.seed(4)
                shuffled <- ann[sample(nrow(ann)), ]
                expect_identical(extract_features_all(shuffled, fmt, col), suppressWarnings(extract_features_do(shuffled, fmt, col)), info = fmt)
                # 'type' as character: the per-locus path
                chr <- ann; chr$type <- as.character(chr$type)
                expect_identical(extract_features_all(chr, fmt, col), suppressWarnings(extract_features_do(chr, fmt, col)), info = fmt)
        }
        expect_error(extract_features_all(ann, "bed", "Name"), "format = 'gtf' or format = 'gff'")
})

test_that("longest_per_group() is identical to the per-gene dplyr::do() result (first of equally long isoforms)", {
        p <- tibble::tibble(transcript_id = paste0("t", 1:12),
                            width = c(5L, 9L, 9L, 3L, 7L, 7L, 1L, 2L, 8L, 8L, 4L, 6L),
                            gene_id = rep(c("g2", "g1", NA, "g1"), 3))
        expect_identical(longest_per_group(p, "gene_id"), longest_per_group_do(p, "gene_id"))
        expect_identical(longest_per_group(p[0, ], "gene_id"), longest_per_group_do(p[0, ], "gene_id"))
})

test_that("retrieve_longest_isoforms() keeps the longest isoform of every gene (gtf and gff)", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        for (fmt in c("gtf", "gff")) {
                out <- tempfile(fileext = ".faa")
                suppressMessages(retrieve_longest_isoforms(x$query$proteome, x$query[[fmt]], out, annotation_format = fmt))
                ids <- sub(" .*", "", names(Biostrings::readAAStringSet(out)))
                expect_length(ids, 20)
                expect_setequal(ids, x$query$ids)       # the '_iso2' isoforms are the shorter ones
        }
})

# ---------------------------------------------------------------- deepclust

test_that("deepclust() with several protein files clusters the table earlier versions built (sorted by ID)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        files <- c(thal_aa(), lyra_aa())
        d <- new_test_dir()
        # as before: sequence by sequence with seqinr into one file, which is read again
        old_dt <- protein_table(suppressMessages(merge_fasta_files(files, seq_type = "protein", dir = d)), seq_type = "protein")
        new_dt <- suppressMessages(deepclust_protein_table(files, seq_type = "protein", dir = d))
        expect_equal(new_dt, old_dt)
        expect_identical(new_dt$geneids, sort(new_dt$geneids, method = "radix"))
        res <- suppressMessages(deepclust(files))
        expect_equal(nrow(res), 40)
        expect_setequal(res$member_id, new_dt$geneids)
})

test_that("deepclust() reads IDs with '|' and ':' (the cluster table is tab-separated)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        files <- system.file(file.path("seqs", c("example_athal_aa.faa", "example_alyra_aa.faa", "example_brapa_aa.faa",
                                                 "example_thalo_aa.faa")), package = "orthologr")
        res <- suppressMessages(deepclust(files))
        expect_identical(names(res), c("representative_id", "member_id"))
        expect_equal(nrow(res), 80)
        expect_true(any(grepl("|PACid:", res$member_id, fixed = TRUE)))
})

test_that("deepclust(direct_input = TRUE) hands the files to DIAMOND unparsed: same IDs, every sequence clustered", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        files <- c(thal_aa(), lyra_aa())
        def <- suppressMessages(deepclust(files))
        dir <- suppressMessages(deepclust(files, direct_input = TRUE))
        expect_setequal(dir$member_id, def$member_id)
        expect_equal(nrow(dir), nrow(def))
        expect_true(all(dir$representative_id %in% dir$member_id))
        one <- suppressMessages(deepclust(thal_aa(), direct_input = TRUE))
        expect_equal(nrow(one), 20)
        expect_error(suppressMessages(deepclust(thal_cds(), seq_type = "cds", direct_input = TRUE)), "needs protein sequences")
        # gzip-compressed input
        d <- new_test_dir()
        gz <- file.path(d, "thal.faa.gz")
        con <- gzfile(gz, "w"); writeLines(readLines(thal_aa()), con); close(con)
        expect_setequal(suppressMessages(deepclust(gz, direct_input = TRUE))$member_id, one$member_id)
})

# ---------------------------------------------------------------- parallel backend

test_that("compute_dnds() gives the same table with 1 core, forked processes and PSOCK processes", {
        skip_on_cran()
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        skip_if(parallel::detectCores() < 2, "one core only")
        one <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        old <- options(orthologr.cluster_type = "PSOCK")
        on.exit(options(old), add = TRUE)
        psock <- suppressMessages(dNdS(thal_cds(), lyra_cds(), comp_cores = 2, print_citation = FALSE))
        options(old)
        expect_equal(as.data.frame(psock), as.data.frame(one), ignore_attr = TRUE)
        expect_equal(nrow(attr(psock, "failed_pairs")), 0)
        skip_on_os("windows")
        fork <- suppressMessages(dNdS(thal_cds(), lyra_cds(), comp_cores = 2, print_citation = FALSE))
        expect_equal(as.data.frame(fork), as.data.frame(one), ignore_attr = TRUE)
})

test_that("dnds_cluster() rejects an unknown cluster type and the package no longer needs foreach / doParallel", {
        old <- options(orthologr.cluster_type = "MPI")
        on.exit(options(old), add = TRUE)
        expect_error(dnds_cluster(2), "FORK.*PSOCK")
        desc <- utils::packageDescription("orthologr")$Imports
        expect_false(grepl("foreach|doParallel", desc))
})
