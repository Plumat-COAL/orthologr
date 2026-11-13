# Best hits are ranked by bit score by default (owner's decision for 0.5.0); best_hit_rule = "evalue" and
# options(orthologr.best_hit_rule = "evalue") give the rule of orthologr <= 0.4.3 exactly. (Second audit: N01.)

thal_1000 <- function() system.file("seqs/ortho_thal_cds_1000.fasta", package = "orthologr")
lyra_1000 <- function() system.file("seqs/ortho_lyra_cds_1000.fasta", package = "orthologr")

hit <- function(q, s, bits, e, len, pid) {
        tibble::tibble(query_id = q, subject_id = s, perc_identity = pid, alig_length = as.integer(len), evalue = e, bit_score = bits)
}

# every tie-break level once (the first row of each block is the hit that must win under "bitscore")
tie_fixture <- function() {
        dplyr::bind_rows(
                # q1: both e-values are 0 (strong hits of a long protein); the paralog has the longer alignment
                hit("q1", "ortholog", 900, 0, 500, 80), hit("q1", "paralog", 700, 0, 800, 55),
                # q2: same bit score: the smaller e-value
                hit("q2", "s_e_small", 300, 1e-80, 200, 60), hit("q2", "a_e_large", 300, 1e-70, 400, 90),
                # q3: same bit score and e-value: the longer alignment
                hit("q3", "s_long", 300, 1e-80, 400, 60), hit("q3", "a_short", 300, 1e-80, 200, 90),
                # q4: ... and same length: the higher percent identity
                hit("q4", "s_ident", 300, 1e-80, 400, 90), hit("q4", "a_less", 300, 1e-80, 400, 60),
                # q5: everything equal: the alphabetically first subject
                hit("q5", "aaa", 300, 1e-80, 400, 90), hit("q5", "bbb", 300, 1e-80, 400, 90),
                # q6: one hit
                hit("q6", "only", 50, 1e-6, 40, 30),
                # q7: the higher bit score wins although its e-value is not the smallest (both rules disagree)
                hit("q7", "high_bits", 500, 1e-100, 300, 70), hit("q7", "low_bits", 450, 0, 300, 70)
        )
}

test_that("the bit-score rule: bit score, then e-value, then alignment length, then identity, then subject_id ...", {
        h <- tie_fixture()
        res <- best_hits(h, "query_id", rule = "bitscore")
        expect_s3_class(res, "grouped_df")
        expect_identical(res$query_id, paste0("q", 1:7))
        expect_identical(res$subject_id, c("ortholog", "s_e_small", "s_long", "s_ident", "aaa", "only", "high_bits"))
        # deterministic: the order of the rows of the hit table does not matter
        set.seed(5)
        for (k in 1:5) {
                shuffled <- h[sample(nrow(h)), ]
                expect_identical(as.data.frame(best_hits(shuffled, "query_id", rule = "bitscore")), as.data.frame(res))
        }
        expect_identical(as.data.frame(best_hits(h[rev(seq_len(nrow(h))), ], rule = "bitscore")), as.data.frame(res))
})

test_that("the e-value rule is the rule of orthologr <= 0.4.3: smallest e-value, longest alignment, first row ...", {
        h <- tie_fixture()
        old <- best_hits(h, "query_id", rule = "evalue")
        expect_identical(old, best_hits_do(h, "query_id"))
        expect_identical(old$subject_id, c("paralog", "s_e_small", "s_long", "s_ident", "aaa", "only", "low_bits"))
        # ... which depends on the order of the rows when everything ties
        expect_identical(best_hits(h[rev(seq_len(nrow(h))), ], rule = "evalue")$subject_id[5], "bbb")
})

test_that("the default rule is 'bitscore'; the option and the argument choose the old rule; other values are errors ...", {
        h <- tie_fixture()
        old_opt <- options(orthologr.best_hit_rule = NULL)
        on.exit(options(old_opt), add = TRUE)
        expect_identical(orthologr_best_hit_rule(), "bitscore")
        expect_identical(best_hits(h)$subject_id[1], "ortholog")
        options(orthologr.best_hit_rule = "evalue")
        expect_identical(orthologr_best_hit_rule(), "evalue")
        expect_identical(best_hits(h)$subject_id[1], "paralog")
        expect_identical(best_hits(h, rule = "bitscore")$subject_id[1], "ortholog")       # the argument wins
        options(orthologr.best_hit_rule = "score")
        expect_error(best_hits(h), "'score' is not supported")
        options(orthologr.best_hit_rule = NULL)
        expect_error(orthologr_best_hit_rule(NA), "must be")
        expect_error(best_hits(h[ , c("query_id", "subject_id", "evalue", "alig_length")], rule = "bitscore"), "needs the columns")
        expect_equal(nrow(best_hits(h[0, ], rule = "bitscore")), 0)
        # the exported functions have the argument, at the end of their signatures, with the option as default
        for (f in list(diamond_best, diamond_rec, blast_best, blast_rec, orthologs, dNdS, divergence_stratigraphy, select_orthologs)) {
                expect_identical(utils::tail(names(formals(f)), 1), "best_hit_rule")
        }
        options(orthologr.best_hit_rule = "evalue")
        expect_identical(eval(formals(dNdS)$best_hit_rule), "evalue")
})

test_that("1,000-gene example: best_hit_rule = 'evalue' reproduces the 0.4.3 tables exactly; 'bitscore' differs only where it should ...", {
        skip_on_cran()
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        fwd <- suppressMessages(diamond(thal_1000(), lyra_1000()))
        rev <- suppressMessages(diamond(lyra_1000(), thal_1000()))
        
        # old rule: identical to the reference implementation (one filter_best_hits() call per query)
        old_bh <- suppressMessages(diamond_best(thal_1000(), lyra_1000(), best_hit_rule = "evalue"))
        expect_identical(old_bh, best_hits_do(fwd))
        a <- best_hits_do(fwd); b <- best_hits_do(rev)
        colnames(b)[1:2] <- c("subject_id", "query_id")
        old_rbh <- suppressMessages(diamond_rec(thal_1000(), lyra_1000(), best_hit_rule = "evalue"))
        expect_identical(old_rbh, dplyr::semi_join(a, b, by = c("query_id", "subject_id")))
        # ... also through the option
        opt <- options(orthologr.best_hit_rule = "evalue")
        expect_identical(suppressMessages(diamond_rec(thal_1000(), lyra_1000())), old_rbh)
        expect_equal(as.data.frame(suppressMessages(orthologs(thal_1000(), lyra_1000(), seq_type = "cds", delete_corrupt_cds = TRUE))),
                     as.data.frame(old_rbh), ignore_attr = TRUE)
        options(opt)
        
        # new rule: every best hit has the highest bit score of its query
        new_bh <- suppressMessages(diamond_best(thal_1000(), lyra_1000()))
        expect_identical(new_bh$query_id, old_bh$query_id)
        top <- tapply(fwd$bit_score, fwd$query_id, max)
        expect_equal(new_bh$bit_score, as.vector(top[new_bh$query_id]))
        expect_true(all(old_bh$bit_score <= new_bh$bit_score))
        # where the two rules choose another subject, the old choice has a lower bit score or the tie was broken
        # by alignment length / row order
        differ <- new_bh$subject_id != old_bh$subject_id
        expect_true(all(old_bh$bit_score[differ] <= new_bh$bit_score[differ]))
        
        # reciprocal best hits under the new rule are reciprocal: the same pairs from both directions
        new_rbh <- suppressMessages(diamond_rec(thal_1000(), lyra_1000()))
        back <- suppressMessages(diamond_rec(lyra_1000(), thal_1000()))
        expect_setequal(paste(new_rbh$query_id, new_rbh$subject_id), paste(back$subject_id, back$query_id))
        expect_false(anyDuplicated(new_rbh$query_id) > 0)
        expect_false(anyDuplicated(new_rbh$subject_id) > 0)
})

test_that("dNdS() and divergence_stratigraphy() pass best_hit_rule on ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        msg <- testthat::capture_messages(new <- dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_match(paste(msg, collapse = ""), "best_hit_rule 'bitscore'", fixed = TRUE)
        msg <- testthat::capture_messages(old <- dNdS(thal_cds(), lyra_cds(), print_citation = FALSE, best_hit_rule = "evalue"))
        expect_match(paste(msg, collapse = ""), "best_hit_rule 'evalue'", fixed = TRUE)
        rbh_old <- suppressMessages(diamond_rec(thal_cds(), lyra_cds(), best_hit_rule = "evalue"))
        rbh_new <- suppressMessages(diamond_rec(thal_cds(), lyra_cds()))
        expect_setequal(paste(old$query_id, old$subject_id), paste(rbh_old$query_id, rbh_old$subject_id))
        expect_setequal(paste(new$query_id, new$subject_id), paste(rbh_new$query_id, rbh_new$subject_id))
        expect_error(dNdS(thal_cds(), lyra_cds(), best_hit_rule = "best"), "'best' is not supported")
        expect_error(divergence_stratigraphy(thal_cds(), lyra_cds(), best_hit_rule = "best"), "'best' is not supported")
        bh <- suppressMessages(dNdS(thal_cds(), lyra_cds(), ortho_detection = "BH", print_citation = FALSE))
        expect_equal(nrow(bh), 20)
})

test_that("blast_best() / blast_rec() use the same rule ...", {
        skip_if_not(blast_available, "BLAST+ is not installed")
        pars <- "-best_hit_score_edge 0.05 -best_hit_overhang 0.25"
        h <- suppressMessages(blast(thal_cds(), lyra_cds(), blast_params = pars))
        expect_identical(suppressMessages(blast_best(thal_cds(), lyra_cds(), best_hit_rule = "evalue")), best_hits_do(h))
        new <- suppressMessages(blast_best(thal_cds(), lyra_cds()))
        expect_equal(new$bit_score, as.vector(tapply(h$bit_score, h$query_id, max)[new$query_id]))
        expect_equal(nrow(suppressMessages(blast_rec(thal_cds(), lyra_cds()))), 20)
})

test_that("select_orthologs() chooses the splice variant of a locus by the same rule ...", {
        x <- make_test_ortholog_input()
        old <- suppressMessages(select_orthologs(x$dnds, x$query$gtf, x$subject$gtf, format = c("gtf", "gtf"), best_hit_rule = "evalue"))
        new <- suppressMessages(select_orthologs(x$dnds, x$query$gtf, x$subject$gtf, format = c("gtf", "gtf")))
        expect_false(anyDuplicated(new$query_gene_locus_id) > 0)
        # per locus: the variant with the highest bit score / the smallest e-value
        d <- suppressMessages(select_orthologs(x$dnds, x$query$gtf, x$subject$gtf, collapse_by = "splice_variant", format = c("gtf", "gtf")))
        for (locus in new$query_gene_locus_id[1:5])
                expect_equal(new$bit_score[new$query_gene_locus_id == locus], max(d$bit_score[d$query_gene_locus_id == locus]))
        for (locus in old$query_gene_locus_id[1:5])
                expect_equal(old$evalue[old$query_gene_locus_id == locus], min(d$evalue[d$query_gene_locus_id == locus]))
})
