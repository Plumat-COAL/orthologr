# Coverage of orthologs()' dispatch to the BLAST searches and its "pairwise only" checks (the searches are replaced
# by stand-ins, so neither BLAST nor DIAMOND is needed), and of the argument checks and error paths of
# promotor_divergence_of_orthologous_genes() (they all come before any distance is computed).

cov_hits <- function() data.table::data.table(query_id = "q1", subject_id = "s1", evalue = 1e-50, bit_score = 300)

test_that("orthologs() hands BH and RBH to the BLAST searches, and cleans their folder on request", {
        seen <- new.env()
        testthat::local_mocked_bindings(
                blast_best = function(...) { seen$bh <- list(...); cov_hits() },
                blast_rec  = function(...) { seen$rbh <- list(...); cov_hits() },
                .package = "orthologr")
        q <- thal_aa()
        s <- lyra_aa()
        res <- suppressMessages(suppressWarnings(orthologs(q, s, ortho_detection = "BH", eval = "1E-10",
                                                           add_params = "-seg no", clean_folders = TRUE)))
        expect_identical(res, cov_hits())
        expect_identical(seen$bh$eval, "1E-10")
        expect_identical(seen$bh$blast_params, "-seg no")
        expect_identical(seen$bh$seq_type, "protein")
        res <- suppressMessages(suppressWarnings(orthologs(q, s, ortho_detection = "RBH", clean_folders = TRUE)))
        expect_identical(res, cov_hits())
        expect_identical(seen$rbh$subject_file, s)
})

test_that("orthologs(): every method is only defined for one subject file", {
        q <- thal_aa()
        two <- c(lyra_aa(), thal_aa())
        for (m in c("DIAMOND_RBH", "DIAMOND_BH", "BH", "RBH"))
                expect_error(suppressMessages(suppressWarnings(orthologs(q, two, ortho_detection = m))),
                             "only defined for pairwise comparisons")
})

# ---------------------------------------------------------------------------------------------------------------

# A promotor folder and an ortholog table folder as promotor_divergence_of_orthologous_genes() reads them.
#   qry / sbj : names of the promotor sequences of the query (Qsp) and the subject (species `sbj_species`)
cov_promotor_input <- function(qry = paste0("G", 1:3, "_prom"), sbj = paste0("S", 1:3, "_prom"),
                               sbj_species = "Ssp") {
        d <- new_test_dir("cov_prom_")
        dir.create(file.path(d, "prom"))
        dir.create(file.path(d, "orth"))
        seqs <- function(ids) Biostrings::DNAStringSet(stats::setNames(rep("ACGTACGTACGTACGTACGT", length(ids)), ids))
        Biostrings::writeXStringSet(seqs(qry), file.path(d, "prom", "Qsp_promotors.fa"))
        Biostrings::writeXStringSet(seqs(sbj), file.path(d, "prom", paste0(sbj_species, "_promotors.fa")))
        ot <- data.frame(query_species = "Qsp", subject_species = "Ssp", query_id = paste0("G", 1:3, ".1"),
                         query_gene_locus_id = paste0("G", 1:3), subject_id = paste0("S", 1:3, ".1"),
                         subject_gene_locus_id = paste0("S", 1:3), dN = 0.01, dS = 0.1, dNdS = 0.1, evalue = 0,
                         bit_score = 100, perc_identity = 90, num_ident_matches = 90, alig_length = 100L,
                         mismatches = 1L, gap_openings = 0L, n_gaps = 0, pos_match = 90, ppos = 90, q_start = 1L,
                         q_end = 100L, q_len = 100, qcov = 100, qcovhsp = 100, s_start = 1L, s_end = 100L,
                         s_len = 100, score_raw = 100)
        utils::write.table(ot, file.path(d, "orth", "Qsp_vs_Ssp.csv"), sep = ";", row.names = FALSE, quote = FALSE)
        list(prom = file.path(d, "prom"), orth = file.path(d, "orth"), out = file.path(d, "pout"), dir = d)
}

cov_promotor <- function(x, ...) {
        suppressMessages(promotor_divergence_of_orthologous_genes(x$prom, x$orth, ...))
}

test_that("promotor_divergence_of_orthologous_genes() checks its folders", {
        x <- cov_promotor_input()
        expect_error(promotor_divergence_of_orthologous_genes(file.path(x$dir, "nope"), x$orth), "valid 'promotor_folder'")
        expect_error(promotor_divergence_of_orthologous_genes(x$prom, file.path(x$dir, "nope")),
                     "valid 'ortholog_tables_folder'")
})

test_that("promotor_divergence_of_orthologous_genes(): species of the two folders must match", {
        x <- cov_promotor_input(sbj_species = "Xsp")
        expect_error(cov_promotor(x, ortholog_promotor_seq_output = x$out), "do not match")
        # the output folder is created before (by default in the working directory)
        expect_true(dir.exists(x$out))
        old <- setwd(x$dir)
        on.exit(setwd(old), add = TRUE)
        expect_error(cov_promotor(x), "do not match")
        expect_true(dir.exists(file.path(x$dir, "ortholog_genes_promotor_seqs")))
})

test_that("promotor_divergence_of_orthologous_genes(): promotors that belong to no ortholog pair", {
        x <- cov_promotor_input(qry = c("Z1_prom", "Z2_prom"))
        expect_error(cov_promotor(x, ortholog_promotor_seq_output = x$out), "Qsp_promotors.fa")
        x <- cov_promotor_input(sbj = c("Z1_prom", "Z2_prom"))
        expect_error(cov_promotor(x, ortholog_promotor_seq_output = x$out), "Ssp_promotors.fa")
        # G1 has a promotor in the query only, S2 in the subject only
        x <- cov_promotor_input(qry = "G1_prom", sbj = "S2_prom")
        expect_error(cov_promotor(x, ortholog_promotor_seq_output = x$out), "has a promotor sequence in both species")
})
