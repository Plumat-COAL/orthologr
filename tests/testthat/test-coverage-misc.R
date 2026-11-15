# Coverage of small functions no other test calls: deprecated wrappers, translate_cds_to_protein_all() and
# plot_diverse_homology_thresholds_core_orthologs() (with the
# helpers testCoreOrthoParamsGeneLocus() / testCoreOrthoParamsSpliceVariant()). No external program is needed.

test_that("the deprecated lnc wrappers call their successors", {
        tbl <- tibble::tibble(query_id = c("g1", "g1", "g2"), subject_species = c("s1", "s2", "s1"),
                              q_len = 100, alig_length = 100L)
        expect_warning(r <- suppressMessages(lnc_map_core_orthologs(tbl, c("s1", "s2"))), "retrieve_core_orthologs")
        expect_identical(r, suppressMessages(retrieve_core_orthologs(tbl, c("s1", "s2"))))
        expect_identical(r$query_id, c("g1", "g1"))
})

test_that("translate_cds_to_protein_all() translates every file of a folder but documentation", {
        input <- new_test_dir("cov_cds_in_")
        writeLines(c(">g1", "ATGAAATTTTAA", ">g2", "ATGGGGCCC"), file.path(input, "sp1.fa"))
        writeLines(c(">h1", "ATGTGGTAG"), file.path(input, "sp2.fa"))
        writeLines("not a fasta file", file.path(input, "documentation.txt"))
        output <- file.path(new_test_dir("cov_aa_out_"), "proteomes")
        msgs <- testthat::capture_messages(translate_cds_to_protein_all(input, output))
        expect_true(any(grepl("Creating new output folder", msgs)))
        expect_setequal(list.files(output), c("sp1.fa", "sp2.fa"))
        aa <- Biostrings::readAAStringSet(file.path(output, "sp1.fa"))
        expect_identical(as.character(aa), c(g1 = "MKF*", g2 = "MGP"))
        # an existing output folder is used as it is
        msgs <- testthat::capture_messages(translate_cds_to_protein_all(input, output))
        expect_false(any(grepl("Creating new output folder", msgs)))
        expect_true(any(grepl("successfully translated", msgs)))

        expect_error(translate_cds_to_protein_all(file.path(tempdir(), "no-input-folder"), output), "seems not to exist")
        only_doc <- new_test_dir("cov_doc_")
        writeLines("x", file.path(only_doc, "documentation.md"))
        expect_error(suppressMessages(translate_cds_to_protein_all(only_doc, output)), "No files were found")
})

test_that("plot_diverse_homology_thresholds_core_orthologs() counts core orthologs per parameter set", {
        skip_if_not_installed("ggplot2")
        skip_if_not_installed("scales")
        # g1, g2: hits in both species that pass every threshold; g3: only in s1 (never a core ortholog)
        tbl <- tibble::tibble(
                query_id = c("g1.1", "g1.1", "g2.1", "g2.1", "g3.1"),
                query_gene_locus_id = c("g1", "g1", "g2", "g2", "g3"),
                subject_species = c("s1", "s2", "s1", "s2", "s1"),
                q_len = 100, alig_length = 100L, qcovhsp = 95, perc_identity = 95)
        expect_error(plot_diverse_homology_thresholds_core_orthologs(tbl, c("s1", "s2"), type = "isoform"),
                     "valid 'type' argument")
        expect_identical(testCoreOrthoParamsGeneLocus(tbl, "p", c("s1", "s2")) |> suppressMessages(),
                         tibble::tibble(parameter = "p", n_core_genes = 2L))
        expect_identical(testCoreOrthoParamsSpliceVariant(tbl, "p", c("s1", "s2"), 2) |> suppressMessages(),
                         tibble::tibble(parameter = "p", n_core_genes = 2L))
        p1 <- suppressWarnings(suppressMessages(plot_diverse_homology_thresholds_core_orthologs(tbl, c("s1", "s2"),
                                                                                              type = "gene_locus")))
        expect_s3_class(p1, "ggplot")
        expect_identical(as.integer(p1$data$n_core_genes), rep(2L, 7))
        expect_identical(levels(p1$data$parameter)[1], "qcovhsp_50_perc_identity_30")
        p2 <- suppressWarnings(suppressMessages(plot_diverse_homology_thresholds_core_orthologs(tbl, c("s1", "s2"),
                                                                                              title = "core")))
        expect_s3_class(p2, "ggplot")
        expect_length(p2$layers, 6)
        expect_identical(p2$labels$title, "core")
})
