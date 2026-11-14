# Coverage of small functions no other test calls: path helpers, deprecated wrappers, readers of ortholog tables,
# translate_cds_to_protein_all(), check_ram() and plot_diverse_homology_thresholds_core_orthologs() (with the
# helpers testCoreOrthoParamsGeneLocus() / testCoreOrthoParamsSpliceVariant()). No external program is needed.

test_that("get_filename(), set_path() and ws_wrap()", {
        expect_identical(get_filename(file.path("a", "b", "genome.v2.fasta")), "genome.v2")
        expect_identical(get_filename("x.fa"), "x")

        sep <- .Platform$file.sep
        wd <- getwd()
        expect_identical(set_path(file.path("some", "dir", "file.txt")), paste0(wd, "file.txt"))
        expect_identical(set_path("file.txt", add.folder = "out"), paste0(wd, sep, "out", "file.txt"))
        # a folder of the working directory with a space is quoted
        d <- file.path(new_test_dir("cov path "), "inner")
        dir.create(d)
        old <- setwd(d)
        on.exit(setwd(old), add = TRUE)
        p <- set_path("f.txt")
        expect_match(p, "'cov path ", fixed = TRUE)
        expect_true(endsWith(p, "innerf.txt"))

        expect_identical(ws_wrap(paste("a", "b c", "d", sep = sep)), paste("a", "'b c'", "d", sep = sep))
        expect_identical(ws_wrap(c(paste("x", "y", sep = sep), "no space")), c(paste("x", "y", sep = sep), "'no space'"))
})

test_that("the deprecated lnc wrappers call their successors", {
        tbl <- tibble::tibble(query_id = c("g1", "g1", "g2"), subject_species = c("s1", "s2", "s1"),
                              q_len = 100, alig_length = 100L)
        expect_warning(r <- filter_core_set_lnc(tbl, c("s1", "s2")), "filter_core_set")
        expect_identical(r, filter_core_set(tbl, c("s1", "s2")))
        expect_warning(r <- suppressMessages(lnc_map_core_orthologs(tbl, c("s1", "s2"))), "retrieve_core_orthologs")
        expect_identical(r, suppressMessages(retrieve_core_orthologs(tbl, c("s1", "s2"))))
        expect_identical(r$query_id, c("g1", "g1"))
})

test_that("OF2Orthologs() and read.ortholog.table() read their files", {
        expect_error(OF2Orthologs(file.path(tempdir(), "no-orthologs.tsv")), "valid path to the orthologs file")
        f <- tempfile(fileext = ".tsv")
        writeLines(c("Orthogroup\tA\tB", "OG0000001\ta1\tb1, b2", "OG0000002\ta2\tb3"), f)
        o <- OF2Orthologs(f)
        expect_identical(names(o), c("Orthogroup", "A", "B"))
        expect_identical(o$B, c("b1, b2", "b3"))

        cols <- c("query_species", "subject_species", "query_id", "query_gene_locus_id", "subject_id",
                  "subject_gene_locus_id", "dN", "dS", "dNdS", "evalue", "bit_score", "perc_identity", "alig_length",
                  "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end")
        f <- tempfile(fileext = ".csv")
        writeLines(c(paste(cols, collapse = ";"),
                     "Ath;Aly;0012;AT1;1e5;AL1;0.1;0.5;0.2;1e-50;300;95.5;100;4;0;1;100;1;100"), f)
        x <- read.ortholog.table(f)
        expect_identical(names(x), cols)
        # IDs stay text
        expect_identical(x$query_id, "0012")
        expect_identical(x$subject_id, "1e5")
        expect_identical(x$alig_length, 100L)
        expect_equal(x$dNdS, 0.2)
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

test_that("check_ram() returns the memory in GB, or NA", {
        skip_if_not_installed("benchmarkme")
        testthat::local_mocked_bindings(get_ram = function() structure(16e9, class = "ram"), .package = "benchmarkme")
        expect_identical(check_ram(), c(ram_in_gb = 16))
        testthat::local_mocked_bindings(get_ram = function() NA, .package = "benchmarkme")
        expect_identical(check_ram(), NA)
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
