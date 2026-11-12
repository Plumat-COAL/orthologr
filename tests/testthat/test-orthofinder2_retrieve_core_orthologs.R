test_that("orthofinder2_retrieve_core_orthologs() keeps the single-copy orthogroups", {
        skip_if_not_installed("readr")
        d <- new_test_dir()
        og <- file.path(d, "Orthogroups.tsv"); sc <- file.path(d, "Orthogroups_SingleCopyOrthologues.txt")
        writeLines(c("Orthogroup\tA\tB\tC", "OG0000000\ta1, a2\tb1\tc1", "OG0000001\ta3\tb3\tc3", "OG0000002\ta4\tb4\tc4",
                     "OG0000003\ta5\t\tc5"), og)
        writeLines(c("OG0000001", "OG0000002"), sc)
        res <- orthofinder2_retrieve_core_orthologs(og, sc)
        expect_identical(names(res), c("orthogroup", "A", "B", "C"))
        expect_identical(res$orthogroup, c("OG0000001", "OG0000002"))
        expect_identical(res$A, c("a3", "a4"))
        expect_error(orthofinder2_retrieve_core_orthologs(file.path(d, "nope.tsv"), sc), "Orthogroups")
        expect_error(orthofinder2_retrieve_core_orthologs(og, file.path(d, "nope.txt")), "SingleCopyOrthogroups")
})
