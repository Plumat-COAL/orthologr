# orthofinder2() is soft-deprecated: PATH lookup, checked run, import of both OrthoFinder layouts

test_that("orthofinder2() is soft-deprecated, looks for orthofinder on the PATH and gives clear errors", {
        d <- new_test_dir()
        expect_warning(expect_error(orthofinder2(file.path(d, "nowhere")), "valid path to the proteome folder"),
                       "soft-deprecated")
        expect_error(suppressWarnings(orthofinder2(d, import_type = "orthogroups_pairwise")), "not implemented")
        expect_error(suppressWarnings(orthofinder2(d, import_type = "something")), "not supported")
        expect_error(suppressWarnings(orthofinder2(d, use_existing_output = TRUE)), "results folder .* does not exist")
        expect_error(suppressWarnings(orthofinder2(d, of_path = file.path(d, "no_bin"))),
                     "don't have OrthoFinder installed")
        expect_false(any(grepl("/opt/miniconda3/bin/orthofinder\"", deparse(body(orthofinder2)), fixed = TRUE)))
})

test_that("orthofinder2(use_existing_output = TRUE) imports single-copy orthogroups of both layouts", {
        skip_if_not_installed("readr")
        for (tbl in c("Orthogroups.tsv", "MCL_Orthogroups.tsv")) {
                d <- new_test_dir()
                og <- file.path(d, "OrthoFinder", paste0("Results_", basename(d)), "Orthogroups")
                dir.create(og, recursive = TRUE)
                writeLines(c("Orthogroup\tA\tB", "OG0000000\ta1, a2\tb1", "OG0000001\ta3\tb3", "OG0000002\ta4\tb4"),
                           file.path(og, tbl))
                expect_error(suppressWarnings(orthofinder2(d, use_existing_output = TRUE, import_type = "orthogroups_core")),
                             "only written by a full OrthoFinder run")
                writeLines(c("OG0000001", "OG0000002"), file.path(og, "Orthogroups_SingleCopyOrthologues.txt"))
                res <- suppressWarnings(orthofinder2(d, use_existing_output = TRUE, import_type = "orthogroups_core"))
                expect_identical(res$orthogroup, c("OG0000001", "OG0000002"))
                expect_identical(suppressWarnings(orthofinder2(d, use_existing_output = TRUE)),
                                 file.path(normalizePath(d), "OrthoFinder", paste0("Results_", basename(d))))
        }
})
