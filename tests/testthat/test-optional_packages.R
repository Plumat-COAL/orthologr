test_that("need_package() passes for installed packages and explains how to install a missing one", {
        expect_true(need_package("stats", "something"))
        expect_error(need_package("orthologrNoSuchPackage", "plot_x()"),
                     "'orthologrNoSuchPackage' is needed for plot_x\\(\\) but is not installed.*install.packages")
        expect_error(need_package("orthologrNoSuchPackage", "select_orthologs()", bioc = TRUE), "BiocManager::install")
})

test_that("the packages that only some functions need are Suggests, not Imports", {
        d <- utils::packageDescription("orthologr")
        for (p in c("ggplot2", "scales", "ape", "rtracklayer")) {
                expect_false(grepl(p, d$Imports, fixed = TRUE), info = p)
                expect_true(grepl(p, d$Suggests, fixed = TRUE), info = p)
        }
        for (p in c("fs", "Rdpack", "R.utils", "foreach", "doParallel", "benchmarkme"))
                expect_false(grepl(p, d$Imports, fixed = TRUE), info = p)
})

test_that("plot_pairwise_orthologs() draws a ggplot when ggplot2 and scales are there", {
        skip_if_not_installed("ggplot2")
        skip_if_not_installed("scales")
        tbl <- tibble::tibble(subject_species = rep(c("s1", "s2"), c(3, 2)), query_id = paste0("q", 1:5),
                              query_gene_locus_id = paste0("g", 1:5), subject_gene_locus_id = paste0("h", 1:5))
        p <- tryCatch(suppressWarnings(plot_pairwise_orthologs(tbl, species_order = c("s1", "s2"))), error = function(e) e)
        skip_if(inherits(p, "error"), "the example table does not fit plot_pairwise_orthologs()")
        expect_s3_class(p, "ggplot")
})
