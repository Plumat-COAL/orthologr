test_that("extract_features() returns gene locus and transcript ids of one locus (gtf and gff)", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        gtf <- tibble::as_tibble(rtracklayer::import(x$query$gtf))
        one <- gtf[gtf$gene_id == "AT1G01010", ]
        res <- extract_features(one, format = "gtf")
        expect_identical(names(res), c("gene_locus_id", "query_id"))
        expect_identical(unique(res$gene_locus_id), "AT1G01010")
        expect_setequal(res$query_id, c("AT1G01010.1", "AT1G01010.1_iso2"))
        gff <- tibble::as_tibble(rtracklayer::import(x$query$gff))
        one <- gff[gff$Name == "AT1G01010", ]
        res <- extract_features(one, format = "gff")
        expect_identical(unique(res$gene_locus_id), "gene-AT1G01010")
        expect_setequal(res$query_id, c("AT1G01010.1", "AT1G01010.1_iso2"))
        expect_error(extract_features(one, format = "bed"), "format = 'gtf' or format = 'gff'")
})
