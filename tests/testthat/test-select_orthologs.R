test_that("select_orthologs() checks its arguments", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        expect_error(select_orthologs(x$dnds, x$query$gtf, x$subject$gtf, collapse_by = "gene", format = c("gtf", "gtf")), "collapse_by")
        expect_error(select_orthologs(x$dnds, x$query$gtf, x$subject$gtf, format = c("gtf", "bed")), "gtf.*gff")
        expect_error(select_orthologs(x$dnds, file.path(x$dir, "nope.gtf"), x$subject$gtf, format = c("gtf", "gtf")), "does not seem to exist")
})

for (fmt in c("gtf", "gff")) {
        test_that(paste0("select_orthologs() keeps one splice variant per gene locus (", fmt, ")"), {
                skip_if_not_installed("rtracklayer")
                x <- make_test_ortholog_input()
                res <- suppressMessages(select_orthologs(x$dnds, x$query[[fmt]], x$subject[[fmt]], collapse_by = "gene_locus", format = c(fmt, fmt)))
                # 25 isoform pairs of 20 gene loci
                expect_equal(nrow(res), 20)
                expect_false(anyDuplicated(res$query_gene_locus_id) > 0)
                expect_false(anyDuplicated(res$subject_gene_locus_id) > 0)
                expect_identical(names(res)[1:9], c("query_species", "subject_species", "query_id", "query_gene_locus_id",
                                                    "subject_id", "subject_gene_locus_id", "dN", "dS", "dNdS"))
                expect_identical(unique(res$query_species), "Athaliana")
                expect_identical(unique(res$subject_species), "Alyrata")
                # of the two isoform pairs of a locus the one with the smaller e-value is kept
                for (locus in res$query_gene_locus_id[1:5]) {
                        kept <- res[res$query_gene_locus_id == locus, ]
                        candidates <- x$dnds[sub("_iso2$", "", x$dnds$query_id) == sub("_iso2$", "", kept$query_id), ]
                        expect_equal(kept$evalue, min(candidates$evalue))
                }
                all_variants <- suppressMessages(select_orthologs(x$dnds, x$query[[fmt]], x$subject[[fmt]], collapse_by = "splice_variant", format = c(fmt, fmt)))
                expect_equal(nrow(all_variants), 25)
                expect_setequal(all_variants$query_id, x$dnds$query_id)
        })
}
