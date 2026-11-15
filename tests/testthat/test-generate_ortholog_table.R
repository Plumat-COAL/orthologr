test_that("generate_ortholog_tables() writes the gene-locus and the splice-variant table", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        dnds_file <- file.path(x$dir, "Athaliana_vs_Alyrata.csv")
        utils::write.table(x$dnds, dnds_file, sep = ";", quote = FALSE, row.names = FALSE)
        expect_error(generate_ortholog_tables(dnds_file, x$query$gtf, x$subject$gtf, output_folder = file.path(x$dir, "nope"),
                                              format = c("gtf", "gtf")), "output_folder")
        out <- file.path(x$dir, "tables")
        dir.create(out)
        suppressMessages(suppressWarnings(generate_ortholog_tables(dnds_file, x$query$gtf, x$subject$gtf, output_folder = out, format = c("gtf", "gtf"))))
        # one sub-folder per kind of table
        by_gene_file <- file.path(out, "orthologs_by_gene_locus_qry_athaliana", "orthologs_by_gene_Athaliana_vs_Alyrata.csv")
        by_sv_file   <- file.path(out, "orthologs_by_splice_variant_qry_athaliana", "orthologs_by_splice_variant_Athaliana_vs_Alyrata.csv")
        expect_true(file.exists(by_gene_file))
        expect_true(file.exists(by_sv_file))
        by_gene <- utils::read.csv(by_gene_file, sep = ";")
        expect_equal(nrow(by_gene), 20)
        expect_equal(nrow(utils::read.csv(by_sv_file, sep = ";")), 25)
})
