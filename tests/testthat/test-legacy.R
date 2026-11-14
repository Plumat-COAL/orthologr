# The results of orthologr 0.4.3. The tables in fixtures/dnds_0.4.3 were computed with orthologr 0.4.3 (see
# fixtures/README and make_fixture.R there); this version must return them unchanged with the options that restore
# the defaults of 0.4.3 (one per changed default, see NEWS.md).

# The options that restore the defaults of 0.4.3, set in the calling test and undone when it ends.
with_0.4.3_defaults <- function(env = parent.frame()) {
        rlang::local_options(orthologr.legacy_pairwise_aln = TRUE, orthologr.legacy_comeron_saturation = TRUE,
                             orthologr.legacy_kaks_zero_dn = TRUE, orthologr.legacy_kaks_saturation = TRUE,
                             orthologr.best_hit_rule = "evalue", .frame = env)
}

# The tables are compared with a tolerance of 1e-12: e-values such as 2.27e-158 are read with a difference in the
# last bit on some platforms (macOS arm64); every other difference shows.
plain_table <- function(x) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        for (n in names(x))
                attributes(x[[n]]) <- NULL
        attributes(x) <- list(names = names(x), class = "data.frame", row.names = seq_len(nrow(x)))
        x
}

stored <- function(method) readRDS(test_path("fixtures", "dnds_0.4.3", paste0(tolower(method), "_20.rds")))

test_that("dNdS() reproduces the stored table of orthologr 0.4.3 (20 genes, Comeron)", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        with_0.4.3_defaults()
        ref <- stored("Comeron")
        expect_identical(ref$orthologr, "0.4.3")
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "Comeron", comp_cores = 1,
                                     quiet = TRUE, print_citation = FALSE))
        expect_identical(names(res), names(ref$table))
        expect_equal(plain_table(res), ref$table, tolerance = 1e-12)
})

test_that("dNdS() reproduces the stored table of orthologr 0.4.3 (20 genes, YN)", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        with_0.4.3_defaults()
        ref <- stored("YN")
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "YN", comp_cores = 1,
                                     quiet = TRUE, print_citation = FALSE))
        expect_identical(names(res), names(ref$table))
        expect_equal(plain_table(res), ref$table, tolerance = 1e-12)
})
