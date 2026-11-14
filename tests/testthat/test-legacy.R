# One switch for the results of orthologr 0.4.3: orthologr_legacy_options() / options(orthologr.legacy = TRUE).
# The tables in fixtures/dnds_0.4.3 were computed with orthologr 0.4.3 (see fixtures/README and make_fixture.R
# there); with the legacy settings this version must return them unchanged.

legacy_names <- c("orthologr.best_hit_rule", "orthologr.legacy_pairwise_aln", "orthologr.legacy_comeron_saturation",
                  "orthologr.legacy_kaks_zero_dn", "orthologr.legacy_kaks_saturation", "orthologr.dS_zero",
                  "orthologr.status_column")

# No orthologr option is set inside the calling test; what was set before, and what the test sets, is undone when
# the test ends.
without_orthologr_options <- function(env = parent.frame()) {
        nm <- names(options())
        nm <- unique(c(nm[startsWith(nm, "orthologr.")], legacy_names, "orthologr.legacy"))
        nm <- setdiff(nm, "orthologr.tmpdir")
        rlang::local_options(!!!stats::setNames(vector("list", length(nm)), nm), .frame = env)
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

test_that("orthologr_legacy_options(set = FALSE) lists the options and sets nothing", {
        without_orthologr_options()
        l <- orthologr_legacy_options(set = FALSE)
        expect_identical(names(l), legacy_names)
        expect_identical(l$orthologr.best_hit_rule, "evalue")
        expect_identical(l$orthologr.status_column, FALSE)
        expect_null(getOption("orthologr.best_hit_rule"))
        expect_error(orthologr_legacy_options(set = NA), "set")
})

test_that("orthologr_legacy_options() sets every option and returns the previous values", {
        without_orthologr_options()
        options(orthologr.best_hit_rule = "bitscore")
        expect_message(old <- orthologr_legacy_options(), "settings of version 0.4.3 are in force")
        expect_identical(getOption("orthologr.best_hit_rule"), "evalue")
        expect_true(getOption("orthologr.legacy_pairwise_aln"))
        expect_true(getOption("orthologr.legacy_comeron_saturation"))
        expect_true(getOption("orthologr.legacy_kaks_zero_dn"))
        expect_true(getOption("orthologr.legacy_kaks_saturation"))
        expect_false(getOption("orthologr.status_column"))
        expect_identical(old$orthologr.best_hit_rule, "bitscore")
        options(old)
        expect_identical(getOption("orthologr.best_hit_rule"), "bitscore")
        expect_null(getOption("orthologr.legacy_pairwise_aln"))
})

test_that("options(orthologr.legacy = TRUE) is a shorthand, and an option that is set itself wins", {
        without_orthologr_options()
        expect_identical(orthologr:::orthologr_best_hit_rule(NULL), "bitscore")
        expect_false(orthologr:::comeron_legacy_saturation())
        expect_true(orthologr:::dnds_status_column())
        options(orthologr.legacy = TRUE)
        expect_identical(orthologr:::orthologr_best_hit_rule(NULL), "evalue")
        expect_true(orthologr:::comeron_legacy_saturation())
        expect_true(orthologr:::kaks_legacy_zero_dn())
        expect_false(orthologr:::dnds_status_column())
        expect_identical(orthologr:::orthologr_option("orthologr.legacy_pairwise_aln", FALSE), TRUE)
        inforce <- orthologr:::orthologr_options_in_force()
        expect_identical(inforce$orthologr.best_hit_rule, "evalue")
        expect_identical(inforce$orthologr.legacy, TRUE)
        # set on its own: wins over the shorthand
        options(orthologr.best_hit_rule = "bitscore", orthologr.status_column = TRUE)
        expect_identical(orthologr:::orthologr_best_hit_rule(NULL), "bitscore")
        expect_true(orthologr:::dnds_status_column())
        # an argument wins over both
        expect_identical(orthologr:::orthologr_best_hit_rule("evalue"), "evalue")
})

test_that("the default of best_hit_rule follows the legacy shorthand in the exported functions", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        without_orthologr_options()
        new <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_identical(attr(new, "orthologr_run")$arguments$best_hit_rule, "bitscore")
        options(orthologr.legacy = TRUE)
        old <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_identical(attr(old, "orthologr_run")$arguments$best_hit_rule, "evalue")
        expect_identical(attr(old, "orthologr_run")$options$orthologr.legacy_pairwise_aln, TRUE)
        expect_false(is.element("status", names(old)))
})

for (how in c("orthologr_legacy_options()", "options(orthologr.legacy = TRUE)")) {
        test_that(paste0(how, " reproduces the stored table of orthologr 0.4.3 (20 genes, Comeron)"), {
                skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
                without_orthologr_options()
                if (startsWith(how, "options")) options(orthologr.legacy = TRUE) else suppressMessages(orthologr_legacy_options())
                ref <- stored("Comeron")
                expect_identical(ref$orthologr, "0.4.3")
                res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "Comeron", print_citation = FALSE))
                expect_identical(names(res), names(ref$table))
                expect_equal(plain_table(res), ref$table, tolerance = 1e-12)
        })

        test_that(paste0(how, " reproduces the stored table of orthologr 0.4.3 (20 genes, YN)"), {
                skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
                skip_if_not(kaks_available, "KaKs_Calculator is not installed")
                without_orthologr_options()
                if (startsWith(how, "options")) options(orthologr.legacy = TRUE) else suppressMessages(orthologr_legacy_options())
                ref <- stored("YN")
                res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "YN", print_citation = FALSE))
                expect_identical(names(res), names(ref$table))
                expect_equal(plain_table(res), ref$table, tolerance = 1e-12)
        })
}

test_that("without the legacy settings the default result differs from 0.4.3 only as documented", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        without_orthologr_options()
        ref <- stored("Comeron")$table
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_identical(setdiff(names(res), names(ref)), "status")
        expect_identical(sort(res$query_id), sort(ref$query_id))
})
