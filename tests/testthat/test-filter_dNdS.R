
# Build a minimal mock dNdS table matching the columns filter_dNdS() uses
mock_dnds <- tibble::tibble(
        query_id   = paste0("gene", 1:6),
        subject_id = paste0("subj", 1:6),
        dN         = c(0.01, 0.02, NA,   0.05, 0.10, 0.00),
        dS         = c(0.10, NA,   0.30, 0.20, 0.50, 0.00),
        dNdS       = c(0.10, NA,   NA,   0.25, 0.20, NA)
)

test_that("filter_dNdS() returns a data frame / tibble", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        expect_true(is.data.frame(result))
})

test_that("filter_dNdS() removes rows with NA in dN", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        expect_false(any(is.na(result[["dN"]])))
})

test_that("filter_dNdS() removes rows with NA in dS", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        expect_false(any(is.na(result[["dS"]])))
})

test_that("filter_dNdS() removes rows where dNdS > threshold", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 0.15)
        expect_true(all(result[["dNdS"]] <= 0.15))
})

test_that("filter_dNdS() retains rows that pass all filters", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        # gene1: dN=0.01, dS=0.10, dNdS=0.10 -> should be retained
        expect_true("gene1" %in% result[["query_id"]])
})

test_that("filter_dNdS() with a very low threshold removes high dNdS rows", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 0.05)
        expect_true(nrow(result) < nrow(mock_dnds))
})

# gene6 has dN = dS = 0: dN/dS is not defined. (This test used to assert that such a row is dropped: it encoded
# the silent loss of gene pairs with dS = 0 as the expected behaviour.)
test_that("filter_dNdS() keeps gene pairs with dS = 0 by default, flagged, and drops them on request", {
        result <- suppressMessages(filter_dNdS(mock_dnds, dnds.threshold = Inf))
        expect_identical(result$query_id, c("gene1", "gene4", "gene5", "gene6"))
        expect_identical(result$dS_zero, c(FALSE, FALSE, FALSE, TRUE))
        expect_true(is.na(result$dNdS[4]))
        expect_message(filter_dNdS(mock_dnds), "1 gene pairs have dS = 0.*1 with dN = 0 as well.*and 0 with dN > 0.*kept with dNdS = NA")
        
        # the behaviour of orthologr <= 0.4.3: the argument, or the option
        old <- dplyr::filter(mock_dnds, !is.na(dN), !is.na(dS), !is.na(dNdS))
        expect_identical(suppressMessages(filter_dNdS(mock_dnds, dnds.threshold = Inf, dS_zero = "drop")), old)
        opt <- options(orthologr.dS_zero = "drop")
        expect_identical(suppressMessages(filter_dNdS(mock_dnds, dnds.threshold = Inf)), old)
        options(opt)
        expect_message(filter_dNdS(mock_dnds, dS_zero = "drop"), "They are left out")
        expect_error(filter_dNdS(mock_dnds, dS_zero = "zero"), "\"keep\" or \"drop\"")
})

test_that("filter_dNdS() tells identical sequences from dN > 0 with dS = 0, and keeps both (N09)", {
        x <- tibble::tibble(query_id = c("identical", "nonsyn_only", "normal", "high", "dS_missing"),
                            subject_id = letters[1:5],
                            dN = c(0, 0.01, 0.1, 0.9, 0.1), dS = c(0, 0, 0.5, 0.3, NA), dNdS = c(NA, NA, 0.2, 3, NA))
        expect_message(res <- filter_dNdS(x), "2 gene pairs have dS = 0.*1 with dN = 0 as well.*and 1 with dN > 0")
        expect_identical(res$query_id, c("identical", "nonsyn_only", "normal"))
        expect_identical(res$dS_zero, c(TRUE, TRUE, FALSE))
        expect_identical(suppressMessages(filter_dNdS(x, dS_zero = "drop"))$query_id, "normal")
        # no dS = 0: the same rows as before, plus the flag column
        y <- x[3:5, ]
        expect_identical(suppressMessages(filter_dNdS(y))[ , names(y)], suppressMessages(filter_dNdS(y, dS_zero = "drop")))
})
