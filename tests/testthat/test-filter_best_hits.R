test_that("filter_best_hits() keeps the hit with the smallest e-value", {
        x <- tibble::tibble(subject_id = c("a", "b", "c"), evalue = c(1e-5, 1e-20, 1e-10), alig_length = c(300L, 100L, 200L))
        expect_identical(filter_best_hits(x)$subject_id, "b")
})

test_that("filter_best_hits() breaks ties by alignment length, then by position", {
        x <- tibble::tibble(subject_id = c("a", "b", "c", "d"), evalue = c(0, 0, 0, 1e-3), alig_length = c(100L, 300L, 300L, 900L))
        expect_identical(filter_best_hits(x)$subject_id, "b")
        expect_equal(nrow(filter_best_hits(x)), 1)
})

test_that("filter_best_hits() stops for an empty table", {
        x <- tibble::tibble(subject_id = character(0), evalue = numeric(0), alig_length = integer(0))
        expect_error(filter_best_hits(x), "empty")
})
