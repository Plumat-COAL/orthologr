
protein_file  <- system.file("seqs/ortho_thal_aa.fasta",  package = "orthologr")
protein_file2 <- system.file("seqs/ortho_lyra_aa.fasta",  package = "orthologr")
cds_file      <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")

expected_cols <- c("representative_id", "member_id")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

# --- Output structure --------------------------------------------------------

test_that("deepclust() returns a tibble", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_true(tibble::is_tibble(result))
})

test_that("deepclust() output has exactly 2 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_equal(colnames(result), expected_cols)
})

test_that("deepclust() output columns are character type", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_type(result[["representative_id"]], "character")
        expect_type(result[["member_id"]],         "character")
})

test_that("deepclust() output is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_gt(nrow(result), 0L)
})

# --- Cluster membership invariants -------------------------------------------

test_that("every input sequence appears as a cluster member", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)

        input_ids <- names(seqinr::read.fasta(protein_file, seqtype = "AA",
                                               as.string = TRUE))
        expect_true(all(input_ids %in% result[["member_id"]]))
})

test_that("all representative_ids are themselves cluster members", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)

        expect_true(all(result[["representative_id"]] %in% result[["member_id"]]))
})

# --- seq_type = "cds" --------------------------------------------------------

test_that("deepclust() runs with seq_type = 'cds'", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = cds_file, seq_type = "cds")
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

# --- Multi-file input --------------------------------------------------------

test_that("deepclust() accepts a vector of input files and returns a tibble", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = c(protein_file, protein_file2))
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

test_that("deepclust() with multiple files contains all sequences from all files", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = c(protein_file, protein_file2))

        ids1 <- names(seqinr::read.fasta(protein_file,  seqtype = "AA", as.string = TRUE))
        ids2 <- names(seqinr::read.fasta(protein_file2, seqtype = "AA", as.string = TRUE))

        expect_true(all(c(ids1, ids2) %in% result[["member_id"]]))
})

test_that("deepclust() errors informatively when a file does not exist", {
        expect_error(
                deepclust(input_file = c(protein_file, "nonexistent_file.fasta")),
                regexp = "not found"
        )
})

# --- approx_id threshold -----------------------------------------------------

# (This test used to compare the numbers of rows, "stricter approx_id yields at least as many rows". deepclust()
# returns one row per sequence whatever the threshold, so it compared 20 with 20 and could not fail. What a stricter
# threshold changes is the number of clusters; it is tested on sequences that do cluster: every protein and a
# copy of it with a few changed residues.)
test_that("deepclust() returns one row per sequence; a stricter approx_id cannot give fewer clusters", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        
        aa <- Biostrings::readAAStringSet(protein_file)
        names(aa) <- paste0("p", seq_along(aa))
        variant <- as.character(aa)
        # every 12th residue of the copy becomes alanine (or glycine): about 92% identity
        variant <- vapply(variant, function(x) {
                ch <- strsplit(x, "")[[1]]
                i <- seq(6, length(ch) - 1, by = 12)
                ch[i] <- ifelse(ch[i] == "A", "G", "A")
                paste(ch, collapse = "")
        }, "")
        both <- c(aa, Biostrings::AAStringSet(stats::setNames(variant, paste0("v", seq_along(aa)))))
        f <- tempfile(fileext = ".faa")
        Biostrings::writeXStringSet(both, f)
        
        lenient <- suppressMessages(deepclust(input_file = f, approx_id = 50))
        strict  <- suppressMessages(deepclust(input_file = f, approx_id = 99))
        n <- length(both)
        expect_equal(nrow(lenient), n)
        expect_equal(nrow(strict), n)
        expect_setequal(lenient$member_id, names(both))
        n_lenient <- length(unique(lenient$representative_id))
        n_strict  <- length(unique(strict$representative_id))
        # the copies cluster with their originals at 50% identity ...
        expect_lt(n_lenient, n)
        expect_lte(n_lenient, length(aa) + 2)
        # ... and not at 99%
        expect_gt(n_strict, n_lenient)
})

# --- Duplicate ID warning (multi-file merge) ---------------------------------

test_that("deepclust() warns when duplicate sequence IDs exist across input files", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")

        # Copy protein_file so all IDs are guaranteed to overlap
        tmp_dup <- tempfile(fileext = ".fasta")
        file.copy(protein_file, tmp_dup)
        on.exit(unlink(tmp_dup), add = TRUE)

        expect_warning(
                deepclust(input_file = c(protein_file, tmp_dup)),
                regexp = "Duplicate sequence IDs"
        )
})

test_that("deepclust() duplicate warning names affected IDs", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")

        tmp_dup <- tempfile(fileext = ".fasta")
        file.copy(protein_file, tmp_dup)
        on.exit(unlink(tmp_dup), add = TRUE)

        first_id <- names(seqinr::read.fasta(protein_file, seqtype = "AA",
                                              as.string = TRUE))[1]
        expect_warning(
                deepclust(input_file = c(protein_file, tmp_dup)),
                regexp = first_id,
                fixed  = TRUE
        )
})

test_that("deepclust() does not warn about duplicate IDs when all IDs are unique", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")

        # thal and lyra have completely distinct IDs
        warned <- FALSE
        withCallingHandlers(
                deepclust(input_file = c(protein_file, protein_file2)),
                warning = function(w) {
                        if (grepl("Duplicate sequence IDs", conditionMessage(w)))
                                warned <<- TRUE
                        invokeRestart("muffleWarning")
                }
        )
        expect_false(warned)
})

# --- save.output -------------------------------------------------------------

test_that("deepclust() writes a TSV when save.output is specified", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        # a folder of its own that holds a stale file of the same name. (The test used to write into the shared
        # tempdir(), where the file of an earlier test existed already: it passed although file.copy() never
        # overwrote, i.e. also when nothing was written.)
        out_dir  <- new_test_dir("deepclust_out_")
        filename <- paste0("deepclust_", basename(protein_file), ".tsv")
        writeLines("STALE", file.path(out_dir, filename))

        res <- suppressMessages(deepclust(input_file = protein_file, save.output = out_dir))

        expect_true(file.exists(file.path(out_dir, filename)))
        written <- utils::read.delim(file.path(out_dir, filename), header = FALSE, stringsAsFactors = FALSE)
        expect_equal(nrow(written), nrow(res))
        expect_setequal(written[[2]], res$member_id)
        # a file path instead of a folder
        one <- file.path(out_dir, "clusters.tsv")
        suppressMessages(deepclust(input_file = protein_file, save.output = one))
        expect_equal(length(readLines(one)), nrow(res))
})
