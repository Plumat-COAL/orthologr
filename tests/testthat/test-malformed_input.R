# Malformed input files: what read.cds() and dNdS() do with an empty file, a truncated record, lower case, a
# protein file, duplicated IDs, Windows line endings, gzip and an N in a sequence.

malformed_files <- function() {
        q <- readLines(thal_cds())
        hdr <- grep("^>", q)
        d <- new_test_dir("malformed_")
        w <- function(name, lines, eol = "\n") {
                f <- file.path(d, name)
                con <- file(f, "wb")
                writeBin(charToRaw(paste0(paste(lines, collapse = eol), if (length(lines)) eol else "")), con)
                close(con)
                f
        }
        with_n <- q
        substr(with_n[hdr[1] + 1], 10, 10) <- "N"
        gz <- file.path(d, "gzip.fa.gz")
        con <- gzfile(gz, "w"); writeLines(q, con); close(con)
        list(
                empty       = w("empty.fa", character(0)),
                no_header   = w("no_header.fa", q[!startsWith(q, ">")]),
                truncated   = w("truncated.fa", c(q, ">AT9G99999.1 a record without sequence")),
                lowercase   = w("lowercase.fa", ifelse(startsWith(q, ">"), q, tolower(q))),
                protein     = thal_aa(),
                duplicate   = w("duplicate.fa", c(q, q[hdr[1]:(hdr[2] - 1)])),
                crlf        = w("crlf.fa", q, eol = "\r\n"),
                gzip        = gz,
                with_n      = w("with_n.fa", with_n)
        )
}

test_that("read.cds(): unreadable and protein files are classed errors; other variants are read", {
        f <- malformed_files()
        expect_error(read.cds(f$empty, format = "fasta"), class = "orthologr_error_bad_input")
        expect_error(read.cds(f$no_header, format = "fasta"), class = "orthologr_error_bad_input")
        expect_error(read.cds(f$protein, format = "fasta"), class = "orthologr_error_bad_input")
        ref <- read.cds(thal_cds(), format = "fasta")
        for (v in c("lowercase", "crlf", "gzip")) {
                x <- read.cds(f[[v]], format = "fasta")
                expect_identical(x$geneids, ref$geneids, info = v)
                expect_identical(toupper(x$seqs), toupper(ref$seqs), info = v)
        }
        # a header without a sequence is read as a record with an empty sequence (dNdS() leaves it out, see below)
        x <- read.cds(f$truncated, format = "fasta")
        expect_identical(nrow(x), nrow(ref) + 1L)
        expect_identical(x$seqs[x$geneids == "AT9G99999.1"], "")
})

test_that("dNdS(): every malformed variant either gives the 20 pairs (dropped records reported) or a classed error", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        f <- malformed_files()
        run <- function(file) suppressWarnings(suppressMessages(dNdS(file, lyra_cds(), print_citation = FALSE)))
        ref <- run(thal_cds())
        for (v in c("lowercase", "crlf", "gzip")) {
                x <- run(f[[v]])
                expect_identical(nrow(x), nrow(ref), info = v)
                expect_equal(x$dNdS, ref$dNdS, info = v)
        }
        x <- run(f$truncated)
        expect_identical(nrow(x), nrow(ref))
        expect_identical(attr(x, "dropped_sequences")$id, "AT9G99999.1")
        expect_match(attr(x, "dropped_sequences")$reason, "no sequence")
        # a CDS with an N: kept by default (0.5.0), left out with options(orthologr.ambiguous_cds = "drop")
        x <- run(f$with_n)
        expect_identical(nrow(x), nrow(ref))
        expect_identical(nrow(attr(x, "ambiguous_sequences")), 1L)
        old <- options(orthologr.ambiguous_cds = "drop")
        x <- run(f$with_n)
        options(old)
        expect_identical(nrow(x), nrow(ref) - 1L)
        expect_match(attr(x, "dropped_sequences")$reason, "other than A, C, G, T")
        expect_error(run(f$duplicate), "more than once", class = "orthologr_error_bad_input")
        expect_error(run(f$empty), class = "orthologr_error_bad_input")
        expect_error(run(f$protein), class = "orthologr_error_bad_input")
})

test_that("options(orthologr.ambiguous_cds = 'keep') (the default) keeps a CDS with an N; its pair gets values close to the N-free ones", {
        expect_identical(translate_cds(c("ATGNCTAAA", "ATGGCTAAA", "atgrctaaa")), c("MXK", "MAK", "MXK"))
        withr_opts <- options(orthologr.ambiguous_cds = "nonsense")
        on.exit(options(withr_opts), add = TRUE)
        expect_error(ambiguous_cds_mode(), class = "orthologr_error_bad_argument")
        options(orthologr.ambiguous_cds = NULL)
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        f <- malformed_files()
        ref <- suppressWarnings(suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE)))
        # the default: no option needed
        expect_identical(ambiguous_cds_mode(), "keep")
        msgs <- character(0)
        x <- withCallingHandlers(suppressWarnings(dNdS(f$with_n, lyra_cds(), print_citation = FALSE)),
                                 message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
        expect_true(any(grepl("ambiguous bases", msgs)))
        expect_identical(nrow(x), nrow(ref))
        expect_identical(nrow(attr(x, "dropped_sequences")), 0L)
        # the sequence, how many of its codons hold such a base, and the summary line
        amb <- attr(x, "ambiguous_sequences")
        expect_identical(nrow(amb), 1L)
        expect_true(amb$ambiguous_codons >= 1L)
        expect_true(any(grepl("1 row has a coding sequence with ambiguous bases", msgs)))
        expect_identical(attr(x, "orthologr_run")$options$orthologr.ambiguous_cds, "keep")
        m <- merge(as.data.frame(x), as.data.frame(ref), by = c("query_id", "subject_id"))
        expect_equal(m$dS.x, m$dS.y, tolerance = 0.05)
        # quiet = TRUE: the message about the sequence is left to the summary
        msgs <- character(0)
        x <- withCallingHandlers(suppressWarnings(dNdS(f$with_n, lyra_cds(), print_citation = FALSE, quiet = TRUE)),
                                 message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
        expect_false(any(grepl("are kept|is kept", msgs)))
        expect_true(any(grepl("ambiguous bases", msgs)))
})
