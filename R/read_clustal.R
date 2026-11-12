# Read an alignment in CLUSTAL format as a named character vector (one aligned sequence per name, in file order).
# Biostrings::readAAMultipleAlignment() and seqinr::read.alignment() do not read every CLUSTAL file the programs
# used here write (T-Coffee ends blocks with a line of blanks; PAL2NAL writes no conservation line): this reader
# is used when they fail.
read_clustal <- function(file) {
        x <- readLines(file, warn = FALSE)
        if (length(x) == 0 || !grepl("^(CLUSTAL|MUSCLE)", x[1]))
                stop("'", file, "' is not a file in CLUSTAL format.", call. = FALSE)
        x <- x[-1]
        # sequence lines start with the name; conservation lines with a blank
        x <- x[nzchar(x) & !grepl("^\\s", x)]
        parts <- strsplit(x, "\\s+")
        if (length(parts) == 0 || any(lengths(parts) < 2))
                stop("'", file, "' is not a file in CLUSTAL format.", call. = FALSE)
        nam <- vapply(parts, `[`, "", 1)
        seq <- vapply(parts, `[`, "", 2)
        ids <- unique(nam)
        res <- vapply(ids, function(i) paste(seq[nam == i], collapse = ""), "")
        if (length(unique(nchar(res))) != 1)
                stop("The sequences of the alignment in '", file, "' do not have the same length.", call. = FALSE)
        res
}
