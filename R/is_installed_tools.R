# Internal checks that an external program can be started. Each program is started once per R session (and per
# `path`); see tool_available(). All return TRUE or stop with a message that says what is missing.

not_installed_message <- function(name) {
        paste0("It seems like you don't have ", name, " installed locally on your machine or the PATH variable to the ",
               name, " program is not set correctly.")
}

is_installed_clustalw <- function(path = NULL) {
        tool_available("clustalw2", version_args = "-help", pattern = "CLUSTAL", path = path,
                       not_found = not_installed_message("ClustalW"))
        return(TRUE)
}

is_installed_clustalo <- function(path = NULL) {
        tool_available("clustalo", version_args = "--help", pattern = "Clustal", path = path,
                       not_found = not_installed_message("Clustal-Omega"))
        return(TRUE)
}

is_installed_mafft <- function(path = NULL) {
        # 'mafft --version' prints e.g. "v7.526 (2024/Apr/26)" ('mafft -v' is not an option: exit status 1)
        tool_available("mafft", version_args = "--version", pattern = "^v[0-9]|MAFFT", path = path,
                       not_found = not_installed_message("MAFFT"))
        return(TRUE)
}

is_installed_muscle <- function(path = NULL) {
        muscle_version_line(path = path)
        return(TRUE)
}

# "MUSCLE v3.8.31 by Robert C. Edgar" or "muscle 5.3.linux64 []"
muscle_version_line <- function(path = NULL) {
        tool_available("muscle", version_args = "-version", pattern = "muscle", path = path,
                       not_found = not_installed_message("MUSCLE"))
}

# The major version of MUSCLE (3 or 5): MUSCLE 5 has another command line (-align / -output instead of -in / -out)
# and writes aligned FASTA only (no CLUSTAL format), with the sequences in an order of its own.
muscle_major_version <- function(path = NULL) {
        line <- muscle_version_line(path = path)
        major <- suppressWarnings(as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1", line)))
        if (is.na(major))
                major <- 3L
        major
}

# Put the sequences of an aligned FASTA file into the order of the ids in `ids` (MUSCLE 5 reorders its input).
reorder_fasta <- function(file, ids) {
        aln <- Biostrings::readAAStringSet(file)
        aln_ids <- sub("\\s.*$", "", names(aln))
        idx <- match(ids, aln_ids)
        if (length(aln) != length(ids) || anyNA(idx))
                stop("The alignment in '", file, "' does not hold the sequences of the input file.", call. = FALSE)
        Biostrings::writeXStringSet(aln[idx], file)
        invisible(file)
}

is_installed_tcoffee <- function(path = NULL) {
        tool_available("t_coffee", version_args = "-version", pattern = "T-COFFEE", path = path,
                       not_found = not_installed_message("T-COFFEE"))
        return(TRUE)
}

is_installed_kaks_calculator <- function(path = NULL) {
        tool_available("KaKs_Calculator", version_args = "-h", pattern = "KaKs_Calculator", path = path,
                       not_found = paste0(not_installed_message("KaKs_Calculator"),
                                          " Please read the INSTALLATION vignette for details."))
        return(TRUE)
}
