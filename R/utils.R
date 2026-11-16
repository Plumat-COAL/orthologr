# Suppress R CMD CHECK NOTEs for data.table's `.()` alias and other NSE symbols
utils::globalVariables(".")


is.dnaSequence <- function(seq){
        
        # only A, C, G, T (upper or lower case). The former pattern "[^(A|C|G|T|a|c|g|t)]" also let '(', ')' and
        # '|' pass.
        # (PCRE on bytes: the same answer as the default regex engine, which took about 50 us per kb)
        return(!grepl("[^ACGTacgt]", seq, perl = TRUE, useBytes = TRUE))
}


# This functions are internal functions to proof if the used tool is provided. 

is.ortho_detection_method <- function(method = NULL){
        
        provided <- c("BH","RBH")
        
        if(is.null(method)){
                message("The following methods are provided: ", paste(provided, collapse = ", "))
                return(FALSE)
        }
        
        return(is.element(method, provided))
}


is.multiple_aln_tool <- function(tool = NULL){
        
        provided <- c("clustalw", "clustalo","muscle", "t_coffee", "mafft")
        
        if(is.null(tool)){
                message("The following methods are provided: ", paste(provided, collapse = ", "))
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.pairwise_aln_tool <- function(tool = NULL){
        
        provided <- c("NW")
        
        if(is.null(tool)){
                message("The following methods are provided: ", paste(provided, collapse = ", "))
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.codon_aln_tool <- function(tool = NULL){
        
        provided <-  c("pal2nal", "coal")
        
        if(is.null(tool)){
                message("The following methods are provided: ", paste(provided, collapse = ", "))
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

#' @title Translate CDS file to Amino Acids file
#' @description This function takes a \code{data.table} object as input and
#' translates the cds sequences stored as \code{seqs} column into the corresponding amino acid
#' sequence.
#' @param dt a \code{data.table} object storing cds sequences in a column named \code{seqs}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
cds2aa <- function(dt, delete_corrupt_cds = TRUE){
        
        if (!is.data.table(dt))
                stop("Your CDS file was not corretly transformed into a data.table object.")
        
        # define visible bindings for global variables
        seqs <- aa <- geneids <- NULL
        
        # the same ID more than once: the sequences cannot be told apart later on (hit tables, joins)
        check_duplicate_ids(dt[ , geneids], what = "coding sequences")
        
        n_total <- nrow(dt)
        
        # omit empty sequences - and say which
        empty <- dt[ , is.na(seqs) | seqs == ""]
        if (any(empty)) {
                report_dropped_cds(dt[empty, geneids], reason = "empty", total = n_total)
                dt <- dt[!empty]
        }
        
        if (delete_corrupt_cds) {
                # omit sequences that are not multiples of 3 (read.cds() has removed and reported them already
                # when the table comes from there)
                not3 <- dt[ , nchar(seqs) %% 3 != 0]
                if (any(not3)) {
                        report_dropped_cds(dt[not3, geneids], reason = "length", total = n_total)
                        dt <- dt[!not3]
                }
        }
        
        # omit sequences consisting of others than ACGT - and say which. With options(orthologr.ambiguous_cds = "keep")
        # (the default since 0.5.0) sequences with IUPAC ambiguity codes (N, R, Y, ...) are kept: a codon with such a
        # letter is translated to X and left to the alignment and estimation programs, which do not compare it (see
        # ?orthologr_options); only other letters are left out. "drop" (0.4.3) leaves them all out.
        not_acgt <- !is.dnaSequence(as.character(dt[ , seqs]))
        if (any(not_acgt) && identical(ambiguous_cds_mode(), "keep")) {
                iupac <- grepl("^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]*\\z", dt[ , seqs], perl = TRUE, useBytes = TRUE)
                kept <- not_acgt & iupac
                if (any(kept)) {
                        ambiguous_record(dt[kept, geneids], dt[kept, seqs], .orthologr_state$current_file)
                        if (!sequence_messages_quiet())
                                message(sum(kept), " coding sequence", if (sum(kept) > 1) "s" else "", " with ambiguous bases (e.g. N) ",
                                        if (sum(kept) > 1) "are" else "is", " kept (options(orthologr.ambiguous_cds = \"keep\"), the default; ",
                                        "\"drop\" leaves them out as 0.4.3 did); codons with such a base are translated to X and not compared.")
                }
                not_acgt <- not_acgt & !iupac
        }
        if (any(not_acgt)) {
                report_dropped_cds(dt[not_acgt, geneids], reason = "not_acgt", total = n_total,
                                   seqs = dt[not_acgt, seqs])
                dt <- dt[!not_acgt]
        }
        
        # translate cds to protein sequences
        tryCatch({
                dt[ , aa := translate_cds(seqs)]
        }, error = function(e) {
                stop(
                        "The input coding sequences could not be translated properly to amino acid sequences: ",
                        conditionMessage(e),
                        "\n",
                        " Please check whether the input file stores valid coding sequences.",
                        call. = FALSE
                )
        })
        
        # inside a dNdS() / dNdS_coal() call: the sequences with in-frame stop codons before the last codon (L6)
        if (!is.null(.orthologr_state$stops) && nrow(dt) > 0)
                stop_codon_record(dt[ , geneids], dt[ , aa], .orthologr_state$current_file, cds = dt[ , seqs])
        
        return(dt)
}

# Stop when a sequence ID occurs more than once. The ID of a sequence is its FASTA header up to the first white
# space; two sequences with the same ID cannot be told apart in the hit tables and joins that follow (the
# translation used to fail with "Wrong argument type in s2c()" / "could not be translated properly").
check_duplicate_ids <- function(ids, what = "sequences", file = NULL) {
        if (is.null(file))
                file <- .orthologr_state$current_file
        dup <- unique(ids[duplicated(ids)])
        if (length(dup) == 0)
                return(invisible(TRUE))
        shown <- utils::head(dup, 10)
        orthologr_abort(
                class = "bad_input",
                length(dup), " sequence ID", if (length(dup) > 1) "s occur" else " occurs", " more than once among the ",
                what, if (!is.null(file)) paste0(" of '", file, "'") else "", ": ",
                paste0("'", shown, "'", collapse = ", "),
                if (length(dup) > length(shown)) paste0(", ... (", length(dup) - length(shown), " more)") else "",
                ". The ID of a sequence is its FASTA header up to the first white space and must be unique within a file. ",
                "Please rename or remove the duplicated sequences.",
                call. = FALSE
        )
}

# Sequences that are left out are never left out silently: a message names them, with their number, the number
# of sequences read and the reason, and inside dNdS() they are collected for attr(<result>, "dropped_sequences").
#   reason : "not_acgt" - characters other than A, C, G, T (cannot be translated)
#            "empty"    - no sequence
#            "length"   - length not a multiple of 3 (delete_corrupt_cds = TRUE)
#   seqs   : the sequences (optional, reason "not_acgt"): used to say how many of them hold N only
report_dropped_cds <- function(ids, reason = "not_acgt", total = NULL, seqs = NULL, file = NULL) {
        if (length(ids) == 0)
                return(invisible(NULL))
        if (is.null(file))
                file <- .orthologr_state$current_file
        reason_text <- switch(
                reason,
                not_acgt = "characters other than A, C, G, T (e.g. N or IUPAC ambiguity codes); they cannot be translated",
                empty    = "no sequence",
                length   = "a length that is not a multiple of 3 (delete_corrupt_cds = TRUE)",
                reason
        )
        dropped_record(ids, reason_text, file)
        if (sequence_messages_quiet())
                return(invisible(ids))
        run_dir <- getOption("orthologr.tmpdir", NULL)
        if (!is.null(run_dir)) {
                # inside one dNdS() call: say it once
                key <- paste0("dropped_cds\r", run_dir, "\r", reason, "\r", paste(ids, collapse = "\r"))
                if (!is.null(.orthologr_state[[key]]))
                        return(invisible(NULL))
                .orthologr_state[[key]] <- TRUE
        }
        shown <- utils::head(ids, 10)
        only_n <- if (!is.null(seqs)) sum(!grepl("[^ACGTNacgtn]", seqs)) else NA
        message(
                length(ids), if (!is.null(total)) paste0(" of ", total), " coding sequence",
                if (length(ids) > 1 || !is.null(total)) "s" else "",
                if (!is.null(file)) paste0(" of '", basename(file), "'") else "",
                if (length(ids) > 1) " have " else " has ", reason_text,
                if (!is.na(only_n) && only_n > 0)
                        paste0(" (", only_n, " of them with N as the only such character, e.g. a 5'-incomplete CDS padded with N)")
                else "",
                " and ", if (length(ids) > 1) "are" else "is", " left out: ",
                paste(shown, collapse = ", "),
                if (length(ids) > length(shown)) paste0(", ... (", length(ids) - length(shown), " more)") else "",
                "."
        )
        invisible(ids)
}

# Inside a dNdS() / dNdS_coal() call with quiet = TRUE: the messages about single input sequences are left to the
# summary (see dropped_start()).
sequence_messages_quiet <- function() {
        !is.null(.orthologr_state$dropped) && isTRUE(.orthologr_state$sequence_messages_quiet)
}

# Collect the coding sequences kept with ambiguous bases (options(orthologr.ambiguous_cds = "keep")) while a dNdS()
# call runs: their IDs and the number of their codons that hold such a base (codons that are not compared).
ambiguous_record <- function(ids, seqs, file) {
        if (is.null(.orthologr_state$ambiguous))
                return(invisible(NULL))
        seqs <- as.character(seqs)
        n_codons <- vapply(seqs, function(x) {
                cod <- substring(x, seq(1, nchar(x) - 2, by = 3), seq(3, nchar(x), by = 3))
                sum(grepl("[^ACGTacgt]", cod))
        }, integer(1), USE.NAMES = FALSE)
        .orthologr_state$ambiguous[[length(.orthologr_state$ambiguous) + 1L]] <- data.frame(
                file             = if (is.null(file)) NA_character_ else as.character(file),
                id               = as.character(ids),
                ambiguous_codons = as.integer(n_codons),
                stringsAsFactors = FALSE
        )
        invisible(NULL)
}

# The table of the sequences kept with ambiguous bases so far in this call (file, id, ambiguous_codons), or NULL.
ambiguous_table <- function() {
        res <- .orthologr_state$ambiguous
        if (length(res) == 0)
                return(NULL)
        res <- unique(do.call(rbind, res))
        rownames(res) <- NULL
        res
}

# Collect the coding sequences with in-frame stop codons before their last codon while a dNdS() call runs (L6 of
# REVIEW-OL-1): their IDs and the number of such codons, counted in the translation ('*' before the last residue).
# A stop codon written with an ambiguity code - TAR (TAA / TAG) or TRA (TAA / TGA), every reading of which is a stop
# codon - is translated to X (options(orthologr.ambiguous_cds = "keep")) and is counted from the coding sequence
# 'cds' (when given). Nothing is changed: the sequences are used as before (COAL and PAL2NAL do not compare the
# codon; see ?dNdS).
stop_codon_record <- function(ids, aa, file, cds = NULL) {
        aa <- as.character(aa)
        inner <- substr(aa, 1L, pmax(nchar(aa) - 1L, 0L))
        n <- nchar(inner) - nchar(gsub("*", "", inner, fixed = TRUE))
        if (!is.null(cds))
                n <- n + ambiguous_stop_codons(cds)
        has <- !is.na(n) & n > 0
        if (!any(has))
                return(invisible(NULL))
        .orthologr_state$stops[[length(.orthologr_state$stops) + 1L]] <- data.frame(
                file        = if (is.null(file)) NA_character_ else as.character(file),
                id          = as.character(ids[has]),
                stop_codons = as.integer(n[has]),
                stringsAsFactors = FALSE
        )
        invisible(NULL)
}

# The number of codons before the last codon of each coding sequence that are stop codons written with an ambiguity
# code (TAR, TRA; any case) in the standard genetic code. Codons are read from the first base.
ambiguous_stop_codons <- function(cds) {
        cds <- as.character(cds)
        n <- integer(length(cds))
        has <- !is.na(cds) & grepl("[Rr]", cds, perl = TRUE, useBytes = TRUE)
        for (i in which(has)) {
                len <- nchar(cds[i]) %/% 3L
                if (len < 2L)
                        next
                cod <- toupper(substring(cds[i], seq(1L, by = 3L, length.out = len - 1L), seq(3L, by = 3L, length.out = len - 1L)))
                n[i] <- sum(cod == "TAR" | cod == "TRA")
        }
        n
}

# The table of the sequences with in-frame stop codons so far in this call (file, id, stop_codons), or NULL.
stop_codon_table <- function() {
        res <- .orthologr_state$stops
        if (length(res) == 0)
                return(NULL)
        res <- unique(do.call(rbind, res))
        rownames(res) <- NULL
        res
}

# Collect the coding sequences that read.cds() trimmed to a multiple of 3 (delete_corrupt_cds = FALSE) while a
# dNdS() call runs; `corrupted_file` holds them as they were read.
trimmed_record <- function(ids, removed, file, corrupted_file) {
        if (is.null(.orthologr_state$trimmed))
                return(invisible(NULL))
        .orthologr_state$trimmed[[length(.orthologr_state$trimmed) + 1L]] <- data.frame(
                file           = as.character(file),
                id             = as.character(ids),
                removed_nt     = as.integer(removed),
                corrupted_file = as.character(corrupted_file),
                stringsAsFactors = FALSE
        )
        invisible(NULL)
}

# The table of the sequences trimmed so far in this call (file, id, removed_nt, corrupted_file; no row twice), or
# NULL when there are none.
trimmed_table <- function() {
        res <- .orthologr_state$trimmed
        if (length(res) == 0)
                return(NULL)
        res <- unique(do.call(rbind, res))
        rownames(res) <- NULL
        res
}

# Collect the sequences left out while a dNdS() call runs (see dropped_start()).
dropped_record <- function(ids, reason, file = NULL) {
        if (is.null(.orthologr_state$dropped))
                return(invisible(NULL))
        .orthologr_state$dropped[[length(.orthologr_state$dropped) + 1L]] <- data.frame(
                file   = if (is.null(file)) NA_character_ else as.character(file),
                id     = as.character(ids),
                reason = reason,
                stringsAsFactors = FALSE
        )
        invisible(NULL)
}

# Start collecting; returns a function that ends it and returns the table (file, id, reason; no row twice).
#   quiet : the run's quiet (dNdS(), dNdS_coal()). With quiet = TRUE the messages about single sequences (left
#           out, trimmed to a multiple of 3) are not shown: they are counted in the summary at the end, and the
#           sequences are listed in attr(<result>, "dropped_sequences") / attr(<result>, "trimmed_sequences").
dropped_start <- function(quiet = FALSE) {
        if (!is.null(.orthologr_state$dropped))
                return(function() NULL)
        .orthologr_state$dropped <- list()
        .orthologr_state$trimmed <- list()
        .orthologr_state$ambiguous <- list()
        .orthologr_state$stops <- list()
        .orthologr_state$sequence_messages_quiet <- isTRUE(quiet)
        function() {
                res <- .orthologr_state$dropped
                .orthologr_state$dropped <- NULL
                .orthologr_state$trimmed <- NULL
                .orthologr_state$ambiguous <- NULL
                .orthologr_state$stops <- NULL
                .orthologr_state$sequence_messages_quiet <- NULL
                if (length(res) == 0)
                        return(data.frame(file = character(0), id = character(0), reason = character(0),
                                          stringsAsFactors = FALSE))
                res <- unique(do.call(rbind, res))
                rownames(res) <- NULL
                res
        }
}

# The reports about input sequences for orthogroups() and dNdS_orthogroups() (M4 of REVIEW-OL-3): as in dNdS(),
# the messages read.cds() and cds2aa() give for every file (sequences trimmed to a multiple of 3, left out, kept
# with ambiguous bases) are collected while the files are read - with quiet = TRUE they are not shown - and the
# tables are kept in the result. Returns a function that ends the collecting and returns
# list(owner, dropped, trimmed, ambiguous, stops): `owner` is FALSE inside a call that already collects (orthogroups()
# called by dNdS_orthogroups()), which then reports them itself; the tables are NULL when empty.
input_report_start <- function(quiet = TRUE) {
        owner <- is.null(.orthologr_state$dropped)
        end <- dropped_start(quiet = quiet)
        function() {
                trimmed <- trimmed_table()
                ambiguous <- ambiguous_table()
                stops <- stop_codon_table()
                dropped <- if (owner) end() else NULL
                if (!is.null(dropped) && nrow(dropped) == 0)
                        dropped <- NULL
                list(owner = owner, dropped = dropped, trimmed = trimmed, ambiguous = ambiguous, stops = stops)
        }
}

# The tables of input_report_start() as the attributes dNdS() uses ("dropped_sequences", "trimmed_sequences",
# "ambiguous_sequences", "stop_codon_sequences"), set on `res` (only those with rows).
input_report_attributes <- function(res, report) {
        tabs <- list(dropped_sequences = report$dropped, trimmed_sequences = report$trimmed,
                     ambiguous_sequences = report$ambiguous, stop_codon_sequences = report$stops)
        for (n in names(tabs)) {
                if (!is.null(tabs[[n]]) && NROW(tabs[[n]]) > 0)
                        attr(res, n) <- tabs[[n]]
        }
        res
}

# One line about the input sequences of a report of input_report_start() for `fun` ("orthogroups()"); NULL when there
# is nothing to say.
input_report_text <- function(report, fun) {
        n <- function(x) if (is.null(x)) 0L else length(unique(paste(x$file, x$id, sep = "\r")))
        parts <- c(
                if (n(report$dropped) > 0) paste0(n(report$dropped), " left out (attr(<result>, \"dropped_sequences\"))"),
                if (n(report$trimmed) > 0) paste0(n(report$trimmed), " with a length that is not a multiple of 3, trimmed at the 3' end ",
                                                  "(delete_corrupt_cds = FALSE; attr(<result>, \"trimmed_sequences\"))"),
                if (n(report$ambiguous) > 0) paste0(n(report$ambiguous), " kept with ambiguous bases, whose codons are not compared ",
                                                    "(attr(<result>, \"ambiguous_sequences\"))"),
                if (n(report$stops) > 0) paste0(n(report$stops), " with in-frame stop codons, used as they are ",
                                                "(attr(<result>, \"stop_codon_sequences\"))")
        )
        if (length(parts) == 0)
                return(NULL)
        paste0(fun, ": input coding sequences - ", paste(parts, collapse = "; "), ".")
}

# Translate coding sequences (character vector) with the standard genetic code; the result is what transl()
# (seqinr::translate()) returns for each of them. Sequences of A, C, G, T whose length is a multiple of 3 are
# translated in one call to Biostrings::translate() (no.init.codon = TRUE: the first codon is translated like any
# other, as seqinr does). All other sequences (ambiguity codes, incomplete last codon, empty) take the old path,
# so that their result cannot differ.
translate_cds <- function(seqs) {
        seqs <- as.character(seqs)
        res <- character(length(seqs))
        if (length(seqs) == 0)
                return(res)
        n <- nchar(seqs)
        # (the regular expressions with PCRE on bytes: the same answers, about 50 times faster on coding sequences)
        fast <- !is.na(seqs) & n > 0 & n %% 3 == 0 & !grepl("[^ACGTacgt]", seqs, perl = TRUE, useBytes = TRUE)
        # ambiguous bases (options(orthologr.ambiguous_cds = "keep")): a codon with such a base becomes X
        fuzzy <- !fast & !is.na(seqs) & n > 0 & n %% 3 == 0 &
                grepl("^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]*\\z", seqs, perl = TRUE, useBytes = TRUE)
        if (any(fuzzy))
                res[fuzzy] <- as.character(Biostrings::translate(
                        Biostrings::DNAStringSet(seqs[fuzzy]),
                        genetic.code   = Biostrings::GENETIC_CODE,
                        no.init.codon  = TRUE,
                        if.fuzzy.codon = "X"
                ))
        if (any(fast))
                res[fast] <- as.character(Biostrings::translate(
                        Biostrings::DNAStringSet(seqs[fast]),
                        genetic.code  = Biostrings::GENETIC_CODE,
                        no.init.codon = TRUE
                ))
        rest <- !fast & !fuzzy
        if (any(rest))
                res[rest] <- vapply(seqs[rest], transl, character(1), USE.NAMES = FALSE)
        res
}

# "keep" (default since 0.5.0) or "drop" (0.4.3; set by orthologr_legacy_options()): what happens to coding sequences
# with IUPAC ambiguity codes (N, R, Y, ...); options(orthologr.ambiguous_cds = ).
ambiguous_cds_mode <- function() {
        m <- orthologr_option("orthologr.ambiguous_cds", "keep")
        if (!is.character(m) || length(m) != 1 || !is.element(m, c("drop", "keep")))
                orthologr_abort("options(orthologr.ambiguous_cds = ) must be \"drop\" or \"keep\".", class = "bad_argument")
        m
}

# The part of each string in 'x' before the first 'sep' (a single character, taken literally); the whole string when
# it holds no 'sep'. What unlist(stringr::str_split(x, sep))[1] returned for one string, also "" for "" and NA for
# NA (strsplit() gives character(0) for "").
first_field <- function(x, sep) {
        sub(paste0("\\Q", sep, "\\E.*$"), "", x, perl = TRUE)
}
