# Suppress R CMD CHECK NOTEs for data.table's `.()` alias and other NSE symbols
utils::globalVariables(".")

test <- function(x){ print(paste0("Test ",x," passed.","\n"))}


is.dnaSequence <- function(seq){
        
        # only A, C, G, T (upper or lower case). The former pattern "[^(A|C|G|T|a|c|g|t)]" also let '(', ')' and
        # '|' pass.
        return(!grepl("[^ACGTacgt]", seq))
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

is.dnds_est_method <- function(method = NULL){
        
        # dNdS estimation methods provided by the KaKs_Calculator 1.2 program
        kaks_calc_methods <- c("MA","MS","NG","LWL","LPB","MLWL","MLPB","YN","MYN","GY","GMYN","ALL","kaks_calc")
        
        provided <- c("Comeron","Li",kaks_calc_methods)
        
        if (is.null(method)){
                message("The following methods are provided: ", paste(provided, collapse = ", "))
                return(FALSE)
        }
        
        return(is.element(method, provided))
}

is.alignment_search_tool <- function(tool = NULL){
        
        provided <-  c("ggsearch","ssearch")
        
        if (is.null(tool)){
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
        # sequences with IUPAC ambiguity codes (N, R, Y, ...) are kept: a codon with such a letter is translated to X
        # and left to the alignment and estimation programs (see ?orthologr_options); only other letters are left out.
        not_acgt <- !vapply(dt[ , seqs], is.dnaSequence, logical(1), USE.NAMES = FALSE)
        if (any(not_acgt) && identical(ambiguous_cds_mode(), "keep")) {
                iupac <- grepl("^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]*$", dt[ , seqs])
                kept <- not_acgt & iupac
                if (any(kept))
                        message(sum(kept), " coding sequence", if (sum(kept) > 1) "s" else "", " with ambiguous bases (e.g. N) ",
                                if (sum(kept) > 1) "are" else "is", " kept (options(orthologr.ambiguous_cds = \"keep\")); ",
                                "codons with such a base are translated to X.")
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
dropped_start <- function() {
        if (!is.null(.orthologr_state$dropped))
                return(function() NULL)
        .orthologr_state$dropped <- list()
        function() {
                res <- .orthologr_state$dropped
                .orthologr_state$dropped <- NULL
                if (length(res) == 0)
                        return(data.frame(file = character(0), id = character(0), reason = character(0),
                                          stringsAsFactors = FALSE))
                res <- unique(do.call(rbind, res))
                rownames(res) <- NULL
                res
        }
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
        fast <- !is.na(seqs) & n > 0 & n %% 3 == 0 & !grepl("[^ACGTacgt]", seqs)
        # ambiguous bases (options(orthologr.ambiguous_cds = "keep")): a codon with such a base becomes X
        fuzzy <- !fast & !is.na(seqs) & n > 0 & n %% 3 == 0 & grepl("^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]*$", seqs)
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

# "drop" (default) or "keep": what happens to coding sequences with IUPAC ambiguity codes (N, R, Y, ...);
# options(orthologr.ambiguous_cds = ).
ambiguous_cds_mode <- function() {
        m <- getOption("orthologr.ambiguous_cds", "drop")
        if (!is.character(m) || length(m) != 1 || !is.element(m, c("drop", "keep")))
                orthologr_abort("options(orthologr.ambiguous_cds = ) must be \"drop\" or \"keep\".", class = "bad_argument")
        m
}
