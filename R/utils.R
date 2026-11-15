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
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(method, provided))
}


is.multiple_aln_tool <- function(tool = NULL){
        
        provided <- c("clustalw", "clustalo","muscle", "t_coffee", "mafft")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.pairwise_aln_tool <- function(tool = NULL){
        
        provided <- c("NW")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.codon_aln_tool <- function(tool = NULL){
        
        provided <-  c("pal2nal")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.dnds_est_method <- function(method = NULL){
        
        # dNdS estimation methods provided by the KaKs_Calculator program (one list for the whole package)
        provided <- c("Comeron", "Li", kaks_calculator_methods())
        
        if (is.null(method)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(method, provided))
}

is.alignment_search_tool <- function(tool = NULL){
        
        provided <-  c("ggsearch","ssearch")
        
        if (is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
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
#' @import data.table
cds2aa <- function(dt, delete_corrupt_cds = TRUE){
        
        if (!is.data.table(dt))
                stop("Your CDS file was not corretly transformed into a data.table object.")
        
        # define visible bindings for global variables
        seqs <- aa <- geneids <- NULL
        
        # omit empty sequences
        dt <- dt[!(is.na(seqs) | seqs == "")]
        
        if (delete_corrupt_cds) {
                # omit sequences that are not multiples of 3
                dt <- dt[nchar(seqs) %% 3 == 0]
        }
        
        # omit sequences consisting of others than ACGT - and say which
        not_acgt <- !vapply(dt[ , seqs], is.dnaSequence, logical(1), USE.NAMES = FALSE)
        if (any(not_acgt)) {
                report_dropped_cds(dt[not_acgt, geneids])
                dt <- dt[!not_acgt]
        }
        
        # translate cds to protein sequences
        tryCatch({
                if (anyDuplicated(dt[ , geneids]) > 0) {
                        # the same ID more than once: keep what the per-gene translation did
                        dt[ , aa := transl(seqs), by = geneids]
                } else {
                        dt[ , aa := translate_cds(seqs)]
                }
                
        }, error = function(e) {
                stop(
                        "The input coding sequences could not be translated properly to amino acid sequences.",
                        "\n",
                        " Please check whether the input file stores valid coding sequences.",
                        call. = FALSE
                )
        })
        
        return(dt)
}

# Coding sequences with characters other than A, C, G, T (N, IUPAC ambiguity codes, ...) cannot be translated and
# are left out. This used to happen silently; the genes are now named in a message (once per dNdS() call, in
# which the same file is read several times).
report_dropped_cds <- function(ids) {
        if (length(ids) == 0)
                return(invisible(NULL))
        run_dir <- getOption("orthologr.tmpdir", NULL)
        if (!is.null(run_dir)) {
                # inside one dNdS() call: say it once
                key <- paste0("dropped_cds\r", run_dir, "\r", paste(ids, collapse = "\r"))
                if (!is.null(.orthologr_state[[key]]))
                        return(invisible(NULL))
                .orthologr_state[[key]] <- TRUE
        }
        shown <- utils::head(ids, 10)
        message(
                length(ids), " coding sequence", if (length(ids) > 1) "s hold" else " holds",
                " characters other than A, C, G, T (e.g. N or IUPAC ambiguity codes) and ",
                if (length(ids) > 1) "are" else "is", " left out: ",
                paste(shown, collapse = ", "),
                if (length(ids) > length(shown)) paste0(", ... (", length(ids) - length(shown), " more)") else "",
                "."
        )
        invisible(ids)
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
        if (any(fast))
                res[fast] <- as.character(Biostrings::translate(
                        Biostrings::DNAStringSet(seqs[fast]),
                        genetic.code  = Biostrings::GENETIC_CODE,
                        no.init.codon = TRUE
                ))
        if (any(!fast))
                res[!fast] <- vapply(seqs[!fast], transl, character(1), USE.NAMES = FALSE)
        res
}
