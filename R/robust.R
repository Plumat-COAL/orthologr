# Internal helpers that make the dN/dS path robust: one working folder per dNdS() call, exit-status checks of
# external programs, and gene pairs that fail without ending the run.

# The folder in which the dNdS pipeline keeps its intermediate files ('_blast_db', '_alignment', '_calculation').
# dNdS() sets a folder of its own for each call (and removes it afterwards); functions called on their own
# (pairwise_aln(), codon_aln(), diamond(), ...) keep using tempdir() as before.
orthologr_tmpdir <- function() {
        d <- getOption("orthologr.tmpdir", NULL)
        if (is.null(d) || !is.character(d) || length(d) != 1 || is.na(d) || !nzchar(d))
                return(tempdir())
        return(d)
}

# Start a working folder for one pipeline call. Returns a function that restores the previous state and removes
# the folder (unless options(orthologr.keep_tmp = TRUE)); call it in on.exit().
orthologr_tmpdir_start <- function(prefix = "orthologr_") {
        run_dir <- tempfile(prefix)
        dir.create(run_dir, recursive = TRUE)
        old <- options(orthologr.tmpdir = run_dir)
        function() {
                options(old)
                if (isTRUE(getOption("orthologr.keep_tmp", FALSE))) {
                        message("Intermediate files were kept in '", run_dir, "' (options(orthologr.keep_tmp = TRUE)).")
                } else {
                        unlink(run_dir, recursive = TRUE, force = TRUE)
                }
                invisible(run_dir)
        }
}

# The last lines of a log file as one string ("" if there is none).
tail_of_file <- function(file, n = 5) {
        if (is.null(file) || !file.exists(file))
                return("")
        x <- tryCatch(readLines(file, warn = FALSE), error = function(e) character(0))
        x <- x[nzchar(trimws(x))]
        if (length(x) == 0)
                return("")
        # a line that reports an error (PAL2NAL: "#---  ERROR: inconsistency between the following pep and nuc
        # seqs  ---#") and what follows it; otherwise the last lines
        err <- grep("ERROR|[Ee]rror|recursion limit|Cannot|Can't", x)
        if (length(err) > 0)
                x <- x[seq(err[1], min(length(x), err[1] + n - 1))]
        x <- utils::tail(x, n)
        x <- ifelse(nchar(x) > 300, paste0(substr(x, 1, 300), " ..."), x)
        paste(x, collapse = "\n")
}

# Stop with a clear message when an external program failed.
#   status  : what system() / shell() returned
#   tool    : the program's name for the message, e.g. "PAL2NAL"
#   context : what it was run for, e.g. "gene pair 'A' / 'B' (q12)"
#   log     : optional file with the program's stderr
#   output  : optional output file that must exist and be non-empty afterwards
check_external <- function(status, tool, context = NULL, log = NULL, output = NULL) {
        what <- if (is.null(context) || !nzchar(context)) "" else paste0(" for ", context)
        detail <- tail_of_file(log)
        detail <- if (nzchar(detail)) paste0("\n", tool, " said: ", detail) else ""
        # the kind of failure (class of the error): an alignment program, or another program
        kind <- c(if (grepl("PAL2NAL|Clustal|MAFFT|MUSCLE|T-Coffee|t_coffee|COAL", tool, ignore.case = TRUE)) "alignment_failed",
                  "tool_failed")
        if (!is.numeric(status) || length(status) != 1 || is.na(status) || status != 0) {
                missing <- identical(as.integer(status), 127L)
                hint <- if (missing) " (the program was not found: is it installed and on the PATH?)" else ""
                orthologr_abort(tool, " failed", what, ": exit status ", paste(status, collapse = " "), hint, ".", detail,
                                class = if (missing) c("tool_missing", kind) else kind)
        }
        if (!is.null(output) && (!file.exists(output) || is.na(file.size(output)) || file.size(output) == 0)) {
                orthologr_abort(tool, " wrote no output", what, " (exit status 0, but '", output, "' is missing or empty).", detail,
                                class = kind)
        }
        invisible(TRUE)
}

# The sequence IDs of a FASTA file (header up to the first white space), in file order.
fasta_ids <- function(file) {
        x <- tryCatch(readLines(file, warn = FALSE), error = function(e) character(0))
        x <- x[startsWith(x, ">")]
        sub("\\s.*$", "", sub("^>\\s*", "", x))
}

# "gene pair 'A' / 'B' (q12)" for messages, from the two-sequence FASTA file of a pair.
pair_context <- function(file, name = NULL) {
        ids <- fasta_ids(file)
        label <- if (length(ids) >= 2) paste0("gene pair '", ids[1], "' / '", ids[2], "'") else "the alignment"
        if (!is.null(name) && nzchar(name))
                label <- paste0(label, " (", name, ")")
        label
}

# The proteins of a CDS table as read.proteome() returns them (data.table: geneids, seqs; keyed by geneids),
# translated with cds2aa().
cds_protein_table <- function(cds_dt, delete_corrupt_cds = TRUE) {
        geneids <- aa <- NULL
        tmp <- cds2aa(data.table::copy(cds_dt), delete_corrupt_cds = delete_corrupt_cds)
        res <- data.table::data.table(geneids = tmp[ , geneids], seqs = toupper(tmp[ , aa]))
        data.table::setkey(res, geneids)
        return(res)
}

# dN, dS and dN/dS as numbers from seqinr::kaks() values, with what cannot be used set to NA: negative values
# (not computable), values of about 10 (seqinr's marker for saturation), and ratios that are not finite.
li_values <- function(ka, ks) {
        dN <- as.numeric(ka)[1]
        dS <- as.numeric(ks)[1]
        dNdS <- dN / dS
        if (is.na(dN) || dN < 0 || dN > 9.9) {
                dN <- NA_real_
                dNdS <- NA_real_
        }
        if (is.na(dS) || dS < 0 || dS > 9.9) {
                dS <- NA_real_
                dNdS <- NA_real_
        }
        if (!is.finite(dNdS))
                dNdS <- NA_real_
        list(dN = dN, dS = dS, dNdS = dNdS)
}

# 'aa_sub_matrix' must be the name of a substitution matrix, a numeric matrix, or NULL.
check_aa_sub_matrix <- function(aa_sub_matrix) {
        if (!is.null(aa_sub_matrix) &&
            !((is.character(aa_sub_matrix) && length(aa_sub_matrix) == 1 && !is.na(aa_sub_matrix)) ||
              (is.matrix(aa_sub_matrix) && is.numeric(aa_sub_matrix))))
                stop("'aa_sub_matrix' must be the name of a substitution matrix (e.g. \"BLOSUM62\"), a numeric substitution matrix, or NULL.", call. = FALSE)
        invisible(TRUE)
}

# The same table for a CDS file: cds_protein_table(read.cds(file)), taken from what was already read and
# translated in this call (see seq_cache_start()).
cds_protein_table_file <- function(file, delete_corrupt_cds = TRUE) {
        geneids <- aa <- NULL
        tmp <- protein_table(file = file, seq_type = "cds", delete_corrupt_cds = delete_corrupt_cds)
        res <- data.table::data.table(geneids = tmp[ , geneids], seqs = toupper(tmp[ , aa]))
        data.table::setkey(res, geneids)
        return(res)
}
