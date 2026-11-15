# Internal helpers shared by the BLAST and DIAMOND wrappers: read a sequence file as a protein table and build a
# search database from such a table.

# The sequences of `file` as a data.table with the columns 'geneids', 'seqs' (as read) and 'aa' (protein).
#   seq_type = "cds"     : coding sequences are translated (see cds2aa() for what is left out)
#   seq_type = "protein" : read as they are
protein_table <- function(file, seq_type = "cds", delete_corrupt_cds = TRUE, ...) {
        
        # inside dNdS(), diamond_rec() and blast_rec() every file is read and translated once
        if (...length() == 0 && seq_cache_active()) {
                key <- seq_cache_key("protein_table", file, seq_type, delete_corrupt_cds)
                hit <- seq_cache_get(key)
                if (is.null(hit)) {
                        hit <- protein_table_read(file, seq_type = seq_type, delete_corrupt_cds = delete_corrupt_cds)
                        seq_cache_set(key, hit)
                }
                return(data.table::copy(hit))
        }
        
        protein_table_read(file, seq_type = seq_type, delete_corrupt_cds = delete_corrupt_cds, ...)
}

protein_table_read <- function(file, seq_type = "cds", delete_corrupt_cds = TRUE, ...) {
        
        seqs <- aa <- geneids <- NULL
        
        # the file that is being read, for the messages about sequences that are left out
        old_file <- .orthologr_state$current_file
        .orthologr_state$current_file <- file
        on.exit(.orthologr_state$current_file <- old_file, add = TRUE)
        
        if (seq_type == "cds") {
                dt <- if (...length() == 0)
                        read_cds_once(file, delete_corrupt_cds = delete_corrupt_cds)
                else
                        data.table::copy(read.cds(file = file, format = "fasta", delete_corrupt_cds = delete_corrupt_cds, ...))
                dt <- tryCatch(
                        cds2aa(dt, delete_corrupt_cds = delete_corrupt_cds),
                        error = function(e) {
                                # what went wrong is kept (e.g. duplicated sequence IDs are named)
                                if (grepl("more than once", conditionMessage(e), fixed = TRUE))
                                        stop(conditionMessage(e), call. = FALSE)
                                stop(sub("\n Please check whether the input file stores valid coding sequences.", "",
                                         conditionMessage(e), fixed = TRUE),
                                     "\n", " Please check whether ", file,
                                     " stores valid coding sequences.", call. = FALSE)
                        }
                )
        }
        
        if (seq_type == "protein") {
                dt <- data.table::copy(read.proteome(file = file, format = "fasta", ...))
                check_duplicate_ids(dt[ , geneids], what = "protein sequences", file = file)
                data.table::setnames(dt, old = "seqs", new = "aa")
        }
        
        if (seq_type == "dna") {
                dt <- data.table::copy(read.genome(file = file, format = "fasta", ...))
        }
        
        return(dt)
}

# Write the 'aa' column of a protein table as FASTA: '>id' and the sequence on one line, the file
# seqinr::write.fasta(as.list(aa), names = ids) wrote, without a loop over the sequences in R.
write_protein_fasta <- function(dt, file) {
        aa <- geneids <- NULL
        seqs <- as.character(dt[ , aa])
        ids  <- as.character(dt[ , geneids])
        writeLines(as.vector(rbind(paste0(">", ids), seqs)), file)
        invisible(file)
}

# ---- every input file is read and translated once per call ----------------------------------------------------
# dNdS(), diamond_rec() and blast_rec() need the same two files several times (query and subject of both search
# directions, then the coding sequences and proteins of the gene pairs). While one of these functions runs, what
# was read is kept (keyed by path, size and modification time of the file and by the arguments) and handed out as
# a copy.

seq_cache_active <- function() !is.null(.orthologr_state$seq_cache)

# Start keeping; returns the function that ends it (a no-op when an outer call keeps already).
seq_cache_start <- function() {
        if (seq_cache_active())
                return(function() invisible(NULL))
        .orthologr_state$seq_cache <- new.env(parent = emptyenv())
        function() {
                .orthologr_state$seq_cache <- NULL
                invisible(NULL)
        }
}

seq_cache_key <- function(what, file, ...) {
        info <- file.info(file)
        paste(what, normalizePath(file, mustWork = FALSE), info$size, as.numeric(info$mtime), ..., sep = "\r")
}

seq_cache_get <- function(key) {
        if (!seq_cache_active())
                return(NULL)
        .orthologr_state$seq_cache[[key]]
}

seq_cache_set <- function(key, value) {
        if (seq_cache_active())
                assign(key, value, envir = .orthologr_state$seq_cache)
        invisible(value)
}

# read.cds(file, format = "fasta") - read once per call, returned as a copy.
read_cds_once <- function(file, delete_corrupt_cds = TRUE) {
        if (!seq_cache_active())
                return(data.table::copy(read.cds(file = file, format = "fasta", delete_corrupt_cds = delete_corrupt_cds)))
        key <- seq_cache_key("read.cds", file, delete_corrupt_cds)
        hit <- seq_cache_get(key)
        if (is.null(hit)) {
                hit <- read.cds(file = file, format = "fasta", delete_corrupt_cds = delete_corrupt_cds)
                seq_cache_set(key, hit)
        }
        data.table::copy(hit)
}

# All files of a database called `db_file` (the FASTA file and what makeblastdb / diamond makedb put next to it).
db_files <- function(db_file) {
        d <- dirname(db_file)
        if (!dir.exists(d))
                return(character(0))
        f <- list.files(d, full.names = TRUE, all.files = TRUE, no.. = TRUE)
        f[startsWith(basename(f), basename(db_file))]
}

# Build a BLAST database `db_file` (a FASTA file is written there, the index files next to it).
make_blast_db <- function(dt, db_file, makedb_type = "protein", path = NULL) {
        unlink(db_files(db_file))
        if (makedb_type == "protein") {
                write_protein_fasta(dt, db_file)
        } else {
                # a nucleotide database holds the sequences as they were read ('seqs')
                seqs <- geneids <- NULL
                seqinr::write.fasta(sequences = as.list(dt[ , seqs]), names = dt[ , geneids], nbchar = 80,
                                    open = "w", file.out = db_file)
        }
        run_tool(
                "makeblastdb",
                args = list("-in", blast_db_arg(db_file), "-input_type", "fasta",
                            "-dbtype", if (makedb_type == "protein") "prot" else "nucl", "-hash_index"),
                path    = path,
                tool    = "makeblastdb",
                context = paste0("'", basename(db_file), "'")
        )
        if (length(db_files(db_file)) < 2)
                stop("makeblastdb wrote no database for '", db_file, "'.", call. = FALSE)
        invisible(db_file)
}

# BLAST+ splits the value of -db / -in at spaces: a path with spaces needs quotes of its own inside the argument.
blast_db_arg <- function(db_file) {
        if (grepl(" ", db_file, fixed = TRUE))
                return(paste0('"', db_file, '"'))
        db_file
}

# Build a DIAMOND database `<db_file>.dmnd` (the FASTA file is written to `db_file`).
make_diamond_db <- function(dt, db_file, path = NULL, comp_cores = 1, quiet = TRUE) {
        unlink(db_files(db_file))
        write_protein_fasta(dt, db_file)
        diamond_makedb(in_file = db_file, db_file = db_file, path = path, comp_cores = comp_cores, quiet = quiet)
}

diamond_makedb <- function(in_file, db_file, path = NULL, comp_cores = 1, quiet = TRUE) {
        run_tool(
                "diamond",
                args = list("makedb", "--in", in_file, "--db", db_file, "--threads", comp_cores,
                            if (quiet) "--quiet"),
                path    = path,
                tool    = "diamond makedb",
                context = paste0("'", basename(in_file), "'"),
                output  = paste0(db_file, ".dmnd"),
                echo    = !quiet
        )
        invisible(db_file)
}

# "x min." / "x sec." for the progress messages.
format_elapsed <- function(seconds) {
        if (seconds >= 60) paste0(round(seconds / 60, 2), " min.") else paste0(round(seconds, 2), " sec.")
}
