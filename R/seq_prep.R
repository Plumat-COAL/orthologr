# Internal helpers shared by the BLAST and DIAMOND wrappers: read a sequence file as a protein table and build a
# search database from such a table.

# The sequences of `file` as a data.table with the columns 'geneids', 'seqs' (as read) and 'aa' (protein).
#   seq_type = "cds"     : coding sequences are translated (see cds2aa() for what is left out)
#   seq_type = "protein" : read as they are
protein_table <- function(file, seq_type = "cds", delete_corrupt_cds = TRUE, ...) {
        
        seqs <- aa <- geneids <- NULL
        
        if (seq_type == "cds") {
                dt <- data.table::copy(read.cds(file = file, format = "fasta", delete_corrupt_cds = delete_corrupt_cds, ...))
                
                if (!is.data.table(dt))
                        stop("Your CDS file was not corretly transformed into a data.table object.")
                
                # omit empty sequences
                dt <- dt[ , .SD[sapply(seqs, function(x) {return(!(is.na(x) || x == ""))})]]
                
                if (delete_corrupt_cds) {
                        # omit sequences that are not multiples of 3
                        dt <- dt[ , .SD[sapply(seqs, function(x) {return(nchar(x) %% 3 == 0)})]]
                }
                
                # omit sequences consisting of others than ACGT
                dt <- dt[ , .SD[sapply(seqs, is.dnaSequence)]]
                
                # translate cds to protein sequences
                tryCatch({
                        dt[ , aa := transl(seqs), by = geneids]
                }, error = function(e) {
                        stop("The input coding sequences could not be translated properly to amino acid sequences.",
                             "\n", " Please check whether ", file, " stores valid coding sequences.", call. = FALSE)
                })
        }
        
        if (seq_type == "protein") {
                dt <- data.table::copy(read.proteome(file = file, format = "fasta", ...))
                data.table::setnames(dt, old = "seqs", new = "aa")
        }
        
        if (seq_type == "dna") {
                dt <- data.table::copy(read.genome(file = file, format = "fasta", ...))
        }
        
        return(dt)
}

# Write the 'aa' column of a protein table as FASTA.
write_protein_fasta <- function(dt, file) {
        aa <- geneids <- NULL
        seqinr::write.fasta(
                sequences = as.list(dt[ , aa]),
                names     = dt[ , geneids],
                nbchar    = 80,
                open      = "w",
                file.out  = file
        )
        invisible(file)
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
        write_protein_fasta(dt, db_file)
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
