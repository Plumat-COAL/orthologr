#' @title Read the genome of a given organism
#' @description This function reads an organism specific genome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readDNAStringSet}} function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details The \code{read.genome} function takes a string specifying the path to the genome file
#' of interest as first argument.
#'
#' It is possible to read in different genome file standards such as \emph{fasta} or \emph{fastq}.
#' Genomes stored in fasta files can be downloaded from https://ensemblgenomes.org, https://www.ncbi.nlm.nih.gov/genome, or https://www.ebi.ac.uk, etc.
#'
#' @examples
#' # reading a genome stored in a fasta file
#' Ath.genome <- read.genome(system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                            format = "fasta")
#' head(Ath.genome, 2)
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

read.genome <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","fastq")))
                orthologr_abort(class = "bad_argument", "Please choose a file format that is supported by this function.")
        
        if (!file.exists(file))
                orthologr_abort(class = "bad_input", "The file path you specified does not seem to exist: '", file,"'.", call. = FALSE)
        
        geneids <- seqs <- NULL
        
        # fasta and fastq are read in the same way (Biostrings::readDNAStringSet(format = ...))
        tryCatch({
                
                read_path <- readable_seq_file(file)
                if (!identical(read_path, file))
                        on.exit(unlink(read_path), add = TRUE)
                genome <- Biostrings::readDNAStringSet(filepath = read_path, format = format, ...)
                
                if (length(genome) == 0)
                        orthologr_abort(class = "bad_input", "The file '", file,"' seems to be empty and does not contain any sequences.")
                
                genome_names <- fasta_header_id(names(genome))
                genome.dt <- data.table::data.table(geneids = genome_names,
                                                    seqs = toupper(as.character(genome)))
                
                data.table::setkey(genome.dt,geneids)
                
        }, error = function(e){ orthologr_abort(class = "bad_input", paste0("File ",file, " could not be read properly. \n",
                                            "Please make sure that ",file," contains only DNA sequences and is in ",format," format."))}
        )
        
        return(genome.dt)
}


#' @title Read the proteome of a given organism
#' @description This function reads an organism specific proteome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the proteome.
#' @param format a character string specifying the file format used to store the proteome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readAAStringSet}} function.
#' @author Hajk-Georg Drost
#' @details The \code{read.proteome} function takes a string specifying the path to the proteome file
#' of interest as first argument.
#'
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{fastq}.
#'
#' Proteomes stored in fasta files can be downloaded from https://www.ebi.ac.uk/reference_proteomes, etc.
#'
#' A FASTA file that holds nucleotide sequences rather than proteins (90 \% or more of the first residues are
#' nucleotide codes - A, C, G, T, U and the IUPAC codes N, R, Y, K, M, S, W, B, D, H, V - and fewer than 1 \% are
#' E, F, I, L, P or Q, letters no nucleotide code uses) is rejected
#' with an error of class \code{orthologr_error_bad_input}: DIAMOND, which searches the proteins, can hang on such
#' input instead of stopping with an error. Coding sequences are read with \code{\link{read.cds}}, or with
#' \code{seq_type = "cds"} in the functions that search them. A file of proteins that is rejected (short
#' low-complexity proteins made of A, G, T, N, C and S, for instance) is read with
#' \code{options(orthologr.check_protein = FALSE)} (see \code{\link{orthologr_options}}).
#'
#' @examples
#' # reading a proteome stored in a fasta file
#' Ath.proteome <- read.proteome(system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'                                format = "fasta")
#' head(Ath.proteome, 2)
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

read.proteome <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","fastq")))
                orthologr_abort(class = "bad_argument", "Please choose a file format that is supported by this function.")
        
        if (!file.exists(file))
                orthologr_abort(class = "bad_input", "The file path you specified does not seem to exist: '", file,"'.", call. = FALSE)
        
        geneids <- seqs <- NULL 
        
        read_path <- readable_seq_file(file)
        if (!identical(read_path, file))
                on.exit(unlink(read_path), add = TRUE)
        
        # nucleotide sequences given as proteins: DIAMOND 2.2 hangs on them (makedb, blastp, deepclust) instead of
        # stopping with an error
        if (format == "fasta")
                check_protein_file(read_path, name = file)
        
                tryCatch({
                        
                        proteome <- Biostrings::readAAStringSet(filepath = read_path, format = format, ...)
                        
                        if (length(proteome) == 0)
                                orthologr_abort(class = "bad_input", "The file '", file,"' seems to be empty and does not contain any sequences.")
                        
                        proteome_names <- fasta_header_id(names(proteome))
                        proteome.dt <- data.table::data.table(geneids = proteome_names,
                                                              seqs = toupper(as.character(proteome)))
                        
                        data.table::setkey(proteome.dt, geneids)
                        
                }, error = function(e){ orthologr_abort(class = "bad_input", "File ",file, " could not be read properly.", "\n",
                                                    "Please make sure that ",file," contains only amino acid sequences and is in ",format," format.")}
                )
        
        return(proteome.dt)
}


#' @title Read the CDS of a given organism
#' @description This function reads an organism specific CDS stored in a defined file format.
#' @param file a character string specifying the path to the file storing the CDS.
#' @param format a character string specifying the file format used to store the CDS, e.g. "fasta", "fastq".
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readDNAStringSet}} function.
#' @author Hajk-Georg Drost
#' @details The \code{read.cds} function takes a string specifying the path to the cds file
#' of interest as first argument.
#'
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{fastq}.
#'
#' CDS stored in fasta files can be downloaded from https://www.ensembl.org/info/data/ftp/index.html, etc.
#'
#' @examples
#' # reading a cds file stored in fasta format
#' Ath.cds <- read.cds(system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                     format = "fasta")
#' head(Ath.cds, 2)
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

read.cds <- function(file, format, delete_corrupt_cds = TRUE, ...) {
        if (!is.element(format, c("fasta", "fastq")))
                orthologr_abort(class = "bad_argument", "Please choose a file format that is supported by this function.")
        
        if (!file.exists(file))
                orthologr_abort(class = "bad_input", "The file path you specified does not seem to exist: '", file,"'.", call. = FALSE)
        
        geneids <- seqs <- NULL
        
        # a byte-order mark or old-Mac line ends (CR only): read a normalised copy (messages name 'file')
        read_path <- readable_seq_file(file)
        if (!identical(read_path, file))
                on.exit(unlink(read_path), add = TRUE)
        
        # a protein file given as coding sequences: Biostrings would drop the letters that are no nucleotides
        # with a warning, and the search would fail later with a message that does not name the cause
        if (format == "fasta" && fasta_looks_like_protein(read_path))
                orthologr_abort(class = "bad_input", "'", file, "' does not look like coding DNA: more than 10 % of the first residues are not ",
                     "nucleotide letters. Is this a protein file? dN/dS and the functions that read coding ",
                     "sequences (seq_type = \"cds\") need the coding sequences; for proteins use seq_type = \"protein\" ",
                     "where a function has it.", call. = FALSE)
        
        tryCatch({
                cds_file <-
                        Biostrings::readDNAStringSet(filepath = read_path, format = format, ...)
                
                if (length(cds_file) == 0)
                        orthologr_abort(class = "bad_input", "The file '", file,"' seems to be empty and does not contain any sequences.")
                
                cds_names <-
                        fasta_header_id(names(cds_file))
                
                cds.dt <-
                        data.table::data.table(geneids = cds_names ,
                                               seqs = tolower(as.character(cds_file)))
                
                
                data.table::setkey(cds.dt, geneids)
                
                mod3 <-
                        function(x) {
                                return((nchar(x) %% 3) == 0)
                        }
                
                all_triplets <- as.logical(cds.dt[ , mod3(seqs)])
                
                n_seqs <- nrow(cds.dt)
                
        }, error = function(e) {
                orthologr_abort(
                        class = "bad_input",
                        "File ",
                        file,
                        " could not be read properly.",
                        "\n",
                        "Please make sure that ",
                        file,
                        " contains only CDS sequences and is in ",
                        format,
                        " format."
                )
        })
        
        if (!all(all_triplets)) {
                # inside dNdS() / dNdS_coal() with quiet = TRUE these messages are left to the summary of the run
                # (the sequences are listed in attr(<result>, "dropped_sequences") / "trimmed_sequences")
                quiet_run <- sequence_messages_quiet()
                if (!quiet_run) message(
                        "There seem to be ",
                        length(which(!all_triplets)),
                        " coding sequences in your input dataset which cannot be properly divided in base triplets, because their sequence length cannot be divided by 3."
                )
                # The file with the corrupted sequences is written to the folder given in
                # options(orthologr.corrupted_cds_dir = ), by default to tempdir(). Its name holds the name of the
                # input file and a code of its whole path (corrupted_cds_file()).
                corrupted_dir <- getOption("orthologr.corrupted_cds_dir", NULL)
                if (is.null(corrupted_dir))
                        corrupted_dir <- tempdir()
                if (!dir.exists(corrupted_dir))
                        dir.create(corrupted_dir, recursive = TRUE)
                corrupted_file <- corrupted_cds_file(file, corrupted_dir)
                
                if (!quiet_run) {
                        message(
                                "A fasta file storing all corrupted coding sequences for inspection was generated and stored at '",
                                corrupted_file,
                                "'."
                        )
                        message("\n")
                }
                corrupted_seqs <- as.data.frame(cds.dt[which(!all_triplets)])
                seq_vector <- corrupted_seqs$seqs
                names(seq_vector) <- corrupted_seqs$geneids
                corrupted_seqs_biostrings <- Biostrings::DNAStringSet(seq_vector, use.names = TRUE)
                # (written under a name of its own and then moved: two calls never write into the same file)
                tmp_corrupted <- tempfile("corrupted_", tmpdir = corrupted_dir, fileext = ".fasta")
                Biostrings::writeXStringSet(corrupted_seqs_biostrings, filepath = tmp_corrupted)
                publish_file(tmp_corrupted, corrupted_file)
                
                if (delete_corrupt_cds) {
                        if (!quiet_run) {
                                message(
                                        "You chose option 'delete_corrupt_cds = TRUE', thus corrupted coding sequences were removed.",
                                        "If after consulting the file '",
                                        corrupted_file,
                                        "' you still wish to retain all coding sequences please specify the argument 'delete_corrupt_cds = FALSE'."
                                )
                                message("\n")
                        }
                        dropped_record(cds.dt[which(!all_triplets), geneids],
                                       "a length that is not a multiple of 3 (delete_corrupt_cds = TRUE)", file)
                        return(cds.dt[-which(!all_triplets) , list(geneids, seqs)])
                }
                
                if (!delete_corrupt_cds) {
                        if (!quiet_run) {
                                message(
                                        "You chose option 'delete_corrupt_cds = FALSE', thus corrupted coding sequences were retained for subsequent analyses.")
                                message(
                                        "The following modifications were made to the CDS sequences that were not divisible by 3:")
                                message(
                                        "- If the sequence had 1 residue nucleotide then the last nucleotide of the sequence was removed.")
                                message(
                                        "- If the sequence had 2 residue nucleotides then the last two nucleotides of the sequence were removed.")
                                message(
                                        "Note: the sequences are trimmed at their 3' end. A coding sequence that is incomplete at its 5' end ",
                                        "(e.g. Ensembl 'cds_start_NF' transcripts that do not start in phase 0) is then translated in the wrong ",
                                        "reading frame; such sequences are best removed (delete_corrupt_cds = TRUE) or put into frame before.")
                                message(
                                        "If after consulting the file '",
                                        corrupted_file,
                                        "' you wish to remove all corrupted coding sequences please specify the argument 'delete_corrupt_cds = TRUE'."
                                )
                        }
                        
                        mod3_residue_1 <-
                                function(x) {
                                        return((nchar(x) %% 3) == 1)
                                }
                        
                        mod3_residue_2 <-
                                function(x) {
                                        return((nchar(x) %% 3) == 2)
                                }
                        
                        residue_1 <- cds.dt[ , mod3_residue_1(seqs)]
                        residue_2 <- cds.dt[ , mod3_residue_2(seqs)]
                        trimmed_ids <- cds.dt[which(!all_triplets), geneids]
                        trimmed_nt <- nchar(cds.dt[which(!all_triplets), seqs]) %% 3
                        
                        residue_1_seqs <- as.character(cds.dt[which(residue_1) , seqs])
                        residue_2_seqs <- as.character(cds.dt[which(residue_2) , seqs])
                        
                        residue_1_seqs_vec <- as.character(sapply(residue_1_seqs, function(x) {
                                substr(x, 1, nchar(x) - 1)
                        }))
                
                        residue_2_seqs_vec <- as.character(sapply(residue_2_seqs, function(x) {
                                substr(x, 1, nchar(x) - 2)
                        }))
                        
                        cds.dt[which(residue_1) , seqs := residue_1_seqs_vec]
                        cds.dt[which(residue_2) , seqs := residue_2_seqs_vec]
                        
                        all_triplets_new <- cds.dt[ , mod3(seqs)]
                        
                        if (any(!all_triplets_new)) {
                                orthologr_abort(class = "bad_input", "Not all sequences of '", file, "' could be trimmed to a multiple of 3.")
                        } else if (!quiet_run) {
                                message("All corrupted CDS were trimmed.")
                        }
                        trimmed_record(trimmed_ids, trimmed_nt, file, corrupted_file)
                        
                        n_seqs_new <- nrow(cds.dt)
                        
                        if(!(n_seqs == n_seqs_new))
                                orthologr_abort(class = "bad_input", "After trimming the coding sequences of '", file, "' to a multiple of 3, ",
                                                n_seqs - n_seqs_new, " sequences were missing.")
                        
                        return(cds.dt)
                }
                
        } else {
                return(cds.dt)
        }
        
}






# The ID of a sequence: its FASTA header up to the first white space (space or tab).
fasta_header_id <- function(headers) {
        if (is.null(headers))
                return(character(0))
        sub("[[:space:]].*$", "", as.character(headers))
}

# How often each byte occurs in the sequence lines of the first 2,000 lines of a FASTA file (upper case; a vector
# of 255 counts, indexed by the byte value).
fasta_residue_counts <- function(file) {
        x <- tryCatch(readLines(file, n = 2000, warn = FALSE), error = function(e) character(0))
        x <- x[!startsWith(x, ">")]
        x <- toupper(paste(x, collapse = ""))
        tabulate(as.integer(charToRaw(x)), 255L)
}

# Does a FASTA file hold proteins rather than nucleotide sequences? Judged from its first 2,000 lines: TRUE when
# more than 10 % of the residues (the letters A-Z and '*') are letters that are no IUPAC nucleotide code (E, F, I,
# L, P, Q, ...).
fasta_looks_like_protein <- function(file) {
        n <- fasta_residue_counts(file)
        total <- sum(n[c(65:90, 42L)])
        if (total < 30)
                return(FALSE)
        nucleotide <- sum(n[as.integer(charToRaw("ACGTUNRYKMSWBDHV"))])
        (total - nucleotide) / total > 0.1
}

# The share of A, C, G, T, U and N among the residues (the letters A-Z and '*') of the first 2,000 lines of a FASTA
# file; NA for fewer than 30 residues. Proteins have about 25 % of these letters, nucleotide sequences nearly 100 %.
fasta_nucleotide_share <- function(file) {
        n <- fasta_residue_counts(file)
        total <- sum(n[c(65:90, 42L)])
        if (total < 30)
                return(NA_real_)
        sum(n[as.integer(charToRaw("ACGTUN"))]) / total
}

# Stop with an orthologr_error_bad_input when the FASTA file `file` looks like nucleotide sequences rather than
# proteins: 90 % or more of its first residues are nucleotide codes (A, C, G, T, U and the IUPAC codes N, R, Y, K,
# M, S, W, B, D, H, V: coding sequences with ambiguous bases count as nucleotides), and fewer than 1 %
# are letters that only proteins have (E, F, I, L, P, Q: no IUPAC nucleotide code, so a nucleotide file has none of
# them; a short low-complexity protein rich in A, G, T, N and C can reach 90 %, but one with a few Q, P, L, ... is
# not taken for nucleotides). DIAMOND 2.2 hangs on nucleotide sequences given as proteins (makedb, blastp,
# deepclust) instead of stopping with an error. Used by read.proteome() (every function that reads proteins) and
# deepclust(direct_input = TRUE). options(orthologr.check_protein = FALSE) switches the check off. `name`: the file
# named in the message.
check_protein_file <- function(file, name = file) {
        if (!protein_check_on())
                return(invisible(NA))
        n <- fasta_residue_counts(file)
        total <- sum(n[c(65:90, 42L)])
        if (total < 30)
                return(invisible(NA))
        nt <- sum(n[as.integer(charToRaw("ACGTUNRYKMSWBDHV"))]) / total
        protein_only <- sum(n[as.integer(charToRaw("EFILPQ"))]) / total
        if (nt >= 0.9 && protein_only < 0.01)
                orthologr_abort(class = "bad_input", "'", name, "' does not look like protein sequences: ",
                                floor(100 * nt), " % of the first residues are nucleotide codes (A, C, G, T, U, N or ",
                                "another IUPAC code) and less than 1 % are E, F, I, L, P or Q. ",
                                "Is this a file of coding sequences? Use seq_type = \"cds\" (the coding sequences are ",
                                "translated) where a function has it, or give the protein file. If these are proteins ",
                                "(e.g. short low-complexity ones), options(orthologr.check_protein = FALSE) switches this ",
                                "check off.", call. = FALSE)
        invisible(nt)
}

# Is the check of check_protein_file() on? options(orthologr.check_protein): TRUE (default) or FALSE.
protein_check_on <- function() {
        opt <- getOption("orthologr.check_protein", TRUE)
        if (!is.logical(opt) || length(opt) != 1 || is.na(opt))
                orthologr_abort(class = "bad_argument", "options(orthologr.check_protein = ) must be TRUE or FALSE.", call. = FALSE)
        opt
}

# The file read.cds() writes the coding sequences of `file` that are not a multiple of 3 long to: in `dir`, named
# after the file and a code of its whole path, so that input files of the same name in different folders get
# files of their own, and the same input file always the same name.
corrupted_cds_file <- function(file, dir) {
        key <- tempfile("path_")
        on.exit(unlink(key), add = TRUE)
        writeBin(charToRaw(enc2utf8(normalizePath(file, winslash = "/", mustWork = FALSE))), key)
        code <- substr(unname(tools::md5sum(key)), 1, 8)
        file.path(dir, paste0(basename(file), "_", code, "_corrupted_cds_seqs.fasta"))
}

# A sequence file (FASTA / FASTQ, plain or gzip) as Biostrings can read it. A UTF-8 byte-order mark at the start of
# the file (written by some Windows editors and by Python's 'utf-8-sig') made the file unreadable ("could not be
# read properly"), and old-Mac line ends (CR without LF) made the whole file one line, which reached DIAMOND as one
# sequence. For such a file the mark is removed and every CR becomes LF, in a copy inside orthologr_tmpdir(), whose
# path is returned (the caller removes it); any other file is returned as it is. The sequences read are those of
# the same file without the mark and with LF line ends, as COAL reads them.
readable_seq_file <- function(file) {
        con <- gzfile(file, "rb")
        head <- tryCatch(readBin(con, "raw", n = 1048576), error = function(e) raw(0), finally = close(con))
        bom <- length(head) >= 3 && identical(head[1:3], as.raw(c(0xef, 0xbb, 0xbf)))
        cr_only <- any(head == as.raw(0x0d)) && !any(head == as.raw(0x0a))
        if (!bom && !cr_only)
                return(file)
        con <- gzfile(file, "rb")
        chunks <- list()
        tryCatch(repeat {
                b <- readBin(con, "raw", n = 16777216)
                if (length(b) == 0)
                        break
                chunks[[length(chunks) + 1]] <- b
        }, finally = close(con))
        bytes <- do.call(c, chunks)
        if (bom)
                bytes <- bytes[-(1:3)]
        if (cr_only)
                bytes[bytes == as.raw(0x0d)] <- as.raw(0x0a)
        out <- tempfile("input_", tmpdir = orthologr_tmpdir(), fileext = ".txt")
        writeBin(bytes, out)
        out
}
