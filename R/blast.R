#' @title Perform a BLAST+ search
#' @description This function performs a BLAST+ search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}.
#' Default is \code{format} = \code{"fasta"}.
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, option is: \code{blast_algorithm} = \code{"blastp"}
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param remote a boolean value specifying whether a remote BLAST search shall be performed.
#' In case \code{remote} = \code{TRUE}, please specify the \code{db} argument. This feature is very experimental,
#' since a query of only a few genes against NCBI nr database, can consume a lot of time and might cause
#' response delay crashes in R.
#' @param db a character string specifying the NCBI data base that shall be queried using remote BLAST.
#' This parameter must be specified when \code{remote} = \code{TRUE} and is \code{NULL} by default.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run BLAST searches.
#' @param blast_params a character string listing the input parameters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param clean_folders a boolean value specifying whether all internal folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the BLAST output shall be stored. E.g. \code{save.output} = \code{getwd()}
#' to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @details This function provides a fast communication between R and BLAST+. It is mainly used as internal functions
#' such as \code{\link{blast_best}} and \code{\link{blast_rec}} but can also be used to perform simple BLAST computations.
#'
#' When using \code{remote} = \code{TRUE}, make sure you specify the \code{db} argument.
#' The following databases can be chosen:
#'
#' \code{db}
#' \itemize{
#' \item "nr"
#' \item "plaza"
#' }
#'
#' Note: When working with remote BLAST, make sure you don't submit too large jobs due to the
#' BLAST query conventions! All in all this functionality is still very experimental and can cause
#' problems due to time out errors when submitting too large queries!
#'
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @references
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#'
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schaeffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' https://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
#'
#' https://blast.ncbi.nlm.nih.gov/Blast.cgi
#' @examples \dontrun{
#' # performing a BLAST search using blastp (default)
#' blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'
#' # performing a BLAST search using blastp (default) using amino acid sequences as input file
#' blast(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'       seq_type     = "protein")
#'
#'
#' # save the BLAST output table in your current working directory
#' blast(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'       seq_type     = "protein",
#'       save.output  = getwd())
#'
#' # in case you are working with a multicore machine, you can also run parallel
#' # BLAST computations using the comp_cores parameter: here with 2 cores
#' blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       comp_cores   = 2)
#'
#'
#'  # running blastp using additional parameters
#'  blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'        blast_params = "-max_target_seqs 1")
#'
#'
#' # running blastp using additional parameters and an external blastp path
#'  blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'        blast_params = "-max_target_seqs 1", path = "path/to/blastp/")
#' }
#'
#' @return A data.table storing the BLAST hit table returned by BLAST.
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{set_blast}}
#' @export
blast <- function(query_file,
                  subject_file,
                  seq_type        = "cds",
                  format          = "fasta",
                  blast_algorithm = "blastp",
                  eval            = "1E-5",
                  max.target.seqs = 10000,
                  delete_corrupt_cds = TRUE,
                  remote          = FALSE,
                  db              = NULL,
                  path            = NULL,
                  comp_cores      = 1,
                  blast_params    = NULL,
                  clean_folders   = FALSE,
                  save.output     = NULL) {
        
        if (!is.element(blast_algorithm, c("blastp")))
                stop(
                        "Please choose a valid BLAST mode. Only 'blastp' is available for this function.",
                        call. = FALSE
                )
        
        if (remote & is.null(db))
                stop(
                        "To use the remote option of blast() please specify the 'db' argument, e.g. db = 'nr'",
                        call. = FALSE
                )
        
        if (!is.null(db)) {
                if (!is.element(db, c("nr", "plaza")))
                        stop("Please choose a database that is supported by remote BLAST.",
                             call. = FALSE)
        }
        
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        
        message("Running ", blast_version(path = path), " ...")
        
        # in case one tries to use more cores than are available
        check_comp_cores(comp_cores)
        
        # every call works in a folder of its own, which is removed when the call ends: a table left by an
        # earlier or failed call cannot be returned as the result of this one
        work <- call_dir_start("blast_")
        on.exit(call_dir_end(work), add = TRUE)
        
        # initialize the BLAST search
        query.dt <- set_blast(
                file     = query_file,
                seq_type = seq_type,
                format   = format,
                delete_corrupt_cds = delete_corrupt_cds
        )[[1]]
        
        filename <- basename(query_file)
        input    <- file.path(work, "query.fasta")
        output   <- file.path(work, "blastresult.tsv")
        
        # where the hit table is kept after the call (as before); a table of an earlier call is removed now
        db_dir <- file.path(orthologr_tmpdir(), "_blast_db")
        kept   <- file.path(db_dir, paste0("blastresult_", filename, ".csv"))
        if (file.exists(kept))
                unlink(kept)
        context  <- paste0("query '", filename, "'")
        
        tryCatch({
                write_protein_fasta(query.dt, input)
        }, error = function(e) {
                stop("The query sequences of '", query_file, "' could not be written to '", input, "': ",
                     conditionMessage(e), call. = FALSE)
        })
        
        if (remote) {
                run_tool(
                        "blastp",
                        args = list("-db", db, "-query", input, "-remote", "-evalue", eval,
                                    "-out", output, "-outfmt", "6"),
                        path = path, tool = blast_algorithm, context = context
                )
        } else {
                # make a BLASTable database of the subject
                subject.dt <- set_blast(
                        file     = subject_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds
                )[[1]]
                database <- make_blast_db(subject.dt, db_file = file.path(work, "subject_db.fasta"), path = path)
                
                # additional blast parameters can be found here:
                # https://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
                run_tool(
                        "blastp",
                        args = list(
                                "-db", blast_db_arg(database),
                                "-query", input,
                                "-evalue", eval,
                                "-max_target_seqs", format(max.target.seqs, scientific = FALSE),
                                "-out", output,
                                "-num_threads", comp_cores,
                                if (!is.null(blast_params)) I(blast_params),
                                "-outfmt",
                                "6 qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp sstart send slen evalue bitscore score"
                        ),
                        path = path, tool = blast_algorithm, context = context
                )
        }
        
        # blastp must have written its table (an empty file when there is no hit)
        if (!file.exists(output))
                stop(blast_algorithm, " wrote no hit table for ", context, ".", call. = FALSE)
        
        blast_table_names <-
                c(
                        "query_id",
                        "subject_id",
                        "perc_identity",
                        "num_ident_matches",
                        "alig_length",
                        "mismatches",
                        "gap_openings",
                        "n_gaps",
                        "pos_match",
                        "ppos",
                        "q_start",
                        "q_end",
                        "q_len",
                        "qcov",
                        "qcovhsp",
                        "s_start",
                        "s_end",
                        "s_len",
                        "evalue",
                        "bit_score",
                        "score_raw"
                )
        
        hit_table <- read_hit_table(output, blast_table_names, tool = "BLAST")
        
        # save the BLAST output file to path save.output
        save_output_copy(output, save.output, paste0("blastresult_", filename, ".csv"))
        
        if (clean_folders) {
                clean_all_folders(db_dir)
        } else {
                if (!file.exists(db_dir))
                        dir.create(db_dir, recursive = TRUE)
                publish_file(output, kept)
        }
        
        return(hit_table)
}

# "blastp: 2.17.0+" - the BLAST+ version line; stops when blastp cannot be started. Asked once per session.
blast_version <- function(path = NULL) {
        tool_available(
                "blastp", version_args = "-version", pattern = "blast", path = path,
                not_found = "It seems like you don't have BLAST installed locally on your machine or the PATH variable to the BLAST program is not set correctly."
        )
}

# Read a tabular hit table (BLAST / DIAMOND output format 6 with the 21 columns used in this package).
read_hit_table <- function(file, col_names, tool) {
        # define the colClasses for faster file streaming
        col_Classes <- c(rep("character", 2),
                         "double",
                         rep("integer", 6),
                         "double",
                         rep("integer", 3),
                         rep("double", 2),
                         rep("integer", 3),
                         rep("double", 3))
        
        tryCatch({
                hit_table <- data.table::fread(
                        input      = file,
                        sep        = "\t",
                        header     = FALSE,
                        colClasses = col_Classes
                )
                
                data.table::setnames(
                        x   = hit_table,
                        old = paste0("V", 1:length(col_names)),
                        new = col_names
                )
                
                data.table::setkeyv(hit_table, c("query_id", "subject_id"))
                
                tibble::as_tibble(hit_table)
                
        }, error = function(e) {
                stop(
                        "The hit table written by ", tool, " could not be read",
                        if (isTRUE(file.size(file) == 0)) " (it is empty: no hits were found)" else "",
                        ". Please check the input files and the parameters handed to ", tool, ".",
                        "\n",
                        "Error: ", conditionMessage(e),
                        call. = FALSE
                )
        })
}
