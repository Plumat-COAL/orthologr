# Integration of COAL (https://github.com/Plumat-COAL/COAL): codon alignments and Ka/Ks estimation without Perl
# or KaKs_Calculator, and a batched dN/dS path that hands every gene pair to COAL in one call.
#
# COAL is optional. It is used in process through the R package 'coalr' when that is installed, otherwise through
# the 'coal-cpp' program (part of the Python package 'coal-align').

# The Ka/Ks methods COAL computes (KaKs_Calculator's names).
coal_methods <- function() {
        c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN")
}

coal_install_message <- function() {
        paste0(
                "COAL is not available. Install its R package (recommended; runs inside R, no other program needed):\n",
                "    remotes::install_github(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\")\n",
                "or put COAL's 'coal-cpp' program on the PATH (pip install coal-align; ",
                "or set options(orthologr.coal_path = \"/path/to/coal-cpp\"))."
        )
}

# Which COAL to use: "coalr" (the R package, in process) or "binary" (the coal-cpp program).
# options(orthologr.coal_backend = "coalr" | "binary") forces one; 'need' = "coalr" asks for the package.
# Returns list(type = , path = ) or stops with the installation message.
coal_backend <- function(need = NULL) {
        forced <- getOption("orthologr.coal_backend", NULL)
        if (!is.null(forced) && !is.element(forced, c("coalr", "binary", "none")))
                stop("options(orthologr.coal_backend) must be \"coalr\" or \"binary\".", call. = FALSE)
        if (identical(forced, "none"))
                stop(coal_install_message(), call. = FALSE)

        if (!identical(forced, "binary") && requireNamespace("coalr", quietly = TRUE))
                return(list(type = "coalr", path = NULL))

        if (identical(forced, "coalr") || identical(need, "coalr"))
                stop(
                        "This needs the R package 'coalr' (COAL's C++ engine inside R). Install it with:\n",
                        "    remotes::install_github(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\")",
                        call. = FALSE
                )

        exe <- coal_binary()
        if (is.null(exe))
                stop(coal_install_message(), call. = FALSE)
        return(list(type = "binary", path = exe))
}

# The coal-cpp program: options(orthologr.coal_path), the environment variable COAL_CPP, then the PATH.
# NULL when there is none.
coal_binary <- function() {
        candidates <- c(getOption("orthologr.coal_path", ""), Sys.getenv("COAL_CPP", ""),
                        unname(Sys.which("coal-cpp")))
        for (exe in candidates) {
                if (is.character(exe) && length(exe) == 1 && nzchar(exe) && file.exists(exe) && !dir.exists(exe))
                        return(normalizePath(exe))
        }
        return(NULL)
}

# Is COAL available (as 'coalr' or as the coal-cpp program)?
is_installed_coal <- function() {
        isTRUE(tryCatch({
                coal_backend()
                TRUE
        }, error = function(e) FALSE))
}

# Run coal-cpp with 'args' (a list; quoted by run_tool()), stdout to 'stdout' (a file) when given. Stops when it
# fails or when 'output' was not written. Started through run_tool() like every other external program.
coal_run_binary <- function(exe, args, context, output, stdout = NULL) {
        run_tool(
                basename(exe),
                args    = args,
                path    = dirname(exe),
                stdout  = stdout,
                tool    = "COAL (coal-cpp)",
                context = context,
                output  = output
        )
        invisible(TRUE)
}

# Read a FASTA file as a named character vector (ids up to the first white space; sequences as written).
read_fasta_strings <- function(file) {
        x <- seqinr::read.fasta(file, as.string = TRUE, set.attributes = FALSE, forceDNAtolower = FALSE,
                                seqtype = "AA", strip.desc = TRUE)
        res <- vapply(x, function(s) as.character(s)[1], character(1))
        names(res) <- sub("\\s.*$", "", names(x))
        return(res)
}

# Write a named character vector as FASTA (one line per sequence) or in CLUSTAL format.
write_alignment <- function(seqs, file, format = "fasta") {
        if (format == "fasta") {
                writeLines(as.vector(rbind(paste0(">", names(seqs)), unname(seqs))), file)
                return(invisible(file))
        }
        width <- 60
        n <- max(nchar(seqs))
        pad <- max(nchar(names(seqs))) + 8   # as PAL2NAL writes it
        lines <- c("CLUSTAL W multiple sequence alignment", "")
        for (start in seq(1, max(n, 1), by = width)) {
                lines <- c(lines,
                           paste0(formatC(names(seqs), width = -pad), substr(seqs, start, start + width - 1)),
                           "")
        }
        writeLines(lines, file)
        invisible(file)
}

# The codon alignment of a protein alignment file and a CDS file with COAL, written to 'file.out'.
coal_codon_aln_file <- function(file_aln, file_nuc, file.out, format = "fasta", aln_format = NULL, context = NULL) {
        backend <- coal_backend()
        if (is.null(context))
                context <- pair_context(file_nuc)
        # as for PAL2NAL: an alignment left by an earlier call is removed, the new one is written under a name of
        # its own and becomes 'file.out' only when it is complete
        if (file.exists(file.out))
                unlink(file.out)
        tmp.out <- tempfile("aln_", tmpdir = dirname(file.out), fileext = ".tmp")
        on.exit(unlink(tmp.out), add = TRUE)

        if (backend$type == "binary") {
                # coal-cpp takes PAL2NAL's command line
                coal_run_binary(backend$path,
                                list(file_aln, file_nuc, "-output", format, "-o", tmp.out),
                                context = context, output = tmp.out)
                publish_file(tmp.out, file.out)
                return(invisible(file.out))
        }

        first <- readLines(file_aln, n = 1, warn = FALSE)
        if (length(first) > 0 && startsWith(first, ">")) {
                rows <- read_fasta_strings(file_aln)
        } else {
                aln <- seqinr::read.alignment(file_aln, format = "clustal", forceToLower = FALSE)
                rows <- structure(toupper(unlist(aln$seq)), names = sub("\\s.*$", "", aln$nam))
        }
        cds <- read_fasta_strings(file_nuc)
        # PAL2NAL pairs the two files by position; so does orthologr (both files hold the pair in the same order)
        if (length(rows) != length(cds))
                stop("COAL: the protein alignment holds ", length(rows), " sequences and the CDS file ", length(cds),
                     " for ", context, ".", call. = FALSE)
        codon <- tryCatch(
                coalr::coal_codon_align(unname(rows), unname(cds)),
                error = function(e) stop("COAL could not build the codon alignment for ", context, ": ",
                                         conditionMessage(e), call. = FALSE)
        )
        names(codon) <- names(cds)
        write_alignment(codon, tmp.out, format = format)
        publish_file(tmp.out, file.out)
        invisible(file.out)
}

# Ka/Ks of a two-sequence codon alignment file with COAL: the table substitutionrate() returns for the
# KaKs_Calculator methods (query_id, subject_id, dN, dS, dNdS, method).
coal_kaks_file <- function(file, est.method, subst_name = NULL) {
        query_id <- NULL
        if (!is.element(est.method, coal_methods()))
                stop("COAL computes the methods ", paste(coal_methods(), collapse = ", "), "; '", est.method,
                     "' needs kaks_tool = \"kaks_calculator\".", call. = FALSE)
        backend <- coal_backend()
        context <- pair_context(file, subst_name)
        seqs <- read_fasta_strings(file)
        if (length(seqs) != 2)
                stop("The codon alignment of ", context, " must hold two sequences, but '", file, "' holds ",
                     length(seqs), ".", call. = FALSE)

        if (backend$type == "coalr") {
                k <- coalr::coal_kaks(seqs[[1]], seqs[[2]], method = est.method)
                if (identical(k$status[1], "error"))
                        stop("COAL could not estimate dN/dS for ", context, ": ", k$note[1], call. = FALSE)
                res <- data.table::data.table(query_id = names(seqs)[1], subject_id = names(seqs)[2],
                                              dN = k$dN[1], dS = k$dS[1], dNdS = k$dNdS[1], method = k$method[1])
        } else {
                calc_dir <- file.path(orthologr_tmpdir(), "_calculation")
                if (!file.exists(calc_dir))
                        dir.create(calc_dir, recursive = TRUE)
                # input and output of this call have names of their own and are removed afterwards
                axt_file <- tempfile(paste0(safe_file_name(basename(file)), "_"), tmpdir = calc_dir, fileext = ".axt")
                out_file <- paste0(axt_file, ".coal.tsv")
                on.exit(unlink(c(axt_file, out_file)), add = TRUE)
                writeLines(c(if (is.null(subst_name)) "pair" else subst_name, seqs[[1]], seqs[[2]], ""), axt_file)
                coal_run_binary(backend$path,
                                list("dnds", axt_file, "-o", out_file, "-m", est.method, "-j", "1", "-q",
                                     "--digits", "17"),
                                context = context, output = out_file)
                k <- utils::read.delim(out_file, sep = "\t", quote = "", na.strings = "NA", stringsAsFactors = FALSE)
                if (nrow(k) < 1)
                        stop("COAL returned no result for ", context, ".", call. = FALSE)
                if (identical(k$status[1], "error"))
                        stop("COAL could not estimate dN/dS for ", context, ": ", k$note[1], call. = FALSE)
                res <- data.table::data.table(query_id = names(seqs)[1], subject_id = names(seqs)[2],
                                              dN = as.numeric(k$dN[1]), dS = as.numeric(k$dS[1]),
                                              dNdS = as.numeric(k$dNdS[1]), method = as.character(k$method[1]))
        }
        data.table::setkey(res, query_id)
        return(res)
}

# A CDS table of read.cds() as the named character vector COAL takes.
cds_vector <- function(cds_dt) {
        geneids <- seqs <- NULL
        structure(cds_dt[ , seqs], names = cds_dt[ , geneids])
}

#' @title Compute dNdS values for two organisms with the COAL engine (batched)
#' @description The fast path of \code{\link{dNdS}}: orthologs are detected with \pkg{orthologr}'s own DIAMOND or
#' BLAST best (reciprocal) hit search, and all orthologous gene pairs are then handed to COAL
#' (\url{https://github.com/Plumat-COAL/COAL}) in one call. COAL translates the coding sequences, aligns each pair
#' of proteins globally (Needleman-Wunsch, BLOSUM62), threads the codons onto the protein alignment and estimates
#' dN, dS and dN/dS, inside the R process and on \code{comp_cores} threads. No file is written and no program is
#' started per gene pair (no Perl / PAL2NAL, no KaKs_Calculator).
#'
#' The result is the table \code{\link{dNdS}} returns: the same columns and types, joined with the hit table.
#' \code{dNdS(..., engine = "coal")} calls this function.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param aligner a character string specifying the sequence aligner used to detect orthologs:
#' \code{"diamond"} (Default) or \code{"blast"}.
#' @param sensitivity_mode the DIAMOND sensitivity mode, see \code{\link{dNdS}}.
#' @param aligner_path a character string specifying the path to the DIAMOND or BLAST program (in case you don't use the default path).
#' @param seq_type a character string specifying the sequence type stored in the input files. Only \code{"cds"} is supported.
#' @param format a character string specifying the file format of the sequence files. Default is \code{format = "fasta"}.
#' @param ortho_detection the orthology inference method: \code{"RBH"} (Default; reciprocal best hit) or \code{"BH"} (best hit).
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input files.
#' @param aligner_params a character string specifying additional parameters that shall be passed to DIAMOND or BLAST.
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND or BLAST hit detection.
#' @param dnds_est.method the dN/dS estimation method(s) computed by COAL: \code{"YN"} (Default; Yang and Nielsen 2000),
#' \code{"NG"}, \code{"LWL"}, \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"} or \code{"MYN"}. The values agree with
#' KaKs_Calculator 2.0's for the same codon alignment. Several methods can be given; the result then has one row
#' per gene pair and method.
#' @param aa_aln_scoring how the pairwise protein alignment is scored:
#' \itemize{
#' \item \code{"pwalign"} (Default): the scoring of \code{\link{pairwise_aln}} (BLOSUM62, gap opening 10, gap
#' extension 4, end gaps penalised), so that the alignments are those \code{dNdS()} computes itself;
#' \item \code{"needle"}: the defaults of EMBOSS needle and of COAL (BLOSUM62, gap opening 10, gap extension 0.5,
#' end gaps free).
#' }
#' @param aa_sub_matrix the substitution matrix: \code{"BLOSUM62"} (Default), \code{"BLOSUM45"}, \code{"BLOSUM80"} or \code{"PAM250"}.
#' @param comp_cores the number of threads COAL (and the DIAMOND or BLAST search) uses.
#' @param quiet a logical value specifying whether the output of the search tool shall be printed to the console.
#' @param print_citation a logical value indicating whether or not the citation message shall be printed.
#' @details
#' Differences from \code{dNdS()} with \code{engine = "orthologr"}:
#' \itemize{
#' \item Comeron's and Li's methods are not available (COAL implements the KaKs_Calculator methods listed above).
#' \item A dN or dS of exactly 0 is returned as 0 (KaKs_Calculator prints values below 1e-6 as \code{NA}); a
#' saturated estimate is \code{NA} (KaKs_Calculator prints a capped value).
#' \item Codon columns with a gap, an ambiguous base or an internal stop codon are skipped; a terminal stop codon is
#' not scored.
#' \item A gene pair that cannot be computed is returned with \code{NA} values and is listed in
#' \code{attr(result, "failed_pairs")}, as in \code{dNdS()}.
#' }
#' COAL's complete table for the gene pairs (standard errors, numbers of sites and substitutions, codons compared,
#' protein identity, status and note) is attached as \code{attr(result, "coal")}.
#'
#' This function needs the R package \pkg{coalr}:
#' \code{remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")}.
#' @author Chloe Newton
#' @return A tibble with the columns of \code{\link{dNdS}}: \code{query_id}, \code{subject_id}, \code{dN}, \code{dS},
#' \code{dNdS}, \code{method} and the columns of the hit table.
#' @examples \dontrun{
#' dNdS_coal(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           dnds_est.method = "YN", comp_cores = 2)
#'
#' # the same through dNdS()
#' dNdS(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      dnds_est.method = "YN", engine = "coal")
#' }
#' @seealso \code{\link{dNdS}}, \code{\link{dNdS_clusters}}, \code{\link{diamond_rec}}
#' @export
dNdS_coal <- function(query_file,
                      subject_file,
                      aligner            = "diamond",
                      sensitivity_mode   = "fast",
                      aligner_path       = NULL,
                      seq_type           = "cds",
                      format             = "fasta",
                      ortho_detection    = "RBH",
                      delete_corrupt_cds = FALSE,
                      aligner_params     = NULL,
                      eval               = "1E-5",
                      dnds_est.method    = "YN",
                      aa_aln_scoring     = "pwalign",
                      aa_sub_matrix      = "BLOSUM62",
                      comp_cores         = 1,
                      quiet              = TRUE,
                      print_citation     = TRUE) {

        if (!is.element(aligner, c("diamond", "blast")))
                stop("engine = \"coal\" detects orthologs with aligner = \"diamond\" or \"blast\".", call. = FALSE)
        if (!is.element(ortho_detection, c("BH", "RBH")))
                stop("engine = \"coal\" supports ortho_detection = \"RBH\" or \"BH\".", call. = FALSE)
        if (!identical(seq_type, "cds"))
                stop("dN/dS estimation needs coding sequences: seq_type = \"cds\".", call. = FALSE)
        check_choice(format, "fasta", "format", hint = "dNdS_coal() reads coding sequences from FASTA files.")
        if (aligner == "diamond")
                check_sensitivity_mode(sensitivity_mode)
        comp_cores <- check_comp_cores(comp_cores)
        if (!is.character(dnds_est.method) || length(dnds_est.method) < 1 || !all(is.element(dnds_est.method, coal_methods())))
                stop(
                        "engine = \"coal\" computes dnds_est.method = ", paste0("\"", coal_methods(), "\"", collapse = ", "),
                        " (KaKs_Calculator's methods). \"",
                        paste(setdiff(dnds_est.method, coal_methods()), collapse = "\", \""),
                        "\" is only available with engine = \"orthologr\".",
                        call. = FALSE
                )
        if (!is.element(aa_aln_scoring, c("pwalign", "needle")))
                stop("aa_aln_scoring must be \"pwalign\" or \"needle\".", call. = FALSE)
        if (!is.character(aa_sub_matrix) || length(aa_sub_matrix) != 1)
                stop("With engine = \"coal\", aa_sub_matrix must be the name of a matrix: \"BLOSUM62\", \"BLOSUM45\", \"BLOSUM80\" or \"PAM250\".", call. = FALSE)
        coal_backend(need = "coalr")

        message("\n")
        message("Starting orthology inference (", ortho_detection, ") and dNdS estimation (",
                paste(dnds_est.method, collapse = ", "), ") with the COAL engine using the follwing parameters:")
        message("query = '", basename(query_file), "'")
        message("subject = '", basename(subject_file), "'")
        message("aligner = ", "'", basename(aligner), "'")
        if (aligner == "diamond")
                message("sensitivity_mode = ", "'", basename(sensitivity_mode), "'")
        message("e-value: ", eval)
        message("aa_aln_scoring = '", aa_aln_scoring, "' (", aa_sub_matrix, ")")
        message("comp_cores = '", comp_cores, "'")
        message("\n")

        # the search step works in a folder of its own, removed when the call ends
        end_tmpdir <- orthologr_tmpdir_start("orthologr_dNdS_")
        on.exit(end_tmpdir(), add = TRUE)

        # the two input files are read and translated once for all steps of this call
        end_seq_cache <- seq_cache_start()
        on.exit(end_seq_cache(), add = TRUE)

        message("Starting Orthology Inference ...")

        if (aligner == "diamond" & ortho_detection == "BH")
                align <- function(...)
                        diamond_best(..., diamond_params = aligner_params, sensitivity_mode = sensitivity_mode)
        if (aligner == "diamond" & ortho_detection == "RBH")
                align <- function(...)
                        diamond_rec(..., diamond_params = aligner_params, sensitivity_mode = sensitivity_mode)
        if (aligner == "blast" & ortho_detection == "BH")
                align <- function(...)
                        blast_best(..., blast_params = aligner_params)
        if (aligner == "blast" & ortho_detection == "RBH")
                align <- function(...)
                        blast_rec(..., blast_params = aligner_params)

        hit.table <-
                align(
                        query_file   = query_file,
                        subject_file = subject_file,
                        delete_corrupt_cds = delete_corrupt_cds,
                        path         = aligner_path,
                        comp_cores   = comp_cores,
                        seq_type     = seq_type,
                        eval         = eval,
                        format       = format
                )

        data.table::setDT(hit.table)
        data.table::setkeyv(hit.table, c("query_id", "subject_id"))

        # (both files were read and translated for the search: they are not read again)
        q_cds <- read_cds_once(query_file, delete_corrupt_cds = delete_corrupt_cds)
        s_cds <- read_cds_once(subject_file, delete_corrupt_cds = delete_corrupt_cds)

        # the gene pairs dNdS() would process: hits whose two genes have a coding sequence that passes the
        # filters of the search step (cds2aa())
        q_ok <- cds_protein_table_file(query_file, delete_corrupt_cds = delete_corrupt_cds)$geneids
        s_ok <- cds_protein_table_file(subject_file, delete_corrupt_cds = delete_corrupt_cds)$geneids
        pairs <- unique(as.data.frame(hit.table)[ , c("query_id", "subject_id")])
        pairs <- pairs[is.element(pairs$query_id, q_ok) & is.element(pairs$subject_id, s_ok), , drop = FALSE]

        if (nrow(pairs) == 0)
                stop("No orthologs could be found! Please check your input files!", call. = FALSE)

        message("Orthology Inference Completed.")
        message("Starting dN/dS Estimation (COAL, ", nrow(pairs), " gene pairs) ...")

        coal_tbl <- coalr::coal_dnds_pairs(
                query_cds   = cds_vector(q_cds),
                subject_cds = cds_vector(s_cds),
                pairs       = pairs,
                methods     = dnds_est.method,
                threads     = comp_cores,
                matrix      = aa_sub_matrix,
                preset      = aa_aln_scoring
        )

        dNdS_tbl <- data.table::data.table(
                query_id   = coal_tbl$seq1,
                subject_id = coal_tbl$seq2,
                dN         = coal_tbl$dN,
                dS         = coal_tbl$dS,
                dNdS       = coal_tbl$dNdS,
                method     = coal_tbl$method
        )
        data.table::setkeyv(dNdS_tbl, c("query_id", "subject_id"))
        dNdS_tbl <- tibble::as_tibble(dNdS_tbl)

        failed <- coal_tbl[coal_tbl$status == "error", c("seq1", "seq2", "note"), drop = FALSE]
        failed <- unique(data.frame(query_id = failed$seq1, subject_id = failed$seq2, reason = failed$note,
                                    stringsAsFactors = FALSE))
        if (nrow(failed) > 0) {
                warning(
                        nrow(failed), " of ", nrow(pairs),
                        " gene pairs could not be processed; their dN, dS and dNdS values are NA. ",
                        "The reasons are in attr(<result>, \"failed_pairs\"). First one: ",
                        "'", failed$query_id[1], "' / '", failed$subject_id[1], "': ", failed$reason[1],
                        call. = FALSE
                )
        }

        subject_id <- NULL
        hit.table_selected <- dplyr::select(hit.table, -subject_id)
        res <- dplyr::inner_join(dNdS_tbl, hit.table_selected, by = "query_id", copy = TRUE)

        attr(res, "failed_pairs") <- failed
        attr(res, "coal") <- coal_tbl

        message("dN/dS Estimation Completed.")

        if (print_citation) {
                message("\n")
                message("Please cite the following paper when using orthologr for your own research:")
                cat("Drost et al. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. 2015. Mol. Biol. Evol. 32 (5): 1221-1231.")
                message("\n")
                message("\n")
        }

        return(res)
}

#' @title Compute dNdS values for all gene pairs within sequence clusters
#' @description Estimates dN, dS and dN/dS for the pairs of genes in each cluster of homologous sequences, e.g. the
#' clusters of \code{\link{deepclust}} or \code{\link{deepclust_realign}}, and summarises each cluster. For every
#' cluster the proteins are aligned (one multiple alignment per cluster with Clustal Omega or MAFFT, or pairwise),
#' the codons are threaded onto the alignment and Ka/Ks is estimated for the pairs of the cluster by COAL
#' (\url{https://github.com/Plumat-COAL/COAL}), in process and on \code{comp_cores} threads.
#' @param cds_file a character vector with the path(s) to the CDS file(s) (fasta) holding the coding sequences of all cluster members.
#' @param clusters the clusters: the tibble returned by \code{\link{deepclust}}, \code{\link{deepclust_realign}} or
#' \code{\link{deepclust_annotate}} (first two columns: \code{representative_id}, \code{member_id}), any data.frame
#' with these two columns first, or the path to such a tab-separated file. If \code{clusters = NULL} (Default),
#' \code{\link{deepclust}} is run on \code{cds_file}.
#' @param dnds_est.method the dN/dS estimation method(s): \code{"YN"} (Default), \code{"NG"}, \code{"LWL"},
#' \code{"MLWL"}, \code{"LPB"}, \code{"MLPB"} or \code{"MYN"}.
#' @param aa_aln_tool how the proteins of clusters with three or more members are aligned:
#' \itemize{
#' \item \code{"clustalo"} (Default): one Clustal Omega multiple alignment per cluster;
#' \item \code{"mafft"}: one MAFFT multiple alignment per cluster;
#' \item \code{"pairwise"}: every pair is aligned on its own (global, BLOSUM62) by COAL; no external program is needed.
#' }
#' Clusters of two members are always aligned pairwise by COAL.
#' @param aa_aln_path a character string specifying the path to the Clustal Omega or MAFFT executable (in case it is not on the PATH).
#' @param max_cluster_size clusters with more usable members than this are not computed for all pairs (the number of
#' pairs grows with the square of the cluster size), but as \code{large_clusters} says. Default is \code{200}.
#' @param large_clusters what to do with clusters above \code{max_cluster_size}: \code{"representative"} (Default; the
#' representative against every other member), \code{"sample"} (\code{pair_sample} pairs drawn reproducibly with
#' \code{seed}) or \code{"skip"}. Such clusters are aligned pairwise.
#' @param pair_sample the number of pairs drawn per large cluster for \code{large_clusters = "sample"}.
#' @param seed the seed of the draw for \code{large_clusters = "sample"}.
#' @param delete_corrupt_cds a logical value indicating whether coding sequences whose length is not a multiple of 3 shall be left out.
#' @param comp_cores the number of threads (COAL, DIAMOND, and the multiple alignment program).
#' @param path a character string specifying the path to the DIAMOND executable (only used when \code{clusters = NULL}).
#' @param ... additional arguments passed to \code{\link{deepclust}} when \code{clusters = NULL}
#' (e.g. \code{approx_id}, \code{mutual_cover}).
#' @details This function needs the R package \pkg{coalr}:
#' \code{remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")}.
#' A gene, pair or cluster that cannot be used never ends the run: see the \code{status} and \code{note} columns and
#' the \code{failed} table.
#' @author Chloe Newton
#' @return A list of three tibbles:
#' \itemize{
#' \item \code{pairs}: one row per gene pair and method: \code{cluster_id}, \code{query_id}, \code{subject_id},
#' \code{method}, \code{dN}, \code{dS}, \code{dNdS}, standard errors, numbers of sites and differences,
#' \code{codons_compared}, \code{prot_aln_length}, \code{prot_identity}, \code{status} and \code{note};
#' \item \code{summary}: one row per cluster and method: \code{cluster_id}, \code{n_members}, \code{n_used},
#' \code{method}, \code{pair_selection}, \code{aligner}, \code{n_pairs}, \code{n_ok}, \code{n_saturated},
#' \code{n_no_data}, \code{n_failed}, and median, mean, quartiles, minimum and maximum of dN, dS and dN/dS;
#' \item \code{failed}: genes and clusters that could not be used, and why.
#' }
#' @examples \dontrun{
#' cds <- system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr')
#' # cluster with DIAMOND DeepClust, then dN/dS within each cluster
#' res <- dNdS_clusters(cds, aa_aln_tool = "pairwise")
#' res$pairs
#' res$summary
#'
#' # with clusters computed before
#' cl <- deepclust(cds, seq_type = "cds")
#' res <- dNdS_clusters(cds, clusters = cl, aa_aln_tool = "clustalo", comp_cores = 2)
#' }
#' @seealso \code{\link{deepclust}}, \code{\link{deepclust_realign}}, \code{\link{dNdS_coal}}, \code{\link{dNdS}}
#' @export
dNdS_clusters <- function(cds_file,
                          clusters           = NULL,
                          dnds_est.method    = "YN",
                          aa_aln_tool        = "clustalo",
                          aa_aln_path        = NULL,
                          max_cluster_size   = 200,
                          large_clusters     = "representative",
                          pair_sample        = 1000,
                          seed               = 1,
                          delete_corrupt_cds = FALSE,
                          comp_cores         = 1,
                          path               = NULL,
                          ...) {

        missing_files <- cds_file[!file.exists(cds_file)]
        if (length(missing_files) > 0)
                stop("The following input file(s) were not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
        if (!is.character(dnds_est.method) || length(dnds_est.method) < 1 || !all(is.element(dnds_est.method, coal_methods())))
                stop("dnds_est.method must be one or more of: ", paste(coal_methods(), collapse = ", "), ".", call. = FALSE)
        if (!is.element(aa_aln_tool, c("clustalo", "mafft", "pairwise")))
                stop("aa_aln_tool must be \"clustalo\", \"mafft\" or \"pairwise\".", call. = FALSE)
        if (!is.element(large_clusters, c("representative", "sample", "skip")))
                stop("large_clusters must be \"representative\", \"sample\" or \"skip\".", call. = FALSE)
        coal_backend(need = "coalr")

        if (is.null(clusters)) {
                message("Clustering the sequences with DIAMOND DeepClust ...")
                clusters <- deepclust(input_file = cds_file, seq_type = "cds", delete_corrupt_cds = delete_corrupt_cds,
                                      path = path, comp_cores = comp_cores, ...)
        } else if (is.character(clusters) && length(clusters) == 1) {
                if (!file.exists(clusters))
                        stop("The clusters file '", clusters, "' was not found.", call. = FALSE)
        } else if (!is.data.frame(clusters) || ncol(clusters) < 2) {
                stop("'clusters' must be the table returned by deepclust() / deepclust_realign() (representative_id, member_id), ",
                     "a file of that table, or NULL.", call. = FALSE)
        }
        if (is.data.frame(clusters))
                clusters <- as.data.frame(clusters)[ , 1:2]

        geneids <- seqs <- NULL
        cds <- do.call(c, lapply(cds_file, function(f) {
                cds_vector(read.cds(file = f, format = "fasta", delete_corrupt_cds = delete_corrupt_cds))
        }))

        msa_aligner <- NULL
        if (aa_aln_tool != "pairwise") {
                exe <- if (is.null(aa_aln_path)) NULL else {
                        if (dir.exists(aa_aln_path)) file.path(aa_aln_path, aa_aln_tool) else aa_aln_path
                }
                msa_aligner <- tryCatch(
                        coalr::coal_msa_aligner(aa_aln_tool, path = exe, threads = 1),
                        error = function(e) stop(conditionMessage(e),
                                                 "\nOr use aa_aln_tool = \"pairwise\" (no external program needed).",
                                                 call. = FALSE)
                )
        }

        message("Estimating dN/dS within ", length(unique(clusters[[1]])), " clusters with COAL ...")
        res <- coalr::coal_dnds_clusters(
                cds              = cds,
                clusters         = clusters,
                aligner          = msa_aligner,
                methods          = dnds_est.method,
                threads          = comp_cores,
                max_cluster_size = max_cluster_size,
                large_clusters   = large_clusters,
                pair_sample      = pair_sample,
                seed             = seed
        )

        pairs <- res$pairs
        names(pairs)[names(pairs) == "seq1"] <- "query_id"
        names(pairs)[names(pairs) == "seq2"] <- "subject_id"

        if (nrow(res$failed) > 0)
                message(nrow(res$failed), " gene(s) or cluster(s) could not be used; see the 'failed' table of the result.")

        list(
                pairs   = tibble::as_tibble(pairs),
                summary = tibble::as_tibble(res$summary),
                failed  = tibble::as_tibble(res$failed)
        )
}
