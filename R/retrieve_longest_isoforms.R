#' @title Retrieve the longest isoforms from a proteome file and save results as fasta file
#' @description Based on a fasta file storing the peptide isoforms of gene loci and
#' an annotation file in \code{gff} (GFF3) or \code{gtf} file format, this function extracts the longest
#' isoform per gene locus and stores the results in a new \code{fasta} file.
#' This procedure enables easier downstream analyses such as orthology inference etc
#' when dealing with proteome \code{fasta} files which usually include isoform peptides.
#' @param proteome_file file path to proteome in \code{fasta} file format.
#' @param annotation_file file path to the corresponding annotation file in \code{gff} (GFF3) or \code{gtf} file format (may be gzip-compressed).
#' @param new_file file path to new file storing only peptide sequences of the longest isoforms.
#' @param annotation_format format of \code{annotation_file}. Options are:
#' \itemize{
#' \item \code{annotation_file = "gff"} (default)
#' \item \code{annotation_file = "gtf"} 
#' }
#' @details
#' \strong{Which gene locus a sequence belongs to} is read from the annotation in the same way for all annotation
#' sources:
#' \itemize{
#' \item \code{gff} (GFF3): the \code{Parent} / \code{ID} relations are followed from every CDS and transcript up
#' to the top-level feature, the gene (CDS -> mRNA -> gene; in NCBI prokaryote annotations CDS -> gene). The gene
#' locus is the \code{ID} of the gene without the prefix \code{gene-} (NCBI; this is the \code{locus_tag}) or
#' \code{gene:} (Ensembl). It is never the gene symbol: genes without a symbol and genes that share one are
#' separate loci.
#' \item \code{gtf}: \code{gene_id} is the gene locus; \code{transcript_id} and \code{protein_id} name the
#' transcript and the protein.
#' }
#' The IDs of the sequence file (the FASTA header up to the first white space) may be protein IDs or transcript IDs.
#' They are compared with the IDs of the annotation as they are (including the spelling with the version that Ensembl
#' gives in a field of its own, e.g. \code{ENST00000373020.9}, and the raw GFF3 IDs such as \code{rna-NM_001.1} or
#' \code{transcript:ENST00000373020}); IDs that are not found are compared without their version suffix, and IDs of
#' NCBI \code{cds_from_genomic} files (\code{lcl|NC_000913.3_cds_NP_414542.1_1}) by their protein accession.
#' Tested with NCBI RefSeq GFF3 + \code{protein.faa}, NCBI GTF, Ensembl GFF3 and Ensembl GTF with \code{pep} and
#' \code{cds} FASTA files. The annotation is read as text in blocks (\code{.gz} files too); it is not imported with
#' \pkg{rtracklayer}.
#' 
#' Of several equally long sequences of a locus the first one in the file is kept. With \code{annotation_format = "gtf"}
#' only genes with \code{gene_biotype "protein_coding"} are kept when the annotation gives biotypes (as in earlier
#' versions; a message says how many sequences of other biotypes are left out); with \code{"gff"} every locus with a
#' sequence is kept. Sequences that are not in the annotation are left out with a message; when no ID fits, the
#' error shows example IDs of both files.
#' @return The path \code{new_file}, invisibly, with the attributes \code{n_loci} and \code{n_sequences}.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # retrieve example data from ENSEMBLGENOMES
#' proteome <- biomartr::getProteome(db = "refseq", organism = "Arabidopsis thaliana")
#' annotation <- biomartr::getGFF(db = "refseq", organism = "Arabidopsis thaliana")
#' # retrieve longest isoforms and store in new file
#' retrieve_longest_isoforms(proteome_file = proteome, 
#'                           annotation_file = annotation, 
#'                           new_file = "Athaliana_pep_longest.fa")
#' # import new file into R session                          
#' Athaliana_pep_longest <- Biostrings::readAAStringSet("Athaliana_pep_longest.fa")
#' }
#' @export
retrieve_longest_isoforms <- function(proteome_file, annotation_file, new_file, annotation_format = "gff") {
        
        if (!file.exists(proteome_file))
                stop("File ", proteome_file, " seems not to exist.", call. = FALSE)
        
        if (!file.exists(annotation_file))
                stop("File ", annotation_file, " seems not to exist.", call. = FALSE)
        
        if (!is.element(annotation_format, c("gtf", "gff")))
                stop("Please specify a valid annotation format. Options are: annotation_file = 'gtf' or annotation_file = 'gff'.", call. = FALSE)
        
        gene_id <- id <- pos <- width <- NULL
        
        message("Extracting longest isoform from ", basename(proteome_file), " ...")
        message("Importing proteome file ", basename(proteome_file), " ...")
        pep_file <- Biostrings::readAAStringSet(proteome_file)
        new_names <- fasta_header_id(names(pep_file))
        check_duplicate_ids(new_names, what = "sequences", file = proteome_file)
        names(pep_file) <- new_names
        
        # check if annotation file is corrupt and remove outliers
        annotation_file <- check_annotation(annotation_file, remove_annotation_outliers = TRUE)
        
        message("Reading the transcript -> gene relations of the ", annotation_format, " file ", basename(annotation_file), " ...")
        tbl <- annotation_table(annotation_file, format = annotation_format)
        message("The annotation holds ", nrow(tbl), " transcripts / coding sequences of ", length(unique(tbl$gene_id)), " gene loci.")
        
        # the gene locus (or loci) of every sequence
        loci <- match_ids_to_genes(new_names, attr(tbl, "keys"))
        if (nrow(loci) == 0)
                annotation_mismatch_error(new_names, tbl, paste0("the sequence file '", basename(proteome_file), "'"), annotation_file)
        
        how <- table(loci$how[!duplicated(loci$id)])
        if (any(names(how) != "exact"))
                message("Sequence IDs were matched to the annotation: ",
                        paste0(how, " ", c(exact = "as they are", without_version = "without their version suffix",
                                           ncbi_cds_id = "by the protein accession in their NCBI 'lcl|..._cds_...' ID")[names(how)],
                               collapse = ", "), ".")
        
        not_found <- setdiff(new_names, loci$id)
        if (length(not_found) > 0)
                message(length(not_found), " of ", length(new_names), " sequences are not in the annotation and are left out, e.g.: ",
                        paste(utils::head(not_found, 5), collapse = ", "), ".")
        
        if (annotation_format == "gtf" && any(!is.na(tbl$gene_biotype))) {
                # as before for GTF files: the genes with gene_biotype 'protein_coding'
                biotype <- tbl$gene_biotype[match(loci$gene_id, tbl$gene_id)]
                other <- !is.na(biotype) & biotype != "protein_coding"
                if (any(other)) {
                        tab <- table(biotype[other & !duplicated(loci$gene_id)])
                        message("Filter for gene_biotype == 'protein_coding': ", length(unique(loci$id[other])), " sequences of ",
                                length(unique(loci$gene_id[other])), " gene loci with another biotype are left out (",
                                paste0(names(tab), ": ", tab, collapse = ", "), ").")
                        loci <- loci[!other]
                }
                if (nrow(loci) == 0)
                        stop("No sequence of '", basename(proteome_file), "' belongs to a gene with gene_biotype 'protein_coding' in '",
                             basename(annotation_file), "'.", call. = FALSE)
        }
        
        # the longest sequence of every gene locus (of several equally long ones the first in the file)
        loci[ , pos := match(id, new_names)]
        loci <- loci[order(pos)]
        loci[ , width := Biostrings::width(pep_file)[pos]]
        message("Retrieving longest isoforms ...")
        res <- longest_per_group(tibble::as_tibble(loci), "gene_id")
        
        # a protein that is the longest of several loci (identical proteins of NCBI prokaryote annotations) is
        # written once
        matched_ids <- unique(res$pos)
        n_loci <- length(unique(res$gene_id))
        
        message("\n")
        message("Writing ", length(matched_ids), " unique longest peptide sequences (from initially ", length(pep_file),
                " isoforms) for ", n_loci, " gene loci to new fasta file ", new_file, " ...",
                if (length(matched_ids) < n_loci)
                        paste0(" (", n_loci - length(matched_ids), " loci share their protein with another locus: identical proteins have one FASTA record.)")
                else "")
        Biostrings::writeXStringSet(pep_file[matched_ids], new_file)
        message("Retrieval finished successfully.")
        
        invisible(structure(new_file, n_loci = n_loci, n_sequences = length(matched_ids)))
}

# The implementation used until 0.4.3 (rtracklayer import; GFF3: CDS rows grouped by the 'gene' attribute, GTF:
# transcript rows with gene_biotype 'protein_coding' matched literally by transcript_id). Kept as the reference
# in the tests for the files it could handle.
retrieve_longest_isoforms_do <- function(proteome_file, annotation_file, new_file, annotation_format = "gff") {
        need_package("rtracklayer", "retrieve_longest_isoforms_do()", bioc = TRUE)
        if (!file.exists(proteome_file))
                stop("File ", proteome_file, " seems not to exist.", call. = FALSE)
        
        if (!file.exists(annotation_file))
                stop("File ", annotation_file, " seems not to exist.", call. = FALSE)
        
        if (!is.element(annotation_format, c("gtf", "gff")))
                stop("Please specify a valid annotation format. Options are: annotation_file = 'gtf' or annotation_file = 'gff'.", call. = FALSE)
        
        '.' <- NULL
        
        message("Extracting longest isoform from ", basename(proteome_file), " ...")
        message("Importing proteome file ", basename(proteome_file), " ...")
        pep_file <- Biostrings::readAAStringSet(proteome_file)
        new_names <- unlist(lapply(pep_file@ranges@NAMES, function(x) {
                unlist(stringr::str_split(x," "))[1]
        }))
        names(pep_file) <- new_names
        
        # check if annotation file is corrupt and remove outliers
        annotation_file <- check_annotation(annotation_file, remove_annotation_outliers = TRUE)
        
        if (annotation_format == "gtf") {
                
                pep_file_tibble <- tibble::tibble(transcript_id = new_names, width = pep_file@ranges@width)
                #print(pep_file_tibble)
                
                message("Importing gtf file ", basename(annotation_file), " ...")
                gtf_import <- tibble::as_tibble(rtracklayer::import(annotation_file))
                message("Filter for gene_biotype == 'protein_coding' AND type == 'transcript' ...")
                gene_biotype <- gene_id <- transcript_id <- type <- width <- NULL
                gtf_import <- dplyr::filter(gtf_import, gene_biotype == "protein_coding", type == "transcript")
                message("After filtering, the gtf file includes ", nrow(gtf_import), " rows.")
                gtf_import <- dplyr::select(gtf_import, gene_id, transcript_id)
                #print(gtf_import)
                # join transcript_ids from annotation file with peptide sequence file
                message("Join transcript_ids from FASTA header with transcript_ids from GTF annotation file ...")
                pep_file_tibble_joined <- dplyr::inner_join(pep_file_tibble, gtf_import, by = "transcript_id")
                message("The joined file contains ", nrow(pep_file_tibble_joined), " rows.")
                message("Removing rows with NA's in columns 'transcript_id', 'width', and 'gene_id'")
                pep_file_tibble_joined <- dplyr::filter(pep_file_tibble_joined, !is.na(transcript_id), !is.na(width), !is.na(gene_id))
                message("After NA removal the file contains ", nrow(pep_file_tibble_joined), " rows.")
                
                #print(pep_file_tibble_joined)
                extract_longest_transcript <- function(x) {
                        pos <- which.max(x$width)[1]
                        return(dplyr::slice(x, pos))
                }
                message("Retrieving longest isoforms ...")
                res <- longest_per_group(pep_file_tibble_joined, "gene_id")
                matched_ids <- stats::na.omit(match(res$transcript_id, pep_file@ranges@NAMES))
                
                if (length(matched_ids) == 0)
                        stop("Gene_ID or Peptide_IDs between the gtf file and the headers of the fasta file did not match!",
                             " Please make sure that the 'transcript_id' column in the gtf file matches the IDs given in the headers of the fasta file.", call. = FALSE)
                if (length(matched_ids) != length(unique(matched_ids)))
                        stop("There should only be unique gene_ids in the file, but non-unique IDs were found! Please check what could have gone wrong. with the input files.", call. = FALSE)
                
                message("\n")
                message("Writing ", length(matched_ids), " unique longest peptide sequences (from initially ", length(pep_file), " isoforms) to new fasta file ", new_file, " ...")
                Biostrings::writeXStringSet(pep_file[matched_ids], new_file)
                message("Retrieval finished successfully.")
        }
        
        if (annotation_format == "gff") {
                
                pep_file_tibble <- tibble::tibble(protein_id = new_names, width = pep_file@ranges@width)
                #print(pep_file_tibble)
                message("Importing gff file ", basename(annotation_file), " ...")
                gff_import <- tibble::as_tibble(rtracklayer::import(annotation_file))
                gene_biotype <- gene <- transcript_id <- type <- width <- NULL
                gff_import <- dplyr::filter(gff_import, type == "CDS")
                protein_id <- NULL
                gff_import <- dplyr::select(gff_import, gene, protein_id)
                #print(gff_import)
                # join transcript_ids from annotation file with peptide sequence file
                pep_file_tibble_joined <- dplyr::inner_join(pep_file_tibble, gff_import, by = "protein_id")   
                pep_file_tibble_joined <- dplyr::filter(pep_file_tibble_joined, !is.na(protein_id), !is.na(width))
                #print(pep_file_tibble_joined)
                extract_longest_transcript <- function(x) {
                        pos <- which.max(x$width)[1]
                        return(dplyr::slice(x, pos))
                }
                message("Retrieving longest isoforms ...")
                res <- longest_per_group(pep_file_tibble_joined, "gene")
                matched_ids <- stats::na.omit(match(res$protein_id, pep_file@ranges@NAMES))
                
                if (length(matched_ids) == 0)
                        stop("Gene_ID or Peptide_IDs between the gtf file and the headers of the fasta file did not match!",
                             " Please make sure that the 'protein_id' column in the gtf file matches the IDs given in the headers of the fasta file.", call. = FALSE)
                
                if (length(matched_ids) != length(unique(matched_ids)))
                        stop("There should only be unique gene_ids in the file, but non-unique IDs were found! Please check what could have gone wrong. with the input files.", call. = FALSE)
                
                message("\n")
                message("Writing ", length(matched_ids), " unique longest peptide sequences (from initially ", length(pep_file), " isoforms) to new fasta file ", new_file, " ...")
                Biostrings::writeXStringSet(pep_file[matched_ids], new_file)
                message("Retrieval finished successfully.")
        }
        
        
}

# For every value of `group_col` the first row with the largest 'width': the rows
# dplyr::do(dplyr::group_by(x, <group_col>), dplyr::slice(., which.max(.$width)[1])) returns, in the same order,
# without calling an R function for every gene.
longest_per_group <- function(x, group_col) {
        if (nrow(x) == 0 || anyNA(x$width))
                return(longest_per_group_do(x, group_col))
        rows <- first_row_per_group(x[[group_col]], list(-x$width))
        res <- x[rows, ]
        dplyr::group_by(res, .data[[group_col]])
}

# The implementation used until 0.4.3; kept as the reference in the tests.
longest_per_group_do <- function(x, group_col) {
        '.' <- NULL
        extract_longest_transcript <- function(x) {
                # (not 'pos': dplyr::slice() looks names up in the table first, and the table of
                # retrieve_longest_isoforms() has a column 'pos')
                longest_row <- which.max(x$width)[1]
                return(x[longest_row, , drop = FALSE])
        }
        dplyr::do(dplyr::group_by(x, .data[[group_col]]), extract_longest_transcript(.))
}
