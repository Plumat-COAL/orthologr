#' @title Select orthologs based on either gene locus or splice variant
#' @description This function selects orthologs based on either gene locus or splice variant.
#' @param dnds_tbl a tibble returned by \code{\link{dNdS}}.
#' @param annotation_file_query file path to query annotation file in either \code{gtf} or \code{gff} format.
#' @param annotation_file_subject file path to query annotation file in either \code{gtf} or \code{gff} format.
#' @param collapse_by locus type by which orthologs should be determined. Options are:
#' \itemize{
#' \item \code{type = "gene_locus"}: the splice variant having the lowest e-value of the BLAST search is chosen to represent the gene locus for the orthology relationship.
#' \item \code{type = "splice_variant"}: all splice variants of homologous gene loci will be returned.
#' }
#' @param format a vector of length 2 storing the annotation file formats of the query annotation file and subject annotation file: either \code{gtf} or \code{gff} format. E.g. \code{format = c("gtf","gtf")}.
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
#' The IDs of the dN/dS table (\code{query_id}, \code{subject_id}) may be protein IDs or transcript IDs.
#' They are compared with the IDs of the annotation as they are (including the spelling with the version that Ensembl
#' gives in a field of its own, e.g. \code{ENST00000373020.9}, and the raw GFF3 IDs such as \code{rna-NM_001.1} or
#' \code{transcript:ENST00000373020}); IDs that are not found are compared without their version suffix, and IDs of
#' NCBI \code{cds_from_genomic} files (\code{lcl|NC_000913.3_cds_NP_414542.1_1}) by their protein accession.
#' Tested with NCBI RefSeq GFF3 + \code{protein.faa}, NCBI GTF, Ensembl GFF3 and Ensembl GTF with \code{pep} and
#' \code{cds} FASTA files. The annotation is read as text in blocks (\code{.gz} files too); it is not imported with
#' \pkg{rtracklayer}.
#' 
#' The columns \code{query_gene_locus_id} and \code{subject_gene_locus_id} hold the gene loci. IDs of the dN/dS table
#' that are not in an annotation are left out with a message; when none is found the error shows example IDs of the
#' table and of the annotation.
#' @author Hajk-Georg Drost
#' @param best_hit_rule how the best hit of a query is chosen among its hits:
#' \itemize{
#' \item \code{"bitscore"} (default since orthologr 0.5.0): the hit with the highest bit score; of several, the one
#' with the smallest e-value, then the longest alignment, then the highest percent identity, then the alphabetically
#' first \code{subject_id}. The choice does not depend on the order of the hit table.
#' \item \code{"evalue"}: the rule of orthologr <= 0.4.3 - the hit with the smallest e-value; of several, the one with
#' the longest alignment; of several, the alphabetically first \code{subject_id}. The bit score is not used.
#' }
#' \strong{This default changed in 0.5.0 and changes which orthologs are returned.} DIAMOND and BLAST report an
#' e-value of 0 for every strong hit of a long protein, so that under the old rule the alignment length decided
#' among them and a paralog could win over the true ortholog. @@BHNUM@@
#' \code{options(orthologr.best_hit_rule = "evalue")} restores the old rule for all functions of a session.
#' @export
select_orthologs <-
        function(dnds_tbl,
                 annotation_file_query,
                 annotation_file_subject,
                 collapse_by = "gene_locus",
                 format = c("gtf", "gtf"),
                 best_hit_rule = getOption("orthologr.best_hit_rule", "bitscore")) {
                
                best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)
                if (!is.element(collapse_by, c("gene_locus", "splice_variant")))
                        stop(
                                "Please specify a collapse method provided by this function: collapse_by = 'gene_locus' or collapse_by = 'splice_variant'.",
                                call. = FALSE
                        )
                
                if (!all(is.element(format, c("gtf", "gff"))))
                        stop(
                                "Please provide a vector of length 2 that contains available formats: either 'gtf' or 'gff'.",
                                call. = FALSE
                        )
                
                if (length(format) != 2)
                        stop(
                                "Please provide a vector of length 2 to specify the query and subject annotation file formats. E.g. format = c('gtf', 'gtf').",
                                call. = FALSE
                        )
                
                if (!file.exists(annotation_file_query))
                        stop(
                                "The file '",
                                annotation_file_query,
                                " does not seem to exist. Please check the path to your annotation file.",
                                call. = FALSE
                        )
                
                if (!file.exists(annotation_file_subject))
                        stop(
                                "The file '",
                                annotation_file_subject,
                                " does not seem to exist. Please check the path to your annotation file.",
                                call. = FALSE
                        )
                
                '.' <- NULL
                
                # (Earlier versions stored the joined table in tempdir() under a name made of the two annotation
                # file names only and reused it in later calls - also for a different dnds_tbl. The table is
                # now computed in every call.)
                {
                        message("Start ortholog selection process by ",
                                collapse_by,
                                " ...")
                        
                        # The gene locus of every query_id / subject_id, from one transcript -> gene table per
                        # annotation (GFF3: Parent / ID relations; GTF: gene_id / transcript_id / protein_id; see
                        # annotation_table()). IDs are matched as they are, then without version suffix.
                        message("Read query annotation file '", annotation_file_query, "' using format ", format[1], " ...")
                        message("Select orthologs in query species ...")
                        extracted_features_qry <- annotation_loci_of_ids(
                                dnds_tbl$query_id, annotation_file_query, format[1], what = "query_id")
                        names(extracted_features_qry) <- c("query_gene_locus_id", "query_id")
                        
                        message("Read subject annotation file '", annotation_file_subject, "' using format ", format[2], " ...")
                        message("Select orthologs in subject species ...")
                        extracted_features_sbj <- annotation_loci_of_ids(
                                dnds_tbl$subject_id, annotation_file_subject, format[2], what = "subject_id")
                        names(extracted_features_sbj) <- c("subject_gene_locus_id", "subject_id")
                        
                        message("Join dnds and annotation tables ...")
                        # join dnds tbl with query annotation
                        joined_dnds_annotation_tbl_qry <-
                                dplyr::inner_join(dnds_tbl,
                                                  extracted_features_qry,
                                                  by = "query_id")
                        
                        if (nrow(joined_dnds_annotation_tbl_qry) == 0)
                                stop("When joining the 'dnds_tbl' and the 'extracted_features_qry' table by 'query_id', no entries were retained. Thus, this join created an empty table. Please check what might have gone wrong during the query feature extraction process.", call. = FALSE)
                        
                        # join dnds tbl with subject annotation
                        joined_dnds_annotation_tbl_sbj <-
                                dplyr::inner_join(dnds_tbl,
                                                  extracted_features_sbj,
                                                  by = "subject_id")
                        
                        if (nrow(joined_dnds_annotation_tbl_sbj) == 0)
                                stop("When joining the 'dnds_tbl' and the 'extracted_features_sbj' table by 'query_id', no entries were retained. Thus, this join created an empty table. Please check what might have gone wrong during the subject feature extraction process.", call. = FALSE)
                        joined_dnds_annotation_tbl_sbj <-
                                dplyr::select(
                                        joined_dnds_annotation_tbl_sbj,
                                        "query_id",
                                        "subject_gene_locus_id",
                                        "subject_id"
                                )
                        
                        # join query and subject tables
                        final_join_dnds_annotation_tbl <-
                                dplyr::inner_join(
                                        joined_dnds_annotation_tbl_qry,
                                        joined_dnds_annotation_tbl_sbj,
                                        by = c("query_id", "subject_id")
                                )
                        
                        if (nrow(final_join_dnds_annotation_tbl) == 0)
                                stop("When joining the 'joined_dnds_annotation_tbl_qry' and the 'joined_dnds_annotation_tbl_sbj' tables by 'query_id' and 'subject_id', no entries were retained. Thus, this join created an empty table. Please check what might have gone wrong during the joining process.", call. = FALSE)
                        
                }
                
                if (collapse_by == "gene_locus") {
                        message(
                                "Select splice variant with smallest e-value for each gene locus of query species ..."
                        )
                        query_gene_locus_id <-
                                subject_gene_locus_id <- NULL
                        res_qry <-
                                best_hits(final_join_dnds_annotation_tbl, "query_gene_locus_id", rule = best_hit_rule)
                        
                        if (nrow(res_qry) == 0)
                                stop("After best splice variant filtering no results were retained. The output was empty ...", call. = FALSE)
                        
                        message(
                                "Select splice variant with smallest e-value for each gene locus of subject species ..."
                        )
                        res_sbj <-
                                dplyr::ungroup(best_hits(final_join_dnds_annotation_tbl, "subject_gene_locus_id", rule = best_hit_rule))
                        
                        if (nrow(res_sbj) == 0)
                                stop("After best splice variant filtering no results were retained. The output was empty ...", call. = FALSE)
                        res_sbj <-
                                dplyr::select(dplyr::ungroup(res_sbj),
                                              "query_id",
                                              "subject_id")
                        
                        res <-
                                dplyr::inner_join(dplyr::ungroup(res_qry),
                                                  res_sbj,
                                                  by = c("query_id", "subject_id"))
                        
                        
                        if (nrow(res) == 0)
                                stop("After joining the best splice variants of the query and subject by 'query_id' and 'subject_id' no overlaps were found. The output was empty ...", call. = FALSE)
                        
                        if (length(unique(res$query_id)) == length(unique(res$query_gene_locus_id))) {
                                message("Good News: ")
                                message(
                                        "Number of unique query_ids matches the number of unique query_gene_locus_ids. Thus, each gene locus has now a representative splice variant. The one having the smallest evalue and longest alignment length if two splice variants have the same evalue."
                                )
                                message("\n")
                        } else {
                                message(
                                        "Please check your input data. It seems that the number of unique query_ids DOES NOT match the number of unique query_gene_locus_ids. Thus, some query gene loci seem to have multiple splice variants as representatives."
                                )
                                message("\n")
                        }
                        
                        if (length(unique(res$subject_id)) == length(unique(res$subject_gene_locus_id))) {
                                message("Good news:")
                                message(
                                        "Number of unique subject_ids matches the number of unique subject_gene_locus_ids. Thus, each gene locus has now a representative splice variant. The one having the smallest evalue and longest alignment length if two splice variants have the same evalue."
                                )
                                message("\n")
                        } else {
                                message(
                                        "Please check your input data. It seems that the number of unique subject_ids DOES NOT match the number of unique subject_gene_locus_ids Thus, some subject gene loci seem to have multiple splice variants as representatives."
                                )
                                message("\n")
                        }
                        
                        if (length(unique(dnds_tbl$query_id)) < length(unique(res$query_id)))
                                stop(
                                        "Something must have gone wrong... The unique number of query_ids in the dNdS_tbl is smaller than the unique number of query_ids in the collapsed output file. Please check your input data.",
                                        call. = FALSE
                                )
                        
                }
                
                if (collapse_by == "splice_variant") {
                        res <- final_join_dnds_annotation_tbl
                        query_id <- subject_id <- NULL
                        res <- dplyr::distinct(
                                res,
                                query_id,
                                query_gene_locus_id,
                                subject_id,
                                subject_gene_locus_id,
                                .keep_all = TRUE
                        )
                }
                
                qry_species_name <-
                        unlist(stringr::str_split(basename(annotation_file_query), "[.]"))[1]
                sbj_species_name <-
                        unlist(stringr::str_split(basename(annotation_file_subject), "[.]"))[1]
                res <-
                        dplyr::mutate(
                                res,
                                query_species = rep(qry_species_name, nrow(res)),
                                subject_species = rep(sbj_species_name, nrow(res))
                        )
                
                evalue <- bit_score <- perc_identity <- NULL
                s_end <-
                        dN <-
                        dS <- query_species <- subject_species <- score_raw <- NULL
                res <- dplyr::select(
                        res,
                        query_species,
                        subject_species,
                        query_id,
                        query_gene_locus_id,
                        subject_id,
                        subject_gene_locus_id,
                        dN,
                        dS,
                        dNdS,
                        evalue,
                        bit_score,
                        perc_identity:score_raw
                )
                
                message("... orthology inference finished successfully!")
                return(res)
        }


# The gene loci of the IDs of a dN/dS table: a tibble (gene_locus_id, id) with one row per ID that is in the
# annotation. Stops with examples of both kinds of IDs when none is.
annotation_loci_of_ids <- function(ids, annotation_file, format, what = "query_id") {
        tbl <- annotation_table(annotation_file, format = format)
        ids <- unique(as.character(ids))
        loci <- match_ids_to_genes(ids, attr(tbl, "keys"))
        if (nrow(loci) == 0)
                annotation_mismatch_error(ids, tbl, paste0("the column '", what, "' of the dN/dS table"), annotation_file)
        # an ID that belongs to several loci (identical proteins): the first locus
        multi <- duplicated(loci$id)
        if (any(multi)) {
                message(length(unique(loci$id[multi])), " ", what, "s belong to more than one gene locus of '",
                        basename(annotation_file), "' (e.g. ", loci$id[multi][1], "); the first locus is used.")
                loci <- loci[!multi]
        }
        how <- table(loci$how)
        if (any(names(how) != "exact"))
                message(sum(how[names(how) != "exact"]), " of ", length(ids), " ", what, "s were matched to the annotation without their version suffix ",
                        "or by the protein accession of an NCBI 'lcl|..._cds_...' ID.")
        not_found <- setdiff(ids, loci$id)
        if (length(not_found) > 0)
                message(length(not_found), " of ", length(ids), " ", what, "s are not in the annotation file '", basename(annotation_file),
                        "' and are left out, e.g.: ", paste(utils::head(not_found, 5), collapse = ", "), ".")
        tibble::tibble(gene_locus_id = loci$gene_id, id = loci$id)
}

# The loci table of the implementation used until 0.4.3 (rtracklayer import, extract_features() per locus; GFF3
# rows grouped by 'Name', GTF: literal transcript_id matches of loci with a 'gene' and CDS rows). Kept as the
# reference in the tests.
annotation_loci_of_ids_do <- function(ids, annotation_file, format) {
        need_package("rtracklayer", "annotation_loci_of_ids_do()", bioc = TRUE)
        imported <- tibble::as_tibble(rtracklayer::import(annotation_file, format = format))
        if (any(imported$source == "ERCC"))
                imported <- imported[imported$source != "ERCC", ]
        if (format == "gff") {
                res <- extract_features_all(imported, format = format, group_col = "Name")
        } else {
                sel <- unique(imported[imported$transcript_id %in% ids, ]$gene_id)
                res <- extract_features_all(imported[imported$gene_id %in% sel, ], format = format, group_col = "gene_id")
        }
        res <- dplyr::ungroup(res)
        tibble::tibble(gene_locus_id = res$gene_locus_id, id = res$query_id)
}
