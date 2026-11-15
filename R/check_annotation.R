#' @title Report over-long lines of an annotation file
#' @description Some annotation files hold lines with more than 65,000 characters (for example CDS and mRNA rows
#' of NCBI RefSeq files with very long \code{Note} attributes). The import function of \pkg{rtracklayer} could not
#' read such files, which is why this function exists. The functions of orthologr that read annotation files
#' (\code{\link{retrieve_longest_isoforms}}, \code{\link{select_orthologs}}, \code{\link{generate_ortholog_tables}},
#' \code{\link{extract_features}}) have their own reader since version 0.5.0, read lines of any length and
#' \strong{no longer call this function}: until then they removed these lines, and with them real features (NCBI
#' TAIR10.1: proteins were reported as not annotated and a shorter isoform was chosen for one gene).
#'
#' \code{check_annotation()} reports how many such lines a file has. A copy without them is only written on
#' request, for other software that cannot read them.
#' @param annotation_file a file path tp the annotation file.
#' @param remove_annotation_outliers shall a copy of \code{annotation_file} without the over-long lines be written?
#' Default is \code{FALSE}: the lines are only reported.
#' If \code{TRUE}, the copy is written to
#' \code{file.path(tempdir(), "orthologr_annotation_without_outliers")} under the same file name (without \code{.gz}) and
#' its path is returned; a warning says how many lines - that is, annotated features - the copy lacks.
#' The \code{annotation_file} itself is not changed (until version 0.4.3 it was overwritten).
#' The removed lines are stored in \code{\link{tempdir}} for further exploration.
#' @return The path of the annotation file: \code{annotation_file} itself when it holds no
#' over-long lines (or when \code{remove_annotation_outliers = FALSE}), otherwise the path of the copy.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # download an example annotation file from NCBI RefSeq
#' Ath_path <- biomartr::getGFF(organism = "Arabidopsis thaliana")
#' # how many over-long lines does it have?
#' orthologr::check_annotation(Ath_path)
#' # orthologr reads the file as it is
#' # a copy without these lines, for software that cannot read them (their features are lost!)
#' Ath_clean <- check_annotation(Ath_path, remove_annotation_outliers = TRUE)
#' }
#' @export

check_annotation <- function(annotation_file, remove_annotation_outliers = FALSE) {
        
        if (!file.exists(annotation_file))
                stop("The file '", annotation_file, "' does not seem to exist. Please provide a valid path to the annotation file.", call. = FALSE)
        check_flag(remove_annotation_outliers, "remove_annotation_outliers")
        
        message("Reading '", annotation_file, "' ...")
        
        # The file is read in blocks: only the numbers of the over-long lines are kept.
        outlier_lines <- annotation_long_lines(annotation_file)
        
        if (length(outlier_lines) == 0)
                return(annotation_file)
        
        if (!remove_annotation_outliers) {
                warning("Your annotation file '", annotation_file, "' contains ", length(outlier_lines),
                        " line", if (length(outlier_lines) != 1) "s", " that have more than 65000 characters. ",
                        "The functions of orthologr read such lines; some other import functions (e.g. of rtracklayer) do not. ",
                        "Only if you need the file for such software: 'remove_annotation_outliers = TRUE' writes a copy ",
                        "without these lines to tempdir() and returns its path (your file is not changed) - the features ",
                        "annotated in these lines are then missing from the copy.", call. = FALSE)
                return(annotation_file)
        }
        
        message("Writing a copy of '", annotation_file, "' without the lines of more than 65000 characters ...")
        
        # The annotation file of the user is left as it is: the other lines are written
        # to a new file in tempdir() that has the same name and extension (without '.gz', the new file
        # is not compressed). The path of the new file is returned.
        out_dir <- file.path(tempdir(), "orthologr_annotation_without_outliers")
        if (!dir.exists(out_dir))
                dir.create(out_dir, recursive = TRUE)
        new_file <- file.path(out_dir, sub("\\.gz$", "", basename(annotation_file)))
        
        outlier_output_file <- file.path(tempdir(), paste0(basename(annotation_file),"_outlier_lines.txt"))
        annotation_copy_lines(annotation_file, new_file, outlier_output_file, outlier_lines)
        
        warning(length(outlier_lines), " line", if (length(outlier_lines) != 1) "s", " of more than 65000 characters ",
                if (length(outlier_lines) != 1) "are" else "is", " not in the copy '", new_file, "': the features annotated there (often CDS or mRNA rows) ",
                "are missing from it. Do not use the copy with orthologr, which reads the original file.", call. = FALSE)
        message("The annotation file without the outlier lines has been stored at '", new_file,
                "' ('", annotation_file, "' was not changed).")
        message("The outlier lines were stored in a separate file to explore at '", outlier_output_file, "'.")
        return(new_file)
}


# Open a text file for reading, also when it is gzip-compressed.
open_text_file <- function(file) {
        if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, open = "r") else file(file, open = "r")
}

# The numbers of the lines of `file` that are longer than `max_chars` characters (read in blocks of `n` lines).
annotation_long_lines <- function(file, n = 200000L, max_chars = 65000L) {
        con <- open_text_file(file)
        on.exit(close(con), add = TRUE)
        res <- integer(0)
        done <- 0L
        repeat {
                x <- readLines(con, n = n, warn = FALSE)
                if (length(x) == 0)
                        break
                long <- which(nchar(x, type = "bytes") > max_chars)
                if (length(long) > 0)
                        res <- c(res, done + long)
                done <- done + length(x)
        }
        res
}

# Copy `file` to `kept_file` without the lines `drop` (line numbers), which go to `dropped_file`.
annotation_copy_lines <- function(file, kept_file, dropped_file, drop, n = 200000L) {
        con <- open_text_file(file)
        on.exit(close(con), add = TRUE)
        kept <- file(kept_file, open = "w")
        on.exit(close(kept), add = TRUE)
        dropped <- file(dropped_file, open = "w")
        on.exit(close(dropped), add = TRUE)
        done <- 0L
        repeat {
                x <- readLines(con, n = n, warn = FALSE)
                if (length(x) == 0)
                        break
                is_drop <- is.element(done + seq_along(x), drop)
                writeLines(x[!is_drop], kept)
                writeLines(x[is_drop], dropped)
                done <- done + length(x)
        }
        invisible(kept_file)
}
