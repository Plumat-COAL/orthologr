#' @title Check whether an annotation file contains outlier lines
#' @description Some annotation files include lines with character lengths greater than 65000. This causes problems when trying to import such annotation files into R using \href{https://rdrr.io/bioc/rtracklayer/man/export.html}{\code{import}}.
#' To overcome this issue, this function screens for such lines
#' in a given annotation file and removes these lines so that
#' \href{https://rdrr.io/bioc/rtracklayer/man/export.html}{\code{import}} can handle the file.
#' @param annotation_file a file path tp the annotation file.
#' @param remove_annotation_outliers shall outlier lines be removed from the input \code{annotation_file}? 
#' If yes, a copy of the annotation file without these lines is written to
#' \code{file.path(tempdir(), "orthologr_annotation_without_outliers")} under the same file name (without \code{.gz}) and
#' its path is returned. The \code{annotation_file} itself is not changed (until version 0.4.3 it was overwritten).
#' The removed outlier lines are stored in \code{\link{tempdir}} for further exploration.
#' @return The path of the annotation file that can be imported: \code{annotation_file} itself when it holds no
#' outlier lines (or when \code{remove_annotation_outliers = FALSE}), otherwise the path of the cleaned copy.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # download an example annotation file from NCBI RefSeq
#' Ath_path <- biomartr::getGFF(organism = "Arabidopsis thaliana")
#' # run annotation file check on the downloaded file
#' orthologr::check_annotation(Ath_path)
#' # several outlier lines were detected, thus we re-run the
#' # function using 'remove_annotation_outliers = TRUE'
#' # to get a copy of the file without the outlier lines
#' Ath_clean <- check_annotation(Ath_path, remove_annotation_outliers = TRUE)
#' }
#' @export

check_annotation <- function(annotation_file, remove_annotation_outliers = FALSE) {
        
        if (!file.exists(annotation_file))
                stop("The file '", annotation_file, "' does not seem to exist. Please provide a valid path to the annotation file.", call. = FALSE)
        
        if (suppressWarnings(requireNamespace("benchmarkme"))) {
                ram <- check_ram()
                available_ram <- ram - ((ram / 100) * 20)
                annotation_file_size <-
                        file.info(annotation_file)$size / 1000000000
                
                if (annotation_file_size > available_ram)
                        stop(
                                "The size of the annotation file '",
                                annotation_file,
                                " exceeds the available RAM (memory) of your computer. Please try running this function with more memory.",
                                call. = FALSE
                        )
        }
        
        message("Importing '", annotation_file, "' ...")
        
        import_annotation <- readr::read_lines(annotation_file)
        outlier_lines <- which(nchar(import_annotation) > 65000)
        
        if (length(outlier_lines) == 0)
                return(annotation_file)
        
        if (length(outlier_lines) > 0 & !remove_annotation_outliers) {
                warning("Your annotation file '", annotation_file, "' contains ", length(outlier_lines), " lines that have more than 65000 characters.",
                        "This will cause the annotation import function to fail importing the annotation file into R.",
                        " If you wish to remove the outlier lines then re-run this function and specify 'remove_annotation_outliers = TRUE': a copy of the annotation file without these lines is then written to tempdir() and its path is returned (your file is not changed). A file with the removed lines will also be stored in the temp directory for you to explore.")
                return(annotation_file)
        }
        
        if (length(outlier_lines) > 0 & remove_annotation_outliers) {
                message("Reading annotation file '", annotation_file, "' and removing all outlier lines with number of characters greater 65000 ...")
                
                # The annotation file of the user is left as it is: the lines that can be imported are written
                # to a new file in tempdir() that has the same name and extension (without '.gz', the new file
                # is not compressed). The path of the new file is returned.
                out_dir <- file.path(tempdir(), "orthologr_annotation_without_outliers")
                if (!dir.exists(out_dir))
                        dir.create(out_dir, recursive = TRUE)
                new_file <- file.path(out_dir, sub("\\.gz$", "", basename(annotation_file)))
                
                readr::write_lines(import_annotation[-outlier_lines], new_file)
                
                message("The annotation file without the outlier lines has been stored at '", new_file,
                        "' ('", annotation_file, "' was not changed).")
                
                outlier_output_file <- file.path(tempdir(), paste0(basename(annotation_file),"_outlier_lines.txt"))
                
                readr::write_lines(import_annotation[outlier_lines], outlier_output_file)
                message("The outlier lines were stored in a separate file to explore at '", outlier_output_file, "'.")
                return(new_file)
        }
        
}
