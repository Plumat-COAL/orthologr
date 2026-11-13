#' @title Import a dNdS table generated with \code{dNdS}
#' @description This function reads a file that stores a dnds table
#' generated with \code{\link{dNdS}}.
#' @param file file path to dNdS table with \code{;} separated columns.
#' @author Hajk-Georg Drost
#' @examples
#' # a small dN/dS table as dNdS() returns it (first columns)
#' dNdS_tbl <- data.frame(query_id   = c("AT1G01010.1", "AT1G01020.1", "AT1G01030.1"),
#'                        subject_id = c("333554|PACid:16033839", "470181|PACid:16064328",
#'                                       "470180|PACid:16054974"),
#'                        dN = c(0.1060, 0.0402, 0.0150), dS = c(0.2540, 0.1340, 0.1160),
#'                        dNdS = c(0.4180, 0.3010, 0.1290))
#'
#' # save the dNdS table as ';' column separated file
#' file <- file.path(tempdir(), "dNdS_tbl.csv")
#' utils::write.table(dNdS_tbl, file, sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)
#'
#' # import the dNdS table into the R session
#' dNdS_tbl_import <- read.dnds.tbl(file)
#' dNdS_tbl_import
#' unlink(file)
#' @export
read.dnds.tbl <- function(file) {
        
        if (!file.exists(file))
                orthologr_abort(class = "bad_input", "The file '",file,"' does not seem to exist. Please specify a valid path to your dnds table.", call. = FALSE)
        
        res <- tibble::as_tibble(readr::read_delim(
                file,
                col_names = TRUE,
                delim = ";",
                col_types = readr::cols("query_id" = readr::col_character(),
                                        "subject_id"= readr::col_character(),
                                        "dN" = readr::col_double(),
                                        "dS"  = readr::col_double(),
                                        "dNdS" = readr::col_double(),
                                        "perc_identity" = readr::col_double(),
                                        "alig_length" = readr::col_integer(),
                                        "mismatches" = readr::col_integer(),
                                        "gap_openings" = readr::col_integer(),
                                        "q_start" = readr::col_integer(),
                                        "q_end" = readr::col_integer(),
                                        "s_start" = readr::col_integer(),
                                        "s_end" = readr::col_integer(),
                                        "evalue" = readr::col_double(),
                                        "bit_score" = readr::col_double()
                                        ))
        )
        
        if (nrow(res) == 0)
                stop("The file '",file,"' seems to be empty. Please provide a non-empty dNdS file.", call. = FALSE)
        return(res)
}