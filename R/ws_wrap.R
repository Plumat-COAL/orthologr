#' @title Wrap whitespace in paths
#' @description Helper function for wrapping white spaces in folder paths.
#' @param paths a string or vector of strings containing paths.
#' @author Hajk-Georg Drost
#' @noRd

ws_wrap <- function(paths){
        
        sep <- .Platform$file.sep
        # (split at every separator, keeping empty parts at the end as stringr::str_split() did)
        paths <- lapply(paths, function(p) {
                parts <- strsplit(paste0(p, sep, "x"), sep, fixed = TRUE)[[1]]
                parts[-length(parts)]
        })
        
        ws.paths <- unlist(lapply(paths, function(x) {
                
                paste(unlist(sapply(x, function(word) {
                        if (grepl(" ", word, fixed = TRUE)) 
                                return(paste0("'", word, "'")) else 
                                        return(word)})), collapse = .Platform$file.sep)
                
        }))
        
        return(ws.paths) 
}
