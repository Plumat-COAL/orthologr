is_installed_blast <- function(path = NULL) {
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        
        # test if a valid BLAST version is installed (the program is started once per session)
        blast_version(path = path)
        return(TRUE)
}
