# from rdiamond (github.com/drostlab/rdiamond)
is_installed_diamond <- function(diamond_exec_path = NULL) {
        if (!is.null(diamond_exec_path))
                diamond_exec_path <- normalizePath(diamond_exec_path, mustWork = FALSE)
        
        # test if a valid DIAMOND version is installed (the program is started once per session)
        diamond_version(path = diamond_exec_path)
        return(TRUE)
}
