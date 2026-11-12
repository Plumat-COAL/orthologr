is_installed_orthofinder <- function(path = NULL) {
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        
        # test if OrthoFinder can be started: 'orthofinder -h' prints its help ("OrthoFinder version ...")
        tool_available(
                "orthofinder", version_args = "-h", pattern = "OrthoFinder", path = path,
                not_found = paste0(
                        "It seems like you don't have OrthoFinder installed locally on your machine or the PATH variable to the OrthoFinder program is not set correctly. ",
                        "Please install it (e.g. 'conda install -c bioconda orthofinder') or use the 'of_path' argument to give the folder that holds the 'orthofinder' executable.")
        )
        return(TRUE)
}
