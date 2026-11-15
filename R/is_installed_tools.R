# Internal checks that an external program can be started. Each program is started once per R session (and per
# `path`); see tool_available(). All return TRUE or stop with a message that says what is missing.

not_installed_message <- function(name) {
        paste0("It seems like you don't have ", name, " installed locally on your machine or the PATH variable to the ",
               name, " program is not set correctly.")
}

is_installed_clustalw <- function(path = NULL) {
        tool_available("clustalw2", version_args = "-help", pattern = "CLUSTAL", path = path,
                       not_found = not_installed_message("ClustalW"))
        return(TRUE)
}

is_installed_clustalo <- function(path = NULL) {
        tool_available("clustalo", version_args = "--help", pattern = "Clustal", path = path,
                       not_found = not_installed_message("Clustal-Omega"))
        return(TRUE)
}

is_installed_mafft <- function(path = NULL) {
        tool_available("mafft", version_args = "-v", pattern = "MAFFT", path = path,
                       not_found = not_installed_message("MAFFT"))
        return(TRUE)
}

is_installed_muscle <- function(path = NULL) {
        tool_available("muscle", version_args = "-version", pattern = "muscle", path = path,
                       not_found = not_installed_message("MUSCLE"))
        return(TRUE)
}

is_installed_tcoffee <- function(path = NULL) {
        tool_available("t_coffee", version_args = "-version", pattern = "T-COFFEE", path = path,
                       not_found = not_installed_message("T-COFFEE"))
        return(TRUE)
}

is_installed_kaks_calculator <- function(path = NULL) {
        tool_available("KaKs_Calculator", version_args = "-h", pattern = "KaKs_Calculator", path = path,
                       not_found = paste0(not_installed_message("KaKs_Calculator"),
                                          " Please read the INSTALLATION vignette for details."))
        return(TRUE)
}
