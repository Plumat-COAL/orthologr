#' @title Count number of orthologous lncRNAs per pairwise species comparison
#' @description This function takes a orthologous lncRNA map generated with \code{\link{map_generator_lnc}}
#' as input and groups the table by species to count the number of orthologous lncRNAs per pairwise species comparison
#' using the reference (query) species specified in \code{\link{map_generator_lnc}}.
#' @param lnc_map a orthologous lncRNA map generated with \code{\link{map_generator_lnc}}.
#' @param species_order  character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
#' The species names must match with the species names present in \code{map_generator_lnc}.
#' @author Hajk-Georg Drost
#' @export  
lnc_map_counts <- function(lnc_map, species_order) {
        
        # lncRNA maps name the species column 'species', the ortholog tables of this package 'subject_species'
        species_col <- intersect(c("species", "subject_species"), names(lnc_map))
        if (length(species_col) == 0)
                stop("'lnc_map' needs a column named 'species' or 'subject_species'.", call. = FALSE)
        all_lncRNAs_df <- dplyr::summarize(dplyr::group_by(lnc_map, .data[[species_col[1]]]), n_orthologs = dplyr::n())
        names(all_lncRNAs_df)[1] <- "subject_species"
        
        all_lncRNAs_df$subject_species <- factor(
                all_lncRNAs_df$subject_species,
                levels = species_order
        )
        
        return(all_lncRNAs_df)
}