#' @title Options and environment variables of orthologr
#' @description Every option (\code{options()}) and environment variable the package reads, with its default.
#' Options are read when a function runs, so they can be set at any time, e.g.
#' \code{options(orthologr.best_hit_rule = "evalue")}.
#' @section Options that restore results of orthologr 0.4.3:
#' Version 0.5.0 changed some defaults and fixed errors that changed numbers. Each of these changes has a switch.
#' \describe{
#' \item{\code{orthologr.best_hit_rule}}{\code{"bitscore"} (default) or \code{"evalue"}: how the best hit of a query
#' is chosen (\code{\link{diamond_best}}, \code{\link{blast_best}}, \code{\link{diamond_rec}}, \code{\link{blast_rec}},
#' \code{\link{orthologs}}, \code{\link{dNdS}}, \code{\link{divergence_stratigraphy}}, \code{\link{orthogroups}};
#' the default of their argument \code{best_hit_rule}). \code{"evalue"} is the rule of 0.4.3.}
#' \item{\code{orthologr.legacy_pairwise_aln}}{\code{FALSE} (default). \code{TRUE}: \code{\link{pairwise_aln}}
#' (and so \code{dNdS()}) aligns two proteins with the code of 0.4.3, which gives the alignments of that version
#' (only in versions with the faster pairwise alignment; see NEWS, "Changes to pairwise alignments").}
#' \item{\code{orthologr.legacy_comeron_saturation}}{\code{FALSE} (default). \code{TRUE}: Comeron's method returns a
#' finite dN / dS for saturated gene pairs as 0.4.3 did (see \code{\link{substitutionrate}}).}
#' \item{\code{orthologr.legacy_kaks_zero_dn}}{\code{FALSE} (default). \code{TRUE}: a gene pair without
#' nonsynonymous change is returned with \code{dN = NA}, \code{dNdS = 0} by the KaKs_Calculator methods, as in 0.4.3
#' (default: \code{dN = 0}; see \code{\link{substitutionrate}}).}
#' \item{\code{orthologr.dS_zero}}{\code{"drop"} (default, as in 0.4.3) or \code{"keep"}: the default of the
#' argument \code{dS_zero} of \code{\link{filter_dNdS}} and \code{\link{divergence_stratigraphy}}.}
#' }
#' @section Other options:
#' \describe{
#' \item{\code{orthologr.cluster_type}}{\code{"FORK"} (default on Linux and macOS) or \code{"PSOCK"} (default and only
#' choice on Windows): the kind of R processes used with \code{comp_cores > 1} in \code{dNdS()} and
#' \code{dNdS_orthogroups()}.}
#' \item{\code{orthologr.keep_tmp}}{\code{FALSE} (default). \code{TRUE}: the folder with the intermediate files of
#' a call (FASTA files of the pairs, alignments, program output) is not removed, and a message names it.}
#' \item{\code{orthologr.tmpdir}}{\code{NULL} (default: \code{tempdir()}): the folder in which functions called on
#' their own (\code{pairwise_aln()}, \code{multi_aln()}, \code{codon_aln()}, \code{substitutionrate()},
#' \code{diamond()}, \code{blast()}, ...) keep the folders \code{_alignment}, \code{_blast_db}, \code{_diamond_db}
#' and \code{_calculation}. \code{dNdS()}, \code{dNdS_orthogroups()} and \code{orthogroups()} set it to a folder of
#' their own for the time of the call and restore it; it is rarely useful to set it by hand.}
#' \item{\code{orthologr.corrupted_cds_dir}}{\code{NULL} (default: the working directory, as in 0.4.3): the folder in
#' which the file \code{<name>_corrupted_cds_seqs.fasta} with the coding sequences whose length is not a multiple
#' of 3 is written.}
#' \item{\code{orthologr.isoform_warning}}{\code{TRUE} (default). \code{FALSE}: \code{dNdS()} and
#' \code{orthologs()} do not check whether the input files hold several isoforms per gene.}
#' \item{\code{orthologr.coal_backend}}{only with the COAL integration (\code{dNdS(engine = "coal")},
#' \code{dNdS_coal()}, \code{dNdS_clusters()}): \code{NULL} (default: the R package \pkg{coalr} when it is
#' installed, else the program \code{coal-cpp}), \code{"coalr"} or \code{"binary"}.}
#' \item{\code{orthologr.coal_path}}{only with the COAL integration: the full path of the program \code{coal-cpp}
#' (default \code{""}: the environment variable \code{COAL_CPP}, then the \code{PATH}).}
#' }
#' The option \code{timeout} of R is raised to at least 600 seconds while \code{diamond_install()} downloads DIAMOND.
#' @section Environment variables:
#' \describe{
#' \item{\code{COAL_CPP}}{only with the COAL integration: the full path of the program \code{coal-cpp}, used when it
#' is not on the \code{PATH} and \code{options(orthologr.coal_path = )} is not set.}
#' \item{\code{PATH}}{the external programs (DIAMOND, BLAST+, perl, KaKs_Calculator, Clustal Omega, ...) are looked
#' up there unless the \code{*_path} argument of a function names their folder.}
#' }
#' @section Parameter strings:
#' The arguments that take further command-line parameters for an external program - \code{params}
#' (\code{\link{multi_aln}}, \code{\link{codon_aln}}), \code{aa_aln_params}, \code{aligner_params},
#' \code{diamond_params}, \code{blast_params}, \code{kaks_calc.params}, \code{deepclust_params},
#' \code{realign_params} - are \strong{lists of arguments, not shell text}. The string is split at white space into
#' arguments; single or double quotes keep an argument with spaces together
#' (\code{"--outfmt 6 --masking 'a b'"}), and a backslash takes the next character literally. Each argument is then
#' handed to the program as it is. Characters that mean something to a shell (\code{;}, \code{|}, \code{&},
#' \code{>}, \code{$( )}, backquotes, \code{*}) are not interpreted: \code{params = "--force; touch x"} gives the
#' program the three arguments \code{--force;}, \code{touch} and \code{x}, which it will reject. Redirections and
#' pipes inside such a string therefore do not work (they did, unintended, up to 0.4.3).
#' @name orthologr_options
#' @aliases orthologr_options orthologr-options
#' @seealso \code{\link{dNdS}}, \code{\link{substitutionrate}}, \code{\link{filter_dNdS}}
NULL
