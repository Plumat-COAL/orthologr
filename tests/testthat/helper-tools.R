# Which external programs are there? (asked once for all test files)
has_program <- function(x) nzchar(Sys.which(x))

blast_available   <- has_program("blastp") && has_program("makeblastdb")
diamond_available <- has_program("diamond")
perl_available    <- has_program("perl")
kaks_available    <- has_program("KaKs_Calculator")

thal_cds <- function() system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
lyra_cds <- function() system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
