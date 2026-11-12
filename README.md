orthologr <img src="man/figures/logo.png" align="right" height="174" width="150" />
=========

## Comparative Genomics with R

### Motivation

The comparative method is a powerful approach in genomics research. Based on our knowledge about the phylogenetic relationships between species, we can study the evolution, diversification, and constraints of biological processes by comparing genomes, genes, and other genomic loci across species. The `orthologr` package aims to provide a framework to perform __large scale__ comparative genomics studies with R. `Orthologr` aims to be as easy to use as possible - from genomic data retrieval to orthology inference and dNdS estimation between several genomes.

In combination with the R package [biomartr](https://github.com/ropensci/biomartr), users can retrieve genomes, proteomes, or coding sequences for several species and use them as input
for orthology inference and dN/dS estimation with `orthologr`. The advantage of using `biomartr` in combination with `orthologr` is that
users can join the new wave of research that promotes and facilitates
[computational reproducibility in genomics studies](https://github.com/ropensci/biomartr#short-package-description) and solve the
issue of comparing genomes with different genome assembly qualities (also referred to as [genome version crisis](https://github.com/ropensci/biomartr#short-package-description)).

You can find a detailed list of all `orthologr` functions here: https://drostlab.github.io/orthologr/reference/index.html

### Citation

**Please cite the following paper in which I introduce `orthologr` when using this package for your own research. This will allow me to continue
working on this software tool and will motivate me to extend its functionality and usability in the next years. Many thanks in advance :)**

> Drost et al. 2015. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ 32 (5): 1221-1231. [doi:10.1093/molbev/msv012](https://doi.org/10.1093/molbev/msv012)

### Short package description

In detail, `orthologr` allows users to perform orthology inference and dN/dS estimation between two genomes or between several genomes. The following methods
to infer orthologous relationships between genes of entire genomes are available in this package:

- __DIAMOND2 reciprocal best hit__ (the default): [`diamond_rec()`](https://drostlab.github.io/orthologr/reference/diamond_rec.html)
- __DIAMOND2 best hit__: [`diamond_best()`](https://drostlab.github.io/orthologr/reference/diamond_best.html)
- [__BLAST reciprocal best hit__](https://drostlab.github.io/orthologr/articles/blast.html#the-blast_rec-function)
- [__BLAST best hit__](https://drostlab.github.io/orthologr/articles/blast.html#the-blast_best-function)
- __DIAMOND DeepClust__ protein clustering: [`deepclust()`](https://drostlab.github.io/orthologr/articles/deepclust.html)

(`orthofinder2()`, a thin wrapper that starts OrthoFinder, is soft-deprecated; see its help page.)

The most useful implementation in `orthologr` is the ability to compute synonymous versus non-synonymous substitution rates (dN/dS)
for all orthologous genes between two __entire genomes__.
Available dN/dS estimation methods are:

- "Comeron": Comeron, J.M. (1995) - the default; built into the package (no external program)
- "Li": Li, W.H. (1993) - computed with `seqinr::kaks()` (no external program)

and, computed with [KaKs_Calculator](https://anaconda.org/bioconda/kakscalculator2) (which has to be installed):

- "NG": Nei, M. and Gojobori, T. (1986)
- "LWL": Li, W.H., et al. (1985)
- "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
- "MLWL" (Modified LWL): Tzeng, Y.H., et al. (2004)
- "YN": Yang, Z. and Nielsen, R. (2000)
- "MYN" (Modified YN): Zhang, Z., et al. (2006)
- "GMYN": Wang, D.P., et al. Biology Direct. (2009)
- "GY": Goldman, N. and Yang, Z. (1994)
- "MS": (Model Selection): based on a set of candidate models, Posada, D. (2003) 
- "MA" (Model Averaging): based on a set of candidate models, Posada, D. (2003)
- "ALL": All models together 

Please find more details [here](https://drostlab.github.io/orthologr/articles/dNdS_estimation.html).

## Install `orthologr`

`orthologr` needs R >= 4.4 and Bioconductor >= 3.19 (for the `pwalign` package).

```r
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install the Bioconductor packages orthologr imports
BiocManager::install(c("Biostrings", "pwalign"))

# install orthologr from GitHub (this also installs its CRAN dependencies)
install.packages("remotes")
remotes::install_github("drostlab/orthologr")
```

Optional packages, needed by some functions only (a message tells you when one is missing):

```r
BiocManager::install("rtracklayer")           # only extract_features() on imported annotations and the *_do() reference code;
                                              # select_orthologs() and retrieve_longest_isoforms() read GFF3 / GTF files themselves
install.packages(c("ggplot2", "scales"))      # the plot_*() functions
install.packages("ape")                       # promotor_divergence_estimation()
remotes::install_github("drostlab/metablastr") # lncRNA orthology: orthologs_lnc(), map_generator_lnc()
```

### External programs

`dNdS()` and the orthology inference functions run programs that have to be installed on your computer:

- __DIAMOND2__ (>= 2.1) for the default sequence search. Either run `orthologr::diamond_install()` once (it
  downloads the DIAMOND release for your operating system into a user-level folder; Linux x86_64, macOS, Windows),
  or install it yourself, e.g. `conda install -c bioconda diamond`.
- __Perl__ for the codon alignments (the PAL2NAL script is part of the package).
- Optional: __BLAST+__ (`aligner = "blast"`), __KaKs_Calculator__ 2.0 (dN/dS methods other than `"Comeron"` and
  `"Li"`), and for `aa_aln_type = "multiple"` one of __ClustalW2__, __Clustal Omega__, __MAFFT__, __MUSCLE__ (3.8 or 5),
  __T-Coffee__.

The file [`environment.yml`](environment.yml) describes a conda environment with R, the R packages and all of these
programs: `conda env create -f environment.yml` (or `mamba` / `micromamba`).
See [Install Prerequisite Tools](https://drostlab.github.io/orthologr/articles/Install.html) for details.

## Use Cases

Learn `orthologr` by reading these tutorials: 

- [Install Prerequisite Tools](https://drostlab.github.io/orthologr/articles/Install.html)
- [Perform Orthology Inference and dNdS Estimation](https://drostlab.github.io/orthologr/articles/dNdS_estimation.html)
- [Perform Only Orthology Inference](https://drostlab.github.io/orthologr/articles/orthology_inference.html)
- [Perform BLAST Searches](https://drostlab.github.io/orthologr/articles/blast.html)
- [Perform Divergence Stratigraphy](https://drostlab.github.io/orthologr/articles/divergence_stratigraphy.html)
- [Perform Sequence Alignments](https://drostlab.github.io/orthologr/articles/sequence_alignments.html)



## Example

### Small example with internal dataset
```r
library(orthologr)

# Detect orthologous genes between a query species and a subject species
# and compute the synonymous versus non-synonymous substitution rates (dN/dS)
# following this paradigm:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) single-core processing (comp_cores = 1)
dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     delete_corrupt_cds = TRUE, # coding sequences that cannot be divided by 3 (triplets) will be removed
     ortho_detection = "RBH", # perform DIAMOND best reciprocal hit orthology inference
     aa_aln_type     = "pairwise", # perform pairwise global alignments of AA seqs 
     aa_aln_tool     = "NW", # using Needleman-Wunsch
     codon_aln_tool  = "pal2nal", # perform codon alignments using the tool Pal2Nal
     dnds_est.method = "Comeron", # use Comeron's method for dN/dS inference
     comp_cores      = 1) # number of compute cores
```

```
# A tibble: 20 x 24
   query_id    subject_id           dN    dS   dNdS perc_identity num_ident_matches alig_length
   <chr>       <chr>             <dbl> <dbl>  <dbl>         <dbl>             <int>       <int>
 1 AT1G01010.1 333554|PACid:1~ 0.108   0.261 0.414           73.2               347         474
 2 AT1G01020.1 470181|PACid:1~ 0.0402  0.104 0.388           91.1               224         246
 3 AT1G01030.1 470180|PACid:1~ 0.0150  0.126 0.118           93.3               335         359
 4 AT1G01040.1 333551|PACid:1~ 0.0168  0.118 0.142           93.4              1840        1969
 5 AT1G01050.1 909874|PACid:1~ 0       0.175 0              100                 213         213
 6 AT1G01060.3 470177|PACid:1~ 0.0442  0.113 0.392           87.5               567         648
 7 AT1G01070.1 918864|PACid:1~ 0.0183  0.106 0.173           92.6               339         366
 8 AT1G01080.1 909871|PACid:1~ 0.0340  0.106 0.322           89.3               268         300
 9 AT1G01090.1 470171|PACid:1~ 0.00984 0.215 0.0458          96.8               420         434
10 AT1G01110.2 333544|PACid:1~ 0.0325  0.122 0.266           87.7               463         528
11 AT1G01120.1 918858|PACid:1~ 0.00307 0.133 0.0232          99.2               525         529
12 AT1G01140.3 470161|PACid:1~ 0.00567 0.131 0.0432          98.5               446         453
13 AT1G01150.1 918855|PACid:1~ 0.132   0.210 0.628           72.6               207         285
14 AT1G01160.1 918854|PACid:1~ 0.119   0.240 0.493           78.8               141         179
15 AT1G01170.2 311317|PACid:1~ 0       0.306 0               85.6                83          97
16 AT1G01180.1 909860|PACid:1~ 0.0370  0.181 0.204           92.6               287         310
17 AT1G01190.1 311315|PACid:1~ 0.0287  0.162 0.177           94.2               502         533
18 AT1G01200.1 470156|PACid:1~ 0.0190  0.168 0.114           95.8               228         238
19 AT1G01210.1 311313|PACid:1~ 0.0207  0.154 0.134           95.3               102         107
20 AT1G01220.1 470155|PACid:1~ 0.0152  0.153 0.0992          96.5              1019        1056
# i 16 more variables: mismatches <int>, gap_openings <int>, n_gaps <int>, pos_match <int>,
#   ppos <dbl>, q_start <int>, q_end <int>, q_len <int>, qcov <dbl>, qcovhsp <dbl>,
#   s_start <int>, s_end <int>, s_len <int>, evalue <dbl>, bit_score <dbl>, score_raw <dbl>
```

When running your own query file, please specify `query_file = "path/to/your/cds.fasta` instead of `system.file(..., package = "orthologr")`. The command `system.file(..., package = "orthologr")` merely references the path to the example file stored in the `orthologr` package itself.


### Example: Computing dN/dS values for all orthologous genes between two genomes

First, users can retrieve all coding sequences from entire genomes using the [biomartr](https://github.com/ropensci/biomartr) package ([see details here](https://docs.ropensci.org/biomartr/articles/Sequence_Retrieval.html#cds-retrieval)).

```r
install.packages("biomartr")
library(biomartr)
# download all coding sequences for Mus musculus
Mmusculus_file <- biomartr::getCDS(organism = "Mus musculus", path = getwd())
# download all coding sequences for Homo sapiens
Hsapiens_file <- biomartr::getCDS(organism = "Homo sapiens", path = getwd())
# compute dN/dS values for Homo sapiens versus Mus musculus
Hs_vs_Mm_dNdS <- 
  dNdS(query_file      = Hsapiens_file,
       subject_file    = Mmusculus_file,
       delete_corrupt_cds = FALSE,
       ortho_detection = "RBH", 
       aa_aln_type     = "pairwise",
       aa_aln_tool     = "NW", 
       codon_aln_tool  = "pal2nal", 
       dnds_est.method = "Comeron", 
       comp_cores      = 1 )
# store result in Excel readable csv file
install.packages("readr")
readr::write_excel_csv(Hs_vs_Mm_dNdS, "Hs_vs_Mm_dNdS.csv")
```

__This way, users can compute dN/dS values for any pairwise genome comparison.__

### On Windows Systems

On Windows the package is built from source, which needs [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
(a separate installer, not an R package). After installing Rtools:

```r
install.packages("remotes")
remotes::install_github("drostlab/orthologr")
```

The article [Setup on Windows](https://drostlab.github.io/orthologr/articles/setup_windows.html) describes how to
get DIAMOND, Perl and the optional programs running on Windows.


## Interfaces implemented in `orthologr`:

The complete list of functions is in the [reference](https://drostlab.github.io/orthologr/reference/index.html).

### dN/dS estimation

* `dNdS()`: Compute dN/dS values for the orthologous genes of two organisms
* `divergence_stratigraphy()`, `divergence_map()`, `filter_dNdS()`: divergence strata from dN/dS values
* `map_generator_dnds()`, `read.dnds.tbl()`: dN/dS tables of one query against many subject organisms
* `substitutionrate()`: dN/dS of one codon alignment

### Sequence searches and orthology inference

* `diamond()`, `diamond_best()`, `diamond_rec()`, `set_diamond()`: DIAMOND2 searches, best hits, reciprocal best hits
* `diamond_install()`, `diamond_path()`: install and locate DIAMOND2
* `blast()`, `blast_best()`, `blast_rec()`, `set_blast()`: the same with BLAST+
* `orthologs()`: orthology inference with one of the methods above
* `deepclust()`, `deepclust_annotate()`, `deepclust_realign()`: protein clustering with DIAMOND DeepClust
* `orthologs_lnc()`, `map_generator_lnc()`, `lnc_map_counts()`: lncRNA homologs (needs `metablastr`)

### Pairwise and multiple sequence alignments

* `multi_aln()`: Compute Multiple Sequence Alignments based on the `clustalw`, `t_coffee`, `muscle`, `clustalo`, and `mafft` programs.
* `pairwise_aln()`: Compute Pairwise Alignments
* `codon_aln()`: Compute a Codon Alignment

### Gene loci, isoforms and ortholog tables

* `retrieve_longest_isoforms()`, `retrieve_longest_isoforms_all()`: keep the longest isoform per gene locus
* `check_annotation()`, `extract_features()`: report over-long lines of an annotation file, features of a locus
* `select_orthologs()`, `generate_ortholog_tables()`, `generate_ortholog_tables_all()`, `import_ortholog_tables_all()`,
  `retrieve_core_orthologs()`: one representative ortholog per gene locus, core orthologs across species
* `plot_pairwise_orthologs()`, `plot_diverse_homology_thresholds()`, `plot_diverse_homology_thresholds_core_orthologs()`

### Read and Write CDS, Genomes, and Proteomes

* `read.cds()`: Read the CDS of a given organism
* `read.genome()`: Read the genome of a given organism
* `read.proteome()`: Read the proteome of a given organism
* `write.proteome()`: Save a proteome in fasta format
* `translate_cds_to_protein()`, `translate_cds_to_protein_all()`: translate CDS files


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/drostlab/orthologr/issues

## Licenses

The `orthologr` package includes source code that has been published under following licenses:

### gestimator

The header of `src/gestimator.cpp` and of the libsequence headers in `src/` reads (quoted as shipped):

```
All files included in `orthologr` that were taken from gestimator are 
also part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.


Modified by Sarah Scharfenberg and Hajk-Georg Drost (2014) to work 
in orthologr without using external libraries from libsequence.

All changes are also free under the terms of GNU General Public License
version 2 of the License, or any later version.

```

### KaKs_Calculator

`orthologr` runs the program KaKs_Calculator when it is installed; none of its code is part of the package. (The
script `parseFastaIntoAXT.pl` of the KaKs_Calculator project, GNU GPL v3, was shipped in `inst/KaKs_Calc_parser`
up to version 0.4.3; it is no longer used and was removed in 0.5.0.)



