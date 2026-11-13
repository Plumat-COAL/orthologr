## `orthologr` 0.5.0

Results change in this version: three defaults were changed and several errors that changed numbers were fixed.
`orthologr_legacy_options()` (or `options(orthologr.legacy = TRUE)`) restores the behaviour of 0.4.3 in one line;
the article **"Upgrading from 0.4.3"** (`vignette("upgrading_from_0.4.3")`) gives the measured effect of every
change. Please report the version of orthologr in the methods of analyses based on dN/dS values.

### New

- **`orthogroups()`**: groups of orthologs of 2 to N species from one all-against-all DIAMOND search (reciprocal
  best hits, connected components) or from DIAMOND DeepClust clusters; `print()`, `summary()`,
  `write_orthogroups()`. **`dNdS_orthogroups()`**: dN/dS for the gene pairs within the groups, with a summary per
  group. OrthoFinder is no longer needed for more than two species; `orthofinder2()` is soft-deprecated.
- **Engine COAL** (optional R package `coalr`): `dNdS(engine = "coal")` / **`dNdS_coal()`** hand all gene pairs to
  [COAL](https://github.com/Plumat-COAL/COAL) in one call (protein alignment, codon alignment and estimation in
  memory, on `comp_cores` threads; 1,000-gene example, 1 core: about 6 s against 50 s); **`dNdS_clusters()`** computes
  dN/dS within DeepClust clusters. `codon_aln(tool = "coal")` and `substitutionrate(kaks_tool = "coal")` replace
  PAL2NAL / KaKs_Calculator for single steps. Methods: NG, LWL, MLWL, LPB, MLPB, YN, MYN, and Comeron and Li with
  a `coalr` that has them.
- **Column `status`** (last column of every dN/dS table): `"ok"`, `"dS_zero"`, `"saturated"`, `"no_data"` or
  `"failed"` for each gene pair - the same words with both engines and in the `pairs` tables of
  `dNdS_clusters()` / `dNdS_orthogroups()` (see `?dNdS`). `options(orthologr.status_column = FALSE)` leaves it out.
- **One summary message at the end of a run**: gene pairs by status, median dS, sequences left out before the
  search, and a note when most pairs are saturated or the median dS is above 1.5 (species too distant for
  pairwise dN/dS). With `quiet = FALSE`, a progress message about every tenth of the gene pairs.
- **Record of the run**: `attr(result, "orthologr_run")` on results of `dNdS()`, `dNdS_orthogroups()` and
  `divergence_stratigraphy()`: package version, arguments, every `orthologr.*` option in force (with its default
  when it is not set, `orthologr.ambiguous_cds` included), engine, program versions, date, size and MD5 of the
  input files, counts by status. **`write_dNdS()`** writes a table as TSV with this record in
  `<file>.run.json`; **`read_dNdS()`** reads both back (needs `jsonlite`). `orthogroups()$params` holds the package version.
- **`orthologr_legacy_options()`** and `options(orthologr.legacy = TRUE)`: all settings of 0.4.3 at once. A test
  compares the result with tables stored from 0.4.3 (20-gene example, Comeron and YN). Only the default engine
  gives the numbers of 0.4.3: the COAL engine cannot make the pairwise alignments of 0.4.3 and with these settings
  warns (`orthologr_warning_legacy_engine`; 20-gene example: 10 of 20 pairs differ from 0.4.3, largest
  difference in dN/dS 0.147, r = 0.984).
- **`?orthologr_options`**: every option and environment variable the package reads.
- **Classed errors and warnings**: the common failures have the class `orthologr_error` plus one for their kind -
  `orthologr_error_tool_missing`, `_package_missing`, `_tool_failed` (with `_alignment_failed` for an alignment
  program), `_no_hits`, `_bad_input`, `_bad_argument` - and the warnings about failed gene pairs and isoforms the
  classes `orthologr_warning_failed_pairs` / `_isoforms` (and `_legacy_engine`, above), so that they can be handled with `tryCatch()` without
  matching message text (`?orthologr_options`, section "Errors and warnings"). Messages are unchanged. The
  argument and input-file checks of `dNdS()`, `dNdS_coal()`, `dNdS_clusters()`, `orthogroups()`, `dNdS_orthogroups()`, `divergence_stratigraphy()`,
  `filter_dNdS()`, the alignment functions, `substitutionrate()` and the readers are classed (17 probes of the
  review: 0 classed before, 17 after); a missing COAL is `orthologr_error_package_missing` (with
  `_tool_missing` when neither `coalr` nor coal-cpp is there); `dNdS()` and `dNdS_coal()` check that their input
  files exist before any program runs.
- `divergence_left_out()`: the gene pairs that are not in a divergence map, with the reason.
- `attr(result, "failed_pairs")` and `attr(result, "dropped_sequences")` on dN/dS tables: gene pairs that could
  not be processed and input sequences that were left out, each with the reason.
- A warning when an input file of `dNdS()`, `divergence_stratigraphy()`, `orthologs()` or `orthogroups()` seems
  to hold several transcripts per gene (`options(orthologr.isoform_warning = FALSE)` switches it off).
- `aa_sub_matrix` (`dNdS()`, `pairwise_aln()`): the substitution matrix of pairwise protein alignments.
- `options(orthologr.ambiguous_cds = "keep")`: coding sequences with IUPAC ambiguity codes (N, R, Y, ...) are kept
  (a codon with such a base is translated to X) instead of being left out (default "drop", as in 0.4.3).
- `best_hit_rule` (`dNdS()`, `orthologs()`, `diamond_best()`, `diamond_rec()`, `blast_best()`, `blast_rec()`,
  `divergence_stratigraphy()`, `select_orthologs()`, `orthogroups()`).
- `dS_zero` (`filter_dNdS()`, `divergence_stratigraphy()`): `"keep"` returns gene pairs with dS = 0 instead of
  leaving them out.
- MUSCLE 5 is supported; `"MLPB"` is an accepted `dnds_est.method`.
- `environment.yml` (conda environment with R, the R packages and the external programs); a pkgdown site
  configuration; articles on orthogroups and on the setup on Windows.

### Changed defaults (and how to get the old behaviour)

| Change | Effect (measured) | Old behaviour |
|---|---|---|
| Best hits are chosen by bit score (was: e-value, then alignment length) | human / mouse: 3.2 % of the best hits change; pairs that contradict Ensembl Compara fall from 2.0 % to 0.4 % | `best_hit_rule = "evalue"` or `options(orthologr.best_hit_rule = "evalue")` |
| Pairwise protein alignments: BLOSUM62, full length, stop codon not aligned | 1,000-gene example: dN or dS changes for a third of the gene pairs, mostly slightly (r of dN/dS old vs new 0.988) | `options(orthologr.legacy_pairwise_aln = TRUE)` |
| Comeron's method: `NA` for a saturated dN or dS (was: a finite, too small value) | 1 dN of 995 pairs (example); 17 dN and 1 dS of 1,500 pairs (human / mouse) | `options(orthologr.legacy_comeron_saturation = TRUE)` |
| KaKs_Calculator methods: `dN = 0` for a pair without nonsynonymous change (was `NA`, with `dNdS = 0`) | 11 of 995 pairs; `filter_dNdS()` keeps them | `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| dN/dS tables end with the column `status` | one more column | `options(orthologr.status_column = FALSE)` |
| T-Coffee runs in its default mode (was `-mode accurate`, which needs further programs) | other alignments with `aa_aln_tool = "t_coffee"` | `aa_aln_params = "-mode accurate"` |
| `aa_aln_params` is passed on to the alignment program (was ignored) | only when it is set | do not set it |
| `divergence_stratigraphy(quiet = TRUE)`, as `dNdS()` (was `FALSE`) | 1,000-gene example: 10 messages instead of 45; same result | `quiet = FALSE` |
| `orthogroups()` and `dNdS_orthogroups()` trim coding sequences whose length is not a multiple of 3 (`delete_corrupt_cds = FALSE`, as `dNdS()` does) | - (new functions) | `delete_corrupt_cds = TRUE` |

All of the first five at once: `orthologr_legacy_options()`.

Further behaviour that differs from 0.4.3 and cannot be switched back:

- A gene pair whose alignment or estimation fails no longer ends the run: it is a row of `NA` with status
  `"failed"` and its reason in `attr(, "failed_pairs")`.
- A missing program (perl, the alignment program, KaKs_Calculator) and unsupported argument values are errors
  before the search starts.
- Parameter strings (`aa_aln_params`, `aligner_params`, `diamond_params`, `blast_params`, `params`,
  `kaks_calc.params`, ...) are lists of arguments, not shell text: `;`, pipes and redirections inside them no
  longer work.
- `dNdS()` writes nothing to `stdout`. With `quiet = TRUE` (default) it gives a line at the start and the
  summary at the end (31 messages before); the citation note is a message shown once per session. Input sequences
  that are left out or trimmed to a multiple of 3 are counted in the summary (listed in
  `attr(, "dropped_sequences")` / `attr(, "trimmed_sequences")`) instead of a block of messages about each kind:
  1,000-gene example, 12 messages before (8 about the 12 trimmed sequences), 2 after (also with
  `engine = "coal"`); `quiet = FALSE` shows the old messages.
- Coding sequences with characters other than A, C, G, T are left out before the search, as before, but no
  longer without a word: they are counted in the summary and listed in `attr(, "dropped_sequences")`.
- Intermediate files are kept in a folder of the call inside `tempdir()` and removed afterwards
  (`options(orthologr.keep_tmp = TRUE)` keeps them); no function changes the working directory or writes into it
  unasked (`read.cds()` wrote `<file>_corrupted_cds_seqs.fasta` there).
- `ggplot2`, `scales`, `ape`, `rtracklayer`, `benchmarkme`, `tidyr`, `biomartr`, `metablastr` moved from Imports
  to Suggests (Imports: 21 -> 11 packages); R >= 4.4 is enough.
- `cdd.path` and `ortho_path` of `dNdS()` are deprecated (never used); `ortho_detection = "Orthofinder2"`, `"PO"`
  and `"DELTA"` are errors with an explanation.

### Fixed

Results:

- `dNdS()` mixed up the two species when query and subject file had the same file name in different folders.
- Sequence IDs are returned as they are in the input files. IDs with a hyphen, with `|`, numeric-looking IDs
  (`007`, `1e5`), IDs of more than 15 characters (MAFFT) and the same ID in both files lost rows, came back
  changed or gave `NA`. The dN/dS rows are joined with the hit table by gene pair (was: by query only).
- Li's method returned text columns and took small values (`9.95e-05`) for the saturation marker.
- `filter_dNdS()` and `divergence_stratigraphy()` say how many gene pairs they leave out and why
  (`divergence_left_out()`); gene pairs with dS = 0 are left out as before unless `dS_zero = "keep"`.
- `divergence_map()`: strata are tie-safe, `n_quantile = 49` or `98` no longer gives `NA` for the largest value,
  a failed binning no longer loses the finished dN/dS run, and the caller's table is not changed.
- Comeron's method (gestimator): input that is not a codon alignment is an R error (a sequence shorter than one
  codon ended the R session; sequences of unequal length gave numbers).
- Best hits no longer depend on the order of the hit table.
- `retrieve_longest_isoforms()` and `select_orthologs()` work on real NCBI and Ensembl GFF3 / GTF files
  (proteins were lost without notice, annotation lines of more than 65,000 characters were deleted, GFF3 loci
  were grouped by name).
- `promotor_divergence_of_orthologous_genes()` pairs promoters by ortholog pair.
- `substitutionrate(kaks_calc.params = )` is added to the method instead of replacing it.

Running programs:

- Stale results can no longer be returned: `blast()`, `diamond()`, `deepclust()` and the functions built on them
  check the exit status and the output of every program and work in files of their own.
- File and folder names with spaces work. Programs are looked up once per session.
- Parallel runs use the `parallel` package only; on Windows (PSOCK) arguments given as variables reach the
  worker processes, and worker processes are ended after an interrupt.
- With `R CMD INSTALL --no-byte-compile`, `devtools::load_all()` or `covr` every gene pair of a default
  `dNdS()` run was `NA` in development versions of 0.5.0 (the native routine `gestimator` was registered with five
  instead of six arguments). The registration is now generated by `Rcpp::compileAttributes()`, and a test
  compares it with the R wrappers.
- `multi_aln()`: T-Coffee and MUSCLE examples run with current versions; CLUSTAL files of T-Coffee are read;
  the MAFFT check uses `mafft --version`.
- `multi_aln(tool = "t_coffee")` and `dNdS(aa_aln_tool = "t_coffee")` returned unaligned rows in development
  versions of 0.5.0 (a bug, no option restores it): T-Coffee 12 was given the sequences with `-infile` and took a
  file whose sequences end in `*` (as the translated coding sequences do) for a finished alignment. They are now
  given with `-seq`. 20-gene example, `aa_aln_type = "multiple"`: before 10 of 20 pairs `"saturated"` (dN/dS
  `NA`), 7 of 20 within 0.1 of ClustalW, r = 0.743 (max difference 0.553); after 0 saturated, 20 of 20 within 0.1
  of ClustalW, r = 0.993 (max difference 0.094). 0.4.3 used `-mode accurate`, which fails without network.
- `pairwise_aln()`, `read.genome(format = "fastq")`, `codon_aln(get_aln = TRUE)`,
  `promotor_divergence_estimation()`, `lnc_map_counts()`, `map_generator_dnds()` with one subject file: errors
  on valid input fixed. `check_annotation(remove_annotation_outliers = TRUE)` no longer overwrites the
  annotation file.
- An unknown `aa_sub_matrix` name (`dNdS()`, `pairwise_aln()`) is an error (`orthologr_error_bad_argument`)
  before the search; it ran the whole search and then failed every gene pair (1,000-gene example,
  `aa_sub_matrix = "FOO"`: before 6.4 s and 995 of 995 pairs `NA` with one warning; after an error at once).
  With `engine = "coal"` the name is checked against the matrices of COAL (BLOSUM62, BLOSUM45, BLOSUM80, PAM250).
- A protein file given where coding sequences are expected is an error that says so (was: a DIAMOND error about
  memory).
- `orthogroups()` runs the check for several transcripts per gene; the check no longer warns for a bacterial
  genome with a handful of genes that have two coding sequences.

### Speed (results unchanged)

- Default `dNdS()` on the 1,000-gene example, 1 core: 744 s (0.4.3) -> 50 s; with 2 cores 27 s. Half of the last
  step comes from building the full-length alignment rows directly instead of with `pwalign::alignedPattern()`
  (97 s -> 50 s); the rest from BLOSUM62 scoring, reading and translating every file once, vectorised
  translation and best-hit selection, and not writing a second copy of every alignment.
- `select_orthologs()`, `retrieve_longest_isoforms()` and `deepclust()` read annotations and sequences in one
  pass.

### Checked against other programs

`dNdS()` (default engine: gestimator for Comeron, KaKs_Calculator 2.0 for YN) against the COAL engine
(`engine = "coal"`), on 300 sampled reciprocal-best-hit pairs of each of three species pairs (NCBI RefSeq coding
sequences; 2026-11-13). Columns: gene pairs in both results, pairs with dS > 0 from both, share of these whose dS
agrees within 1 %, median relative difference of dS, Pearson r of dS and of dN/dS, pairs with dS from one engine
only.

| Species pair | Method (default engine) | pairs | dS both | within 1 % | median rel. diff. | r dS | r dN/dS | one only |
|---|---|---|---|---|---|---|---|---|
| S. cerevisiae / C. glabrata | Comeron (gestimator) | 283 | 73 | 63.0 % | 0.0018 | 0.952 | 0.989 | 3 |
| S. cerevisiae / C. glabrata | YN (KaKs_Calculator) | 283 | 38 | 65.8 % | 2e-06 | 0.883 | 0.965 | 245 |
| S. cerevisiae / K. lactis | Comeron (gestimator) | 265 | 68 | 69.1 % | 0.0023 | 0.977 | 0.999 | 1 |
| S. cerevisiae / K. lactis | YN (KaKs_Calculator) | 265 | 31 | 71.0 % | 2.3e-06 | 0.965 | 0.999 | 234 |
| Salmonella LT2 / Shigella flexneri | Comeron (gestimator) | 294 | 274 | 95.3 % | 0 | 0.999 | 0.998 | 0 |
| Salmonella LT2 / Shigella flexneri | YN (KaKs_Calculator) | 294 | 278 | 92.4 % | 1.2e-06 | 0.972 | 0.987 | 16 |

For the close pair the engines agree (median difference 0, r >= 0.97). The remaining differences come from the
protein alignments (the two engines may place a gap differently where several alignments score the same). The
yeast pairs are saturated at synonymous sites (median dS above 2): there COAL returns `NA` (status "saturated")
for most YN estimates for which KaKs_Calculator still prints a number (dS 4 - 5); where both give a value they
agree. DIAMOND `sensitivity_mode`: "fast" (default) finds 5.0 % (S. cerevisiae / C. glabrata), 8.6 %
(S. cerevisiae / K. lactis) and 1.5 % (Salmonella / Shigella) fewer reciprocal best hits than "default", at about
the same run time on these genomes; the default is unchanged, and `dNdS()` suggests "default" when the median
identity of the best hits is below 70 %.

### Internal

- External programs are started through one helper (`run_tool()`); `compute_dnds()` has one code path for serial
  and parallel runs.
- Tests: 16 of the 20 empty test files now hold tests; pinned dN/dS values of fixed codon alignments, the 20
  known ortholog pairs end to end, a tie-break fixture for both best-hit rules, IDs end to end, stored 0.4.3
  tables for the legacy settings, malformed input files (empty, truncated record, lower case, protein instead of
  CDS, duplicated IDs, Windows line endings, gzip, an N in a sequence). The suite passes without any optional
  program.
- GitHub Actions: R CMD check on Linux, macOS and Windows with the external programs installed, one leg without
  them; triggers for pull requests, the default branch and a manual start. Checked with `--as-cran`.
- Examples that need no external program are run by `R CMD check`; the others are guarded by the programs they
  need.
- `.travis.yml`, `orthologr.yml`, `Read-and-delete-me` and the unused `inst/KaKs_Calc_parser` were removed.

## `orthologr` version 0.4.3

### New Function

- `diamond_install()` downloads the pre-compiled DIAMOND2 binary matching the user's operating system (Linux x86_64, macOS, Windows) into the R user cache. `diamond()`, `diamond_best()`, `diamond_rec()`, `set_diamond()`, `deepclust()` and `deepclust_realign()` use it automatically if `path = NULL` and `diamond` is not on the `PATH`. `diamond_path()` reports which DIAMOND executable is used.
- `deepclust()` enables `diamond deepclust` functionalities for clustering protein sequences. Multiple input `.fasta` files or `.fasta.gz` are supported.
- `deepclust_realign()` enables `diamond deepclust` functionalities for realigning protein sequences from the output of `deepclust()`. Multiple input `.fasta` files or `.fasta.gz` are supported.
- `deepclust_annotate()` enables `diamond deepclust` functionalities for annotating protein sequences from the output of `deepclust()` based on user provided names for each proteome file (or the name of the proteome file if not provided). `deepclust_annotate()` can also be used to annotate the output of `deepclust_realign()`. Multiple input `.fasta` files or `.fasta.gz` are supported.

### Miscellaneous

- Unit tests for `deepclust()`, `deepclust_realign()` and `deepclust_annotate()` have been added to ensure the correct functionality of the new features.
- Unit tests for `dNdS()`, `divergence_stratigraphy()` and `divergence_map()` functions have been created.

## `orthologr` version 0.4.2

### New Function

- `diamond()`, `set_diamond()`, `diamond_best()` and `diamond_rec()` enable a massive speed-up in pairwise sequence alignment functionalities.

### New Features

- a logo is now set for `orthologr`.
- DIAMOND2 is now used by default in `dNdS()` and `divergence_stratigraphy()` unless blast is specified (`aligner = "blast"`).

## `orthologr` version 0.4.1

### New Features

- the `divergence_stratigraphy()` and `divergence_map()` functions now include the parameter `n_quantile`, which enables users to choose the number of quantiles to generate for the `divergence map`. This could allow users to get a higher-resolution `divergence map` if `n_quantile` is greater than 10. Alternatively, this can resolve the issue of empty divergence strata when deciling the dNdS values for closely related organisms with dNdS = 0 for over 10% of the genes.

## `orthologr` version 0.4.0

### New Functions

- new function `check_annotation()` helps to detect corrupt GFF or GTF annotation files and removes such outlier lines

### New Features

- the `generate_ortholog_tables()` and `retrieve_longest_isoforms()` now include the new `check_annotation()` function to capture corrupt GFF or GTF files and fix them 
- adding a new argument `of_path` to `orthofinder2()` to allow users to specify their own path
to their locally installed `orthofinder` executable

- adding new argument `task` to `map_generator_lnc()` and `orthologs_lnc()` to allow users
to use the full `blastn` range provided by [blast_nucleotide_to_nucleotide()](https://drostlab.github.io/metablastr/reference/blast_nucleotide_to_nucleotide.html)
- adding new argument `path` to `map_generator_lnc()` to allow users to specify their local installation path of BLAST

## `orthologr` version  0.3.0

### New Functions

- new function `plot_pairwise_orthologs()` allows users to plot pairwise orthologs tables for multiple pairwise comparisons

- new function `retrieve_core_orthologs()` allows users to retrieve a core set of orthologous gene loci from several pairwise ortholog tables

- new functions `generate_ortholog_tables()` and `generate_ortholog_tables_all()` allow users to generate ortholog tables by gene locus and splice varaint for a set of species

- new function `retrieve_longest_isoforms_all()` allows users to specify folders
and retrieve the longest splice variants for all proteomes stored in a folder

- new functions `translate_cds_to_protein()` and `translate_cds_to_protein_all()` which
translate coding sequences into amino acid sequences for single or multiple files

### New Features

- in `orthologs()` the default value of `delete_corrupt_cds` changed from `delete_corrupt_cds = TRUE` to `delete_corrupt_cds = FALSE` to be consistent with
`dNdS()` and `divergence_stratigraphy()`

- the `divergence_stratigraphy()` function received a new argument `dnds_est.method`
which now allows users to select different dNdS estimation methods when running `divergence_stratigraphy()` (suggested by Momir Futo)

- the `divergence_stratigraphy()` function now allows to change the `eval` argument which wasn't passed down to the `dNdS()` call within the function (Many thanks to Momir Futo) 

- the function `map.generator()` was renamed to `map_generator_dnds()` to be more consistent with the notation of other functions

- the function `map.generator.lnc()` was renamed to `map_generator_lnc()` to be more consistent with the notation of other functions

- the function `DivergenceMap()` was renamed to `divergence_map()` to be more consistent with the notation of other functions

- the function `DivergenceMap()` was renamed to `divergence_map()` to be more consistent with the notation of other functions

- the function `orthologs.lnc` was renamed to `orthologs_lnc` to be more consistent with the notation of other functions

- the function `OF2CoreOrthologs()` was renamed to `orthofinder2_retrieve_core_orthologs()`

### Removed Functions

- the function `advanced_blast()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `advanced_makedb()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `blast.nr()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `delta.blast()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `ProteinOrtho()` is not supported anymore and thus is not available to users anymore


## `orthologr` version  0.2.0

### New Features

- new function `retrieve_longest_isoforms()` which enables retrieval of the longest isoforms from a proteome file and save results as fasta file for downstream analyses
- new function `OF2CoreOrthologs()` to retrieve core orthologs across multiple species from Orthofinder2 output
- new function `extract_features()`: Helper function to extract gene loci and splice variant IDs from GFF files
- new function `filter_best_hits()`: Helper function to select best BLAST hit based on minimum evalue
- new function `generate_ortholog_tables()`: Generate ortholog tables by gene locus and splice varaint

## `orthologr` version  0.1.0

### New Features

- `read.cds()` now trimms corrupted CDS (= CDS not divisible by 3) when `delete_corrupt_cds = FALSE` is specified
- the new default value of argument `delete_corrupt_cds` in `dNdS()` is now `FALSE`.
Thus, given the new trimming feature in `read.cds()`, corrupted CDS equences will be trimmed before being translated.


## `orthologr` version  0.0.5

### Bug fixes

- The default setting of the `BLAST` argument `max_target_seqs 1` was removed from `blast_best()` and `blast_rec()` due to
the misunderstood functionality of the `BLAST` argument (See details [here](https://doi.org/10.1093/bioinformatics/bty833) and #9 ; Many thanks to @armish)

## `orthologr` version  0.0.4

### New Features 

- Users can now control via the new `delete_corrupt_cds` argument in `dNdS()` and related downstream functions
whether or not corrupted input coding sequences shall be removed prior to dN/dS inference.
In case corrupted CDS exist, the `dNdS()` function will now generate a separate fasta file which
stores all corrupted CDS so that they can be investigated. See issue #8 for details.

#### Function updates

- `dNdS()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `read.cds()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `cds2aa()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `set_blast()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast_best()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast_rec()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

## `orthologr` version  0.0.3

- Fixing internal path bug that caused that wrong pal2nal paths were generated when using multiple sequence aligners -> see issue https://github.com/drostlab/orthologr/issues/5 (Many thanks to Dr. Mario López-Pérez)

## `orthologr` version  0.0.2

- Fixing a major bug that caused KaKs_Calculator to not be able to correctly 
parse the kaks computation output (Many thanks to Hongyi Li who spotted the bug and found a solution).
