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
- **Column `status`** (last column of every dN/dS table): `"ok"`, `"dS_zero"`, `"saturated"`, `"no_data"` or
  `"failed"` for each gene pair - the same words with both engines and in the `pairs` tables of
  `dNdS_clusters()` / `dNdS_orthogroups()` (see `?dNdS`). `options(orthologr.status_column = FALSE)` leaves it out.
- **One summary message at the end of a run**: gene pairs by status, median dS, sequences left out before the
  search, and a note when most pairs are saturated or the median dS is above 1.5 (species too distant for
  pairwise dN/dS). With `quiet = FALSE`, a progress message about every tenth of the gene pairs.
- **Record of the run**: `attr(result, "orthologr_run")` on results of `dNdS()`, `dNdS_orthogroups()` and
  `divergence_stratigraphy()`: package version, arguments, every `orthologr.*` option in force (with its default
  when it is not set, `orthologr.ambiguous_cds` included), engine, program versions (KaKs_Calculator's from the
  "Version:" line of its `-h`), date, size and MD5 of the input files, counts by status; the call as written
  (`call`) and, when it is of another function, its name and all its arguments (`called_as`); the DIAMOND /
  BLAST+ command lines of the run (`commands`). **`write_dNdS()`** writes a table as TSV with this record in
  `<file>.run.json`; **`read_dNdS()`** reads both back (needs `jsonlite`). `orthogroups()$params` holds the package version.
- **`orthologr_legacy_options()`** and `options(orthologr.legacy = TRUE)`: all settings of 0.4.3 at once. A test
  compares the result with tables stored from 0.4.3 (20-gene example, Comeron and YN).
- **`?orthologr_options`**: every option and environment variable the package reads.
- **Classed errors and warnings**: the common failures have the class `orthologr_error` plus one for their kind -
  `orthologr_error_tool_missing`, `_package_missing`, `_tool_failed` (with `_alignment_failed` for an alignment
  program), `_no_hits`, `_bad_input`, `_bad_argument` - and the warnings about failed gene pairs and isoforms the
  classes `orthologr_warning_failed_pairs` / `_isoforms`, so that they can be handled with `tryCatch()` without
  matching message text (`?orthologr_options`, section "Errors and warnings"). Messages are unchanged. The
  argument and input-file checks of `dNdS()`, `orthogroups()`, `dNdS_orthogroups()`, `divergence_stratigraphy()`,
  `filter_dNdS()`, the alignment functions, `substitutionrate()` and the readers are classed; `dNdS()` checks
  that its input files exist before any program runs.
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
| KaKs_Calculator methods: dN and dS that are 0 are 0 (KaKs_Calculator prints `NA` below 1e-6): identical sequences have dN = dS = 0 and status `"dS_zero"` (was `NA` and `"saturated"`), a pair with nonsynonymous differences only has dS = 0 | E. coli K-12 / O157, YN (3,663 pairs): 230 pairs `"saturated"` -> `"dS_zero"` (the summary no longer counts them as saturated), 12 dN of about 1e-12 -> 0; a file against itself (20 genes): "20 saturated" and the note that the species are probably too distant -> "20 with dS = 0"; 1,000-gene example: 0 pairs | `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| KaKs_Calculator methods: an estimate that KaKs_Calculator replaced with a capped value because the correction for multiple substitutions is not defined is `NA`, status `"saturated"` (was returned as an estimate, status `"ok"`): Ka or Ks = 99, JC69's cap 0.75 log(sites) of YN / MYN, a class of sites counted with 0 substitutions in LWL / MLWL / LPB / MLPB (`?substitutionrate`; where the printed counts are within rounding of 0, the counts are computed again from the codon alignment with KaKs_Calculator's arithmetic, so that the status is what KaKs_Calculator did) | S. cerevisiae / C. glabrata, YN (4,164 pairs): `"ok"` 4,164 -> 605, `"saturated"` 0 -> 3,559 (172 with dS = 99, whose median dN/dS was 0.0028 against 0.081 overall, and 3,387 at JC69's cap, dS 2.56 - 5.99); `filter_dNdS(dnds.threshold = 2)` kept 4,164 rows (the 172 with dS = 99 included) -> 605; E. coli K-12 / O157, YN: 21 of 3,663 pairs (divergent pairs of 35 - 88 % identity, dS 2.8 - 4.8); 1,000-gene example: YN 1 of 995, NG 0; LWL family, S. cerevisiae / C. glabrata: Ks `NA` for 1,925 (LWL, LPB) and 2,100 (MLWL, MLPB) of 4,164 pairs - every pair whose class KaKs_Calculator found undefined, and no other (checked against a KaKs_Calculator that reports its decision, also on S. cerevisiae / K. lactis MLPB and E. coli; with the printed counts alone an intermediate version took 21 pairs the wrong way, among them 2 MLPB Ks of 0.97 and 2.26 that are underestimates) | `options(orthologr.legacy_kaks_saturation = TRUE)` |
| KaKs_Calculator methods: an estimate printed as `nan` (Ka or Ks) is a failed estimate - `NA`, status `"failed"`, the reason with the method in `attr(result, "failed_pairs")` and the failed-pairs warning (was `NaN`, in one pair with a dN, and status `"saturated"`); `substitutionrate()` lists it in `attr(, "kaks_failed")` | S. cerevisiae / C. glabrata, MYN (4,164 pairs): 50 pairs `"saturated"` -> `"failed"` (all of them iterations that did not converge, not saturation: COAL estimates 43 of them, dS 0.12 - 4.1), the dN of one pair 0.0547 -> `NA`; YN, NG and the LWL family on the same pairs, and E. coli K-12 / O157: no `nan` | `options(orthologr.legacy_kaks_saturation = TRUE)` |
| dN/dS tables end with the column `status` | one more column | `options(orthologr.status_column = FALSE)` |
| T-Coffee runs in its default mode (was `-mode accurate`, which needs further programs) | other alignments with `aa_aln_tool = "t_coffee"` | `aa_aln_params = "-mode accurate"` |
| `aa_aln_params` is passed on to the alignment program (was ignored) | only when it is set | do not set it |
| A pairwise protein alignment of more than 5e8 cells (length x length; pwalign needs about 3 bytes per cell, so about 1.5 GB) is not computed: the gene pair is returned with status `"failed"` and the reason (was: aligned whatever the memory; `options(orthologr.max_aln_cells)`, `?pairwise_aln`) | human 1,000 transcripts / all mouse CDS, 2 cores: 1 of 946 pairs (titin, 35,991 x 35,390 residues, 1.27e9 cells) failed instead of dS 0.50; peak memory 4.3 GB -> 0.79 GB, run time 87 s -> 42 s; 1,000-gene example: 0 pairs | `options(orthologr.max_aln_cells = Inf)` |
| `divergence_stratigraphy(quiet = TRUE)`, as `dNdS()` (was `FALSE`) | 1,000-gene example: 10 messages instead of 45; same result | `quiet = FALSE` |
| `orthogroups()` and `dNdS_orthogroups()` trim coding sequences whose length is not a multiple of 3 (`delete_corrupt_cds = FALSE`, as `dNdS()` does) | - (new functions) | `delete_corrupt_cds = TRUE` |

All of these except the T-Coffee mode, `aa_aln_params`, `divergence_stratigraphy(quiet)` and `delete_corrupt_cds`
at once: `orthologr_legacy_options()`.

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
  1,000-gene example, 12 messages before (8 about the 12 trimmed sequences), 2 after; `quiet = FALSE` shows the
  old messages.
- Coding sequences with characters other than A, C, G, T are left out before the search, as before, but no
  longer without a word: they are counted in the summary and listed in `attr(, "dropped_sequences")`.
- Intermediate files are kept in a folder of the call inside `tempdir()` and removed afterwards
  (`options(orthologr.keep_tmp = TRUE)` keeps them until the end of the R session, which removes `tempdir()`;
  `options(orthologr.keep_tmp = "<folder>")` keeps them in that folder, also after `Rscript` or a batch job, and
  the message says which applies); no function changes the working directory or writes into it unasked
  (`read.cds()` wrote `<file>_corrupted_cds_seqs.fasta` there; the summary says when that file is in `tempdir()`).
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
- `dnds_est.method = "MLPB"` (default engine, KaKs_Calculator) works: `substitutionrate()` returned `NULL` for it
  and `dNdS()` ended after the whole search with an error about a missing column (20-gene example: before an
  error after 4 s; after 20 pairs, the values KaKs_Calculator gives for `-m MLPB`). There is one list of the
  KaKs_Calculator methods for the whole package, and a method that `substitutionrate()` does not handle is an
  error (`orthologr_error_bad_argument`) instead of `NULL`.

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
- A protein file given where coding sequences are expected is an error that says so (was: a DIAMOND error about
  memory).
- `orthogroups()` runs the check for several transcripts per gene; the check no longer warns for a bacterial
  genome with a handful of genes that have two coding sequences.
- FASTA files with a UTF-8 byte-order mark were rejected ("could not be read properly"), and files with old-Mac
  line ends (CR only) reached DIAMOND as one sequence and failed after about 10 s with a memory error. Both are now
  read like the plain file (same result as for the file without the mark / with LF line ends; also gzipped).
- `dnds_est.method = "kaks_calc"` in `dNdS()` / `dNdS_orthogroups()` is an error before the search (it ran the
  search and failed every gene pair: `kaks_calc.params` exists only in `substitutionrate()`); `"GY"`, `"GMYN"`,
  `"MS"`, `"MA"` and `"ALL"` are documented in `?dNdS` (MA and MS take about 6.5 s per gene pair); with `"ALL"`
  the summary counts gene pairs x methods (20-gene example: "340 rows (20 gene pairs x 17 methods)", was "340 gene
  pairs").
- `dNdS_orthogroups()` summary: `n_saturated` counted pairs with dS above `dS_saturation`, while the status
  `"saturated"` means that dN or dS could not be estimated (60 yeast groups: 296 pairs with status "saturated",
  `sum(n_saturated)` 30). `n_saturated` now counts status `"saturated"`; the new columns
  `n_dS_zero`, `n_no_data` and `n_failed` count the other statuses, and `n_dS_above` is the old count.
- `read_dNdS()` returns IDs with a double quote as written (`d"q` came back as `d""q`), and non-ASCII text with
  the bytes and the encoding marks it had when it was written (kept in `<file>.run.json`), so that IDs read back
  equal those of `dNdS()` in every locale. (An intermediate version marked all text as UTF-8: in a C / POSIX
  locale a non-ASCII ID then no longer matched the ID of the run - `==`, `match()` and `merge()` failed; and
  `write_dNdS()` wrote such IDs of `attr(, "failed_pairs")` as escapes like `<c3><bc>` in a C locale.)
- An invalid `options(orthologr.max_aln_cells)` (`-1`, `"5e8"` as text, `NA`, several values) is an error
  (`orthologr_error_bad_argument`, naming the value) before the search in `dNdS()` and `dNdS_orthogroups()`; it ran
  the search and then failed every gene pair with the message about the option (20-gene example: 20 of 20 pairs
  "failed"). `options(orthologr.keep_tmp)` of `NA` or of several values is an error as well (was taken for
  `FALSE`), and a folder that cannot be created is `orthologr_error_bad_argument` (was an unclassed error). A run of
  one gene pair is summarised as "1 gene pair".
- A DIAMOND older than 2.1.0 is used with a warning (`orthologr_warning_tool_version`; it was used without a word);
  a program that is there but crashes when asked for its version is `orthologr_error_tool_failed` ("could not be
  started ... exit status 139"), no longer "not installed"; a failed DIAMOND search is `orthologr_error_tool_failed`
  and `comp_cores` above the number of cores `orthologr_warning_comp_cores` (both were unclassed). The column `qcov`
  of DIAMOND hit tables equals `qcovhsp` (DIAMOND has no per-subject coverage; documented in `?diamond`, kept for
  the column layout of BLAST tables and of 0.4.3).

### Speed (results unchanged)

- Default `dNdS()` on the 1,000-gene example, 1 core: 744 s (0.4.3) -> 50 s; with 2 cores 27 s. Half of the last
  step comes from building the full-length alignment rows directly instead of with `pwalign::alignedPattern()`
  (97 s -> 50 s); the rest from BLOSUM62 scoring, reading and translating every file once, vectorised
  translation and best-hit selection, and not writing a second copy of every alignment.
- `select_orthologs()`, `retrieve_longest_isoforms()` and `deepclust()` read annotations and sequences in one
  pass.

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
