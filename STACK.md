# The upstream-v2 branch series

A series of branches on top of drostlab/orthologr master (06c3098, version 0.4.3). Each branch builds on the
previous one and is meant to pass R CMD check on its own. Each was checked locally (R CMD check --as-cran and the
full test suite with the external programs); GitHub Actions (Linux release / devel / oldrel, Linux without any
optional program, Windows, macOS) passed only on the branches and commits named under "CI evidence" below. The
least controversial changes come first; the branches that change
results come after the ones that do not; COAL appears first in branch 15. Branches 16 - 19 (version 0.6.0 in
development: tests, performance, fixes) were appended after 15 when `integration/all` had moved on; the existing
branches 01 - 15 (and later 01 - 18) were not changed. The tree of `upstream-v2/19-tool-batching-and-fixes` is the
tree of `integration/all` (736fc5c) plus this file.

`tests/testthat/test-legacy.R` (from branch 01 on) compares `dNdS()` on the 20-gene example with tables stored
from 0.4.3 (Comeron and YN). Every branch keeps it passing: a branch that changes a default adds the option that
restores it to that test, so the test shows that the old numbers can still be had.

| Branch | What it changes | Results change? | Old behaviour |
|---|---|---|---|
| `upstream-v2/01-ci-and-tests` | GitHub Actions (with and without the external programs, Linux / macOS / Windows); stored 0.4.3 tables and `test-legacy.R`; no package code | no | - |
| `upstream-v2/02-robustness` | a failing gene pair is an `NA` row with a reason instead of ending the run; per-call working folder; same-named input files; hyphenated IDs; Li's method numeric | only fixes (same-named files, hyphenated IDs, Li's columns) | - |
| `upstream-v2/03-tool-calls` | every program through `run_tool()` (exit status and output checked, no `setwd()`); parameter strings are argument lists, not shell text (R07) | no (`;`, pipes in parameter strings no longer work) | - |
| `upstream-v2/04-options-and-inputs` | arguments checked up front, missing programs found before the search (R09), documented options work, MUSCLE 5 / T-Coffee, input fixes; `dnds_est.method = "MLPB"` works with KaKs_Calculator (one list of its methods; O01) | T-Coffee default mode; `aa_aln_params` now reaches the program; fix: T-Coffee gets the sequences with `-seq` (with `-infile` it returned rows ending in `*` unaligned - 10 of 20 pairs saturated on the 20-gene example) | `aa_aln_params = "-mode accurate"`; do not set `aa_aln_params` |
| `upstream-v2/05-vectorise` | each file read and translated once, vectorised best hits / isoforms, `parallel` instead of foreach / doParallel; PSOCK workers on Windows (R08, R13; working folder handed over) | no | - |
| `upstream-v2/06-deps-and-docs` | version 0.5.0, optional packages to Suggests, README / vignettes / man pages corrected, `environment.yml`, pkgdown, native routines registered by Rcpp | no | - |
| `upstream-v2/07-ids-and-joins` | IDs returned as in the files, results joined by gene pair, duplicate IDs an error, dropped sequences reported, isoform warning, same ID in both files (R06) | only fixes (lost or changed IDs) | `options(orthologr.isoform_warning = FALSE)` silences the new warning |
| `upstream-v2/08-annotation` | longest isoforms and gene loci on real NCBI / Ensembl GFF3 / GTF; long annotation lines kept (R03) | only fixes | - |
| `upstream-v2/09-pairwise-alignment` | pairwise protein alignments: BLOSUM62, full length, no terminal stop; full-length rows built directly (about twice as fast); no alignment above 5e8 cells (about 1.5 GB; O12) | **yes**: dN or dS of about a third of the pairs of the 1,000-gene example, mostly slightly (r = 0.988); human / mouse: 1 of 946 pairs (titin) failed instead of aligned in 4.3 GB | `options(orthologr.legacy_pairwise_aln = TRUE)`, `options(orthologr.max_aln_cells = Inf)` |
| `upstream-v2/10-strata-and-saturation` | tie-safe strata, dS = 0 policy, gestimator input checks, Comeron saturation (R17), dN = 0 from KaKs_Calculator (R14), left-out pairs reported (R23); KaKs_Calculator's capped values of saturated estimates are NA (O02; LWL family at 1 - 2P - Q = 0 as KaKs_Calculator decided, V02), zeros it prints as NA are 0 (identical sequences, nonsynonymous-only pairs; O03); an estimate printed as `nan` is a failed pair (V07) | **yes**: Comeron NA for a saturated dN / dS (1 of 995 pairs); Comeron's two-position pathways instead of libsequence's branches (1,000-gene example: dS of 474 of 995 pairs, median 0 %, p95 2.9 %); a Comeron class with fewer than 3 sites cannot make dN / dS NA; dN = 0 instead of NA (11 of 995); YN S. cerevisiae / C. glabrata 3,559 of 4,164 dS -> NA, E. coli K-12 / O157 230 NA -> 0 and 21 dS -> NA; LWL S. cerevisiae / C. glabrata 1,925 Ks -> NA (the classes KaKs_Calculator found undefined), MYN 50 `nan` pairs failed | `options(orthologr.legacy_comeron_saturation = TRUE)`, `options(orthologr.comeron_twosubs = "libsequence")`, `options(orthologr.comeron_min_sites = 0)`, `options(orthologr.legacy_kaks_zero_dn = TRUE)`, `options(orthologr.legacy_kaks_saturation = TRUE)` |
| `upstream-v2/11-best-hit-rule` | best hits by bit score; order-independent tie-breaks (R20) | **yes**: 3.2 % of human / mouse best hits | `best_hit_rule = "evalue"` or `options(orthologr.best_hit_rule = "evalue")` |
| `upstream-v2/12-golden-values` | pinned dN / dS of fixed codon alignments (Comeron for both two-position pathway orders), known ortholog pairs, tie-break fixture; fragile tests repaired | no (tests only) | - |
| `upstream-v2/13-orthogroups` | new `orthogroups()`, `write_orthogroups()`, `dNdS_orthogroups()` (engine "orthologr"; quartile columns `_q25` / `_q75`); `orthofinder2()` soft-deprecated | no (new functions) | - |
| `upstream-v2/14-run-record-and-legacy` | column `status`, summary message, record of the run, `write_dNdS()` / `read_dNdS()`, `orthologr_legacy_options()`, `?orthologr_options` (R22), upgrade article; classed errors and warnings (`orthologr_error_*`, `orthologr_warning_*`), `options(orthologr.ambiguous_cds = "keep")`, a hint for `sensitivity_mode` with distant species, malformed-input tests; `aa_sub_matrix` checked before the search; sequences trimmed or left out counted in the summary under `quiet = TRUE` (`attr(, "trimmed_sequences")`); fixes of the deep test: `keep_tmp` folder (O04), BOM / CR input (O05), legacy run time documented (O06), `kaks_calc` / `ALL` (O07), status counts in the `dNdS_orthogroups()` summary (O09), quotes in `read_dNdS()` (O10), DIAMOND version and crash, classed DIAMOND error and `comp_cores` warning (O14) | one more column; `divergence_stratigraphy(quiet = TRUE)` by default (progress in interactive sessions); **yes**: DIAMOND `sensitivity_mode = "default"` (was `"fast"`; yeasts +5.3 to +9.4 % reciprocal best hits, human / mouse +0.9 %, same precision), coding sequences with ambiguous bases kept (`ambiguous_cds = "keep"`; 1,000-gene example 2 more pairs) | `options(orthologr.status_column = FALSE)`; `quiet = FALSE`; `options(orthologr.sensitivity_mode = "fast")`; `options(orthologr.ambiguous_cds = "drop")`; all at once: `orthologr_legacy_options()` |
| `upstream-v2/15-coal` | optional COAL engine (`coalr`): `dNdS(engine = "coal")`, `dNdS_coal()`, `dNdS_clusters()`, `codon_aln(tool = "coal")`, `substitutionrate(kaks_tool = "coal")`; a classed warning when the settings of 0.4.3 meet the COAL engine (which cannot reproduce them); with `coalr` installed `dNdS()` builds its codon alignments with COAL (O11); no failed-pairs warning for saturated pairs of `dNdS_orthogroups(engine = "coal")` (O08); cluster status `"saturated"` in `dNdS_clusters()` (O13); the COAL engine with Comeron's two-position pathways and minimum class size (`coalr` with `comeron_twosubs`, `comeron_min_sites`); with `coalr` installed, pairs above `max_aln_cells` aligned by COAL's pwalign preset (the alignment pwalign would return; titin identical); selenocysteine `U` identical with PAL2NAL and COAL; COAL's codon-alignment errors classed; maintainer notes; this file | no for calls that do not ask for COAL (the codon alignments of COAL and PAL2NAL are byte-identical: 12,605 of 12,605 pairs) | do not use `engine = "coal"`; `options(orthologr.codon_aln_tool = "pal2nal")`; `options(orthologr.large_pair_aligner = "none")` |
| `upstream-v2/16-coverage-tests` | tests only: the COAL glue (`R/coal.R`, with a fake `coal-cpp` script), `map_generator_dnds()` / `map.generator.lnc()`, `orthologs_lnc()`, Kimura80, `orthologs()` dispatch, `orthofinder2()`'s run path (fake OrthoFinder), `gestimator()` output modes, readers' error paths, small helpers; two tests of known bugs (B1, B2) are skipped here and run from 18 on | no (tests only) | - |
| `upstream-v2/17-performance` | `dNdS()` about 2 - 2.5 times faster (whole call, 1,000-gene example: 29.2 -> 12.1 CPU s; S. cerevisiae / C. glabrata 106.6 -> 53.9): protein alignments of up to 100 pairs in one `pwalign::pairwiseAlignment()` call, COAL codon alignments from the rows in memory, `gestimator()`'s redundancy and Grantham tables built once per process, no `data.table` sub-assignment per pair, sequence checks with PCRE; `stringr` (and `stringi`) no longer imported; `test-perf.R` (batched against per-pair code); NEWS 0.6.0 with the measured times | no: every table identical to 0.5.0 (all columns, attributes, run record; 20-gene, 1,000-gene, S. cerevisiae / C. glabrata; default, `engine = "coal"`, PAL2NAL, Comeron / Li / YN, 1 and 2 cores, legacy switch) | - |
| `upstream-v2/18-fixes` | Comeron's min-sites rule (D9 of 10) only where the value that results is > 0 (D9 refinement; COAL #76); `dNdS_clusters()` with a named `cds_file` vector (B1) and `codon_aln(tool = "coal")` on CLUSTAL files without a conservation line (B2), both bugs of 15; `dNdS_orthogroups()` summary `_min` / `_max` (K04); errors for an unknown `dnds_est.method` list each engine's methods and suggest the likely one (U09); a note on `library(orthologr)` that results differ from 0.4.3, once per session (U06); README links to the repository's articles (U14); docs: titin timing of the large-pair route (K02), D8's and D9's measured effects (K03) | **yes** (D9 refinement): a Comeron dN or dS that the min-sites rule of 10 left at exactly 0 although a change was seen is `NA` (status `"saturated"`) instead of 0 (`"dS_zero"`); no value > 0 changes; real data (20-gene, 1,000-gene, S. cerevisiae / C. glabrata): no value changes; seeded short pairs of 5 - 20 codons: 2 of 160 dS 0 -> NA. Otherwise new summary columns, new messages | `options(orthologr.comeron_min_sites = 0)` (NA for every undefined small class); 0.4.3's numbers: `orthologr_legacy_options()` (unaffected); `options(orthologr.quiet_upgrade = TRUE)` hides the note |
| `upstream-v2/19-tool-batching-and-fixes` | PAL2NAL and KaKs_Calculator with one process per chunk of up to 100 gene pairs instead of one per pair (`inst/pal2nal/pal2nal_batch.pl` compiles `pal2nal.pl` once and forks per pair; one AXT file per chunk for KaKs_Calculator; a pair the batch does not compute is computed on its own; not under the legacy switch; `options(orthologr.batch_tools)`); KaKs_Calculator output read with its known header (KaKs_Calculator 2.0 may leave out or repeat the header line of the gamma methods), so GMYN is batched too; rows and fields of a chunk checked per record (AUDIT-OLP2 M1, L1); fixes of the review REVIEW-OL-1: a warning for a `coalr` without COAL #76 by its feature flag (H4), Comeron's two flags in the recorded `coal pipeline` command (M1), `kaks_tool = "coal"` chooses the program of the seven KaKs_Calculator methods only and the unreachable COAL branch for Comeron / Li is removed (M2), COAL's engine version in the run record (M3), `engine = "coal"` with several methods and `kaks_tool` (L1), `comp_cores` a whole number >= 1 (L2), coding sequences with in-frame stop codons reported (L6), unused internal functions removed and `benchmarkme` out of Suggests (L4), a test that skipped on any error (L5), docs (L3, L7, L8, M4, M5), Chloe Newton in Authors@R (L10); CI: no `coalr` on the leg without tools, `coalr` pinned to COAL commit 03612da, an R 4.4 leg (M6); faster input-type check, gestimator reader and no unread per-pair files with COAL in memory (P1 - P3) | no for the values of a table except GMYN: PAL2NAL / KaKs_Calculator per chunk give the tables of one process per pair (all columns, attributes, `failed_pairs`, run record; every KaKs_Calculator method, 1 and 2 cores, 20-gene, 1,000-gene and yeast data; NEWS 0.6.0, Performance), and the review fixes gave identical tables on 40 of 40 runs. **Yes** (from NEWS 0.6.0, Fixed / Other changes): GMYN (and `"ALL"`) pairs whose header line KaKs_Calculator left out now have values instead of failing (f64de0f; output with one header line first reads as before); run record of the COAL engine: `engine_settings$coal_command` has `--comeron-twosubs` / `--comeron-min-sites` and `engine_settings` the fields `comeron_twosubs`, `comeron_min_sites`, so the recorded command reproduces the run's Comeron values (8e46d2a, M1), `tools$coal` recorded (ae448f0, M3), COAL neither required nor recorded for Comeron / Li with `kaks_tool = "coal"` (a93711e, M2); Comeron with an old `coalr`: the `orthologr_warning_coal_version` warning where it was silent (f67fb5c, H4; values unchanged); new attribute `stop_codon_sequences` and summary line (a320e69, L6); `comp_cores = 1.5` is an error (b51353e, L2) | `options(orthologr.batch_tools = FALSE)` (one process per pair, the same tables); the GMYN reading and the run-record fields are fixes without a switch |

The fixes of the deep test of 2026-11-14 (O01 - O14) were added in the same way: O01 in 04 (with the test of every
method in 06), O12 in 09, O02 / O03 in 10, the parts that use the classed conditions, the run summary and the
options page (O04 - O07, O09, O10, O14) in 14, the COAL parts (O08, O11, O13) in 15.

The fixes of the re-test of 2026-11-14 (V01 - V07) likewise: V04 (`options(orthologr.max_aln_cells)` checked before
the search) in 09 and, for `dNdS_orthogroups()`, in 13; V02 (LWL-family saturation decided by KaKs_Calculator's own
arithmetic where its printed counts are within rounding of 0) and V07 (an estimate printed as `nan` is a failed
pair, not saturated) in 10, the status parts of V07 in 14; V03 (`read_dNdS()` keeps the bytes and encoding marks of
non-ASCII IDs), the `keep_tmp` checks of V06 and the documentation of non-converged YN / MYN iterations (V01) in 14;
the COAL parts (V01: `kaks_tool = "coal"` reports them; V05: `codon_aln_tool = "coal"` without COAL uses PAL2NAL;
V06: the `max_aln_cells` reason mentions `engine = "coal"`) in 15. Results change only in 10 (LWL family: pairs at
1 - 2P - Q = 0; `nan` estimates), both under `options(orthologr.legacy_kaks_saturation = TRUE)` as before.

The fixes of the integration pass with COAL of 2026-11-14 (C04, C05, C07) likewise: the run record's
KaKs_Calculator version (C04), the user's call, `called_as` and the DIAMOND / BLAST+ `commands` (C05, C07) in 14;
the COAL parts (`engine_settings` with the equivalent `coal pipeline` command line, pwalign's BLOSUM62 for
`engine = "coal"`, the documentation of both directions in `?dNdS_coal`, and the agreement of the two engines
measured with COAL's pwalign tie-break, in NEWS) in 15. No result changes for calls that do not ask for COAL.

The decisions of 2026-11-14 (DECIDE, items 1, 3, 4, 7, 8, 9) likewise: Comeron's two-position pathways and minimum
class size (D8, D9) in 10, with the pinned values of both orders in 12; the quartile names of `dNdS_orthogroups()`
(D7b) in 13; `sensitivity_mode = "default"` (D1), `ambiguous_cds = "keep"` (D7a), the progress line of
`divergence_stratigraphy()` (D7c), the new options in the legacy switch and the "Upgrading from 0.4.3" block of NEWS
in 14; the COAL parts (the engine's Comeron settings, large pairs aligned by COAL (D3), selenocysteine (D4), the
classed COAL error) in 15. The legacy switch still reproduces the stored tables of 0.4.3 (test-legacy.R).

The fixes of the third review (R01 - R23) are in the branches they correct (named in the table and in the
commit messages); there is no separate fix-up branch.

Branches 16 - 19 are different: they were appended after 15 instead of being folded into the branches they touch,
because folding them in would have meant merging forward through every later branch for changes that leave no
earlier branch broken (each of 01 - 15 still passes its own check). Where they correct an earlier branch:

- 18's D9 refinement changes the min-sites rule of 10 (and the COAL engine's copy in 15). A maintainer who takes
  10 without 18 gets dS = 0 (status `"dS_zero"`) for the rare pair whose only synonymous change falls in a class
  of fewer than 3 sites; no real data set measured has such a pair.
- 18's B1 and B2 fix `dNdS_clusters()` (a named vector of CDS files found no CDS) and `codon_aln(tool = "coal")`
  (CLUSTAL files without a conservation line read as empty), both added in 15 and both only with `coalr`. 16 adds
  their tests skipped, 18 makes them run.
- 17 changes no result. Its first commit is the net effect of two commits of `integration/all` (5a35092, which
  removed unused internal functions, and 329e65c, which restored them because 16's tests call them): only the
  `stringr` calls are replaced. Each commit of 16 - 18 names the commit of `integration/all` it was taken from.
- 19 corrects 15 (and 18) in places: the documentation of `kaks_tool = "coal"` and the COAL requirement for
  Comeron / Li (M2), the recorded COAL version (M3) and `coal pipeline` command (M1, of C05), the `kaks_tool` check
  with `engine = "coal"` (L1), and the `coalr` version check (H4, which 18's D9 refinement needs). A maintainer who
  takes 15 without 19 gets correct values but an incomplete run record and documentation that claims COAL computes
  Comeron / Li with `kaks_tool = "coal"`.
- In 19, "KaKs_Calculator per chunk: not for GMYN" (5b36cff) takes GMYN out of the batch and
  f64de0f puts it back once the output is read with its known header; both are kept, because the first message
  explains the KaKs_Calculator 2.0 behaviour that the second handles.
- The commits of 19 are those of `integration/all` from 62545f6 to 90d2863 (the merges of `perf/0.6-kc` and
  `fix/review-ol-1`), in order and with their authorship. Five are one commit of `integration/all` together with a
  later follow-up of it (the NEWS placeholder filled by 2d969ff; test fixes 816420a, 2d5d2db; the CI fix 90d2863);
  their messages name both. The other commits are unchanged.
- NEWS has a 0.6.0 section from 17 on; the `Version` field stays 0.5.0 (as on `integration/all`), to be raised at
  release.

CI evidence (GitHub Actions, all six legs green; seven from 19 on), by branch and commit, as of 2026-11-15:

- `upstream-v2/01-ci-and-tests` at 0445681 (run 43493418485);
- `upstream-v2/04-options-and-inputs` at 2897cd9 (run 43531501018) and `upstream-v2/05-vectorise` at b70fd05 (run
  43531506648) - both before the fixes of the deep test (O01), which were added to 04 and merged into 05 later;
- `upstream-v2/08-annotation` at 447cb11 (run 43493414706);
- `upstream-v2/15-coal` at e0472ce (run 43507540285);
- `upstream-v2/15-coal` at 0bd0b00, its current tip (run 43616613722);
- `integration/all` at 173bb3e (run 43505256122, second attempt: the macOS leg of the first was cancelled) and at
  da9c518 (run 43563306610: the fixes of the deep test, O01 - O14);
- the trees of the tips of 16, 17 and 18 without this file: `integration/all` at 9dbb233 (run 43632907574; tip of
  16), at 4e56d32 (run 43641948334; tip of 17) and `fix/d9-ol` at b3e50fa (runs 43649614517, 43649621781; tip of
  18, the tree that `integration/all` 071f933 has). The run of the tip of 18 itself is on the Actions page; the
  runs of the tips of 16 and 17 were cancelled to save runner time (their trees passed as listed). Locally, with
  the external programs and `coalr`, the full test suite at each tip: 16 571 tests (3,890 expectations), 17 581
  (4,031), 18 585 (4,084), 0 failed (skipped: no `metablastr`, no TAIR files, root user, `coalr` installed, and in
  16 the two tests of B1 / B2); R CMD check --as-cran --no-manual on 18: 0 ERROR, 0 WARNING, 6 NOTEs (`Remotes`,
  the time server, `-march`, an `\examples` line of `?diamond_install` over 100 characters, the `metablastr`
  cross-reference, a 6 s example; none new since 15).
- the tree of the tip of 19 without this file: `fix/review-ol-1` at 90d2863 (run 43711434023, all seven legs: the
  six above and R 4.4, the oldest R that `DESCRIPTION` allows), whose tree `integration/all` 736fc5c has; before the
  review fixes, `integration/all` at 6a764a6 (run 43698427951; the merge of `perf/0.6-kc`). The runs of the tips of
  19 and 20 themselves are on the Actions page (19's was cancelled to save runner time; 20 differs from 19 only in
  install commands, CI triggers and notes).

The run of `upstream-v2/15-coal` at aa43248 (43563309087; the same tree as `integration/all` da9c518) failed on
its Ubuntu release leg in `test-substitutionrate.R`: KaKs_Calculator wrote no result line for one of its methods,
as on two legs of the runs of fdb3730 / b9f978f before (the run of da9c518 passed). Locally every method gives a
result. The test now requires values from the approximate methods (NG, LWL, MLWL, LPB, MLPB, YN, MYN) and accepts
an error naming the method for the model-based ones (GY, GMYN, MA, MS, ALL); the error names the method, so a
further failure shows which one it is.

No other branch has a completed run: the runs of 02, 03, 06, 07 and 09 - 14 were cancelled, because the fork's runners
start about one job at a time and a push of all branches queues more than 15 runs. Their evidence is local:
R CMD check --as-cran --no-manual with DIAMOND, BLAST+, Clustal Omega / ClustalW, MAFFT, MUSCLE, T-Coffee,
KaKs_Calculator and coalr, 0 ERROR and 0 WARNING (NOTEs from the environment) - in the independent review of
2026-11-13 on 02, 05, 10 and 13, and in earlier passes on every branch. The fixes of that review (T01 - T07,
T14, T19) were checked locally on `integration/all`, `upstream-v2/04-options-and-inputs`,
`upstream-v2/14-run-record-and-legacy` and `upstream-v2/15-coal`; the CI runs of the commits that hold them are
listed on the Actions page, not here. The fixes of the deep test of 2026-11-14 (O01 - O14) were checked with
R CMD check --as-cran on `integration/all` (0 ERROR, 0 WARNING, 4 NOTEs from the environment; tests 3,309 passed,
0 failed) and on `upstream-v2/14-run-record-and-legacy` (tests 2,896 passed, 0 failed; its one WARNING, a test
that asked for `coalr` before 15 lists it in Suggests, was fixed in 10 and merged forward, and the static checks
were run again), and with the test files they touch on 04, 06, 09 and 10.
Fixes found later were added as commits to the branch they belong to and
merged forward through the later branches. To get a green run on a branch without one, start the workflow by
hand for that branch (`workflow_dispatch`), one branch at a time.

**Proposal: splitting branch 14 later.** `upstream-v2/14-run-record-and-legacy` holds several independent
changes; the branches are pushed and are not rewritten. A maintainer-facing series could split it without
touching the existing branches, by building new branches from `upstream-v2/13-orthogroups` and taking the
changes over file by file (`git checkout upstream-v2/14-run-record-and-legacy -- <files>` plus hand edits where
a file mixes them), each with its tests, in this order:

1. `14a-status-and-summary`: column `status`, `dnds_status()`, the summary message, `quiet` (incl. the counting of
   trimmed / left-out sequences), `options(orthologr.status_column)`; `R/run_summary.R` (status and summary
   parts), `R/dNdS.R`, `R/utils.R`, `R/io.R`, `test-run_summary.R` (first half).
2. `14b-run-record`: `attr(, "orthologr_run")`, `write_dNdS()` / `read_dNdS()` (jsonlite in Suggests).
3. `14c-legacy-switch`: `orthologr_legacy_options()`, `options(orthologr.legacy)`, `?orthologr_options`, the
   upgrade article (needs 14a for the status column option).
4. `14d-classed-conditions`: `R/conditions.R`, the `orthologr_abort()` calls, `test-conditions.R`,
   malformed-input tests; largely independent of 14a - 14c (the calls inside their code move with them).
5. `14e-inputs`: `options(orthologr.ambiguous_cds = "keep")`, the `sensitivity_mode` hint, the `aa_sub_matrix`
   check.

The tree of the last new branch must equal that of `upstream-v2/14-run-record-and-legacy` (checked with
`git diff --stat`), so that 15 can be merged on top unchanged.

Differences from `integration/all`: this file, and its line in `.Rbuildignore`.
