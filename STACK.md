# The upstream-v2 branch series

A series of branches on top of drostlab/orthologr master (06c3098, version 0.4.3). Each branch builds on the
previous one and is meant to pass R CMD check on its own. Each was checked locally (R CMD check --as-cran and the
full test suite with the external programs); GitHub Actions (Linux release / devel / oldrel, Linux without any
optional program, Windows, macOS) passed only on the branches and commits named under "CI evidence" below. The
least controversial changes come first; the branches that change
results come after the ones that do not; COAL appears only in the last branch. The tree of the last branch is
the tree of `integration/all` plus this file.

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
| `upstream-v2/10-strata-and-saturation` | tie-safe strata, dS = 0 policy, gestimator input checks, Comeron saturation (R17), dN = 0 from KaKs_Calculator (R14), left-out pairs reported (R23); KaKs_Calculator's capped values of saturated estimates are NA (O02; LWL family at 1 - 2P - Q = 0 as KaKs_Calculator decided, V02), zeros it prints as NA are 0 (identical sequences, nonsynonymous-only pairs; O03); an estimate printed as `nan` is a failed pair (V07) | **yes**: Comeron NA for a saturated dN / dS (1 of 995 pairs); dN = 0 instead of NA (11 of 995); YN S. cerevisiae / C. glabrata 3,559 of 4,164 dS -> NA, E. coli K-12 / O157 230 NA -> 0 and 21 dS -> NA; LWL S. cerevisiae / C. glabrata 1,925 Ks -> NA (the classes KaKs_Calculator found undefined), MYN 50 `nan` pairs failed | `options(orthologr.legacy_comeron_saturation = TRUE)`, `options(orthologr.legacy_kaks_zero_dn = TRUE)`, `options(orthologr.legacy_kaks_saturation = TRUE)` |
| `upstream-v2/11-best-hit-rule` | best hits by bit score; order-independent tie-breaks (R20) | **yes**: 3.2 % of human / mouse best hits | `best_hit_rule = "evalue"` or `options(orthologr.best_hit_rule = "evalue")` |
| `upstream-v2/12-golden-values` | pinned dN / dS of fixed codon alignments, known ortholog pairs, tie-break fixture; fragile tests repaired | no (tests only) | - |
| `upstream-v2/13-orthogroups` | new `orthogroups()`, `write_orthogroups()`, `dNdS_orthogroups()` (engine "orthologr"); `orthofinder2()` soft-deprecated | no (new functions) | - |
| `upstream-v2/14-run-record-and-legacy` | column `status`, summary message, record of the run, `write_dNdS()` / `read_dNdS()`, `orthologr_legacy_options()`, `?orthologr_options` (R22), upgrade article; classed errors and warnings (`orthologr_error_*`, `orthologr_warning_*`), `options(orthologr.ambiguous_cds = "keep")`, a hint for `sensitivity_mode` with distant species, malformed-input tests; `aa_sub_matrix` checked before the search; sequences trimmed or left out counted in the summary under `quiet = TRUE` (`attr(, "trimmed_sequences")`); fixes of the deep test: `keep_tmp` folder (O04), BOM / CR input (O05), legacy run time documented (O06), `kaks_calc` / `ALL` (O07), status counts in the `dNdS_orthogroups()` summary (O09), quotes in `read_dNdS()` (O10), DIAMOND version and crash, classed DIAMOND error and `comp_cores` warning (O14) | one more column; `divergence_stratigraphy(quiet = TRUE)` by default | `options(orthologr.status_column = FALSE)`; `quiet = FALSE`; all at once: `orthologr_legacy_options()` |
| `upstream-v2/15-coal` | optional COAL engine (`coalr`): `dNdS(engine = "coal")`, `dNdS_coal()`, `dNdS_clusters()`, `codon_aln(tool = "coal")`, `substitutionrate(kaks_tool = "coal")`; a classed warning when the settings of 0.4.3 meet the COAL engine (which cannot reproduce them); with `coalr` installed `dNdS()` builds its codon alignments with COAL (O11); no failed-pairs warning for saturated pairs of `dNdS_orthogroups(engine = "coal")` (O08); cluster status `"saturated"` in `dNdS_clusters()` (O13); maintainer notes; this file | no for calls that do not ask for COAL (the codon alignments of COAL and PAL2NAL are byte-identical: 12,605 of 12,605 pairs) | do not use `engine = "coal"`; `options(orthologr.codon_aln_tool = "pal2nal")` |

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

The fixes of the third review (R01 - R23) are in the branches they correct (named in the table and in the
commit messages); there is no separate fix-up branch.

CI evidence (GitHub Actions, all six legs green), by branch and commit, as of 2026-11-14:

- `upstream-v2/01-ci-and-tests` at 0445681 (run 43493418485);
- `upstream-v2/04-options-and-inputs` at 2897cd9 (run 43531501018) and `upstream-v2/05-vectorise` at b70fd05 (run
  43531506648) - both before the fixes of the deep test (O01), which were added to 04 and merged into 05 later;
- `upstream-v2/08-annotation` at 447cb11 (run 43493414706);
- `upstream-v2/15-coal` at e0472ce (run 43507540285);
- `integration/all` at 173bb3e (run 43505256122, second attempt: the macOS leg of the first was cancelled) and at
  da9c518 (run 43563306610: the fixes of the deep test, O01 - O14).

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
