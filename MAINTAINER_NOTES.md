# For the maintainer

Notes that need a decision by the maintainer of `drostlab/orthologr`. This file is not part of the built package
(`.Rbuildignore`).

## 1. Licence statements are inconsistent (finding R05 of the third review) - nothing was changed

| Where | What it says |
|---|---|
| `DESCRIPTION:12` | `License: GPL-2` (version 2 **only**) |
| `LICENSE`, `COPYING` | text of the GNU GPL version 2 |
| `src/Comeron95.h:11`, `src/Comparisons.h:11`, `src/Grantham.h:11`, `src/GranthamWeights.h:11`, `src/Kimura80.h:11`, `src/RedundancyCom95.h:11`, `src/Sequence.h:11`, `src/SingleSub.h:11`, `src/Sites.h:11`, `src/TwoSubs.h:11` (also `src/ThreeSubs.h`, `src/gestimator.cpp`) | libsequence: "GNU General Public License [...] either **version 3** of the License, or (at your option) any later version" |
| the same files, line 27 | the orthologr changes to them: "version 2 of the License, or any later version" |
| `inst/pal2nal/pal2nal.v14/` | PAL2NAL v14 (GPL; see its README) |
| `README.md`, section "Licenses" | quotes the libsequence header |

Code under "GPL version 3 or later" cannot be distributed as part of a work that is licensed "GPL-2 only": the
package as a whole can only be GPL (>= 3). The `License` field was **not** changed on this fork, because that is
the copyright holders' decision.

Further fact: the optional package **coalr** (COAL's R package, used by `dNdS(engine = "coal")`,
`dNdS_coal()`, `dNdS_clusters()`, `dNdS_orthogroups(engine = "coal")` on the branch with the COAL integration)
is **GPL-3** (only). It is in `Suggests` and loaded only when installed; orthologr does not contain its code.
A GPL-2-only orthologr that calls a GPL-3 package at run time is a second reason to move to version 3.

**Proposal**

- `DESCRIPTION`: `License: GPL (>= 3)`
- replace the text of `LICENSE` / `COPYING` by the GPL version 3
- add `inst/COPYRIGHTS` with, in substance:

  ```
  orthologr: Copyright (C) Hajk-Georg Drost and contributors. GPL (>= 3).

  src/gestimator.cpp and the headers Comeron95.h, Comparisons.h, Grantham.h, GranthamWeights.h, Kimura80.h,
  RedundancyCom95.h, Sequence.h, SingleSub.h, Sites.h, ThreeSubs.h, TwoSubs.h are taken from libsequence
  (Copyright (C) Kevin Thornton), GNU GPL version 3 or (at your option) any later version; modified by
  Sarah Scharfenberg and Hajk-Georg Drost (2014) and by the orthologr contributors (2026: memory, input checks,
  saturation flags), modifications under GPL version 2 or any later version.

  inst/pal2nal/pal2nal.v14: PAL2NAL v14, Copyright (C) Mikita Suyama, GNU GPL (see inst/pal2nal/pal2nal.v14/README).
  ```

## 2. Defaults and results that differ from orthologr 0.4.3, and the option that restores each

| Change in 0.5.0 | Restored by |
|---|---|
| Best hits are ranked by bit score (was: e-value, then alignment length) | `options(orthologr.best_hit_rule = "evalue")` or the argument `best_hit_rule = "evalue"` |
| Pairwise protein alignments (BLOSUM62; only on the branch with the faster pairwise alignment) | `options(orthologr.legacy_pairwise_aln = TRUE)` |
| Comeron's method: `NA` instead of a finite value when a correction is undefined in a class that has sites | `options(orthologr.legacy_comeron_saturation = TRUE)` |
| KaKs_Calculator methods: `dN = 0` instead of `dN = NA` for a pair without nonsynonymous change | `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| Li's method returns numeric columns; a value below 9.9 that R compared as text is no longer `NA` | none (bug fix) |
| dN/dS tables end with the column `status` (`"ok"`, `"dS_zero"`, `"saturated"`, `"no_data"`, `"failed"`) | `options(orthologr.status_column = FALSE)` |
| T-Coffee default mode (was `-mode accurate`) | `aa_aln_params = "-mode accurate"` / `params = "-mode accurate"` |
| `aa_aln_params` is passed on to the alignment program (was ignored) | do not set it |
| Parameter strings are argument lists, not shell text (no `;`, pipes, redirections) | none (security fix) |
| A missing program (perl, aligner, KaKs_Calculator) is an error before the search (was: all pairs `NA` + warning on this fork; an error at the first pair in 0.4.3) | none |
| Annotation lines of more than 65,000 characters are read (were deleted) | none (bug fix); `check_annotation(remove_annotation_outliers = TRUE)` still writes a copy without them |

**Not changed any more** (a default introduced earlier on this fork was taken back after the third review):
`filter_dNdS()` and `divergence_stratigraphy()` leave gene pairs with dS = 0 out, as 0.4.3 did
(`dS_zero = "drop"`); new are a message with the numbers, `attr(, "left_out")` / `divergence_left_out()`, and
`dS_zero = "keep"` / `options(orthologr.dS_zero = "keep")` to keep them.

**One switch:** `orthologr_legacy_options()` sets the five options of this table (and `orthologr.dS_zero =
"drop"`); `options(orthologr.legacy = TRUE)` is the shorthand. With them the 20-gene and 1,000-gene examples give
the tables of 0.4.3 (master 06c3098) for Comeron and YN. `tests/testthat/test-legacy.R` compares the 20-gene
result with tables stored from 0.4.3 (`tests/testthat/fixtures/dnds_0.4.3/`, made by `make_fixture.R` there with
the 0.4.3 package installed) on every run, for both ways of setting the switch.

## 3. Things that are set for the fork and must change when this is merged

| Where | Now | Change to |
|---|---|---|
| `DESCRIPTION`, `Remotes:` | `coalr=Plumat-COAL/COAL/packaging/r/coalr` - COAL's default branch, which holds the R package coalr since COAL's pull request #54 was merged (syntax `package=owner/repo/subdir[@ref]`; the `coalr=` prefix is needed because the package is not named like its repository - without it pak, which GitHub Actions uses, reports "Can't find package called coalr") | a release tag (`...coalr@v<version>`) once COAL has one with the Comeron / Li methods |
| install command of orthologr in `README.md` (two places), `environment.yml`, `vignettes/orthologr.Rmd` | `remotes::install_github("Plumat-COAL/orthologr@integration/all")`, so that the command installs the version the page describes | `remotes::install_github("drostlab/orthologr")` |
| `.github/workflows/R-CMD-check.yaml`, `on: push: branches` | also `integration/**`, `upstream/**`, `upstream-v2/**`, `fix/**` (the fork's working branches) and `workflow_dispatch` | the fork branches can be removed; `workflow_dispatch` is worth keeping |

GitHub Actions runs on the fork. All six legs (Linux release / devel / oldrel-1, Linux without the optional
programs, Windows, macOS; the programs installed with micromamba) passed on `upstream-v2/15-coal` (run
43507540285) and on `integration/all` (run 43505256122, macOS leg re-run after a cancellation), and on
`upstream-v2/01-ci-and-tests` and `upstream-v2/08-annotation`; the runs of the other upstream-v2 branches were
cancelled because the fork's runners start about one job at a time (see `STACK.md` on `upstream-v2/15-coal`). The
CI installs coalr through the `Remotes` field; if the COAL repository cannot be read by the workflow's token, or the ref
disappears, the dependency step fails on every leg - remove the `Remotes` line (coalr is in `Suggests` and every
use is guarded) rather than leave it dangling.

## 4. Proposals that were deliberately not made in 0.5.0, and what 0.5.0 does about them

- **`sensitivity_mode = "default"` as the default of `dNdS()`, `orthologs()`, `orthogroups()`** (now `"fast"`,
  DIAMOND's mode "for hits of more than 90 % identity"). Measured with `diamond_rec()`, 2 cores:

  | reciprocal best hits | `"fast"` | `"default"` | `"sensitive"` |
  |---|---|---|---|
  | *A. thaliana* / *A. lyrata*, 1,000-gene example | 983 | 985 | 985 |
  | *S. cerevisiae* / *C. glabrata* | 4,164 (7.3 s) | 4,384, +5.3 % (8.2 s) | 4,502, +8.1 % (24.4 s) |

  The pairs found only with `"default"` have a median identity of 37 % and a median length of 298 amino acids
  (not only short proteins); the search is about 10 s of a 200 s `dNdS()` run. It would be the third change of
  the ortholog set in one release, so 0.5.0 only documents it (`?dNdS`, the dN/dS vignette). Suggested for 0.6.0.
- **Keeping coding sequences with `N`** (translate the codon to `X` and let the estimator skip it) instead of
  leaving the sequence out. The loss is 0 to 1 sequence in RefSeq sets but 2.5 % (human chr21) to 4.7 % (mouse
  chr10/16/17) of Ensembl transcripts, typically one leading `N`. The default is unchanged ("drop"); 0.5.0
  documents the behaviour, reports the sequences, and offers it as an opt-in: `options(orthologr.ambiguous_cds =
  "keep")` (recorded in `attr(, "orthologr_run")$options`). Making "keep" the default is a decision for 0.6.0.

**Done in 0.5.0** (were listed here as not made):

- **Classed conditions**: `orthologr_error_*` / `orthologr_warning_*` (`?orthologr_options`, "Errors and
  warnings"). The argument and input-file checks of the main functions (`dNdS()`, `dNdS_coal()`,
  `dNdS_clusters()`, `orthogroups()`, `dNdS_orthogroups()`, `divergence_stratigraphy()`, `filter_dNdS()`, the
  alignment functions, `substitutionrate()`, the readers) and every check for a program or package are classed;
  105 calls are classed; about 190 `stop()` calls, mostly in other, older functions and inside single gene
  pairs, are still plain errors (`?orthologr_options` says which functions are classed).
- **Shared code of `dNdS()` and `dNdS_coal()`**, first part: the search step is one function, `dnds_search()`.
  Still duplicated: the argument checks, the start message and the set-up of the call (working folder, sequence
  cache, collection of left-out sequences) - a second step for 0.6.0.

## 5. Changes of this pass that touch existing behaviour

- `orthogroups()` and `dNdS_orthogroups()` have `delete_corrupt_cds = FALSE` (trim), the default of `dNdS()`;
  they had `TRUE`. `dNdS()` itself is unchanged. (`diamond_rec()`, `diamond_best()`, `deepclust()` called on
  their own keep `TRUE`.) For two species `orthogroups()` therefore returns the reciprocal best hits of
  `dNdS()` / `diamond_rec(delete_corrupt_cds = FALSE)`: 995 instead of 983 pairs on the 1,000-gene example.
- `dNdS()` output: a new last column `status`; one start line and one summary instead of 31 messages (also on
  real data: sequences that are left out or trimmed to a multiple of 3 are counted in the summary and listed in
  `attr(, "dropped_sequences")` / `attr(, "trimmed_sequences")`; 1,000-gene example 2 messages, `quiet = FALSE`
  36); nothing on `stdout`. `divergence_stratigraphy()` has `quiet = TRUE` too (was `FALSE`; 45 -> 10 messages).
  Scripts that parse the messages or take columns by position from the end are affected.
- `dNdS(quiet = FALSE)` no longer prints one line per alignment and estimation call; it reports progress.
- The detailed, review-by-review log that used to be the NEWS entry is `MAINTAINER_CHANGELOG.md`.

The first Actions run on the fork failed before any R code ran: micromamba received `"diamond>=2.1"` with its
double quotes (inside a folded YAML block no shell removes them) and could not parse the version, and pak could
not match the `Remotes` entry to coalr. Both are corrected (`diamond>=2.1` unquoted in the two workflows; the
`coalr=` prefix). The later runs pass on all legs (section 3).
