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
| T-Coffee default mode (was `-mode accurate`) | `aa_aln_params = "-mode accurate"` / `params = "-mode accurate"` |
| `aa_aln_params` is passed on to the alignment program (was ignored) | do not set it |
| Parameter strings are argument lists, not shell text (no `;`, pipes, redirections) | none (security fix) |
| A missing program (perl, aligner, KaKs_Calculator) is an error before the search (was: all pairs `NA` + warning on this fork; an error at the first pair in 0.4.3) | none |
| Annotation lines of more than 65,000 characters are read (were deleted) | none (bug fix); `check_annotation(remove_annotation_outliers = TRUE)` still writes a copy without them |

**Not changed any more** (a default introduced earlier on this fork was taken back after the third review):
`filter_dNdS()` and `divergence_stratigraphy()` leave gene pairs with dS = 0 out, as 0.4.3 did
(`dS_zero = "drop"`); new are a message with the numbers, `attr(, "left_out")` / `divergence_left_out()`, and
`dS_zero = "keep"` / `options(orthologr.dS_zero = "keep")` to keep them.

With all legacy options set, the 20-gene and 1,000-gene examples give the tables of 0.4.3 (master 06c3098)
for Comeron and YN; see NEWS.

## 3. Things that are set for the fork and must change when this is merged

| Where | Now | Change to |
|---|---|---|
| `DESCRIPTION`, `Remotes:` | `Plumat-COAL/COAL/packaging/r/coalr@improve/pass2` (syntax `owner/repo/subdir@ref`). The R package coalr is not on COAL's default branch yet; `improve/pass2` is the branch that has it with the Comeron / Li methods. | `Plumat-COAL/COAL/packaging/r/coalr` (default branch) or `...@<release tag>` once COAL's pull requests are merged |
| the coalr install command in `README.md`, `?codon_aln`, `?dNdS_coal`, `?dNdS_clusters` and in the messages of `R/coal.R` | `remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr", ref = "improve/pass2")` | the same without `ref` (or with the release tag) |
| install command of orthologr in `README.md` (two places), `environment.yml`, `vignettes/orthologr.Rmd` | `remotes::install_github("Plumat-COAL/orthologr@integration/all")`, so that the command installs the version the page describes | `remotes::install_github("drostlab/orthologr")` |
| `.github/workflows/R-CMD-check.yaml`, `on: push: branches` | also `integration/**`, `upstream/**`, `upstream-v2/**`, `fix/**` (the fork's working branches) and `workflow_dispatch` | the fork branches can be removed; `workflow_dispatch` is worth keeping |

GitHub Actions has not run on the fork (Actions are not enabled in the fork's settings): the macOS and Windows
legs of the workflow and the installation of the programs with micromamba are untested there. The CI installs
coalr through the `Remotes` field; if the COAL repository cannot be read by the workflow's token, or the ref
disappears, the dependency step fails on every leg - remove the `Remotes` line (coalr is in `Suggests` and every
use is guarded) rather than leave it dangling.
