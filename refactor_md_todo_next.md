# BrainEnrich Refactor - Next TODO (Execution Handoff)

This file lists only the **remaining work** after the current refactor pass, in execution order.

## 0) Environment & branch
1. Checkout branch: `codex/dev`
2. Ensure R packages installed: `devtools`, `testthat`, `roxygen2`, `ggplot2`, `mockery`
3. Work from repo root: `E:/xhmhc/BrainEnrich`

## 1) Hardening pass (code-level)
1. Re-check there is no runtime DOSE/enrichplot dependency in code paths:
   - Search keywords: `DOSE`, `gseaResult`, `enrichplot`, `getFromNamespace("...","DOSE")`
2. Remove/replace old README install/troubleshooting text that still tells users to install DOSE/enrichplot.
3. Confirm `brainenrich()` and `brainscore.lm_test()` both return `EnrichRes`, with:
   - `analysis_type` set correctly
   - stable `table` schema
4. Confirm `brainscore()` and `brainscore.simulate()` are unchanged as native score/simulation outputs (no new class wrapping).

## 2) Native plotting API (if still missing or incomplete)
Implement/export first-party plotting functions (minimum set):
1. `plot_terms(x, type = c("dot","bar","volcano","ridge"))`
2. `plot_core_genes(x, term_id, mode = c("impact","direction"))`
3. `plot_heatmap_terms(x, ...)`
4. `plot_term_network(x, ...)` (lightweight term-term network)

Notes:
- Input class target: `EnrichRes`
- Keep implementation lightweight; avoid new heavy runtime dependencies.
- Use clear error messages when required columns are missing.

## 3) Parity harness (must add if not present)
Create scripts under `scripts/validation/`:
1. `run_parity_baseline.R`
2. `run_parity_refactor.R`
3. `compare_parity.R`

Required settings:
- fixed seed
- `n_cores = 1`
- `n_perm = 10` for all permutation paths
- gene set source: `SynGO`

Coverage to compare:
1. `brainenrich`
2. `brainscore`
3. `brainscore.lm_test`
4. `brainscore.simulate`

Output artifacts:
1. machine-readable report (`.rds` or `.json`)
2. markdown summary (`reports/parity_summary.md`) with PASS/FAIL + mismatch details

## 4) Tests (deterministic and offline-first)
1. Ensure tests do not install packages at runtime.
2. Ensure default tests do not require network.
3. Keep permutation tests using `n_perm = 10`.
4. Prefer `SynGO` in enrichment workflow tests.
5. Run targeted tests first:
   - `test-brainenrich.R`
   - `test-brainscore.lm_test.R`
   - `test-get_annoData.R`
   - `test-get_geneExp.R`
   - `test-job_cat.R`
   - `test-job_splitter.R`
6. Then run full suite on high-compute machine.

## 5) Documentation/NAMESPACE sync
1. Run:
   - `Rscript -e "devtools::document()"`
2. Verify:
   - `NAMESPACE` exports/imports are correct
   - man pages reflect `EnrichRes`/`EnrichAnno`
   - no DOSE-centric return class text remains

## 6) CI and packaging checks
1. Confirm GitHub workflows trigger on normal push/PR (no commit-message gate).
2. Validate package build:
   - `R CMD build .`
3. Optional but recommended:
   - `R CMD check --as-cran <tarball>`

## 7) Data strategy follow-through
1. Keep large GS payload plan on `BrainEnrichData`.
2. In main package (`BrainEnrich`), keep only small example GS.
3. Confirm loaders prefer:
   - `BrainEnrichData::extdata` if installed
   - local/cache fallback otherwise
4. Add/maintain checksum + version manifest for large resources.

## 8) Final acceptance checklist
All items below must be true before merge:
1. `library(BrainEnrich)` shows no DOSE/enrichplot runtime behavior.
2. No runtime DOSE references remain in `R/`, `NAMESPACE`, `DESCRIPTION`.
3. `brainenrich` and `brainscore.lm_test` return valid `EnrichRes`.
4. `brainscore` and `brainscore.simulate` remain native (no forced class conversion).
5. Parity validation passes (or mismatches documented with concrete reasons).
6. Tests pass under deterministic settings.
7. Docs/README are updated to native ecosystem terminology.

## 9) Suggested execution command sequence
```powershell
git checkout codex/dev

Rscript -e "devtools::load_all('.')"
Rscript -e "testthat::test_file('tests/testthat/test-brainenrich.R', reporter='summary')"
Rscript -e "testthat::test_file('tests/testthat/test-brainscore.lm_test.R', reporter='summary')"
Rscript -e "testthat::test_file('tests/testthat/test-get_annoData.R', reporter='summary')"
Rscript -e "testthat::test_file('tests/testthat/test-get_geneExp.R', reporter='summary')"
Rscript -e "testthat::test_file('tests/testthat/test-job_cat.R', reporter='summary')"
Rscript -e "testthat::test_file('tests/testthat/test-job_splitter.R', reporter='summary')"

Rscript -e "devtools::document()"
Rscript -e "devtools::test()"
```
