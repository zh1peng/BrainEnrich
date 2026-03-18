# BrainEnrich Refactor Spec (Agent-Ready, No DOSE Ecosystem)

## 0) Mission
Refactor BrainEnrich to fully remove `DOSE`/`enrichplot` runtime dependency and replace it with a native BrainEnrich ecosystem.

This is a full cut-over (no backward compatibility, no compatibility adapter).

## 1) Locked Decisions (Do Not Re-negotiate)
- Branch: `codex/dev` (this is the implementation branch for the requested `dev` workstream).
- No backward compatibility to `gseaResult`.
- No usage of DOSE namespace internals.
- `brainenrich` and `brainscore.lm_test` must both return one shared native class: `EnrichRes`.
- `brainscore` and `brainscore.simulate` stay score/simulation-native outputs (no new class requirement).
- Move large GS resource payloads to a separate package `BrainEnrichData`; keep only a small example GS in `BrainEnrich`.
- Validation defaults:
  - `n_perm = 10` whenever permutation-based testing is needed.
  - Use `SynGO` as the default moderate gene-set source in tests/parity runs.
- HCP-scale validation is deferred and run later on a separate machine.

## 2) Scope of Work
In scope:
- Core architecture refactor.
- Output object redesign.
- Annotation pipeline redesign.
- Correctness bug fixes identified in prior review.
- Parallel/reproducibility hardening.
- Performance upgrades for heavy paths.
- CI/test/doc updates.

Out of scope:
- Any compatibility bridge to DOSE objects.
- HCP-scale long-run benchmark execution.

## 3) Target Architecture

### 3.1 Native result contract
- Shared S3 class: `EnrichRes` (used by `brainenrich` and `brainscore.lm_test`).
  - `analysis_type` (`"brainenrich"` or `"brainscore_lm_test"`)
  - `table` (term-level primary output)
  - `gene_sets`
  - `perm_scores`
  - `params`
  - `diagnostics`
  - Optional fields used only when needed by `brainscore.lm_test`:
    - `model_table`
    - `np_table`
    - `core_genes`
- `brainscore` remains score-oriented output (list/data.frame/matrix contract).
- `brainscore.simulate` remains simulation-oriented output (list/data.frame contract).

### 3.2 Native annotation object
- `EnrichAnno` list contract:
  - `term2gene` with columns `term_id`, `gene_id`
  - `term2name` with columns `term_id`, `term_name`
  - `meta` (source/version/time)

### 3.3 Native plotting interface (first release)
- `plot_terms(x, type = c("dot","bar","volcano","ridge"))`
- `plot_core_genes(x, term_id, mode = c("impact","direction"))`
- `plot_heatmap_terms(x, ...)`
- `plot_term_network(x, ...)` (network-lite)

## 4) File-Level Execution Plan

## Phase A - Dependency Decoupling
1. Update `DESCRIPTION`
   - Remove `DOSE` from `Imports`.
   - Remove `enrichplot` from `Suggests` if not used in new docs.
2. Update `NAMESPACE`
   - Remove `import(DOSE)` and `importClassesFrom(DOSE,gseaResult)`.
   - Remove imports that only existed for DOSE bridge.
3. Remove all DOSE calls from code:
   - `getFromNamespace("calculate_qvalue","DOSE")`
   - `getFromNamespace("TERM2NAME","DOSE")`
   - `getFromNamespace("getGeneSet","DOSE")`
   - `getFromNamespace("geneSet_filter","DOSE")`
   - `getFromNamespace("build_Anno","DOSE")`
4. Remove all `new("gseaResult", ...)` paths.

Definition of done:
- No runtime reference to DOSE in `R/*.R`, `NAMESPACE`, or exported outputs.

## Phase B - Output Object Refactor
1. Add object constructors and validators in a new module (for example `R/result_objects.R`).
2. Refactor:
   - `R/brainenrich.R` to return `EnrichRes`.
   - `R/brainscore.lm_test.R` to return `EnrichRes` with `analysis_type = "brainscore_lm_test"`.
   - Keep `R/brainscore.R` and `R/brainscore.simulate.R` as non-class score/simulation outputs (do not introduce new classes).
3. Ensure constructor/validator and print/summary/data.frame helpers exist for `EnrichRes`.


Definition of done:
- `brainenrich` and `brainscore.lm_test` return validated `EnrichRes`.
- `brainscore` and `brainscore.simulate` keep stable non-class output contracts.

## Phase C - Annotation Contract Refactor
1. Refactor `R/get_geneSets.R`:
   - `get_annoData()` returns `EnrichAnno`.
   - `get_geneSetList()` derives list from `EnrichAnno$term2gene`.
   - `get_termDescription()` uses `EnrichAnno$term2name`.
   - `filter_geneSetList()` rewired without DOSE filter helper.
   - `Table2Anno()` and `Anno2Table()` convert to/from `EnrichAnno`.
2. Keep column naming stable and explicit (`term_id`, `gene_id`, `term_name`).

Definition of done:
- Annotation helpers are fully native and deterministic.

## Phase D - Mandatory Correctness Fixes
1. Fix empty-result branch crash in `brainscore.lm_test`.
2. Unify method API in `brainenrich`/`brainscore`:
   - Accept string method names or user-supplied function.
   - Avoid `match.arg` rejection when function provided.
3. Remove `readline <- NULL` in `R/misc.R`.
4. Fix `%%` boolean misuse in `R/hpc_helper.R`.
5. Core-gene directionality:
   - Keep impact score (`abs delta`) and add direction labels:
     - `driver` (same direction as term effect)
     - `buffer` (opposite direction)

Definition of done:
- Prior reviewed correctness issues are covered by tests and closed.

## Phase E - Data Strategy and Packaging
1. Prepare split-package strategy:
   - `BrainEnrich`: algorithms and minimal example data, including a small bundled GS example.
   - `BrainEnrichData`: large GS payloads (full gene-set collections).
2. Remove write-to-install-dir behavior in loaders.
3. Add checksum/version manifest for large resource files.

Definition of done:
- Main package no longer depends on writable install directory behavior.
- Large GS resources are sourced from `BrainEnrichData`, while `BrainEnrich` remains lightweight with a small example GS.

## Phase F - Parallel Safety and Performance
1. Cluster safety:
   - Wrap all clusters with `on.exit(stopCluster(cl), add = TRUE)`.
2. Reproducibility:
   - Add deterministic RNG stream handling in parallel blocks.
3. Performance target 1:
   - Optimize `aggregate_geneSetList` via precomputed indices / matrix path.
4. Performance target 2:
   - Refactor `simple_lm` to matrix batch regression (`model.matrix + lm.fit`).

Definition of done:
- No orphan worker leaks on error and measurable speedup on benchmark fixtures.

## Phase G - CI, Tests, and Docs
1. CI:
   - Remove commit-message-only trigger gates.
   - Ensure PR/Push run checks by default.
2. Tests:
   - Remove runtime package installs in tests.
   - Remove mandatory network fetches in default test path.
   - Add deterministic seeds.
3. Docs:
   - Rewrite examples to `EnrichRes` (for enrichment outputs) and native score/simulation outputs.
   - Fix README encoding and code-block closure issues.

Definition of done:
- CI green by default, tests deterministic offline, docs reflect native API.

## 5) Validation Protocol (Agent Must Follow Exactly)

## 5.1 Deterministic test configuration
- Set fixed seed for all parity/test scripts.
- Force `n_cores = 1` in parity validation.
- Use `n_perm = 10` for permutation workflows.
- Use `SynGO` in validation pipelines as the default moderate gene set.

## 5.2 Parity coverage
Compare baseline (pre-refactor) vs refactor for:
- `brainenrich`
- `brainscore`
- `brainscore.lm_test`
- `brainscore.simulate`

Required checks:
- ID/key/set membership exact match.
- Table dimensions exact match.
- Numeric outputs tolerance match (`all.equal`-style, tolerance `1e-6`).

## 5.3 Utility smoke tests
- `job_splitter` and `job_cat` complete a sample cycle.
- Null-model helpers execute successfully at sample scale (`n_perm=10`).

## 5.4 Reporting artifacts
- Emit parity report:
  - machine-readable (`.rds` or `.json`)
  - markdown summary (`PASS/FAIL`, mismatch table)
- Record exact test settings (`seed`, `n_perm=10`, dataset source, gene set source `SynGO`).

## 6) Acceptance Criteria (Release Gate)
- `library(BrainEnrich)` loads without DOSE/enrichplot runtime behavior.
- No DOSE/gseaResult usage remains in core runtime path.
- Four main workflows pass parity checks with tolerance criteria.
- Mandatory bug fixes are validated by tests.
- CI passes on default triggers.
- `brainenrich` and `brainscore.lm_test` return `EnrichRes` with correct `analysis_type`.
- `brainscore` and `brainscore.simulate` are not forced into new class wrappers.
- Documentation matches `EnrichRes` + native score/simulation outputs and plotting APIs.

## 7) Suggested Task Breakdown for Multi-Agent Execution
- Agent A: dependency/object refactor (Phases A+B)
- Agent B: annotation/data refactor (Phases C+E)
- Agent C: correctness + parallel/performance refactor (Phases D+F)
- Agent D: CI/tests/docs + parity harness (Phase G + validation protocol)

Coordination rule:
- Merge order: A -> C -> B -> D, then final parity run.

## 8) Handoff Notes for High-Compute Machine
- HCP-level utility/stress testing is not required for this implementation pass.
- After merge on `codex/dev`, run HCP validation separately and append report.
