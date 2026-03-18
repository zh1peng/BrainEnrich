# BrainEnrich Full DOSE-Decoupling Refactor Plan

## Summary
- Work on branch `codex/dev`.
- Fully remove runtime dependencies on `DOSE` and `enrichplot`.
- Do not keep backward compatibility with DOSE objects/ecosystem.
- Build a BrainEnrich-native object model, annotation model, and plotting layer.
- Prove output parity on sample data using numeric tolerance checks.

## Execution Steps
1. **Create and prepare branch**
   - Create/switch to `codex/dev`.
   - Freeze a baseline from current package behavior for parity comparison.

2. **Remove DOSE coupling from core**
   - Remove `DOSE` imports/usages from `DESCRIPTION`, `NAMESPACE`, and R code.
   - Remove all `gseaResult` construction and DOSE namespace calls.
   - Refactor:
     - `brainenrich` -> native `be_enrich_result`
     - `brainscore` -> native `be_score_result`
     - `brainscore.lm_test` -> native `be_lm_result`
     - `brainscore.simulate` -> native `be_sim_result`

3. **Replace annotation/data contracts**
   - Introduce native annotation object `be_anno` with:
     - `term2gene(term_id, gene_id)`
     - `term2name(term_id, term_name)`
     - `meta`
   - Refactor `get_geneSetList`, `get_termDescription`, `filter_geneSetList`, `Table2Anno`, `Anno2Table` to native contract.

4. **Integrate correctness fixes during refactor**
   - Fix `brainscore.lm_test` empty-result branch stability.
   - Make `cor_method/aggre_method` consistently accept both strings and user functions.
   - Remove `readline <- NULL` behavior.
   - Fix logical misuse in `job_cat` (`%%` in boolean logic).

5. **Data strategy and packaging**
   - Keep core compute logic in `BrainEnrich`.
   - Move large resource payloads to separate `BrainEnrichData`.
   - Remove write-to-package-install-dir download behavior.
   - Add dataset version/checksum manifest.

6. **Performance and parallel reliability**
   - Add safe cluster lifecycle (`on.exit(stopCluster, add = TRUE)`).
   - Add deterministic parallel RNG stream controls.
   - Optimize `aggregate_geneSetList` backend to reduce repeated gene-set hit recomputation.
   - Replace `simple_lm` nested tidyverse loop with matrix/batch regression path.

7. **Native plotting and docs**
   - Implement first-party plotting APIs for core chart set:
     - dot, bar, volcano, ridge, heatmap, network-lite.
   - Rewrite docs/vignettes to remove DOSE/enrichplot assumptions.
   - Fix README encoding and code-block issues.

8. **CI and test hardening**
   - Remove fragile commit-message-only triggers in CI.
   - Eliminate test-time package installation and mandatory network dependency.
   - Keep deterministic seeds and stable local fixtures for default test runs.

## Validation Protocol
1. **Parity baseline setup**
   - Capture baseline outputs from current implementation on canonical sample data.
   - Use deterministic execution settings for both baseline and refactored runs:
     - fixed seed
     - `n_cores = 1`
     - fixed `n_perm`
     - fixed sample datasets

2. **Workflow parity checks**
   - Compare outputs for:
     - `brainenrich`
     - `brainscore`
     - `brainscore.lm_test`
     - `brainscore.simulate`
   - Enforce:
     - exact match: IDs, dimensions, ordering keys, term membership
     - numeric match: tolerance-based equality (`all.equal` style), target tolerance `1e-6`

3. **Utility and null-model smoke checks**
   - Validate utility paths on sample scale:
     - `job_splitter`
     - `job_cat`
     - key null-model helpers (`spin_brain`, `resample_gene`, `coexp_matched`) for functional execution

4. **Reporting**
   - Save parity artifacts in machine-readable form.
   - Produce a short markdown parity summary with pass/fail status and mismatches.
   - Note that HCP-scale utility validation is intentionally deferred for execution on a higher-compute machine.

## Acceptance Criteria
- Package loads and runs with **no runtime dependency** on `DOSE` or `enrichplot`.
- Core outputs are native BrainEnrich objects (`be_enrich_result`, `be_score_result`, `be_lm_result`, `be_sim_result`).
- All listed correctness fixes are merged and covered by tests.
- Sample-data parity checks pass for all four core workflows under tolerance rules.
- CI runs by default on push/PR and passes.
- Documentation is updated to native ecosystem semantics and README rendering issues are fixed.


修复 brainscore.lm_test 空结果分支的潜在崩溃
理由：nrow(res)==0 时后面仍使用 selected.gs 构建 gseaResult，有未定义风险。
位置：brainscore.lm_test.R (line 233), brainscore.lm_test.R (line 279)

统一“自定义方法”API，修正文档与实现不一致
理由：文档说 cor_method/aggre_method 支持函数，但入口函数先 match.arg() 会直接拒绝函数。
位置：brainenrich.R (line 16), brainenrich.R (line 115), brainscore.R (line 11), brainscore.R (line 54)

修复两个明确逻辑问题
理由：会导致运行时异常或逻辑失效。
位置：readline <- NULL 会破坏交互输入 misc.R (line 25)；%% 被误用为逻辑与 hpc_helper.R (line 196)

改造数据下载/缓存路径策略
理由：现在写到包安装目录 find.package()/extdata，在只读库常失败；建议改为用户缓存目录并加 checksum/version manifest。
位置：get_geneExp.R (line 41), get_geneSets.R (line 52)

并行稳定性增强（on.exit(stopCluster) + RNG 流控制）
理由：当前大量 makeCluster/stopCluster，中途报错会遗留 worker；并且并行随机过程可复现性不足。
位置示例：brainscore.R (line 120), find_core_genes.R (line 49), resample_gene_functions.R (line 91)

性能优化优先做两块
理由：这是你包里最重计算路径。
位置：

aggregate_geneSetList 里重复 gene-set 命中计算，可改成预计算索引/稀疏矩阵批量聚合 aggregate_functions.R
simple_lm 现在是 tidyverse 嵌套逐模型，建议改 model.matrix + lm.fit 批处理 simple_lm.R
CI 和测试体系建议“去脆弱化”
理由：目前 CI 依赖 commit message 才触发，容易漏检；测试里还有运行期安装包/网络依赖。
位置：R-CMD-check.yml (line 15), CodeCov.yml (line 16), tests/testthat.R (line 15)

发布工程化建议（特别是体积）
理由：仓库当前约 203MB，大数据文件多；后续 CRAN/安装体验会受影响，建议拆分“核心包 + 数据包/远程资源”。
相关目录：inst/extdata

文档可读性快速修复
理由：README 有乱码和未闭合代码块，影响首印象和安装指引。
位置：Readme.md (line 19), Readme.md (line 77)

## Out of Scope (Current Cycle)
- No backward compatibility bridge to `gseaResult` or enrichplot ecosystem.
- Full HCP-scale benchmarking/validation (to be executed later on dedicated compute environment).


BrainEnrich 全量去 DOSE 重构计划（无后向兼容，建立自有生态）
1) 摘要
目标：彻底移除 DOSE/enrichplot 依赖，重建结果对象、注释结构、绘图体系与数据分发方式。
策略：一次性全量重构 v1，首版提供核心 4-6 类可视化，数据拆分为独立 BrainEnrichData 包。
结果：BrainEnrich 成为轻量核心计算包；BrainEnrichData 负责大数据；不做任何后向兼容桥接。
2) 关键实现变更（按优先级）
P0-核心架构替换

删除 DOSE（Imports）与 enrichplot（Suggests）依赖，清理 NAMESPACE 中所有 DOSE 引用。
新建原生结果对象（S3）：
be_enrich_result（群体富集）
be_score_result（个体得分/空模型得分）
be_lm_result（个体层统计检验）
be_sim_result（仿真结果）
新建原生注释对象 be_anno：term2gene(term_id,gene_id) + term2name(term_id,term_name) + meta。
brainenrich、brainscore、brainscore.lm_test、brainscore.simulate 保留函数名，但返回上述新对象（不再返回 gseaResult）。
P0-正确性修复（必须并入本次）

修复 brainscore.lm_test 空结果分支未定义对象风险。
统一 cor_method/aggre_method：支持字符与函数两种输入，移除入口 match.arg 导致的函数拒绝。
修复 readline <- NULL 破坏交互与 %% 误用逻辑运算问题。
引入方向性核心基因输出：在现有“影响大小”之外，增加 driver/buffer 方向标签，避免 LOO 解释冲突。
P1-性能与并行稳定性

聚合加速：对 mean/meanabs/meansqr/maxmean 提供矩阵化/稀疏化后端，减少重复命中计算。
simple_lm 重写为 model.matrix + lm.fit 批处理路径（替代 tidyverse 嵌套逐模型）。
并行执行统一封装：on.exit(stopCluster, add=TRUE)、clusterSetRNGStream、统一 n_cores/seed 规则。
coexp_matched 采样流程优化：去重与匹配过程减少高频 list 去重开销。
P1-数据与发布工程化

新建 BrainEnrichData 独立包，迁移 inst/extdata 大文件（表达矩阵、基因集、perm 资源）。
主包仅保留小型示例数据；get_geneExp/get_annoData 改为从 BrainEnrichData 读取（无写入包安装目录逻辑）。
为数据增加版本号与 checksum 清单，确保可重复性和资源一致性。
P2-绘图生态与文档

提供首版原生绘图（核心图先齐）：dot/bar/volcano/ridge/heatmap/network-lite。
重写 vignette，移除 enrichplot/clusterProfiler/DOSE 示例路径。
修复 README 编码乱码与未闭合代码块，统一安装/排错指引。
P2-CI/测试去脆弱化

移除“按 commit message 才触发”的 CI 条件，改为 PR/Push 默认触发。
测试分层：离线单元测试（默认）+ 数据/集成测试（显式触发）。
清除测试中的运行期安装与网络硬依赖；用固定 seed + 本地最小夹具数据保证稳定性。
3) 公共接口与输出变更（决策已锁定）
brainenrich(...) -> be_enrich_result
result_table（term 级统计，含 p_value/p_adj、核心基因字段）
gene_scores、null_scores、gene_sets、params、diagnostics
brainscore(...) -> be_score_result
scores（empirical 或 null 列表/矩阵）、params、diagnostics
brainscore.lm_test(...) -> be_lm_result
model_table、np_table、core_genes、params
brainscore.simulate(...) -> be_sim_result
sim_table、meta
新增通用绘图接口（原生）：
plot_terms(x, type=...)
plot_core_genes(x, term_id, mode=...)
删除/弃用：所有 gseaResult 相关依赖与返回假设（无兼容层）。
4) 测试计划与验收标准
单元测试
注释对象转换、基因集过滤、p 值计算、LOO 核心基因（含方向标签）、聚合后端一致性。
集成测试
四条主流程（群体、个体、lm_test、simulate）在最小数据集下稳定可复现。
三类 null model 与 HPC job_splitter/job_cat 回归测试。
性能回归
aggregate_geneSetList 与 simple_lm 在固定基准数据上的耗时显著下降（目标 >=30%）。
工程验收
library(BrainEnrich) 不再触发 DOSE 相关加载/提示。
DESCRIPTION/NAMESPACE 无 DOSE/enrichplot 运行依赖。
CI 默认触发并稳定通过。
5) 里程碑与默认假设
里程碑
M1（架构与接口定稿）：结果对象/注释对象/API 签名冻结。
M2（核心替换+P0 修复）：主流程可跑且无 DOSE。
M3（性能+并行+数据包）：速度、稳定性、数据拆包落地。
M4（绘图+文档+CI 测试）：对外发布候选版。
默认假设（已确认）
全量重构 v1，一次切换，不做后向兼容。
数据采用独立 BrainEnrichData 包。
首版绘图深度为“核心图先齐”。