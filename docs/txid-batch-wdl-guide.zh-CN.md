# TxID 多 GTF 批处理 WDL 使用手册

本文说明如何使用
[`txid_batch.wdl`](../workflows/txid_batch.wdl)，在同一个参考基因组和参考注释
上下文中导入多个独立产生的 GTF，生成稳定的 TxID、逐输入重写 GTF、映射表、
cohort catalog 和版本化 SQLite registry。

> **适用范围**：这里的“多个 GTF”可以来自不同样本、不同上游工具，或同一样本的
> 不同工具结果。所有 GTF 必须属于同一个参考基因组 assembly。不同 assembly 不能
> 放进同一次工作流运行。

> **功能边界**：TxID 是转录本身份和互操作层，不负责转录本发现、reads 比对或表达量
> 估计。WDL 输入应是 IsoQuant、FLAIR、StringTie、TALON 等工具已经生成的转录本 GTF。

## 1. 工作流概览

工作流使用 WDL 1.1，是不含 `import` 的单文件工作流。执行关系如下：

```text
reference_fasta + reference_gtf + assembly
                    │
                    ▼
                txid init
                    │
                    ▼
                 gtfs
                    │
                    ▼
        生成 import-manifest.tsv
                    │
                    ▼
     txid batch（按确定性顺序逐个导入）
        ┌───────────┼───────────┐
        ▼           ▼           ▼
  rewritten GTF  mapping TSV  SQLite registry
                                │
                     ┌──────────┴──────────┐
                     ▼                     ▼
                txid export          txid validate
                     │                     │
                     ▼                     ▼
              cohort catalog       validation JSON
```

主要阶段为：

1. 读取参考 FASTA，计算 sequence-collection fingerprint，并用参考 GTF 初始化 registry。
2. 从每个 GTF 文件名自动生成样本标识，并写入批处理 manifest。
3. 对输入 GTF 进行完整解析、contig/坐标验证、精确身份计算和参考分类。
4. 依次、原子地写入同一个 SQLite registry，并为每个输入生成重写 GTF 和 mapping TSV。
5. 导出 cohort catalog，并重新计算和校验 registry 中保存的 digest。

工作流故意不对 GTF 使用 WDL `scatter`。WDL task 的文件系统彼此隔离，scatter
不能安全地让多个容器共同修改一个 SQLite registry。单个 task 内的导入顺序是确定的，
输入数组顺序、manifest 行顺序不会改变精确结构 ID。

## 2. 运行前准备

### 2.1 所需文件

提交前准备以下数据：

| 文件 | 要求 |
|---|---|
| 参考 FASTA | 与上游比对和转录本生成使用的 assembly 一致；支持 FASTA 或 `.gz` |
| 参考 GTF | 定义 known transcript/gene 的注释上下文；支持 GTF 或 `.gtf.gz` |
| 输入 GTF | 至少一个，每个文件包含 exon-derived transcript models |
| contig alias 表 | 可选；仅在 GTF 与 FASTA contig 名称不一致时使用 |

参考 FASTA 和参考 GTF 必须满足：

- FASTA contig 名称唯一，序列非空。
- GTF 坐标为 1-based closed。
- 每条 exon 记录包含 `transcript_id`；建议同时包含 `gene_id`。
- 每个 transcript 的 exon 位于同一 contig、同一链，且不能互相重叠。
- GTF 中的最大坐标不能超过 FASTA 对应 contig 长度。
- 所有输入 GTF 必须来自同一个 FASTA assembly。

TxID 不会自动假设 `1` 等于 `chr1`。需要别名时，提供两列、tab 分隔的文件：

```text
alias	primary
1	chr1
2	chr2
MT	chrM
```

`primary` 必须是 FASTA 中真实存在的名称。不要用简单删除 `chr` 前缀的方式批量制造
别名；别名表应来自经过核对的参考元数据。

### 2.2 执行环境

需要满足以下条件：

- WDL 引擎支持 WDL 1.1，例如 MiniWDL 或 Cromwell 92。
- 后端能够运行 Docker，或平台能够把 `runtime.docker` 映射到容器运行时。
- 引擎能够访问输入 JSON 中的所有文件；云平台通常需要先上传并替换成本平台 URI。
- 当前发布镜像是 Linux/amd64。
- 运行磁盘能够容纳本地化输入、SQLite registry、重写 GTF 和 mapping 文件。

默认 Docker 镜像为：

```text
dawang02/txid:0.1.2-jupyter
```

它包含 TxID 0.1.2 和 JupyterLab 4.6.2。可先检查：

```bash
docker pull dawang02/txid:0.1.2-jupyter
docker run --rm dawang02/txid:0.1.2-jupyter txid --version
docker run --rm dawang02/txid:0.1.2-jupyter jupyter lab --version
```

论文或归档分析建议把输入 JSON 中的镜像改为不可变 digest：

```text
dawang02/txid@sha256:889dd70c540fa2895820ec457e67b6a47676ceff6e9402a4cd395b8ecaaffd92
```

## 3. 多个 GTF 如何组织

工作流只需要一个普通的一维文件数组：

```wdl
Array[File] gtfs
```

`sample` 自动取 GTF 的文件名并去掉 `.gtf` 或 `.gtf.gz` 后缀，`tool` 固定记录为
`unspecified`：

| 下标 | `gtfs` | 自动生成的 `sample` | manifest 中的 `tool` |
|---:|---|---|---|
| 0 | `/data/S1.isoquant.gtf` | `S1.isoquant` | `unspecified` |
| 1 | `/data/S1.flair.gtf` | `S1.flair` | `unspecified` |
| 2 | `/data/S2.isoquant.gtf.gz` | `S2.isoquant` | `unspecified` |

对应 JSON 为：

```json
{
  "TxIDBatch.gtfs": [
    "/data/S1.isoquant.gtf",
    "/data/S1.flair.gtf",
    "/data/S2.isoquant.gtf.gz"
  ]
}
```

必须遵守以下规则：

- `gtfs` 至少包含一个文件。
- 数组中不使用 `null`。
- `annotation_name` 不能为空，不能包含 tab 或换行。
- 建议每个 GTF 使用可识别且互不重复的文件名，因为文件名会成为自动生成的 `sample`。
- 同一个文件不应重复列出，除非确实希望记录重复 observation。
- 若必须准确记录真实样本名和上游工具 provenance，请直接使用 `txid batch` CLI 的
  manifest 接口；这个简化 WDL 不要求额外的平行数组。

与 CycloneCell 的 `Array[Array[File?]] fastq_groups` 不同，TxID 的每个元素已经是一个
完整 GTF，不需要内层分组，也不使用数组后的 `+` 非空量词。非空约束由 task 内部
显式检查，以兼容目标 WDL 平台。

## 4. 输入参数

JSON 键使用 `工作流名.参数名` 的完整形式，即 `TxIDBatch.<参数名>`。

### 4.1 数据和身份上下文参数

| 参数 | WDL 类型 | 必填 | 默认值 | 说明 |
|---|---|---:|---|---|
| `gtfs` | `Array[File]` | 是 | — | 需要赋予 TxID 的一个或多个输入 GTF |
| `reference_gtf` | `File` | 是 | — | known transcript/gene 的参考注释 |
| `reference_fasta` | `File` | 是 | — | 身份上下文使用的参考基因组序列 |
| `assembly` | `String` | 是 | — | 用户可读 assembly 名称，例如 `GRCh38` |
| `annotation_name` | `String` | 否 | `reference` | 注释版本名称，例如 `GENCODE-v49` |
| `contig_aliases` | `File?` | 否 | 未设置 | `alias` 到 FASTA `primary` contig 的显式映射 |

`assembly` 名称会记录到 registry；结构身份同时由 FASTA 内容计算出的 assembly
fingerprint 约束。同名但序列内容不同的 FASTA 不会被视为同一身份上下文。

当前 WDL manifest 固定写入 `format=gtf`，因此此 WDL 只接受 GTF。TxID CLI 本身支持
GFF3，但若需要 GFF3，应使用 CLI 或另行扩展 WDL 的逐输入格式参数。

### 4.2 精确模式与 fuzzy 模式

| 参数 | WDL 类型 | 必填 | 默认值 | 说明 |
|---|---|---:|---|---|
| `fuzzy_splice_tolerance_bp` | `Int?` | 否 | 未设置 | splice boundary 容差，单位 bp |
| `fuzzy_end_tolerance_bp` | `Int?` | 否 | 未设置 | TSS/TES 或单 exon 两端容差，单位 bp |

默认是精确模式，不需要任何容差。若启用 fuzzy 模式：

- 两个参数必须同时设置。
- 两个参数必须是大于或等于 0 的整数。
- fuzzy grouping 使用确定性的 complete-linkage 规则。
- fuzzy cluster ID 为 `txid:FC1.*`，不会替代或修改精确 `SC1`、`TF1`、`SE1` ID。
- 可能连接多个 cluster 的 ambiguous bridge 会被显式标记，不会进行不安全的传递合并。

例如：

```json
{
  "TxIDBatch.fuzzy_splice_tolerance_bp": 5,
  "TxIDBatch.fuzzy_end_tolerance_bp": 25
}
```

如果分析只需要可复现的精确身份，建议省略这两个字段。

### 4.3 运行资源参数

| 参数 | WDL 类型 | 必填 | 默认值 | 说明 |
|---|---|---:|---|---|
| `docker_image` | `String` | 否 | `dawang02/txid:0.1.2-jupyter` | task 使用的运行镜像 |
| `cpu` | `Int` | 否 | `1` | WDL 后端为整个 task 申请的 CPU |
| `memory_gb` | `Int` | 否 | `4` | WDL 后端为整个 task 申请的内存，GB |
| `disk_gb` | `Int` | 否 | `20` | WDL 后端本地磁盘请求，GB |

当前 TxID batch 在一个 task 内顺序处理输入，`cpu` 调大不会自动并行导入。参考 FASTA
解析期间会在内存中保存序列内容，因此大型 assembly 或大量 transcript models 可能需要
提高 `memory_gb`。

WDL 内置的 4 GB/20 GB 默认值适合小型或试运行数据。仓库示例针对全基因组使用场景
从 16 GB 内存和 50 GB 磁盘开始；实际生产值应根据参考 FASTA、参考 GTF、输入 GTF
总大小和平台本地化开销调整。

`disks` runtime 属性的解释依赖 WDL 后端。有些本地 MiniWDL 配置会忽略它；此时应在
MiniWDL/Docker 的运行目录或后端配置中保证实际磁盘空间。`disk_gb` 至少应覆盖全部
本地化输入、重写 GTF、mapping、registry 和工作目录开销。

## 5. 完整输入 JSON 示例

仓库提供可直接复制修改的
[`txid_batch.inputs.example.json`](../workflows/txid_batch.inputs.example.json)。

下面是两个 GTF 的精确模式示例：

```json
{
  "TxIDBatch.gtfs": [
    "/data/sample01.gtf",
    "/data/sample02.gtf"
  ],
  "TxIDBatch.reference_gtf": "/reference/gencode.v49.gtf",
  "TxIDBatch.reference_fasta": "/reference/GRCh38.primary_assembly.genome.fa",
  "TxIDBatch.assembly": "GRCh38",
  "TxIDBatch.annotation_name": "GENCODE-v49",
  "TxIDBatch.docker_image": "dawang02/txid:0.1.2-jupyter",
  "TxIDBatch.cpu": 1,
  "TxIDBatch.memory_gb": 16,
  "TxIDBatch.disk_gb": 50
}
```

如果需要 contig alias 和 fuzzy grouping，可增加：

```json
{
  "TxIDBatch.contig_aliases": "/reference/GRCh38.contig-aliases.tsv",
  "TxIDBatch.fuzzy_splice_tolerance_bp": 5,
  "TxIDBatch.fuzzy_end_tolerance_bp": 25
}
```

JSON 不支持注释。提交前可检查语法：

```bash
python -m json.tool my-txid-inputs.json >/dev/null
```

## 6. 提交工作流

以下命令均在 TxID 仓库根目录执行。

### 6.1 MiniWDL

先做静态检查：

```bash
miniwdl check workflows/txid_batch.wdl
```

运行工作流：

```bash
miniwdl run \
  workflows/txid_batch.wdl \
  -i my-txid-inputs.json \
  --dir output/txid-cohort
```

MiniWDL 会返回包含 `dir` 和 `outputs` 的 JSON。实际目录结构依赖 MiniWDL 版本和
配置；下游自动化应读取 workflow output 映射，不要硬编码 `call-RunTxIDBatch` 路径。

### 6.2 Cromwell

使用本地 Cromwell JAR：

```bash
java -jar cromwell-92.jar run \
  workflows/txid_batch.wdl \
  --inputs my-txid-inputs.json
```

在平台页面提交时：

1. 上传 `workflows/txid_batch.wdl`。
2. 选择 `TxIDBatch` workflow。
3. 按参数表上传参考 FASTA、参考 GTF 和输入 GTF。
4. 确认每个输入 GTF 的文件名能够作为自动生成的 sample 标识。
5. 保留默认 Docker Hub 镜像，或替换成平台能够访问的镜像地址。

该 WDL 是单文件自包含工作流，不需要额外的 WDL import ZIP。

## 7. 输出说明

WDL 引擎可能复制、链接或重新本地化最终文件。应以 workflow output 映射返回的路径
为准，不要依赖 task 内部的相对路径。

### 7.1 工作流公开输出

| workflow output | WDL 类型 | 内容 |
|---|---|---|
| `registry` | `File` | `txid.sqlite`，包含身份、参考上下文、observations 和 provenance |
| `catalog` | `File` | `txid.catalog.tsv`，cohort 级精确 form 汇总 |
| `manifest` | `File` | 实际执行使用的 `import-manifest.tsv` |
| `initialization_summary` | `File` | `txid.init.json`，assembly/annotation fingerprint 和参考 transcript 数 |
| `batch_summary` | `File` | `txid.batch.json`，每个输入的 import ID、checksum、输出路径和 transcript 数 |
| `validation_summary` | `File` | `txid.validation.json`，完整性、外键和 digest 校验结果 |
| `rewritten_gtfs` | `Array[File]` | 每个输入对应一个重写后的 `.txid.gtf` |
| `mappings` | `Array[File]` | 每个输入对应一个 `.mapping.tsv` |

输出文件名格式为：

```text
<sample>.<tool>.<annotation_name>.<input-checksum前10位>.txid.gtf
<sample>.<tool>.<annotation_name>.<input-checksum前10位>.mapping.tsv
```

空格和不适合作为文件名的字符会转换为下划线。由于 batch 会确定性排序输入，且数组
输出来自 `glob()`，不要假设输出数组下标仍等于输入数组下标。应通过文件名、mapping
中的 `sample`/`tool` 字段或 `batch_summary` 建立对应关系。

### 7.2 重写 GTF

重写 GTF 使用确定性排序和 1-based closed 坐标。主要行为为：

- exact known transcript 保留参考 `gene_id` 和 `transcript_id`。
- novel transcript 若唯一归属于 known gene，保留该参考 `gene_id`，并使用 TxID form ID。
- genuinely new 或 unresolved locus 使用 registry 管理的 `txid:GL1.*` accession。
- multi-exon transcript 记录 `txid_sc` 和 `txid_form`。
- single-exon transcript 使用 `txid:SE1.*` form，不存在 `txid_sc`。
- 上游属性会尽量保留；与 TxID 保留字段冲突时使用 `txid_upstream_` 命名空间。
- fuzzy 模式下额外记录 `txid_fuzzy_cluster` 和 bridge status，精确 ID 仍保留。

输入 GTF 永远不会被原地修改。

### 7.3 Mapping TSV

每个 mapping 表包含：

| 列 | 说明 |
|---|---|
| `sample`、`tool` | 输入数组中提供的 provenance |
| `original_transcript_id`、`original_gene_id` | 上游原始标识符 |
| `txid_transcript_id`、`txid_gene_id` | 重写 GTF 使用的 transcript/gene ID |
| `txid_sc` | multi-exon splice-chain ID；single-exon 为空 |
| `txid_form` | exact transcript form ID，推荐作为跨样本矩阵 join key |
| `classification` | `known`、`novel_in_known_gene`、`ambiguous_gene` 或 `new_locus` |
| `annotation_name` | 本次参考注释上下文 |
| `gene_candidates` | ambiguous assignment 的候选 gene，逗号分隔 |
| `fuzzy_cluster` | fuzzy 模式的 `FC1` ID；精确模式为空 |
| `fuzzy_bridge_status` | `unambiguous`、`ambiguous_bridge` 或空 |

### 7.4 Cohort catalog

`catalog` 每行对应一个 exact form，包含：

- `form_id` 和可选的 `splice_chain_id`。
- 当前输出 gene assignment。
- 观察到的 classifications。
- observation、sample 和 tool 数量。

构建 sample-by-transcript matrix 时，应使用 `mapping.txid_form` 作为精确列键。不要用
`fuzzy_cluster` 替代精确列，除非分析明确需要、并完整记录了 tolerance。

## 8. 运行成功后的快速检查

首先检查 validation JSON：

```bash
jq . /path/to/txid.validation.json
```

成功运行应满足：

```json
{
  "issues": [],
  "valid": true
}
```

`summary` 还会报告 assembly、annotations、imports、observations、exact forms、
splice chains、loci 和 classification 数量。

继续检查：

```bash
test "$(find /path/to/outputs -name '*.txid.gtf' | wc -l)" -eq 2
test "$(find /path/to/outputs -name '*.mapping.tsv' | wc -l)" -eq 2
head -n 3 /path/to/txid.catalog.tsv
head -n 3 /path/to/sample01.mapping.tsv
```

其中数字 `2` 应替换为输入 GTF 数量。还应确认：

- `batch_summary.imports` 数量等于输入 GTF 数量。
- 每个 import 的 `transcripts` 大于 0。
- `validation_summary.valid` 为 `true` 且 `issues` 为空。
- known transcript 的参考 ID 得到保留。
- 不同样本中的相同结构获得相同 exact TxID。
- 同一 splice chain 的不同 transcript ends 具有相同 `SC1`、不同 `TF1`。

## 9. 常见问题与排错

### 9.1 平台提示不支持数组后的 `+`

当前 WDL 使用普通 `Array[File]`，没有 `Array[File]+` 语法。若平台
仍显示旧声明，重新上传最新的 `txid_batch.wdl`，并确认缓存的 workflow 版本已更新。

如果日志同时出现 `mapfile: not found`、`Bad substitution` 或 `[[: not found`，说明平台
仍在运行早期的 Bash 专用版本。当前 WDL 的 task command 已改为 POSIX `/bin/sh`
兼容写法，不再使用 `set -eu`、`mapfile`、Bash 数组、`[[ ... ]]` 或 `pipefail`；命令
失败由显式退出状态检查处理。重新上传 WDL 后再提交任务。

### 9.2 `gtfs must contain at least one input GTF`

`gtfs` 是空数组。至少上传一个 GTF。

### 9.3 自动生成的 sample 名称不符合预期

WDL 用本地化后的 GTF 文件 basename 自动生成 `sample`，并去掉 `.gtf` 或 `.gtf.gz`
后缀。请在提交前把输入文件改成可识别且互不重复的名称。如果平台本地化时重写文件
名，或需要精确的 sample/tool provenance，请改用 TxID CLI manifest 接口。

### 9.4 参考 contig 不存在

典型错误包含：

```text
contig '1' is absent from reference 'GRCh38'
```

确认 FASTA 和 GTF assembly 一致。若只是经过验证的 contig 命名差异，提供
`contig_aliases`；不要未经验证地批量改名。

### 9.5 transcript 坐标超过 contig 长度

这通常说明 FASTA 与 GTF assembly 不一致、使用了不同 patch/haplotype 集合，或 GTF
损坏。TxID 会拒绝输入，不能通过放宽 tolerance 绕过 assembly 验证。

### 9.6 GTF 没有可解析 transcript

确认文件有 9 个 tab 分隔字段，存在 `exon` feature，exon attributes 中包含
`transcript_id`。当前 WDL 强制按 GTF 解析，不能直接传 GFF3。

### 9.7 `duplicate GTF attribute key 'tag'`

这是旧版 TxID 0.1.0 的限制。GENCODE 等注释会在同一 feature 上合法地重复 `tag` 或
`ont`。TxID 0.1.2 会保留这些多值属性，不需要预处理或删除它们。确认 WDL 使用
`dawang02/txid:0.1.2-jupyter` 或对应 digest。重复的 `gene_id`、`transcript_id` 仍会
被拒绝，因为它们会使 feature 身份或父级关系产生歧义。

### 9.8 fuzzy 参数只设置了一个

工作流会报：

```text
fuzzy mode requires both fuzzy tolerances
```

同时提供 splice 和 end tolerance，或同时删除这两个字段。负数也会被拒绝。

### 9.9 镜像拉取失败或 `txid` 命令不存在

确认计算节点可以访问 Docker Hub、镜像地址拼写正确、节点为 amd64，并且平台没有
覆盖 `docker_image`。可在可访问 Docker 的节点先运行：

```bash
docker run --rm dawang02/txid:0.1.2-jupyter txid --version
```

### 9.10 `database or disk is full`

增加 `disk_gb` 或清理 WDL 后端运行分区。注意某些本地后端会忽略 WDL `disks` 请求，
此时需要调整后端工作目录所在文件系统，而不是只修改 JSON。

### 9.11 内存不足或 exit 137

增加 `memory_gb`。参考 FASTA 在 fingerprint 计算期间需要读取序列内容；全基因组 FASTA
和大型 GTF 比小型测试数据需要更多内存。

### 9.12 某个后期输入失败后已有 import

`txid batch` 对每个输入分别使用原子 transaction。后期输入失败不会破坏之前已提交的
import，但整个 WDL task 会以非零状态结束。WDL 引擎重试通常会使用新的 task 工作目录
和新 registry；如果手工在同一目录重试，应先理解现有 registry 状态，不要覆盖输入。

## 10. Exact identity 与 fuzzy grouping

TxID 的 exact ID 是默认、永久的结构身份：

- `txid:SC1.*`：assembly、contig、strand 和完整 intron chain。
- `txid:TF1.*`：multi-exon splice chain 加 TSS/TES。
- `txid:SE1.*`：single-exon start/end form。

Fuzzy grouping 是独立的可选分析层：

- `txid:FC1.*` 是给定 tolerance 下的 versioned cluster。
- 两个不同 exact ID 不会因为属于同一 fuzzy cluster 而被别名化或覆盖。
- 比较不同运行的 fuzzy 结果时，必须同时比较 splice/end tolerance 和软件版本。

`novel` 不出现在永久结构 ID 中。Novelty 是相对于 `annotation_name` 对应参考注释的
classification；同一结构换一个注释版本可能由 novel 变为 known，但 exact TxID 不变。

## 11. 重跑、增量导入与可复现性

每次 `TxIDBatch` WDL 运行都会创建一个新的 `txid.sqlite`，不会接收已有 registry 作为
输入。因此：

- 一次 cohort 分析应把希望共同汇总的 GTF 全部放入同一个输入 JSON。
- 需要在既有 registry 上持续增量导入时，应使用 `txid add`/`txid batch` CLI，而不是
  重新提交本 WDL。
- 修改参考 FASTA 会改变 assembly fingerprint 和 exact ID 身份上下文。
- 修改参考 GTF 可能改变 known/novel classification 和保留的参考 aliases，但在同一
  FASTA assembly 上不会改变相同结构的 exact TxID。
- 精确 ID 不依赖 sample 名称、tool 名称、输入顺序、数据库行号或时间戳。

建议归档以下内容：

- WDL 文件和完整输入 JSON。
- Docker image digest。
- 参考 FASTA、参考 GTF、contig alias 表及其 checksum。
- `initialization_summary`、`batch_summary` 和 `validation_summary`。
- SQLite registry、mapping、重写 GTF 和 catalog。

## 12. 提交前检查清单

提交前逐项确认：

- [ ] WDL 文件是 1.1 版本，数组声明不含 `+`。
- [ ] `gtfs` 非空，每个输入 GTF 的文件名可识别且尽量唯一。
- [ ] 所有输入都是 GTF/GTF.GZ，不是 GFF3。
- [ ] FASTA、参考 GTF 和全部输入 GTF 属于同一 assembly。
- [ ] contig 名称一致，或提供经过验证的 alias 表。
- [ ] `annotation_name` 能明确标识参考注释版本。
- [ ] fuzzy tolerance 要么全部省略，要么两个都设置为非负整数。
- [ ] Docker 镜像可从计算节点拉取，生产运行优先使用 digest。
- [ ] 内存和磁盘资源与 FASTA、GTF 总规模匹配。
- [ ] 运行结束后检查 `validation_summary.valid == true`。
