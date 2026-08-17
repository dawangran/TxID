# TxID `multi-add` 多 GTF WDL

[`txid_multi_add.wdl`](../workflows/txid_multi_add.wdl) 接收一个参考基因组
FASTA、一个参考注释 GTF 和多个待处理 GTF，在同一个 task 中直接运行
`txid multi-add`。工作流不创建或接收 manifest。

## 输入

```text
reference_fasta + reference_gtf + gtfs[]
                         │
                         ▼
                     txid init
                         │
                         ▼
                  txid multi-add
                         │
            ┌────────────┼────────────┐
            ▼            ▼            ▼
      rewritten GTF  mapping TSV  SQLite registry
                                      │
                               export + validate
```

主要参数：

| 参数 | WDL 类型 | 说明 |
|---|---|---|
| `reference_fasta` | `File` | 用于 assembly fingerprint 和坐标校验的参考 FASTA |
| `reference_gtf` | `File` | known transcript/gene 参考注释 |
| `gtfs` | `Array[File]` | 一个或多个待赋予 TxID 的 GTF |
| `assembly` | `String` | assembly 名称，例如 `GRCh38` |
| `annotation_name` | `String` | 参考注释名称，默认 `reference` |
| `tool` | `String` | 所有输入共用的上游工具名称，默认 `unspecified` |
| `contig_aliases` | `File?` | 可选的两列 contig alias 表 |
| `fuzzy_splice_tolerance_bp` | `Int?` | 可选 splice 容差 |
| `fuzzy_end_tolerance_bp` | `Int?` | 可选 ends 容差，必须与 splice 容差同时设置 |
| `docker_image` | `String` | 运行镜像，默认 `dawang02/txid:0.1.3-jupyter` |

样本名由每个输入 GTF 的 basename 自动推导，例如 `donor01.gtf` 对应
`donor01`。所有输入必须属于同一 assembly、同一参考注释上下文，并共享这里记录的
`tool`。不同工具或不同注释版本应拆分运行，或使用 `txid batch` manifest 接口。

## 运行

先从当前代码构建包含 `multi-add` 的镜像：

```bash
docker build -t dawang02/txid:0.1.3-jupyter .
```

修改示例输入文件
[`txid_multi_add.inputs.example.json`](../workflows/txid_multi_add.inputs.example.json)，
然后运行：

```bash
miniwdl check workflows/txid_multi_add.wdl

miniwdl run \
  workflows/txid_multi_add.wdl \
  -i workflows/txid_multi_add.inputs.example.json
```

在云端 WDL 平台运行时，把 `docker_image` 替换为已经推送到平台可访问 registry
的不可变镜像标签或 digest。

## 输出

- `registry`：版本化 SQLite registry。
- `catalog`：确定性 cohort catalog。
- `rewritten_gtfs`：每个输入对应的 TxID 重写 GTF。
- `mappings`：每个输入对应的 mapping TSV。
- `initialization_summary`、`multi_add_summary`、`export_summary`：各阶段 JSON。
- `validation_summary`：registry 完整性和 digest 校验 JSON。

工作流不使用 scatter，因为多个隔离 task 不能安全地共同修改一个 SQLite registry。
输入会在一个 task 中通过优化后的 multi-file 路径按确定性顺序处理。
