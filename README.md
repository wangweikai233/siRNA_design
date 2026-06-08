# siRNA 自动设计工具

基于 DSIR 预测 + 多标准筛选 + BLAST 验证的自动化 siRNA 设计流程。

## 流程概览

```
① NCBI获取CDS → ② DSIR 19nt评分 → ③ 多标准筛选 → ④ BLAST特异性验证 → ⑤ .docx报告
```

## 依赖

```bash
pip install python-docx
```

Python 标准库即可运行（`re`, `json`, `argparse`, `urllib`），`python-docx` 仅用于生成 `.docx` 报告，缺失时会自动回退为 `.txt`。

## 用法

```bash
# 按基因名自动搜索（推荐）
python sirna_design.py --gene SUN1

# 指定转录本（结果可精确复现）
python sirna_design.py --transcript NM_001130965.3

# 使用自定义CDS序列文件（FASTA格式）
python sirna_design.py --cds my_sequence.fasta

# 指定输出目录
python sirna_design.py --gene SUN1 --output ./results

# 安静模式（只输出最终结果）
python sirna_design.py --gene SUN1 --quiet

# 指定物种
python sirna_design.py --gene ENSA --species "Rattus norvegicus"
```

## 参数

| 参数 | 说明 |
|------|------|
| `--gene` | 基因名（与 `--transcript` / `--cds` 三选一） |
| `--transcript` | RefSeq 转录本 ID（如 NM_001130965.3） |
| `--cds` | FASTA 格式 CDS 文件路径 |
| `--species` | 物种（默认 Homo sapiens） |
| `--output` / `-o` | 输出目录（默认当前目录） |
| `--quiet` / `-q` | 安静模式 |

## 设计标准

严格遵循 2025.7.6 siRNA 设计流程文档：

| 标准 | 要求 |
|------|------|
| DSIR 模式 | 19nt, NA(Nn)NN |
| DSIR Score | ≥ 90（不够可降至 ≥ 50） |
| Corrected Score | ≥ 90 优先（不够放宽） |
| 位置 | > 100aa（从起始密码子算起） |
| 第 2 位碱基 | A |
| GC 含量 | 8-11 / 19 (38-52%) |
| F1 (sense 5'端) | G 或 C |
| R1 (antisense 5'端) | A 或 U |
| GC/CG 连续二核苷酸 | 不允许 ≥ 2 对连续 |
| 连续相同碱基 | 不允许 ≥ 4 个 |
| UAUAUA / AUAUAU | 不允许 |
| BLAST 脱靶 | 人类转录组中无显著脱靶 |
| 两条 siRNA 间距 | > 25bp |

## 输出

- 终端输出完整设计流程和订购序列
- `{基因名}_siRNA设计报告.docx` — 包含序列信息、筛选标准对照、BLAST 结果、订购表
- 如未安装 `python-docx`，回退为 `.txt` 格式

## 注意

- 需要网络连接（NCBI 获取 CDS、DSIR 评分、BLAST 验证均需联网）
- DSIR 网站为 http://biodev.extra.cea.fr/DSIR/DSIR.html，如网站不可用则流程中断
- BLAST 使用 NCBI API，19nt 短序列查询的 E 值偏高属正常现象，程序已通过物种限制和 E 值阈值过滤
- 不同转录本可能产生相同 siRNA 序列（靶向区域在变异体中保守）
